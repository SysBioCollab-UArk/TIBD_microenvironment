from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from Lemaire2004 import model
from sympy import Piecewise, sympify
from util import get_exp_data
import os

# experimental data
exp_data = get_exp_data(os.path.join('DATA', 'Mouse_Data_March2024.csv'))

# ########## Bone ##########
Monomer('Bone')
Parameter('Bone_0', 100)  # percentage

# Make OB bone production (OC bone consumption) rate increase (decrease) with decreasing (increasing) concentration
mean_OB = 1.06E-02 * 1000  # fM, from data
mean_OC = 9.16E-04 * 1000  # fM, from data

# === MODEL #1 ===
# '''
Parameter('k1_B_builds_bone', 0)  # /day 1.25e8
Parameter('k2_B_builds_bone', 100)  # fM
Parameter('nB', 1)  # 1.1
Parameter('k1_C_consumes_bone', 1e6)  # /day
Parameter('k2_C_consumes_bone', 1)  # fM
Parameter('nC', 1)  # 0.8
alias_model_components()
Expression('k_B_builds_bone',
           Piecewise((0, B_obs <= 0),
                     (k1_B_builds_bone * (k2_B_builds_bone ** nB + B_obs ** nB) / B_obs ** nB, True)))
Expression('k_C_consumes_bone',
           Piecewise((0, C_obs <= 0),
                     (k1_C_consumes_bone * (k2_C_consumes_bone ** nC + C_obs ** nC) / C_obs ** nC, True)))

# estimate rate constant for OB -> OB + Bone to get bone density ~ 100
k1C = k1_C_consumes_bone.value
k2C = k2_C_consumes_bone.value
k2B = k2_B_builds_bone.value
k1_B_builds_bone.value = k1C * (k2C ** nC.value + mean_OC ** nC.value) / mean_OC ** (nC.value-1) * Bone_0.value * \
                         mean_OB ** (nB.value - 1) / (k2B ** nB.value + mean_OB ** nB.value) * \
                         1.33  # 1.8
kB_label = r'$k_{1B} \cdot \frac{k_{2B}^{n_B} + [B]^{n_B}}{[B]^{n_B}}$'
kC_label = r'$k_{1C} \cdot \frac{k_{2C}^{n_C} + [C]^{n_C}}{[C]^{n_C}}$'
# '''
# === MODEL #2 ===
'''
Parameter('k1_B_builds_bone', 2.3e6)  # /day
Parameter('k2_B_builds_bone', 0.4)  # fM 0.8
Parameter('nB', 1.7)  # 2
Parameter('kB_star', 1.4e9)  # /day
Parameter('k1_C_consumes_bone', 1e6)  # /day
Parameter('k2_C_consumes_bone', 1)  # fM
Parameter('nC', 1)  # 0.8
Parameter('kC_star', 2e6)  # /day
alias_model_components()
Expression('k_B_builds_bone',
           k1_B_builds_bone + kB_star * k2_B_builds_bone ** nB / (k2_B_builds_bone ** nB + B_obs ** nB))
Expression('k_C_consumes_bone',
           k1_C_consumes_bone + kC_star * k2_C_consumes_bone ** nC / (k2_C_consumes_bone ** nC + C_obs ** nC))

# estimate rate constant for OB -> OB + Bone to get bone density ~ 100
k1C = k1_C_consumes_bone.value
k2C = k2_C_consumes_bone.value
k2B = k2_B_builds_bone.value
k1_B_builds_bone.value = ((k1C + kC_star.value * k2C ** nC.value / (k2C ** nC.value + mean_OC ** nC.value)) * mean_OC * \
                         Bone_0.value / mean_OB - kB_star.value * k2B ** nB.value / (k2B ** nB.value + mean_OB ** nB.value)) * \
                         1.93
kB_label = r'$k_{1B} + k_B^\ast \cdot \frac{k_{2B}^{n_B}}{k_{2B}^{n_B} + [B]^{n_B}}$'
kC_label = r'$k_{1C} + k_C^\ast \cdot \frac{k_{2C}^{n_C}}{k_{2C}^{n_C} + [C]^{n_C}}$'
'''

alias_model_components()
Initial(Bone(), Bone_0)
Observable('Bone_tot', Bone())
Rule('B_builds_bone', B() >> B() + Bone(), k_B_builds_bone)
Rule('C_consumes_bone', C() + Bone() >> C(), k_C_consumes_bone)

# ######### Tumor #########
Monomer('Tumor')
Parameter('Tumor_0', 0)  # fM (1000 fM = 1 pM)

# standard exponential growth
Parameter('k_tumor_div_basal', 1)  # 1/day
Parameter('k_tumor_dth', 0.6)  # 1/day

# carrying capacity (logistic growth)
Parameter('CC_ON', 1)
Parameter('N', 600)  # fM

# allee effect
Parameter('ALLEE_ON', 0)
Parameter('A', 1.5)  # 1/fM

# growth enhancement due to TGF-beta
Parameter('k_tumor_div_PTHrP', 0.1)  # /day
alias_model_components()
Expression('k_tumor_div', k_tumor_div_basal + k_tumor_div_PTHrP * pi_C)

Observable('Tumor_tot', Tumor())
alias_model_components()

Expression('k_tumor_cc', CC_ON * (k_tumor_div - k_tumor_dth) / N)
Expression('k_tumor_allee',
           Piecewise((0, Tumor_tot <= 0),
                     (ALLEE_ON * (k_tumor_div - k_tumor_dth) / (A * Tumor_tot), True)))
alias_model_components()

Initial(Tumor(), Tumor_0)
Rule('tumor_division', Tumor() >> Tumor() + Tumor(), k_tumor_div)
Rule('tumor_death', Tumor() >> None, k_tumor_dth)
Rule('tumor_cc', Tumor() + Tumor() >> Tumor(), k_tumor_cc)
Rule('tumor_allee', Tumor() >> None, k_tumor_allee)

# Make PTH (i.e., PTHrP) synthesis rate a function of the number of tumor cells
Parameter('k_tumor_PTHrP', 5e2)  # 1/day
alias_model_components()
IP.expr = sympify(k_tumor_PTHrP * pi_C * Tumor_tot)  # PTHrP expression depends on the TGF-beta occupancy pi_C

# Hypothesize that tumor cells secrete a factor that promotes OB death
Parameter('k_tumor_OB', 0.002)  # 1/fM-day  # 0.01
alias_model_components()
Rule('AOB_death_tumor', Tumor() + B() >> Tumor(), k_tumor_OB)

# Hypothesize that tumor cells secrete a factor that promotes OC creation
Parameter('k_tumor_OC', 0.004)  # 1/fM-day  # 0.01
alias_model_components()
Rule('AOC_creation_tumor', Tumor() >> Tumor() + C(), k_tumor_OC)

# ######### BISPHOSPHONATES #########

Monomer('Bisphos')
Parameter('Bisphos_0', 0)  # fM
Parameter('k_bisphos_AOC', 1)  # 1/fM-day
alias_model_components()
Initial(Bisphos(), Bisphos_0)
Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)

# run simulation(s)

# start with zero cells and run to equilibrium
R_0.value = 0.0
B_0.value = 0.0
C_0.value = 0.0

# modify some constants
DA.value = 2  # /day
kB.value = 0.022  # /day

# create simulator object
tspans = []
outputs = []
linestyles = []
sim = ScipyOdeSimulator(model, verbose=True)

# equilibration
tspans.append(np.linspace(-500, 0, 5001))  # 100, 1001)
outputs.append(sim.run(tspan=tspans[-1]))
linestyles.append('-')

# Add a small amount of tumor cells
tspans.append(np.linspace(0, 30, 31))  # days
initials = outputs[-1].species[-1]
idx_tumor = [str(sp) for sp in model.species].index('Tumor()')  # get index of Tumor species
initials[idx_tumor] = 1  # fM
outputs.append(sim.run(tspan=tspans[-1], initials=initials))
linestyles.append('-')

# Add bisphosphonate (Zoledronic acid) at day 6
tspans.append(np.linspace(6, 30, 241))  # days
initials = outputs[-1].species[6]
idx_bisphos = [str(sp) for sp in model.species].index('Bisphos()')  # get index of Bisphos species
initials[idx_bisphos] = 1  # fM
outputs.append(sim.run(tspan=tspans[-1], initials=initials))
linestyles.append('--')

# plot time courses
obs_to_plot = [['OB_tot'], ['C_obs'], ['Bone_tot'], ['Tumor_tot']]
refs = [[exp_data[('Tumor', 'OB')]], [exp_data[('Tumor', 'OC')]], [exp_data[('Tumor', 'Bone')]],
        [exp_data[('Tumor', 'Tumor')]]]
leg_labels = [['OB'], ['OC'], ['Bone'], ['Tumor']]
colors = [['b'], ['r'], ['g'], ['orange']]
ylabels = ['concentration (pM)', 'concentration (pM)', 'relative BV/TV', 'concentration (pM)']
scales = [1/1000, 1/1000, 1/Bone_0.value, 1/1000]
for obs, ref, leg_label, color, ylabel, scale in zip(obs_to_plot, refs, leg_labels, colors, ylabels, scales):
    print(obs)
    plt.figure()
    for i in range(len(obs)):
        add_label = True
        for tspan, output, ls in zip(tspans, outputs, linestyles):
            label = leg_label[i] if add_label else None
            plt.plot(tspan/7, output.observables[obs[i]] * scale, lw=2, color=color[i], ls=ls, label=label)
            add_label = False
        if obs[i] == 'Bone_tot':
            plt.ylim(top=1.65, bottom=-0.05)  # top=1.2
        # plot reference (experimental) data
        if ref[i] is not None:
            plt.errorbar(*ref[i], ecolor='k', ls='None', marker='o', ms=8,  capsize=10, color=color[i],
                         label='experiment')
    plt.xlabel('time (week)')
    plt.ylabel(ylabel)
    plt.xlim(left=0, right=4)
    plt.legend(loc=0)
    plt.tight_layout()

# plot OB bone production and OC bone consumption rates over time
plt.figure()
plt.plot(tspans[-2]/7, outputs[-2].expressions['k_B_builds_bone'], lw=2, color='b', label=kB_label)
plt.plot(tspans[-1]/7, outputs[-1].expressions['k_B_builds_bone'], lw=2, color='b', ls='--')
plt.xlabel('time (week)')
plt.ylabel('value (/day)')
plt.legend(loc=0, fontsize=14)
plt.tight_layout()

plt.figure()
plt.plot(tspans[-2]/7, outputs[-2].expressions['k_C_consumes_bone'], lw=2, color='r', label=kC_label)
plt.plot(tspans[-1]/7, outputs[-1].expressions['k_C_consumes_bone'], lw=2, color='r', ls='--')
plt.xlabel('time (week)')
plt.ylabel('value (/fM-day)')
plt.legend(loc=0, fontsize=14)
plt.tight_layout()

plt.show()
