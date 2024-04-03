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
Parameter('k_B_builds_bone', 1.25e8)  # /day
Parameter('k_C_consumes_bone', 1e6)  # /day

alias_model_components()

Initial(Bone(), Bone_0)
Observable('Bone_tot', Bone())
Rule('B_builds_bone', B() >> B() + Bone(), k_B_builds_bone)
Rule('C_consumes_bone', C() + Bone() >> C(), k_C_consumes_bone)

# ######### Tumor #########

Monomer('Tumor')
Parameter('Tumor_0', 0)  # fM (1000 fM = 1 pM)
# standard exponential growth
Parameter('k_tumor_div', 1)  # 1/day
Parameter('k_tumor_dth', 0.6)  # 1/day
# carrying capacity (logistic growth)
Parameter('CC_ON', 1)
Parameter('N', 600)  # fM
# allee effect
Parameter('ALLEE_ON', 1)
Parameter('A', 1.5)  # 1/fM

alias_model_components()

Observable('Tumor_tot', Tumor())

alias_model_components()

Expression('k_tumor_cc', CC_ON * (k_tumor_div - k_tumor_dth) / N)
Expression('k_tumor_allee', Piecewise((0, Tumor_tot <= 0),
                                      (ALLEE_ON * (k_tumor_div - k_tumor_dth) / (A * Tumor_tot), True)))

alias_model_components()

Initial(Tumor(), Tumor_0)
Rule('tumor_division', Tumor() >> Tumor() + Tumor(), k_tumor_div)
Rule('tumor_death', Tumor() >> None, k_tumor_dth)
Rule('tumor_cc', Tumor() + Tumor() >> Tumor(), k_tumor_cc)
Rule('tumor_allee', Tumor() >> None, k_tumor_allee)

# TODO: Add PTH (i.e., PTHrP) synthesis rate as a function of the number of tumor cells

Parameter('k_tumor_PTHrP', 5e2)
alias_model_components()
IP.expr = sympify(k_tumor_PTHrP * Tumor_tot)

# TODO: Add rules for tumor growth enhancement by TGF-beta (osteoclasts).
#  QUESTION: Does TGF-beta enhance the division rate, increase the rate of PTHrP production from tumor cells, or both?

# ...

# #########################

# start with zero cells and run to equilibrium
R_0.value = 0.0
B_0.value = 0.0
C_0.value = 0.0

# modify some constants
DA.value = 2  # /day
kB.value = 0.022  # /day

# estimate rate constant for OB -> OB + Bone to get bone density ~ 100
mean_OB = 1.06E-02 * 1000  # from data
mean_OC = 9.16E-04 * 1000  # from data
k_B_builds_bone.value = k_C_consumes_bone.value * mean_OC / mean_OB * 100 * 2

# run simulation(s)
tspans = []
outputs = []
sim = ScipyOdeSimulator(model, verbose=True)

# equilibration
tspans.append(np.linspace(-500, 0, 5001))  # 100, 1001)
outputs.append(sim.run(tspan=tspans[-1]))

# inject PTH (or PTHrP from tumor cells)
tspans.append(np.linspace(0, 30, 301))
initials = outputs[-1].species[-1]
idx = [str(sp) for sp in model.species].index('Tumor()')
initials[idx] = 1  # fM
outputs.append(sim.run(tspan=tspans[-1], initials=initials))

# plot time courses
obs_to_plot = [['OB_tot'], ['C_obs'], ['Bone_tot'], ['Tumor_tot']]
refs = [[exp_data[('Tumor', 'OB')]], [exp_data[('Tumor', 'OC')]], [exp_data[('Tumor', 'Bone')]],
        [exp_data[('Tumor', 'Tumor')]]]
leg_labels = [['OB'], ['OC'], ['Bone'], ['Tumor']]
colors = [['b'], ['r'], ['g'], ['orange']]
ylabels = ['concentration (pM)', 'concentration (pM)', 'relative BV/TV', 'concentration (pM)']
scales = [1/1000, 1/1000, 1/100, 1/1000]
for obs, ref, leg_label, color, ylabel, scale in zip(obs_to_plot, refs, leg_labels, colors, ylabels, scales):
    print(obs)
    plt.figure()
    for i in range(len(obs)):
        add_label = True
        for tspan, output in zip(tspans, outputs):
            label = leg_label[i] if add_label else None
            plt.plot(tspan/7, output.observables[obs[i]] * scale, lw=2, color=color[i], label=label)
            add_label = False
        # plot reference (experimental) data
        if ref[i] is not None:
            plt.errorbar(*ref[i], ecolor='k', ls='None', marker='o', ms=8,  capsize=10, color=color[i],
                         label='experiment')
    plt.xlabel('time (week)')
    plt.ylabel(ylabel)
    plt.xlim(left=0, right=4)
    plt.legend(loc=0)
    plt.tight_layout()

plt.show()
