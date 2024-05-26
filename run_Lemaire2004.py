from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from util import get_exp_data
import os
from Lemaire2004 import create_model_elements as create_Lemaire_ME
from Harris2024 import create_model_elements as create_Harris_ME

# experimental data
exp_data = get_exp_data(os.path.join('DATA', 'Mouse_Data_May2024.csv'))

Model()

# create the model
create_Lemaire_ME()
OB_OC_model = 1
create_Harris_ME(OB_OC_model)
alias_model_components()

if OB_OC_model == 1:
    kB_label = r'$k_{1B} \cdot \frac{k_{2B}^{n_B} + [B]^{n_B}}{[B]^{n_B}}$'
    kC_label = r'$k_{1C} \cdot \frac{k_{2C}^{n_C} + [C]^{n_C}}{[C]^{n_C}}$'
elif OB_OC_model == 2:
    kB_label = r'$k_{1B} + k_B^\ast \cdot \frac{k_{2B}^{n_B}}{k_{2B}^{n_B} + [B]^{n_B}}$'
    kC_label = r'$k_{1C} + k_C^\ast \cdot \frac{k_{2C}^{n_C}}{k_{2C}^{n_C} + [C]^{n_C}}$'

# run simulation(s)
ADD_BISPHOSPHONATE = True
SAVE_FIGS = False

# start with zero cells and run to equilibrium
R_0.value = 0.0
B_0.value = 0.0
C_0.value = 0.0

# modify some constants
# DA.value = 0.7  # 2  # /day
kB.value = 0.013  # 0.022  # /day

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
if ADD_BISPHOSPHONATE:
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
            top = 2.1 if ADD_BISPHOSPHONATE else 1.2
            plt.ylim(top=top, bottom=-0.05)  # top=1.2 2.1
        # plot reference (experimental) data
        if ref[i] is not None:
            plt.errorbar(*ref[i], ecolor='k', ls='None', marker='o', ms=8,  capsize=10, color=color[i],
                         label='experiment')
    plt.xlabel('time (week)')
    plt.ylabel(ylabel)
    plt.xlim(left=0, right=4)
    plt.legend(loc=0)
    plt.tight_layout()

    # Save plot as PDF
    if SAVE_FIGS:
        outfilename = 'FIG'
        for label in leg_label:
            outfilename += '_%s' % label
        outfilename += '_ZA.pdf' if ADD_BISPHOSPHONATE else '.pdf'
        plt.savefig(outfilename, format='pdf')

# plot OB bone production and OC bone consumption rates over time
plt.figure()
if ADD_BISPHOSPHONATE:
    plt.plot(tspans[-2] / 7, outputs[-2].expressions['k_B_builds_bone'], lw=2, color='b', label=kB_label)
    plt.plot(tspans[-1] / 7, outputs[-1].expressions['k_B_builds_bone'], lw=2, color='b', ls='--')
else:
    plt.plot(tspans[-1] / 7, outputs[-1].expressions['k_B_builds_bone'], lw=2, color='b', label=kB_label)
plt.xlabel('time (week)')
plt.ylabel('value (/day)')
plt.legend(loc=0, fontsize=14)
plt.tight_layout()

plt.figure()
if ADD_BISPHOSPHONATE:
    plt.plot(tspans[-2] / 7, outputs[-2].expressions['k_C_consumes_bone'], lw=2, color='r', label=kC_label)
    plt.plot(tspans[-1] / 7, outputs[-1].expressions['k_C_consumes_bone'], lw=2, color='r', ls='--')
else:
    plt.plot(tspans[-1] / 7, outputs[-1].expressions['k_C_consumes_bone'], lw=2, color='r', label=kC_label)
plt.xlabel('time (week)')
plt.ylabel('value (/fM-day)')
plt.legend(loc=0, fontsize=14)
plt.tight_layout()

plt.show()
