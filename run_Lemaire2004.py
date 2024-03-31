from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from Lemaire2004 import model

alias_model_components()

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
tspans.append(np.linspace(0, 500, 5001))  # 100, 1001)
outputs.append(sim.run(tspan=tspans[-1]))

# inject PTH (or PTHrP from tumor cells)
tspans.append(np.linspace(500, 1000, 5001))
initials = outputs[-1].species[-1]
outputs.append(sim.run(tspan=tspans[-1], initials=initials, param_values={'IP': 1e3 * 1000}))

# plot time courses
obs_to_plot = [['OB_tot', 'C_obs'], ['Bone_tot']]
# ref_vals = [[mean_OB, mean_OC], [None]]
ref_vals = [[None, None], [None]]
leg_labels = [['OB', 'OC'], ['Bone']]
colors = [['b', 'r'], ['g']]
ylabels = ['concentration (pM)', 'relative BV/TV']
scales = [1/1000, 1]
for obs, ref, leg_label, color, ylabel, scale in zip(obs_to_plot, ref_vals, leg_labels, colors, ylabels, scales):
    plt.figure()
    for i in range(len(obs)):
        add_label = True
        for tspan, output in zip(tspans, outputs):
            label = leg_label[i] if add_label else None
            plt.plot(tspan, output.observables[obs[i]] * scale, lw=2, color=color[i], label=label)
            add_label = False
        if ref[i] is not None:
            plt.plot(tspan[-1], ref[i], '*', ms=12, color=color[i])
    plt.xlabel('time (d)')
    plt.ylabel(ylabel)
    plt.legend(loc=0)
    plt.tight_layout()

plt.show()
