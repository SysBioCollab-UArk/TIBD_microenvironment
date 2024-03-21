from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from Lemaire2004_reduced import model

alias_model_components()

# start with zero cells and run to equilibrium
R_0.value = 0.0
B_0.value = 0.0
C_0.value = 0.0

# estimate rate constant for OB -> OB + Bone to get bone density ~ 100
mean_OB = 1.06E-02  # from data
mean_OC = 9.16E-04  # from data
k_B_builds_bone.value = k_C_consumes_bone.value * mean_OC / mean_OB * 100 * 2

# run simulation
tspan = np.linspace(0, 500, 5001)  # 100, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

# plot time courses
obs_to_plot = [['OB_tot', 'C_obs'], ['Bone_tot']]
ref_vals = [[mean_OB, mean_OC], [None]]
leg_labels = [['OB', 'OC'], ['Bone']]
colors = [['b', 'r'], ['g']]
ylabels = ['concentration (pM)', 'relative BV/TV']
for obs, ref, leg_label, color, ylabel in zip(obs_to_plot, ref_vals, leg_labels, colors, ylabels):
    plt.figure()
    for i in range(len(obs)):
        plt.plot(tspan, output.observables[obs[i]], lw=2, color=color[i], label=leg_label[i])
        if ref[i] is not None:
            plt.plot(tspan[-1], ref[i], '*', ms=12, color=color[i])
    plt.xlabel('time (d)')
    plt.ylabel(ylabel)
    plt.legend(loc=0)
    plt.tight_layout()

plt.show()
