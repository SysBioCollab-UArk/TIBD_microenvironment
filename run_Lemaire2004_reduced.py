from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
import numpy as np
import matplotlib.pyplot as plt
from Lemaire2004_reduced import model

alias_model_components()

R_0.value = 0.0
B_0.value = 0.0
C_0.value = 0.0

Observable('OB_tot', R()+B())

tspan = np.linspace(0, 100, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

obs_to_plot = ['OB_tot', 'C_obs']
colors = ['b', 'r']  # , 'g']
labels = ['OB', 'OC']  # ['ROB', 'AOB', 'AOC']
for obs, color, label in zip(obs_to_plot, colors, labels):
    plt.plot(tspan, output.observables[obs], lw=2, color=color, label=label)
plt.xlabel('time (d)')
plt.ylabel('concentration (pM)')
plt.legend(loc=0)
plt.tight_layout()

plt.show()
