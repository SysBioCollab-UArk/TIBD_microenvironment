from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from util import get_exp_data
import os

Model()

Monomer('Tumor')
Parameter('Tumor_0', 1)  # fM (1000 fM = 1 pM)
Initial(Tumor(), Tumor_0)
Observable('Tumor_tot', Tumor())

# standard exponential growth
Parameter('k_tumor_div', 1)  # 1/day
Parameter('k_tumor_dth', 0.6)  # 1/day
Rule('tumor_division', Tumor() >> Tumor() + Tumor(), k_tumor_div)
Rule('tumor_death', Tumor() >> None, k_tumor_dth)

# carrying capacity (logistic growth)
Parameter('CC_ON', 1)
Parameter('N', 600)  # fM
Expression('k_tumor_cc', CC_ON * (k_tumor_div - k_tumor_dth) / N)
Rule('tumor_cc', Tumor() + Tumor() >> Tumor(), k_tumor_cc)

# allee effect
Parameter('ALLEE_ON', 1)
Parameter('A', 1.5)  # 1/fM
Expression('k_tumor_allee', ALLEE_ON * (k_tumor_div - k_tumor_dth) / (A * Tumor_tot))
Rule('tumor_allee', Tumor() >> None, k_tumor_allee)

# simulations
sim = ScipyOdeSimulator(model, verbose=True)

plt.figure(constrained_layout=True)

# exponential growth
tspan = np.linspace(0, 20, 201)
output = sim.run(tspan=tspan, param_values={'CC_ON': 0, 'ALLEE_ON': 0})
plt.plot(tspan/7, output.observables['Tumor_tot']/1000, lw=2, label='exponential')  # convert fM to pM for plot

# logistic growth
tspan = np.linspace(0, 30, 301)
output = sim.run(tspan=tspan, param_values={'ALLEE_ON': 0})
plt.plot(tspan/7, output.observables['Tumor_tot']/1000, lw=2, label='logistic')  # convert fM to pM for plot

# allee-logistic growth
tspan = np.linspace(0, 30, 301)
output = sim.run(tspan=tspan)
plt.plot(tspan/7, output.observables['Tumor_tot']/1000, lw=2, label='allee-logistic')  # convert fM to pM for plot

# experimental data
exp_data = get_exp_data(os.path.join('DATA', 'Mouse_Data_March2024.csv'))
plt.errorbar(*exp_data[('Tumor', 'Tumor')], ls='None', marker='o', ms=8, capsize=10, color='k', label='experiment')

plt.ylim(bottom=-0.05, top=1.05)
plt.legend(loc=0)
plt.xlabel('time (week)')
plt.ylabel('concentration (pM)')

plt.savefig('allee_logistic_tumor.pdf', format='pdf')

plt.show()
