from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# This is a reproduction of the computational model presented in Lemaire et al., J. Theor. Biol. 229, 293-309 (2004).
# https://doi.org/10.1016/j.jtbi.2004.03.023

Model()

Monomer('R')
Monomer('B')
Monomer('C')

Parameter('R_0', 0.7734e-3)
Parameter('B_0', 0.7282e-3)
Parameter('C_0', 0.9127e-3)

Initial(R(), R_0)
Initial(B(), B_0)
Initial(C(), C_0)

Observable('R_obs', R())
Observable('B_obs', B())
Observable('C_obs', C())

# constants
Cs = 5e-3  # pM
DA = 0.7  # /day
dB = 0.7  # /day
DC = 2.1e-3  # pM/day
DR = 7e-4  # pM/day
f0 = 0.05  # unitless
Parameter('IL', 0)  # pM/day (range: 0-1e6)
Parameter('IO', 0)  # pM/day (range: 0-1e6)
Parameter('IP', 0)  # pM/day (range: 0-1e6)
K = 10  # pM
k1 = 1e-2  # /pM-day
k2 = 10  # /day
k3 = 5.8e-4  # /pM-day
k4 = 1.7e-2  # /day
k5 = 0.02  # /pM-day
k6 = 3  # /day
Parameter('kB', 0.189)  # /day
KLP = 3e6  # unitless
kO = 0.35  # /day
KOP = 2e5  # /day
kP = 86  # /day
rL = 1e3  # pM/day
SP = 250  # pM/day

# compound constants
# Pbar = IP/kP
P0 = SP/kP
Ps = k6/k5

# ratios
Expression('pi_P', (IP/kP + P0) / (IP/kP + Ps))
Expression('pi_C', (C_obs + f0*Cs) / (C_obs + Cs))
Expression('pi_L', k3/k4 * KLP*pi_P*B_obs * (1 + IL/rL) / (1 + k3*K/k4 + k1/k2/kO*(KOP/pi_P*R_obs + IO)))

# cell dynamics rules
Expression('DR_pi_C', DR*pi_C)
Expression('DB_over_pi_C', f0*dB/pi_C)
Expression('DC_pi_L', DC*pi_L)
Expression('DA_pi_C', DA*pi_C)

Rule('ROB_creation', None >> R(), DR_pi_C)
Rule('ROB_to_AOB', R() >> B(), DB_over_pi_C)
Rule('AOB_death', B() >> None, kB)
Rule('AOC_creation_death', None | C(), DC_pi_L, DA_pi_C)

# cell injections
Parameter('kf_AOB', 0)  # pM/day
Parameter('kf_AOC', 0)  # pM/day
Parameter('kf_ROB', 0)  # pM/day
Parameter('kr_AOB', 0)  # pM/day
Parameter('kr_AOC', 0)  # pM/day
Parameter('kr_ROB', 0)  # pM/day
Rule('AOB_injection', None | B(), kf_AOB, kr_AOB)
Rule('AOC_injection', None | C(), kf_AOC, kr_AOC)
Rule('ROB_injection', None | R(), kf_ROB, kr_ROB)

# simulations
sim = ScipyOdeSimulator(model, verbose=True)

tspan1 = np.linspace(0, 20, 201)
tspan2 = np.linspace(20, 80, 601)
tspan3 = np.linspace(80, 140, 601)

# perturbations
perturb = [
    [{'kf_AOB': 0}, {'kf_AOB': 1e-4}, {'kf_AOB': 0}],
    [{'kf_AOC': 0}, {'kf_AOC': 1e-4}, {'kf_AOC': 0}],
    [{'kf_ROB': 0}, {'kf_ROB': 1e-4}, {'kf_ROB': 0}],
    [{'kr_AOB': 0}, {'kr_AOB': 8.3e-2}, {'kr_AOB': 0}],  # something wrong here
    [{'kr_AOC': 0}, {'kr_AOC': 2.9e-1}, {'kr_AOC': 0}],  # something wrong here
    [{'kr_ROB': 0}, {'kr_ROB': 1.2e-1}, {'kr_ROB': 0}],  # something wrong here
    [{'IP': 0}, {'IP': 1e3}, {'IP': 0}],
    [{'IO': 0}, {'IO': 2e5}, {'IO': 0}],
    [{'IL': 0, 'IO': 0}, {'IL': 1e4, 'IO': 0}, {'IL': 1e4, 'IO': 9e4}]  # something wrong with IL here
]

colors = ['b', 'r', 'g']
labels = ['ROB', 'AOB', 'AOC']
fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, figsize=(10,8))
row = 0
col = 0
for p in perturb:
    print(p)
    output1 = sim.run(param_values=p[0], tspan=tspan1)
    output2 = sim.run(param_values=p[1], tspan=tspan2, initials=output1.species[-1])
    output3 = sim.run(param_values=p[2], tspan=tspan3, initials=output2.species[-1])
    for obs, color, label in zip(model.observables, colors, labels):
        axs[row, col].plot(tspan1, output1.observables[obs.name], lw=2, color=color, label=label)
        axs[row, col].plot(tspan2, output2.observables[obs.name], lw=2, color=color)
        axs[row, col].plot(tspan3, output3.observables[obs.name], lw=2, color=color)
    if row == 2:
        axs[row, col].set_xlabel('time (d)')
    if col == 0:
        axs[row, col].set_ylabel('concentration (pM)')
    axs[row, col].set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
    axs[row, col].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axs[row, col].legend(loc=0)
    if (row+1) % 3 == 0:
        col += 1
        row = 0
    else:
        row += 1

# for i, ode in enumerate(model.odes):
#     print('%s: %s' % (model.species[i], ode))

plt.tight_layout()
plt.show()
