from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
# Parameters
# Initials
# Rules
# Observables
# simulation commands

# Monomers
Monomer('P', ['pr'])
Monomer('Pr', ['p'])
Monomer('O', ['l'])
Monomer('L', ['k'])
Monomer('K', ['l'])

Observable('P_free', P(pr=None))
Observable('Pr_free', Pr(p=None))
Observable('Pr_P_bound', Pr(p=1) % P(pr=1))
Observable('O_free', O(l=None))
Observable('L_free', L(k=None))
Observable('K_free', K(l=None))
Observable('O_L_bound', O(l=1) % L(k=1))
Observable('K_L_bound', K(l=1) % L(k=1))

Expression('piP_actual', Pr_P_bound / (Pr_free + Pr_P_bound))
Expression('piP_theory', P_free / (P_free + k6/k5))

R_0 = 0.0007734  # pM Initial Number of ROB's
B_0 = 0.0007282  # pM Initial Number of AOB's
Kop = 2e5  # Minimal rate of production of OPG per cell
Klp = 3e6  # pmol/pmol cells, Maximum number of RANKL attached to each cell surface
rL = 1e3  # pM/day, Rate of RANKL production and elimination

#Parameter for P and Pr
Parameter('Sp', 250)  # pM/day
Parameter('Ip', 0)  # pM/day range from 0 to 10^6
Parameter('kp', 86)  # 1/day
Parameter('k5', 0.02)  # 1/pM-day
Parameter('k6', 3)  # 1/day (Rate of PTH binding with receptor)

# Parameters for O
Expression('So',Kop/piP_actual*R_0)
Parameter('Io',0 )  # pM/day range from 0 to 10^6
Parameter('ko', 0.35)  # 1/day (Rate of elimination of OPG)
Parameter('k1', 1e-2)  # 1/pMday ( Rate of OPG-RANKL binding)
Parameter('k2', 10)  # 1/day (Rate of OPG-RANKL unbinding)

# Parameters for L
Expression('Sl', rL * (1-(L_free + O_L_bound + K_L_bound)/(Klp * piP_actual * B_0)))
Parameter('Il', 0)  # pM/day range from 0 to 10^6

# Parameters for K
Parameter('k3', 5.8e-4)  # 1/pMday Rate of RANK-RANKL binding
Parameter('k4', 1.7e-2)  # 1/day Rate of RANK-RANKL unbinding



# Model Parameter values
# R_0 = 0.0007734#pM Initial Number of ROB's
# B_0 = 0.0007282#pM Initial Number of AOB's
# K_o_p = 2e5  # Minimal rate of production of OPG per cell
# K_p_L = 3e6  # pmol/pmol cells, Maximum number of RANKL attached to each cell surface
# K_ = 10  # pM Fixed concentration of RANK
# r_l = 1e3  # pM/day, Rate of RANKL production and elimination
I_o = 0  # (0-10^6) pM/day, Rate of administration of OPG
I_L = 0  # (0-10^6) pM/day, Rate of administration of RANKL


#  Pbar, concentration of P at steady state
Pbar = (Sp.value + Ip.value) / kp.value
pi_P = Pbar / (Pbar + k6.value/k5.value)
RTP = 100  # Total Number of PTH receptors
po_ = (Pbar / pi_P) * R_0 + Ip  # Production rate of OPG being down-regulated by PTH
O_ = (K_o_p / k_o * pi_P) * R_0 + (I_o / k_o)
L_ = (K_p_L * pi_P * B_0) / 1 + k3 * K_ / k4 + k1 * O_ / k2 * (1 + I_L / r_l)



Parameter('P_init', 100)
Parameter('PrP_init', 0)  # pi_P * RTP)
Parameter('Pr_init', RTP)  # RTP - PrP_init.value)
Parameter('O_init', 100)
Parameter('L_init', 100)
Parameter('K_init', 100)


#Initials
Initial(P(pr=None), P_init)
Initial(Pr(p=1) % P(pr=1), PrP_init)
Initial(Pr(p=None), Pr_init)
Initial(O(l=None), O_init)
Initial(L(k=None), L_init)
Initial(K(l=None), K_init)

'''
In the the paper there is number of receptors per cell
Represented by R(p,t), stands for number of receptor per cell, it a constant
(B + R), it stands for concentration of

'''

#Parameters

# Reaction schemes of binding of OPG-RANKL and RANK-RANKL

#Rules
Rule('P_creation_basal', None >> P(pr=None), Sp)
Rule('P_creation_external', None >> P(pr=None), Ip)
Rule('P_destruction', P(pr=None) >> None, kp)
Rule('P_binds_PR', P(pr=None) + Pr(p=None) | P(pr=1) % Pr(p=1), k5, k6)
# Rules for O + L
Rule('O_creation_basal', None >> O(l=None), So)
Rule('O_creation_external', None >> O(l=None), Io)
Rule('O_destruction', O(l=None) >> None, ko)
Rule('O_bind_L', O(l=None) + L(k=None) | O(l=1) % L(k=1), k1, k2)
Rule('L_creation_basal', None >> L(k=None), Sl)
Rule('L_creation_external', None >> L(k=None), Il)
# Rule('L_destruction', L(k=None) >> None, dl)
#Rules for L + K
Rule('K_binds_to_L', K(l=None) +L(k=None) | K(l=1) % L(k=1), k3, k4)



tspan =np.linspace(0,.2,101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
# for i,ODE in enumerate(model.odes):
#     print(model.species[i], ODE)
# quit()

result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)


plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc='best')

plt.figure()
for exp in model.expressions:
    plt.plot(tspan, result.expressions[exp.name], lw=2, label=exp.name)
    # print(exp)
plt.xlabel('time')
plt.ylabel('value')
plt.legend(loc='best')

plt.show()
