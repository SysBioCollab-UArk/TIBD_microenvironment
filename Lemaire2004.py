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
Monomer('P',['pr'])
Monomer('Pr',['p'])
Monomer('O',['l'])
Monomer('L',['k'])
Monomer('K',['l'])


R_0 = 0.0007734 #pM Initial Number of ROB's
B_0 = 0.0007282 #pM Initial Number of AOB's

Parameter('Sp', 250) #pM/day
Parameter('Ip', 0)# pM/day range from 0 to 10^6
Parameter('kp', 86) # 1/day
Parameter('k5', 0.02)# 1/pM-day
Parameter('k6', 3) #1/day

Pbar = (Sp.value + Ip.value) / kp.value
pi_P = Pbar / (Pbar + k6.value/k5.value)
RTP = 100 # Total Number of PTH receptors


Parameter('P_init', 100)
Parameter('PrP_init', 0) #pi_P * RTP)
Parameter('Pr_init', RTP) # RTP- PrP_init.value)
Parameter('O_init', 100)
Parameter('L_init', 100)
Parameter('K_init', 100)

#Initials
Initial(P(pr=None), P_init)
Initial(Pr(p=1) % P(pr=1), PrP_init)
Initial(Pr(p=None), Pr_init)
Initial(O(l=None), O_init)
Initial(L(k=None), L_init)
Initial(K(l=None),K_init)


'''
In the the paper there is number of receptors per cell
Represented by
R(p,t), stands for number of receptor per cell, it a constant
(B + R), it stands for concentration of

'''

#Parameters

# Reaction schemes of binding of OPG-RANKL and RANK-RANKL
Parameter('po', 1)
Parameter('do', 1)
Parameter('k1', 1)
Parameter('k2', 1)
Parameter('pl', 1)
Parameter('dl', 1)



#Rules
Rule('P_creation_basal',None >>P(pr=None), Sp)
Rule('P_creation_external',None>>P(pr=None), Ip)

Rule('P_destruction',P(pr=None) >> None, kp)
Rule('P_binds_PR',P(pr=None) + Pr(p=None) | P(pr=1) % Pr(p=1),k5,k6)
# Rules for O + L
# Rule('O_creation', None >> O(l=None),po)
# Rule('O_destruction', O(l=None) >> None,do)
# Rule('O_bind_L',O(l=None) + L(k=None) | O(l=1) % L(k=1),k1,k2)
# Rule('L_creation',None >> L(k=None),pl)
# Rule('L_destruction',L(k=None) >>None,dl)
#Rules for L + K
# where does K comes from?
'''
k is not considered a model variable

'''
Parameter('k3', 1)
Parameter('k4', 1)
# Rule('K_binds_to_L',K(l=None) +L(k=None) | K(l=1) % L(k=1),k3,k4)


#Observables

# Do we want to see Overall PTH production(ONLY?)
# Do we want to see the

Observable('P_free', P(pr=None))
Observable('Pr_free', Pr(p=None))
Observable('Pr_P_bound', Pr(p=1) % P(pr=1))
# Observable('O_free', O(l=None))
# Observable('L_free', L(k=None))
# Observable('K_free', K(l=None))


Expression('piP_actual', Pr_P_bound / (Pr_free + Pr_P_bound))
Expression('piP_theory', P_free / (P_free + k6/k5))


tspan=np.linspace(0,.2, 101)
sim=ScipyOdeSimulator(model,tspan,verbose=True)
result=sim.run()

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



