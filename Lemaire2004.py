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
Monomer('P'['pr'],)
Monomer('Pr',['p'])
Monomer('O',['l'])
Monomer('L',['o'],['k'])
Monomer('K',['l'])


'''
In the the paper there is number of receptors per cell
Represented by
R(p,t), stands for number of receptor per cell, it a constant
(B + R), it stands for concentration of


'''

#Parameters
Parameter('pp_P_creation', 1)
Parameter('dp_destruction', 1)
Parameter('k5, k5_P_binds_Pp', 1)
# Reaction schemes of binding of OPG-RANKL and RANK-RANKL
Parameter('po_O_creation', 1)
Parameter('do_O_destruction', 1)
Parameter('k1_k2_O_bind_L', 1)
Parameter('k3_k4_K_binds to L', 1)
Parameter('pl_L_creation', 1)
Parameter('dl_L_destruction', 1)

#Rules
Rule('P_creation',P(None)>>P(1), pp)
Rule('P_destruction',P(1) >> P(None), dp)
Rule('P_binds_Pp',P([p=None]) + Pr([p=None]) | P([pr=1]), Pr([p=1]),k5,k6)
# Rules for O + L
Rule('O_creation',O(None) >> O(1),po)
Rule('O_destruction', o(1) >> O(None),do)
Rule('O_bind_L',O(l=None) + L(o=None) | O(l=1) % L(o=1),k1,k2)
Rule('L_creation',L(None) >> l(1),pl)
Rule('L_destruction',L(1) >>L(None),dl)
#Rules for L + K
# where does K comes from?
'''
k is not considered a model variable

'''
Rule('K_binds to L',K(l=None) +L(k=None) | K(l=1) +L(k=1),k3,k4)

#Initials
Initial(P(None), P_init)
Initial(Pp(p=None), Pp_Init)
Initial(O(l=None), 0_Init), O_Init)
Initial(L(0=None, (k=None)), L_init)

#Observables

# Do we want to see Overall PTH production(ONLY?)

Observable('P_creation', Pr(loc='creation'))
Observable('P_destruction', Pr(loc='destruction'))
Observable('O_creation', L(loc='creation'))
Observable('O_destruction', K(loc='destruction'))


span=np.linspace(0,100,101)
sim=ScipyOdeSimulator(model,tspan,verbose=True)
result=sim.run()




plt.plot(tspan,result.observables['P_creation'], lw=2,label='creation')
plt.plot(tspan,result.observables['P_destruction'], lw=2,label='destruction')
plt.plot(tspan,result.observables['O_creation'], lw=2,label='creation')
plt.plot(tspan,result.observables['O_destruction'],lw=2,label='destruction')

plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc='best')

plt.show()



