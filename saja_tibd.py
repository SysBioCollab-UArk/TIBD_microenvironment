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





Monomer('T') # Tumor cell
Monomer('P') # (P) PTHrP monomer
Monomer('B', ['state'], {'state': ['x','L']}) # (B) osteoblast monomer
Monomer('C') # (C) osteoclast monomer
Monomer('Beta') # (Beta) TGF-Beta monomer

Parameter('T_init',100)
Initial(T(), T_init)
Parameter('B_init', 100)
Initial(B(state='x'), B_init)
Parameter('C_init', 100)
Initial(C(state='cc'), C_init)

Parameter('ktp_tumor', 0.05)
Rule('PTHrP_production' , T() >> T() + P() , ktp_tumor )

Parameter('k_BX_BL', 0.05)
Rule('Ob_express_RANKL' , P() + B(state='x') >> P() + B(state='L') , k_BX_BL)

Parameter('k_C_cc', 1)
Rule('Oc_doupled' , B(state='L') + C() >> B(state='L') + C() + C(), k_C_cc)

Parameter('k_T_div', 1)
Rule('T_divides', T() >> T() + T(), k_T_div)

Parameter('k_T_death', 1)
Rule('T_dies', T() >> None, k_T_death)

Observable('Tumor_cells',T())
Observable('PTHrP', P())
Observable('BL' , B(state='L'))
Observable('BX', B(state='x'))
Observable('Cc', C(state='c'))
Observable('C2c', C(state='cc'))


tspan = np.linspace(0, 10, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

obs = result.observables
print(model.observables)
for o in model.observables:
    plt.plot(tspan, obs[o.name] , lw=2, label=o.name)

plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()






#Rule('P_binds_B',