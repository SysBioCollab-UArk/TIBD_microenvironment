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
Monomer('Beta', ['state'], {'state': ['a','i']}) # (Beta) TGF-Beta monomer

Parameter('T_init',100)
Initial(T(), T_init)
Parameter('B_init', 100)
Initial(B(state='x'), B_init)
Parameter('C_init', 100)
Initial(C(), C_init)
Parameter('betaa_init', 100)
Initial(Beta(state='a') , betaa_init)
Parameter('Beta_init', 100)
Initial(Beta(state='i'), Betai_init)
Parameter('B_init', 100)
Initial(B(), B_init)


Parameter('ktp_tumor', 0.05)
Rule('PTHrP_production' , T() >> T() + P() , ktp_tumor )

Parameter('k_BX_BL', 0.05)
Rule('Ob_express_RANKL' , P() + B(state='x') >> P() + B(state='L') , k_BX_BL)

Parameter('k_C_cc', 1)
Rule('Oc_interact_Ob' , B(state='L') + C() >> B(state='L') + C() + C(), k_C_cc)

Parameter('k_Oc_Beta', 1)
Rule('Oc_express_TGFB', C() >> C() + Beta(state='a') , k_Oc_Beta)


Parameter('K_Beta_tumor' , 1)
Rule('Beta_binds_T', Beta(state='a') + T() >> Beta(state='a') + T() + T() , K_Beta_tumor)

Parameter('k_T_div', 1)
Rule('T_divides', T() >> T() + T(), k_T_div)

Parameter('k_OC_div', 1)
Rule('T_divides', C() >> C() + C(), k_OC_div)

Parameter('k_OB_div', 1)
Rule('T_divides', B() >> B() + B(), k_OB_div)

Parameter('k_T_death', 1)
Rule('T_dies', T() >> None, k_T_death)

Parameter('k_OC_death', 1)
Rule('OC_dies', C() >> None, K_OC_death)

Parameter('k_OB_death', 1)
Rule('OB_dies', B() >> None, k_OB_death)



Observable('Tumor_cells',T())
Observable('PTHrP', P())
Observable('BL' , B(state='L'))
Observable('BX', B(state='x'))
Observable('OC', C())
Observable('OB', B())
Observable('TGF_Beta_A',Beta(state='a'))
Observable('TGF_Beta_I', Beta(state='i'))




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






