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

# TODO: Include bone in the model
Monomer('T')  # Tumor cell
Monomer('P')  # (P) PTHrP monomer
Monomer('B', ['state'], {'state': ['x', 'L']})  # (B) osteoblast monomer
Monomer('C')  # (C) osteoclast monomer
Monomer('Beta')  # (Beta) TGF-Beta monomer
Monomer('Bone')

Parameter('T_init', 100)
Initial(T(), T_init)
Parameter('B_init', 100)
Initial(B(state='x'), B_init)
Parameter('C_init', 100)
Initial(C(), C_init)
Parameter('beta_init', 100)
Initial(Beta(), beta_init)
Parameter('Bone_init', 10000)
Initial(Bone(), Bone_init)

Parameter('ktp_tumor', 1)
Rule('PTHrP_production', T() >> T() + P(), ktp_tumor)

Parameter('k_BX_BL', 1)
Rule('Ob_express_RANKL', P() + B(state='x') >> P() + B(state='L'), k_BX_BL)

Parameter('k_C_cc', 0.1)
Rule('Oc_interact_Ob' , B(state='L') + C() >> B(state='L') + C() + C(), k_C_cc)

Parameter('k_Oc_Beta', 0.0001)
# TGF-beta being produced by OC consumption of bone.
Rule('Oc_express_TGFB', C() + Bone() >> C() + Beta(), k_Oc_Beta)

Parameter('K_Beta_tumor', 0.1)
Rule('T_divides_beta', Beta() + T() >> Beta() + T() + T(), K_Beta_tumor)

Parameter('k_T_div', 1)
Rule('T_divides', T() >> T() + T(), k_T_div)

Parameter('k_OC_div', 1)
Rule('C_divides', C() >> C() + C(), k_OC_div)

Parameter('k_OB_div', 1)
Rule('Bx_divides', B(state='x') >> B(state='x') + B(state='x'), k_OB_div)
Rule('BL_divides', B(state='L') >> B(state='L') + B(state='L'), k_OB_div)

Parameter('k_Beta_synth', 0.1)
Rule('OC_create_beta' , C() >> C() + Beta(), k_Beta_synth)
Rule('OB_create_beta' , B() >> B() + Beta(), k_Beta_synth)
Rule('T_create_beta' , T() >> T() + Beta(), k_Beta_synth)


Parameter('k_T_death', 10)
Rule('T_dies', T() >> None, k_T_death)

Parameter('k_OC_death', 10)
Rule('OC_dies', C() >> None, k_OC_death)

Parameter('k_OB_death', 1)
Rule('OB_dies', B() >> None, k_OB_death)

Parameter('k_PTHrP_deg', 1)
Rule('PTHrP_degrades', P() >> None, k_PTHrP_deg)

Parameter('k_TGF_deg', 1)
Rule('TGF_degrades', Beta() >> None, k_TGF_deg)

#Observable('Tumor_cells', T())
Observable('PTHrP', P())
Observable('BL', B(state='L'))
Observable('BX', B(state='x'))
Observable('OC', C())
Observable('OB', B())
Observable('TGF_Beta', Beta())
Observable('bone_frac', Bone())


tspan = np.linspace(0, 1, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

obs = result.observables
obs['bone_frac'] = obs['bone_frac'] / Bone_init.value * 100  # convert into a value between 0 and 100
#print(model.observables)
for o in model.observables:
    plt.plot(tspan, obs[o.name], lw=2, label=o.name)

plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc=0)

plt.show()






