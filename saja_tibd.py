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
Monomer('P')  # PTHrP
Monomer('B', ['state'], {'state': ['x', 'L']})  # osteoblast
Monomer('C')  # osteoclast
Monomer('Beta')  # TGF-Beta
Monomer('Bone')

Parameter('T_init', 100)  # 100
Parameter('B_init', 100)
Parameter('C_init', 100)
Parameter('Beta_init', 0)
Parameter('Bone_init', 10000)

Initial(T(), T_init)
Initial(B(state='x'), B_init)
Initial(C(), C_init)
Initial(Beta(), Beta_init)
Initial(Bone(), Bone_init)

Parameter('ktp_tumor', 1)
Rule('PTHrP_production', T() >> T() + P(), ktp_tumor)

Parameter('kf_BX_BL', 1)
Rule('Ob_express_RANKL', P() + B(state='x') >> P() + B(state='L'), kf_BX_BL)

Parameter('kr_BX_BL', 100)  # 100
Rule('Ob_unexpress_RANKL', B(state='L') >> B(state='x'), kr_BX_BL)

Parameter('k_C_cc', 0.1)
Rule('Oc_interact_Ob', B(state='L') + C() >> B(state='L') + C() + C(), k_C_cc)

Parameter('k_Oc_Bone', 0.001)  # 0.001
# TGF-beta being produced by OC consumption of bone.
Rule('Oc_consumes_bone', C() + Bone() >> C() + Beta(), k_Oc_Bone)

Parameter('k_Ob_Bone',  0.1)
Rule('Ob_produces_bone', B() + Beta() >> B() + Bone(), k_Ob_Bone)

Parameter('k_Beta_tumor', 0.1)  # 0.14 # 0.01
Rule('T_divides_beta', Beta() + T() >> Beta() + T() + T(), k_Beta_tumor)

Parameter('k_T_div', 1)
Parameter('k_T_div_cc', k_T_div.value/2000)
Rule('T_divides', T() | T() + T(), k_T_div, k_T_div_cc)

Parameter('k_OC_div', 0.1)  # 0.1
Parameter('k_OC_div_cc', k_OC_div.value/2000)
Rule('C_divides', C() | C() + C(), k_OC_div, k_OC_div_cc)

Parameter('k_OB_div', 1)
Parameter('k_OB_div_cc', k_OB_div.value/2000)
Rule('Bx_divides', B(state='x') | B(state='x') + B(state='x'), k_OB_div, k_OB_div_cc)
Rule('BL_divides', B(state='L') | B(state='L') + B(state='L'), k_OB_div, k_OB_div_cc)

Parameter('k_Beta_synth', 0.01)
Rule('OC_create_beta', C() >> C() + Beta(), k_Beta_synth)
Rule('OB_create_beta', B() >> B() + Beta(), k_Beta_synth)
Rule('T_create_beta', T() >> T() + Beta(), k_Beta_synth)

Parameter('k_T_death', 10)  # 10
Rule('T_dies', T() >> None, k_T_death)

Parameter('k_OC_death', 1)
Rule('OC_dies', C() >> None, k_OC_death)

Parameter('k_OB_death', 1)
Rule('OB_dies', B() >> None, k_OB_death)

Parameter('k_PTHrP_deg', 1)
Rule('PTHrP_degrades', P() >> None, k_PTHrP_deg)

Parameter('k_TGF_deg', 10)  # 1
Rule('TGF_degrades', Beta() >> None, k_TGF_deg)

Observable('Tumor_cells', T())
# Observable('osteoblasts', B())
Observable('osteoblasts_RANKL', B(state='L'))
Observable('osteoblasts_no_lig', B(state='x'))
Observable('osteoclasts', C())
Observable('PTHrP', P())
Observable('TGF_Beta', Beta())
Observable('bone_frac', Bone())

tspan = np.linspace(0, 2.5, 1001)  # 100
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

obs = result.observables
obs['bone_frac'] = obs['bone_frac'] / Bone_init.value * 100  # convert into a value between 0 and 100
# print(model.observables)
for o in model.observables:
    plt.plot(tspan, obs[o.name], lw=2, label=o.name)

plt.xlabel('time')
plt.ylabel('concentration')
# plt.yscale('log')
plt.legend(loc=0)

plt.show()
