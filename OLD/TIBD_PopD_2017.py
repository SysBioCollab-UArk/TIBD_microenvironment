from pysb import *
from pysb.simulator import ScipyOdeSimulator
# from pysb.logging import NAMED_LOG_LEVELS
from pysb.logging import EXTENDED_DEBUG
import numpy as np
import matplotlib.pyplot as plt
import sympy

Model()

Monomer('C')  # osteoclasts
Monomer('B')  # osteoblasts
Monomer('T')  # tumor cells
Monomer('Z')  # bone density
Monomer('Drug')

Parameter('C_0', 13)  # 15)
Initial(C(), C_0)

Parameter('B_0', 300)  # 316)
Initial(B(), B_0)

Parameter('T_0', 1)
Initial(T(), T_0)

Parameter('Z_0', 100)
Initial(Z(), Z_0)

Parameter('Drug_0', 0)
Initial(Drug(), Drug_0)

Observable('Obs_C', C())
Observable('Obs_B', B())
Observable('Obs_T', T())
Observable('Obs_Z', Z())
Observable('Obs_Drug', Drug())

Parameter('LT', 100)  # max tumor size

# Osteoclasts
Parameter('alpha1', 3)
Parameter('g11', 1.1, nonnegative=False)
Parameter('r11', 0.005, nonnegative=False) 
Parameter('g21', -0.5, nonnegative=False) 
Parameter('r21', 0.0, nonnegative=False) 
Expression('kC_growth', alpha1 * Obs_C**(g11*(1+r11*Obs_T/LT)) * Obs_B**(g21*(1+r21*Obs_T/LT)))
# #### Below is the erroneous equation from the Matlab code #####
# Expression('kC_growth', alpha1 * Obs_C**(g11*(1+r11*Obs_T/LT)) * Obs_B**(g21*(1-r21*Obs_T/LT)))
################################################################
Parameter('beta1', 0.2)
# Rule('Osteoclast_growth', None | C(), kC_growth, beta1)
Rule('Osteoclast_growth', None >> C(), kC_growth)
Rule('Osteoclast_death', C() >> None, beta1)

# Osteoblasts
Parameter('alpha2', 4)
Parameter('g12', 1, nonnegative=False)
Parameter('r12', 0, nonnegative=False)
Parameter('g22', 0, nonnegative=False)
Parameter('r22', 0.2, nonnegative=False)
Expression('kB_growth', alpha2 * Obs_C**(g12/(1+r12*Obs_T/LT)) * Obs_B**(g22-r22*Obs_T/LT))
Parameter('beta2', 0.02)
Parameter('kdrug_B', 0.001)
Expression('kB_death', beta2-kdrug_B*Obs_Drug)
# Rule('Osteoblast_growth', None | B(), kB_growth, beta2)
Rule('Osteoblast_growth', None >> B(), kB_growth)
Rule('Osteoblast_death', B() >> None, kB_death)

# Tumor mass
Parameter('gammaT', 0.005) 
Parameter('kdrug_T', 0.008)
Expression('kT_growth', (gammaT-kdrug_T*Obs_Drug)*sympy.log(LT/Obs_T))
# #### Below is the erroneous equation from the Matlab code #####
# Expression('kT_growth', (gammaT-kdrug_T*Obs_Drug)*sympy.log(LT/Obs_T, 10))
################################################################
Rule('Tumor_growth', T() >> T() + T(), kT_growth)

# Bone density
Expression('Gamma', (g12/(1+r12)) * (g21*(1+r21)) - (1-g11*(1+r11))*(1-g22+r22))
Expression('Cbar', (beta1/alpha1)**((1-(g22-r22))/Gamma) * 
           (beta2/alpha2)**(g21*(1+r21)/Gamma))
# #### Below is the erroneous equation from the Matlab code #####
# Expression('Cbar', (beta1/alpha1)**((1-(g22*(1-r22)))/Gamma) * 
#            (beta2/alpha2)**(g21*(1+r21)/Gamma))
################################################################
Expression('Bbar', (beta1/alpha1)**(g12/(1+r12)/Gamma) * 
           (beta2/alpha2)**((1-g11*(1+r11))/Gamma))

Parameter('k1', 0.0748) 
# Expression('kZ_deg', k1*((Obs_C-Cbar)/Obs_Z)*(Obs_C > Cbar))
Expression('kZ_deg', sympy.Piecewise((0, Obs_C < Cbar), (0, Obs_Z <= 0), 
                                     (k1*(Obs_C-Cbar)/Obs_Z, True)))
Rule('Bone_degradation', Z() >> None, kZ_deg)

Parameter('k2', 0.0006395) 
# Expression('kZ_syn', k2*(Obs_B-Bbar)*(Obs_B > Bbar))
Expression('kZ_syn', sympy.Piecewise((0, Obs_B < Bbar), (0, Obs_Z <= 0),
                                     (k2*(Obs_B-Bbar), True)))
Rule('Bone_synthesis', None >> Z(), kZ_syn)

# '''

tspan = np.linspace(0, 1200, 12001)  # 600, 6001)

sim = ScipyOdeSimulator(model, tspan, verbose=True)  
#                         , integrator_options = {'atol' : 1e-2, 'rtol' : 1e-2}) #, cleanup=False)
x = sim.run()
x_obs = x.observables

tspan2 = np.linspace(int(tspan[-1]), 2500, (2500-int(tspan[-1]))*10+1)

drug_idx = [str(sp) for sp in model.species].index('Drug()')
initials = x.species[-1]
initials[drug_idx] = 1

# x_obs = np.append(x_obs, sim.run(tspan=tspan2, initials=initials).observables)

names = ['Osteoclasts', 'Osteoblasts', 'Tumor', 'Bone']
steady = [Cbar, Bbar, LT.value, 0]
for obs, name, ss in zip(model.observables,names,steady):
    print(obs, name, ss)
    plt.figure()
#     np.append(tspan,tspan2)
    plt.plot(tspan, x_obs[obs.name], lw=2, label=name)
#     plt.annotate('%s(%.0f) = %.1f' % (obs.name, tspan[-1], x_obs[obs.name][-1]), 
#                  (0.65,0.85), xycoords='axes fraction', fontweight='bold')
#     plt.annotate('Theory = %.1f' % ss, 
#                  (0.65,0.8), xycoords='axes fraction', fontweight='bold')
#     plt.plot([0,tspan[-1]], [ss, ss], 'r--', lw=2)
#     plt.plot([600,600], [0,plt.ylim()[1]], 'r--', lw=2, label='Drug addition')
    plt.legend(loc=0)
    plt.xlabel('Time (days)')
    plt.ylabel('Cell density')
#     plt.ylim(ymin=0)

# plt.figure()
# for e in ['kZ_deg', 'kZ_syn']:
#     plt.plot(tspan, x.expressions[e], lw=2, label=e)
# plt.plot(tspan, x.expressions['kZ_syn']-x.expressions['kZ_deg'], lw=2, label='net growth')
# plt.legend(loc=0)
# 
# plt.figure()
# plt.plot(tspan, x.observables['Obs_C']-Cbar, lw=2, label='C(t)-Cbar')
# plt.plot(tspan, x.observables['Obs_B']-Bbar, lw=2, label='B(t)-Bbar')
# plt.legend(loc=0)

plt.show()

# '''


