from pysb import *
from pysb.util import alias_model_components
from MODULES.Lemaire2004 import create_model_elements as create_Lemaire_ME
from MODULES.Harris2024 import create_model_elements as create_Harris_ME

Model()

create_Lemaire_ME()
create_Harris_ME(OB_OC_BONE_MODEL=1)

alias_model_components()

# modify a few parameter values
R_0.value = 0.0  # fM
B_0.value = 0.0  # fM
C_0.value = 0.0  # fM
kB.value = 0.013  # /day

# print(model.parameters)

#PYDREAM_IT no-sample R_0
#PYDREAM_IT no-sample B_0
#PYDREAM_IT no-sample C_0
# PYDREAM_IT prior Cs uniform 1
# PYDREAM_IT prior DA uniform 1
# PYDREAM_IT prior dB uniform 1
# PYDREAM_IT prior DC uniform 1
# PYDREAM_IT prior DR uniform 1
#PYDREAM_IT no-sample f0
#PYDREAM_IT no-sample IL
#PYDREAM_IT no-sample IO
#PYDREAM_IT no-sample IP_const
# PYDREAM_IT prior K uniform 1
# PYDREAM_IT prior k1 uniform 1
# PYDREAM_IT prior k2 uniform 1
# PYDREAM_IT prior k3 uniform 1
# PYDREAM_IT prior k4 uniform 1
# PYDREAM_IT prior k5 uniform 1
# PYDREAM_IT prior k6 uniform 1
# PYDREAM_IT prior kB uniform 1
# PYDREAM_IT prior KLP uniform 1
# PYDREAM_IT prior kO uniform 1
# PYDREAM_IT prior KOP uniform 1
# PYDREAM_IT prior kP uniform 1
# PYDREAM_IT prior rL uniform 1
# PYDREAM_IT prior SP uniform 1
#PYDREAM_IT no-sample Bone_0
# PYDREAM_IT prior k1_B_builds_bone uniform 1
# PYDREAM_IT prior k2_B_builds_bone uniform 1
#PYDREAM_IT no-sample nB
# PYDREAM_IT prior k1_C_consumes_bone uniform 1
# PYDREAM_IT prior k2_C_consumes_bone uniform 1
#PYDREAM_IT no-sample nC
#PYDREAM_IT no-sample Tumor_0
# PYDREAM_IT prior k_tumor_div_basal uniform 1
# PYDREAM_IT prior k_tumor_dth uniform 1
#PYDREAM_IT no-sample CC_ON
#PYDREAM_IT prior N uniform 0.3
#PYDREAM_IT no-sample ALLEE_ON
#PYDREAM_IT no-sample A
# PYDREAM_IT prior k_tumor_div_TGFb uniform 1
# PYDREAM_IT prior k_tumor_PTHrP uniform 1
# PYDREAM_IT prior k_tumor_OB uniform 1
# PYDREAM_IT prior k_tumor_OC uniform 1
#PYDREAM_IT no-sample Bisphos_0
#PYDREAM_IT no-sample k_bisphos_AOC

# def likelihood(position):
#     y = np.copy(position)
#     logp_data = [0] * n_experiments
#     for n in range(n_experiments):
#         param_values[rates_mask] = 10 ** y
#         # equilibration
#         equil = solver.run(tspan=np.linspace(-500, 0, 2), param_values=param_values)
#         # add tumor cells
#         initials = equil.species[-1]
#         idx_tumor = [str(sp) for sp in model.species].index('Tumor()')  # get index of Tumor species
#         initials[idx_tumor] = 1  # fM
#         sim = solver.run(tspan=[0] + tspan[n], param_values=param_values, initials=initials).all
#         sim = sim[1:]  # don't have data for time 0
#         # calculate log-likelihood
#         for sp in like_data[n].keys():
#             logp_data[n] += np.sum(like_data[n][sp].logpdf(sim[sp][tspan_mask[n][sp]]))
#         if np.isnan(logp_data[n]):
#             logp_data[n] = -np.inf
#     return sum(logp_data)
