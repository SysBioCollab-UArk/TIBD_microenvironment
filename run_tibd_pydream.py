from MODELS.TIBD_PopD_v1 import model
from MODULES.perturbations import add_bisphosphonate_components
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *

add_bisphosphonate_components()

solver = ScipyOdeSimulator(model)

# Experiment A
tumor_injection = SequentialInjections(solver, t_equil=500, perturb_time_value=('Tumor()', 0, 1))

# Experiment B
# tumor_bisphos_injection = SequentialInjections(solver, t_equil=500,
#                                                perturb_time_amount={'Tumor()': (0, 1), 'Bisphos()': (6, 1)})
perturb_time_value = [('Tumor()', 0, 1), [('Tumor()', 0, 1), ('Bisphos()', 6, 1)]]
scale_by_eidx_time = {'Tumor_tot': (0, 14)}  # scale by output at t=14 in expt 0
multi_exp_injection = ParallelExperiments(solver, t_equil=500, perturb_time_value=perturb_time_value,
                                          scale_by_eidx_time=scale_by_eidx_time)

# Experiment C
bisphos_injection = SequentialInjections(solver, t_equil=500, perturb_time_value=('Bisphos()', 6, 1))


custom_priors = {'N': ('uniform', 0.3)}  # , 'nB': ('norm', 1), 'nC': ('norm', 1)}
no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'Tumor_0', 'CC_ON', 'nB', 'nC',
             'ALLEE_ON', 'A', 'Bisphos_0']
obs_labels = {'Bone_tot': 'bone density', 'C_obs': 'osteoclasts', 'OB_tot': 'osteoblasts',
              'Tumor_tot': 'tumor cells'}

exp_data_file = os.path.join('DATA', 'TIBD_PopD_Data.csv')

param_expts_map = {'k_tumor_div_basal': [['A'], ['B', 'C']],  #, 'D']],
                   'k_tumor_dth': [['A'], ['B', 'C']],  #, 'D']],
                   'k_tumor_div_TGFb': [['A'], ['B', 'C']],  #, 'D']],
                   'k_tumor_PTHrP': [['A'], ['B', 'C']],  #, 'D']],
                   'k_tumor_OB': [['A'], ['B', 'C']],  #, 'D']],
                   'k_tumor_OC': [['A'], ['B', 'C']]}  #, 'D']]}

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      # [tumor_injection] * 2 + [tumor_bisphos_injection, bisphos_injection],
                                      [tumor_injection, multi_exp_injection, bisphos_injection],
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True})
