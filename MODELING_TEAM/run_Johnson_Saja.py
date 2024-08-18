from MODELS.TIBD_PopD_v1 import model
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import SequentialInjections
from param_calibration import *

solver = ScipyOdeSimulator(model)
tumor_injection = SequentialInjections(solver, t_equil=500, perturb_time_amount={'Tumor()': (0, 1)})

custom_priors = {'N': ('uniform', 0.3)}

no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'nB', 'nC', 'Tumor_0', 'CC_ON',
             'ALLEE_ON', 'A']
obs_labels = {'Bone_tot': 'bone density', 'C_obs': 'osteoclasts', 'OB_tot': 'osteoblasts', 'Tumor_tot': 'tumor cells'}

exp_data_file = '../DATA/TIBD_John2011_PBS_Tumor'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      tumor_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True)
