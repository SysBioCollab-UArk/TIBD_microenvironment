from pysb import *
from MODULES.Lemaire2004 import create_model_elements as create_Lemaire_ME
from MODULES.Harris2024 import create_model_elements as create_Harris_ME
from MODULES.perturbations import add_bisphosphonate_components

Model()

create_Lemaire_ME()
create_Harris_ME(OB_OC_BONE_MODEL=1)
# add_bisphosphonate_components()

# modify a few parameter values
R_0.value = 0.0  # fM
B_0.value = 0.0  # fM
C_0.value = 0.0  # fM
kB.value = 0.013  # /day

#PYDREAM_IT no-sample R_0
#PYDREAM_IT no-sample B_0
#PYDREAM_IT no-sample C_0
#PYDREAM_IT no-sample f0
#PYDREAM_IT no-sample IL
#PYDREAM_IT no-sample IO
#PYDREAM_IT no-sample IP_const
#PYDREAM_IT no-sample Bone_0
#PYDREAM_IT no-sample nB
#PYDREAM_IT no-sample nC
#PYDREAM_IT no-sample Tumor_0
#PYDREAM_IT no-sample CC_ON
#PYDREAM_IT prior N uniform 0.3
#PYDREAM_IT no-sample ALLEE_ON
#PYDREAM_IT no-sample A
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
#         sim = solver.run(tspan[n], param_values=param_values, initials=initials).all
#         # calculate log-likelihood
#         for sp in like_data[n].keys():
#             logp_data[n] += np.sum(like_data[n][sp].logpdf(sim[sp][tspan_mask[n][sp]]))
#         if np.isnan(logp_data[n]):
#             logp_data[n] = -np.inf
#     return sum(logp_data)

if __name__ == '__main__':

    from param_calibration import *
    from sim_protocols import tumor_injection, tumor_injection_2
    import glob

    custom_priors = {'N': ('uniform', 0.3)}
    no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'nB', 'nC', 'Tumor_0', 'CC_ON',
                 'ALLEE_ON', 'A', 'Bisphos_0', 'k_bisphos_AOC']

    # TODO: Still need to figure out how to define the sample_idxs lists for all experiments

    calibrator = ParameterCalibration(model,
                                      'TIBD_PopD_Mouse_Data.csv',
                                      tumor_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5)

    # print('Plotting log-likelihoods')
    # logps_files = glob.glob('dreamzs*logps*')
    # calibrator.plot_log_likelihood(logps_files, cutoff=2)

    # print('Plotting parameter distributions')
    # samples_files = glob.glob('dreamzs*params*')
    # param_samples = calibrator.plot_param_dist(samples_files,
    #                                            labels=[calibrator.model.parameters[i].name
    #                                                    for i in calibrator.parameter_idxs],
    #                                            cutoff=2)

    # groups = [[i for i in range(0, 12)],
    #           [i for i in range(12, 30)]]
    # labels = [[model.parameters[parameters_idxs[i]].name for i in group] for group in groups]
    # plot_param_dist(samples, labels, groups=groups)

    # print('Plotting time courses')
    # tspan = np.linspace(0, 30, 31)  # time points after tumor injection
    # calibrator.plot_timecourses(calibrator.model, calibrator.raw_data, calibrator.sim_protocols, tspan, param_samples,
    #                             calibrator.parameter_idxs,
    #                             xlabel='time (day)',
    #                             ylabels=['relative BV/TV'] + ['concentration (fM)'] * 3,
    #                             locs=['upper right', 'upper left', 'upper right', 'upper left'])

    # samples_idxs = [[i for i in range(30)], [i for i in range(30)]]  # different params for different experiments TODO

    # plt.show()
