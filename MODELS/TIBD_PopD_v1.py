from pysb import *
from MODULES.Lemaire2004 import create_model_elements as create_Lemaire_ME
from MODULES.Harris2024 import create_model_elements as create_Harris_ME

Model()

create_Lemaire_ME()  # TODO: consider passing in a scaling factor to adjust concentration units (e.g., pM to cells/mm)
create_Harris_ME(OB_OC_BONE_MODEL=1)
# add_bisphosphonate_components()

# modify a few parameter values
R_0.value = 0.0  # fM
B_0.value = 0.0  # fM
C_0.value = 0.0  # fM
kB.value = 0.013  # /day

if __name__ == '__main__':

    from param_calibration import *
    from SIM_PROTOCOLS.sim_protocols import tumor_injection
    import glob

    exp_data_file = '../DATA/TIBD_PopD_Mouse_Data_2expts.csv'
    sim_protocols = [tumor_injection] * 2
    custom_priors = {'N': ('uniform', 0.3)}
    no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'nB', 'nC', 'Tumor_0', 'CC_ON',
                 'ALLEE_ON', 'A', 'Bisphos_0', 'k_bisphos_AOC']

    # If parameter values can vary between experiments, define groups of experiments where they are the same.
    # If a parameter value is unique to an experiment, the group will be length 1. The default is for parameter values
    # to be common across all experiments, in which case there will be one group containing indices for all experiments.
    # This is handled internally, so only need to define cases where param values are NOT common across all expts.
    param_expts_map = {'Cs': [['A', 'B']],  # this is just here to show what the default is
                       'k_tumor_div_basal': [['A'], ['B']],  # this and the following params can vary btwn the 2 expts
                       'k_tumor_dth': [['A'], ['B']],
                       'k_tumor_div_TGFb': [['A'], ['B']],
                       'k_tumor_PTHrP': [['A'], ['B']],
                       'k_tumor_OB': [['A'], ['B']],
                       'k_tumor_OC': [['A'], ['B']]}

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocols,
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)

    calibrator.run(niterations=50000, nchains=5, plot_results=True, plot_tc_args={'separate_plots': False})

    # can plot the results yourself using the code below
    # logps_files = glob.glob('dreamzs*logps*')
    # samples_files = glob.glob('dreamzs*params*')
    # calibrator.create_figures(logps_files, samples_files, show_plots=True, plot_tc_args={'separate_plots': False})

    '''
    print('Plotting log-likelihoods')
    logps_files = glob.glob('dreamzs*logps*')
    calibrator.plot_log_likelihood(logps_files, cutoff=2)

    print('Plotting parameter distributions')
    samples_files = glob.glob('dreamzs*params*')
    # get groups of parameters common for different sets of experiments
    param_groups = []
    expt_idxs = [n for n in range(calibrator.n_experiments)]
    for n in expt_idxs:
        for combo in combinations(expt_idxs, n + 1):  # all possible combinations of expt groups
            param_list_1 = [set(calibrator.samples_idxs[i]) for i in combo]
            param_list_2 = [set(calibrator.samples_idxs[i]) for i in expt_idxs if i not in combo]
            if len(param_list_2) == 0:
                param_list_2 = [set()]  # empty set
            # get parameters unique to this group of expts (might be none)
            param_diff = set.intersection(*param_list_1) - set.intersection(*param_list_2)
            if len(param_diff) > 0:
                param_groups.append(sorted(list(param_diff)))
    # get the parameter names to label the plots in the histogram
    param_labels = []
    for group in param_groups:
        param_labels.append([])
        for i in group:
            for s_idx_list in calibrator.samples_idxs:
                if i in s_idx_list:
                    param_labels[-1].append(
                        calibrator.model.parameters[calibrator.parameter_idxs[s_idx_list.index(i)]].name)

    param_samples = calibrator.plot_param_dist(samples_files, labels=param_labels, groups=param_groups, cutoff=2)

    print('Plotting time courses')
    tspans = [np.linspace(
        calibrator.tdata[n][0], calibrator.tdata[n][-1],
        int((calibrator.tdata[n][-1] - calibrator.tdata[n][0]) * 10 + 1))
        for n in range(calibrator.n_experiments)]

    calibrator.plot_timecourses(calibrator.model, calibrator.raw_data, calibrator.sim_protocols, tspans, param_samples[:100],
                                calibrator.parameter_idxs, calibrator.samples_idxs, separate_plots=False,
                                xlabel='time (day)',
                                ylabels=[['relative BV/TV'] + ['concentration (fM)'] * 3] * calibrator.n_experiments,
                                locs=[['upper right', 'upper left', 'upper right',
                                      'upper left']] * calibrator.n_experiments)

    plt.show()
    '''
