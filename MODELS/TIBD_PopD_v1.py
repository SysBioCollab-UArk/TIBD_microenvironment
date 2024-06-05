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

if __name__ == '__main__':

    from param_calibration import *
    from sim_protocols import tumor_injection
    import glob

    exp_data_file = 'TIBD_PopD_Mouse_Data_2expts.csv'
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

    calibrator.run(niterations=50000, nchains=5)

    '''
    print('Plotting log-likelihoods')
    logps_files = glob.glob('dreamzs*logps*')
    calibrator.plot_log_likelihood(logps_files, cutoff=2)

    print('Plotting parameter distributions')
    samples_files = glob.glob('dreamzs*params*')
    labels = [calibrator.model.parameters[i].name for i in calibrator.parameter_idxs]
    param_samples = calibrator.plot_param_dist(samples_files, labels=labels, cutoff=2)

    # code to test plotting parameter distributions using groups
    # groups = [[i for i in range(0, 12)], [i for i in range(12, 30)]]
    # labels = [[calibrator.model.parameters[calibrator.parameter_idxs[i]].name for i in group] for group in groups]
    # filenames = ['fig_PyDREAM_histograms_PARAMS_%d_%d' % (i, j) for (i, j) in [(0, 11), (12, 29)]]
    # param_samples = calibrator.plot_param_dist(samples_files, labels=labels, groups=groups, cutoff=2,
    #                                            save_plot=filenames)

    print('Plotting time courses')
    tspan = np.linspace(0, 30, 31)  # time points after tumor injection
    calibrator.plot_timecourses(calibrator.model, calibrator.raw_data, calibrator.sim_protocols, tspan, param_samples,
                                calibrator.parameter_idxs,
                                xlabel='time (day)',
                                ylabels=['relative BV/TV'] + ['concentration (fM)'] * 3,
                                locs=['upper right', 'upper left', 'upper right', 'upper left'])

    plt.show()
    '''
