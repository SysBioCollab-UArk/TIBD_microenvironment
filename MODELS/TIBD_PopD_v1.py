from pysb import *
from MODULES.Lemaire2004 import create_model_elements as create_Lemaire_ME
from MODULES.Harris2024 import create_model_elements as create_Harris_ME

Model()

create_Lemaire_ME()  # TODO: consider passing in a scaling factor to adjust concentration units (e.g., pM to cells/mm)
create_Harris_ME(OB_OC_BONE_MODEL=1)

# modify a few parameter values
R_0.value = 0.0  # fM
B_0.value = 0.0  # fM
C_0.value = 0.0  # fM
kB.value = 0.013  # /day


if __name__ == '__main__':
    from param_calibration import *
    from pysb.simulator import ScipyOdeSimulator
    from SIM_PROTOCOLS.sim_protocols import SequentialInjections
    from MODULES.perturbations import *
    import pandas as pd
    from matplotlib.lines import Line2D

    add_bisphosphonate_components()
    add_RANKLi_components()
    add_tumor_RANK_components(tumorK_migration=True)
    add_tumor_migration_components()

    solver = ScipyOdeSimulator(model)

    tumor_tibia_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ("Tumor()", 1)})
    tumor_bisphos_tibia_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ('Tumor()', 1), 6: ('Bisphos()', 1)})
    tumor_RANKLi_tibia_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ('Tumor()', 1), 6: ('RANKLi()', 10)})

    tumor_cardiac_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ("Tumor_out()", 1)})
    tumor_bisphos_cardiac_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ('Tumor_out()', 1), 6: ('Bisphos()', 1)})
    tumor_RANKLi_cardiac_injection = \
        SequentialInjections(solver, t_equil=500, time_perturb_value={0: ('Tumor_out()', 1), 6: ('RANKLi()', 10)})

    custom_priors = {'N': ('uniform', 0.3)}  # , 'nB': ('norm', 1), 'nC': ('norm', 1)}
    no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'Tumor_WT_0', 'CC_ON', 'nB', 'nC',
                 'ALLEE_ON', 'A']
    obs_labels = {'Bone_tot': 'Bone Density', 'C_obs': 'Osteoclasts', 'OB_tot': 'Osteoblasts',
                  'Tumor_tot': 'Tumor Cells'}

    # Run a simulation
    tspan = np.linspace(0, 28, 281)
    params_df = pd.read_csv('fit_params_PAR.csv')  # parental parameter values
    params_dict = dict(zip(params_df['parameter'], params_df['value']))
    params_df = pd.read_csv('fit_params_BONE.csv', index_col='parameter')  # bone clone parameter values
    new_params = {
        'k_tumorKL_div_basal': 'k_tumor_div_basal',
        'k_tumorKL_div_TGFb': 'k_tumor_div_TGFb',
        'k_tumorKL_dth': 'k_tumor_dth',
        'k_tumorKL_PTHrP': 'k_tumor_PTHrP',
        'k_tumorKL_OB': 'k_tumor_OB'
    }
    for p in new_params.keys():
        params_dict[p] = params_df.loc[new_params[p], 'value']
    param_values = [params_dict.get(p.name, p.value) for p in model.parameters]

    results = [
        ('Ctrl tibia',      tumor_tibia_injection.run(tspan, param_values)),
        ('+ZA tibia'    ,   tumor_bisphos_tibia_injection.run(tspan, param_values)),
        ('+RANKLi tibia',   tumor_RANKLi_tibia_injection.run(tspan, param_values)),
        ('Ctrl cardiac',    tumor_cardiac_injection.run(tspan, param_values)),
        ('+ZA cardiac',     tumor_bisphos_cardiac_injection.run(tspan, param_values)),
        ('+RANKLi cardiac', tumor_RANKLi_cardiac_injection.run(tspan, param_values))
    ]

    # plot simulation results
    fig_all, axs_all = plt.subplots(2, 2, sharex=True, constrained_layout=True)
    axs_flat_all = axs_all.flatten()
    fig_tumor, axs_tumor = plt.subplots(len(results) // 2, 2, sharex=True, sharey='row', constrained_layout=True,
                                        figsize=(1 * 6.4, 1.25 * 4.8))
    axs_flat_tumor = axs_tumor.flatten(order='F')
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, (label, result) in enumerate(results):
        # tumor plot
        axs_flat_tumor[i].plot(tspan, result['Tumor_tot'], lw=3, ls='-', color=colors[i % (len(results) // 2)],
                               label='Tumor cells')
        axs_flat_tumor[i].plot(tspan, result['Tumor_obs'], lw=2, ls='--', color='0.5', label=r'RANK$^-$')
        axs_flat_tumor[i].plot(tspan, result['TumorK_tot'], lw=2, ls=':', color='0.5', label=r'RANK$^+$')
        axs_flat_tumor[i].plot(tspan, result['TumorK_u'], lw=1, ls='--', color='r', label=r'K$^+$ unbound')
        axs_flat_tumor[i].plot(tspan, result['TumorK_L'], lw=1, ls=':', color='r', label=r'K$^+$-RANKL')
        if i == len(results) - 1:
            axs_flat_tumor[i].set_xlabel('Time (day)')
        axs_flat_tumor[i].set_ylabel('Amount')
        axs_flat_tumor[i].set_xticks(range(0, int(tspan[-1]) + 1, 7))
        axs_flat_tumor[i].annotate(label, xy=(0.1, 0.85), xycoords='axes fraction', color=colors[i % (len(results) // 2)],
                                   fontsize=12, fontweight='bold')
        # all observables plot
        for obs, ax in zip(obs_labels.keys(), axs_flat_all):
            ls, lw = ('-', 3) if i < len(results) // 2 else ('--', 2)
            ax.plot(tspan, result[obs], lw=lw, ls=ls, color=colors[i % (len(results) // 2)], label=label)
            if obs in ['OB_tot', 'Tumor_tot']:
                ax.set_xlabel('Time (day)')
            if obs in ['Bone_tot', 'OB_tot']:
                ax.set_ylabel('Amount')
            ax.set_xticks(range(0, int(tspan[-1]) + 1, 7))
            ax.set_title(obs_labels[obs])
    # # Leave space on the right for the legend
    fig_all.subplots_adjust(right=0.8)
    fig_tumor.subplots_adjust(right=0.8)
    # Add figure-level legends
    handles, labels = axs_flat_all[0].get_legend_handles_labels()
    fig_all.legend(handles, labels, frameon=False, loc='outside right upper')
    handles, labels = axs_flat_tumor[0].get_legend_handles_labels()
    handles[0] = Line2D([0], [0], color='k', lw=3, linestyle='-')
    fig_tumor.legend(handles, labels, frameon=False, loc='outside right upper')

    plt.show()

    # EXAMPLE 1
    # Simple example using a ParameterCalibration object to run PyDREAM for a data set with only one experiment
    '''
    exp_data_file = '../DATA/TIBD_PopD_Mouse_Data.csv'
    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      tumor_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)
    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True)
    '''

    # EXAMPLE 2
    # A more complex example involving two experiments where some parameter values can vary between experiments
    '''
    exp_data_file = '../DATA/TIBD_PopD_Mouse_Data_2expts.csv'

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
                                      [tumor_injection] * 2,
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)
    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True, 
                   plot_tc_args={'separate_plots': False})
    '''

    # EXAMPLE 3
    # Create plots using output files already generated by PyDREAM using the calibrator.create_figures() function
    '''
    import glob
    import os
    # get the existing PyDREAM output files
    path = os.getcwd()  # path to where PyDREAM generated files are - default is current working directory
    logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))  # need to pass ALL 'logps' files
    samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))  # need to pass ALL 'params' files
    # create the ParameterCalibration object
    exp_data_file = os.path.join(path, 'TIBD_PopD_Data.csv')
    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      tumor_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)
    # call the 'create_figures' function
    calibrator.create_figures(logps_files, samples_files, obs_labels=obs_labels, show_plots=True,
                              plot_tc_args={'separate_plots': False})
    '''

    # EXAMPLE 4
    # Run simulations for the experimental condition and for predicting the effect of a perturbation, e.g., treatment
    # with zoledronic acid, using output files already generated by PyDREAM. To do this, we need to use the
    # 'ParameterCalibration.plot_timecourses' function directly. This function allows for an arbitrary number of
    # simulation protocols. Experimental data is optional and will only be plotted for sim_protocol[i] if i < n_expts.
    '''
    import glob
    import os
    # get the existing PyDREAM output files
    path = os.getcwd()  # path to where PyDREAM generated files are - default is current working directory
    logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))  # need to pass ALL 'logps' files
    samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))  # need to pass ALL 'params' files
    # create the ParameterCalibration object (not strictly necessary since 'plot_timecourses' is a class method)
    exp_data_file = '../DATA/TIBD_PopD_Mouse_Data.csv'
    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      tumor_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)
    # get various arguments needed for 'plot_timecourses'
    param_labels = [model.parameters[i].name for i in calibrator.parameter_idxs]
    param_samples = calibrator.plot_param_dist(samples_files, param_labels, show_plot=True)
    tumor_bisphos_injection = SequentialInjections(solver, t_equil=500, 
                                                   perturb_time_amount={'Tumor()': (0, 1), 'Bisphos()': (6, 1)})
    sim_protocols = [tumor_injection, tumor_bisphos_injection]  # 1 experiment, 1 perturbation
    tspans = [np.linspace(calibrator.tdata[0][0], calibrator.tdata[0][-1],
                          int((calibrator.tdata[0][-1] - calibrator.tdata[0][0]) * 10 + 1))] * len(sim_protocols)
    xlabel = 'time (day)'
    ylabels = [['amount (relative % BV/TV)'] + ['amount (fM)'] * 3] * len(sim_protocols)
    leg_labels = [['bone density', 'osteoclasts', 'osteoblasts', 'tumor cells'],
                  ['bone density (+ZA)', 'osteoclasts (+ZA)', 'osteoblasts (+ZA)', 'tumor cells (+ZA)']]
    calibrator.plot_timecourses(model, tspans, sim_protocols, param_samples, calibrator.parameter_idxs,
                                observables=calibrator.observables * len(sim_protocols), exp_data=calibrator.raw_data,
                                show_plot=True, separate_plots=False, xlabel=xlabel, ylabels=ylabels,
                                leg_labels=leg_labels)
    '''
