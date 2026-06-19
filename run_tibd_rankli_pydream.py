from MODELS.TIBD_PopD_v1 import model
from MODULES.perturbations import *
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *
from pysb.util import alias_model_components
import os
import pandas as pd
import sympy

this_dir = os.path.dirname(__file__)

expt_data_file = os.path.join(this_dir, 'DATA', 'Canon_OPGFc_2008.csv')
expt_data = pd.read_csv(expt_data_file)

# add_bisphosphonate_components()
add_RANKLi_components()
add_tumor_RANK_components(tumorK_migration=True)
add_tumor_migration_components()

solver = ScipyOdeSimulator(model)

alias_model_components()

Expression('Tumor_BLI', sympy.log(Tumor_tot, 10))  # 'Tumor_tot' is all tumor cells in the bone
Expression('OC_score', C_obs)  # For fitting to OC Score data from Fig. 2A of Canon et al. (2008)

alias_model_components()

# Canon2008_A: PBS-Tumor Control
canon_2008_A = SequentialInjections(
    solver, t_equil=5000, time_perturb_value={0: ('TumorK_out()', 1)})
expt_data_A = expt_data.loc[expt_data['expt_id'] == 'Canon2008_A', ['observable', 'time', 'average']]
protocol_A = ScaleBkProtocol(canon_2008_A, ['Tumor_BLI'], expt_data_A)

# Canon2008_B: OPG-Fc 3.0 mg/kg Preventative
canon_2008_B = SequentialInjections(
    solver, t_equil=5000, time_perturb_value={0: [('TumorK_out()', 1), ('RANKLi()', 10)]})
expt_data_B = expt_data.loc[expt_data['expt_id'] == 'Canon2008_B', ['observable', 'time', 'average']]
protocol_B = ScaleBkProtocol(canon_2008_B, ['Tumor_BLI'], expt_data_B)

# Canon2008_C: OPG-Fc 0.3 mg/kg Therapeutic
canon_2008_C = SequentialInjections(
    solver, t_equil=5000, time_perturb_value={0: ('TumorK_out()', 1), 7: ('RANKLi()', 1)})
expt_data_C = expt_data.loc[expt_data['expt_id'] == 'Canon2008_C', ['observable', 'time', 'average']]
protocol_C = ScaleBkProtocol(canon_2008_C, ['Tumor_BLI'], expt_data_C)

# Canon2008_D: OPG-Fc 3.0 mg/kg Therapeutic
canon_2008_D = SequentialInjections(
    solver, t_equil=5000, time_perturb_value={0: ('TumorK_out()', 1), 7: ('RANKLi()', 10)})
expt_data_D = expt_data.loc[expt_data['expt_id'] == 'Canon2008_D', ['observable', 'time', 'average']]
protocol_D = ScaleBkProtocol(canon_2008_D, ['Tumor_BLI'], expt_data_D)

# Canon2008_OCScore:
# 0: PBS-Tumor Control, 1: OPG-Fc 0.3 mg/kg Preventative, 2: OPG-Fc 3.0 mg/kg Preventative
tpv = [
    {0: ('TumorK_out()', 1)},
    {0: [('TumorK_out()', 1), ('RANKLi()', 1)]},  # set 1 = 0.3 mg/kg
    {0: [('TumorK_out()', 1), ('RANKLi()', 10)]}  # set 10 = 3.0 mg/kg
]
protocol_E = Canon2008Protocol(solver, t_equil=5000, time_perturb_value=tpv)

# Additional PyDREAM parameters
custom_priors = {'N': ('uniform', 0.3)}
no_sample = ['R_0', 'B_0', 'C_0', 'f0', 'IL', 'IO', 'IP_const', 'Bone_0', 'Tumor_0', 'CC_ON', 'nB', 'nC', 'ALLEE_ON',
             'A', 'Bisphos_0', 'RANKLi_0', 'TumorK_0', 'k_T_adapt', 'TumorK_out_0', 'Tumor_out_0']
obs_labels = {'Bone_tot': 'bone density', 'C_obs': 'osteoclasts', 'OB_tot': 'osteoblasts', 'Tumor_tot': 'tumor cells'}
param_expts_map = None
sim_protocols = [protocol_A, protocol_B, protocol_C, protocol_D, protocol_E]


if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      expt_data_file,
                                      sim_protocols,
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True}, restart=False)

    quit()

    import matplotlib.pyplot as plt

    # simulation results
    tspan = np.linspace(0, 25, 251)
    param_dir = os.path.join(this_dir, 'MODELS')
    params_df = pd.read_csv(os.path.join(param_dir, 'fit_params_PAR.csv'))  # parental parameter values
    params_dict = dict(zip(params_df['parameter'], params_df['value']))
    params_df = pd.read_csv(os.path.join(param_dir, 'fit_params_BONE.csv'), index_col='parameter')  # bone clone parameter values
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
    results = {'Canon2008_A': protocol_A.run(tspan, param_values).reshape(-1, len(tspan)),
               'Canon2008_B': protocol_B.run(tspan, param_values).reshape(-1, len(tspan)),
               'Canon2008_C': protocol_C.run(tspan, param_values).reshape(-1, len(tspan)),
               'Canon2008_D': protocol_D.run(tspan, param_values).reshape(-1, len(tspan)),
               'Canon2008_OCScore': protocol_E.run(tspan, param_values).reshape(-1, len(tspan))}

    labels_dict = {
        'OC_score': 'Osteoclasts score',
        'Tumor_BLI': 'Bioluminescence',
        'Bone_tot': 'Bone Density'
    }

    expt_ids = expt_data['expt_id'].unique()
    for expt_id in expt_ids:
        print('expt_id:', expt_id)
        data_expt_id = expt_data[expt_data['expt_id'] == expt_id]
        observables = data_expt_id['observable'].unique()
        for obs in observables:
            print('  obs:', obs)
            data_obs = data_expt_id[data_expt_id['observable'] == obs]
            plt.figure(obs, constrained_layout=True)

            # plot experimental data
            alt_expt_ids = data_obs['alt_expt_id'].unique()
            time_units = data_obs['time_units'].unique()[0]
            amount_units = data_obs['amount_units'].unique()
            colors = []
            for i, alt_expt_id in enumerate(alt_expt_ids):
                data_alt_expt_id = data_obs[data_obs['alt_expt_id'] == alt_expt_id]
                time = data_alt_expt_id['time']
                average = data_alt_expt_id['average']
                stderr = data_alt_expt_id['stderr']
                label = alt_expt_id
                if len(amount_units) > 1:
                    label += ' (%s)' % amount_units[i]
                p = plt.errorbar(time, average, yerr=stderr, fmt='o', capsize=6, label=label)
                colors.append(p.lines[0].get_color())

            # plot simulation results
            for i, result in enumerate(results[expt_id][obs]):
                plt.plot(tspan, result, color=colors[i], lw=3)
            # finish the plot
            plt.xlabel('Time (%s)' % time_units, fontsize=14)
            ylabel = labels_dict.get(obs, obs)
            if len(amount_units) == 1:
                ylabel += ' (%s)' % amount_units[i]
            plt.ylabel(ylabel, fontsize=14)
            plt.tick_params(labelsize=14)
            plt.xlim(left=-1)
            plt.legend(loc='best', fontsize=12)

            # plot C_obs() simulation results
            if expt_id == 'Canon2008_OCScore':
                plt.figure(constrained_layout=True)
                for i, result in enumerate(results[expt_id]['C_obs']):
                    plt.plot(tspan, result, color=colors[i], lw=3, label=alt_expt_ids[i])
                plt.xlabel('Time (%s)' % time_units, fontsize=14)
                plt.ylabel('Osteoclasts (fM)', fontsize=14)
                plt.tick_params(labelsize=14)
                plt.legend(loc='best', fontsize=12)


    plt.show()
