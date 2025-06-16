import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import os
from math import isclose
import importlib
import pandas as pd


def is_in_array(x, arr, rel_tol=1e-9):
    for a in arr:
        if isclose(x, a, rel_tol=rel_tol):
            return True
    return False


def find_exact_index(tspan, target, rel_tol=1e-9):
    """Find the index of an exact floating-point match in a list.

    Raises an Exception if no exact match is found.
    """
    for i, t in enumerate(tspan):
        if isclose(t, target, rel_tol=rel_tol):
            return i
    raise ValueError(f"No exact match found for target {target} in tspan.")


def find_closest_index(tspan, target):
    best_idx = 0
    min_diff = abs(tspan[0] - target)

    for i in range(1, len(tspan)):  # Start from index 1 to avoid redundant checks
        diff = abs(tspan[i] - target)
        if diff < min_diff:
            best_idx, min_diff = i, diff
        else:
            break  # Exit early if difference starts increasing

    return best_idx


def get_exp_data(filepath, show_plots=False, save_plots=False):

    data = np.genfromtxt(filepath, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    # print(data.dtype.names)
    cond_cellTypes = np.unique([(d['Condition'], d['CellType']) for d in data], axis=0)
    # print(cond_cellTypes)
    return_data = {}

    for cc in cond_cellTypes:
        print(cc)

        # plot all data points
        time = [d['Week'] for d in data if d['Condition'] == cc[0] and d['CellType'] == cc[1]]
        conc = [d['Concentration'] for d in data if d['Condition'] == cc[0] and d['CellType'] == cc[1]]
        weeks = np.unique(time)

        # plot averages and standard errors
        avg_conc = [np.mean([conc[i] for i in range(len(conc)) if time[i] == w]) for w in weeks]
        se_conc = [[conc[i] for i in range(len(conc)) if time[i] == w] for w in weeks]
        se_conc = [np.std(x, ddof=1) / np.sqrt(len(x)) for x in se_conc]

        return_data[(cc[0], cc[1])] = (weeks, avg_conc, se_conc)

        if any([show_plots, save_plots]):
            plt.figure(cc[1], constrained_layout=True)
            p = plt.plot(time, conc, 'o', mfc='none')
            plt.errorbar(weeks, avg_conc, se_conc, marker='o', ms=10, capsize=10, lw=3, color=p[0].get_color(),
                         label='%s (%s)' % (cc[1], cc[0]))
            plt.xlabel('Week')
            ylabel = 'Relative BV/TV' if cc[1] == 'Bone' else 'Concentration (pM)'
            plt.ylabel(ylabel)
            plt.xticks(ticks=weeks, labels=weeks)
            # Remove error bars from legend
            handles, labels = plt.gca().get_legend_handles_labels()
            handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
            plt.legend(handles, labels, loc=0)

    if save_plots:  # save_plots can be False, True, or a string, which is the output filename
        filename = os.path.split(filepath)[-1].split('.')[0]+'.pdf' if save_plots is True else save_plots
        plt.savefig(filename, format='pdf')

    if show_plots:
        plt.show()

    return return_data


# Helper function for getting a good x-axis upper limit for dose-response plots
def round_up_nice(x):
    if x == 0:
        return 1
    exponent = np.floor(np.log10(x))
    fraction = x / 10 ** exponent
    nice_fraction = np.select(
        [fraction <= 1, fraction <= 2, fraction <= 5],
        [1, 2, 5],
        default=10
    )
    return nice_fraction * 10 ** exponent


# Helper function for getting a good x-axis lower limit for dose-response plots
def round_down_nice(x):
    if x <= 0:
        return 0
    exponent = np.floor(np.log10(x))
    fraction = x / 10 ** exponent
    nice_fraction = np.select(
        [fraction < 2, fraction < 5, fraction < 10],
        [1, 2, 5],
        default=10
    )
    return nice_fraction * 10 ** exponent


def get_sim_and_expt_data(path, run_pydream_filename):

    # import everything from run_<...>_pydream.py file that's in the path
    run_pydream_file = os.path.join(path, run_pydream_filename)
    import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
    module = importlib.import_module(str(import_string))  # import the module

    # get the path to the experimental data file referenced in the run_<...>_pydream.py file that's in the path
    exp_data_file = os.path.normpath(module.exp_data_file) if os.path.isabs(module.exp_data_file) else \
        os.path.normpath(os.path.join(path, module.exp_data_file))
    expt_data = pd.read_csv(exp_data_file)
    print('expt_data:', list(expt_data.columns))

    sim_data_file = os.path.join(path, 'SIM_DATA.csv')
    sim_data = pd.read_csv(sim_data_file)
    print('sim_data:', list(sim_data.columns))

    return sim_data, expt_data


def add_plot_to_fig(sim_data, expt_data, expt_id, xvals=None, label_dict=None, **kwargs):

    if xvals is None:
        xvals_sim = sim_data.loc[sim_data['sim_id'] == expt_id, 'time'].to_numpy()
        xvals_expt = expt_data.loc[expt_data['expt_id'] == expt_id, 'time'].to_numpy()
    else:
        xvals_sim = xvals_expt = xvals

    # simulation data
    legend_label = '%s (sim)' % label_dict.get(expt_id, expt_id)
    yval_min = sim_data.loc[sim_data['sim_id'] == expt_id, 'yval_min'].to_numpy()
    yval_max = sim_data.loc[sim_data['sim_id'] == expt_id, 'yval_max'].to_numpy()
    p = plt.plot(xvals_sim, (yval_min + yval_max) / 2, ls='--', label=legend_label)
    plt.fill_between(xvals_sim, yval_min, yval_max, alpha=0.25, color=p[0].get_color(), label='x')

    # experimental data
    legend_label = '%s (expt)' % label_dict.get(expt_id, expt_id)
    avg = expt_data.loc[expt_data['expt_id'] == expt_id, 'average'].to_numpy()
    stderr = expt_data.loc[expt_data['expt_id'] == expt_id, 'stderr'].to_numpy()
    plt.errorbar(xvals_expt, avg, yerr=stderr, fmt='o', ms=8, capsize=6, color=p[0].get_color(),
                 label=legend_label)

    plt.title(kwargs.get('title', None), fontsize=16, fontweight='bold')
    plt.xlabel(kwargs.get('xlabel'), fontsize=16)
    plt.ylabel(kwargs.get('ylabel'), fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best')


def merge_legend():
    # merge line and fill_between legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    n_sims = len([label for label in labels if 'sim' in label])
    handles = ([(handles[n], handles[n + 1]) for n in range(0, n_sims * 2 - 1, 2)] +
               list(handles[n_sims * 2:]))
    labels = [labels[n] for n in range(0, n_sims * 2 - 1, 2)] + list(labels[n_sims * 2:])
    # reorder handles and labels
    new_handles = []
    new_labels = []
    for n in range(n_sims):
        new_handles += [handles[:n_sims][n], handles[n_sims:][n]]
        new_labels += [labels[:n_sims][n], labels[n_sims:][n]]
    plt.legend(new_handles, new_labels, loc='best', fontsize=12)


def plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict=None, show_plot=True):

    if label_dict is None:
        label_dict = {}

    if isinstance(directories, str):
        directories = [directories]

    if len(directories) != len(expt_doses):
        raise ValueError("'expt_doses' and 'directories' must have the same length")

    for i, (directory, expt_doses_list) in enumerate(zip(directories, expt_doses)):
        print('Directory:', directory)

        if isinstance(expt_doses_list, dict):
            expt_doses_list = [expt_doses_list]

        # set the 'path' variable to the directory where the SIM_DATA.csv, run_<...>_pydream.py, and expt data file are
        path = os.path.join(basepath, directory)

        # get simulation and experimental data
        sim_data, expt_data = get_sim_and_expt_data(path, run_pydream_filename)

        for j, expt_doses_dict in enumerate(expt_doses_list):
            experiments = [key for key in expt_doses_dict.keys() if key not in ['xlabel', 'title']]
            all_observables = expt_data.loc[expt_data['expt_id'].isin(experiments), 'observable'].unique()
            suffix = {obs: [] for obs in all_observables}  # for the filename
            min_conc = np.inf
            max_conc = 0
            for expt in [key for key in expt_doses_dict.keys() if key not in ['xlabel', 'title']]:
                # get dose-response concentrations
                conc = expt_doses_dict[expt]
                if min(conc) < min_conc:
                    min_conc = min(conc)
                if max(conc) > max_conc:
                    max_conc = max(conc)
                # loop over observables for this experiment (could be different observables for different experiments)
                these_observables = expt_data.loc[expt_data['expt_id'] == expt, 'observable'].unique()
                for obs in these_observables:
                    # for each experiment, create a figure for each observable
                    plt.figure('%d_%d_%s' % (i, j, obs), constrained_layout=True, figsize=(6.4 * 0.7, 4.8 * 0.9))
                    # get units for the observable and make the y-axis label
                    y_units = expt_data.loc[
                        (expt_data['expt_id'] == expt) & (expt_data['observable'] == obs), 'amount_units'].unique()
                    if len(y_units) > 1:
                        raise Exception("Multiple units for observable '%s' in experiment '%s'" % (obs, expt))
                    y_label = '%s (%s)' % (label_dict.get(obs, obs), label_dict.get(y_units[0], y_units[0]))
                    # plot simulation and experimental data
                    add_plot_to_fig(sim_data, expt_data, expt, xvals=conc, label_dict=label_dict,
                                    title=expt_doses_dict.get('title', None), xlabel=expt_doses_dict.get('xlabel'),
                                    ylabel=y_label)
                    suffix[obs].append(expt.split('_'))  # for the filename

            # complete and save the figures for each observable
            for obs in all_observables:
                plt.figure('%d_%d_%s' % (i, j, obs))

                # adjust y-limits
                plt.ylim(bottom=0)

                # adjust x-limits
                xmin = round_down_nice(min_conc)
                xmax = round_up_nice(max_conc)
                buffer = 0.05 * (xmax - xmin)  # 5% buffer
                plt.xlim(left=xmin - buffer, right=xmax + buffer)

                # merge line and fill_between legend handles
                merge_legend()

                # save figure to file
                filename = 'fig_PyDREAM_DRC'
                for k in range(np.array(suffix[obs]).shape[1]):
                    filename += '_' + '_'.join(np.unique([x[k] for x in suffix[obs]]))
                filename += '_%s' % obs
                plt.savefig(os.path.join(str(path), filename))

    if show_plot:
        plt.show()


def plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, label_dict=None, show_plot=True):

    if label_dict is None:
        label_dict = {}

    if isinstance(directories, str):
        directories = [directories]

    if isinstance(tc_ids, str):
        tc_ids = [tc_ids]

    for directory in directories:
        print('Directory:', directory)

        # set the 'path' variable to the directory where the SIM_DATA.csv, run_<...>_pydream.py, and expt data file are
        path = os.path.join(basepath, directory)

        # get simulation and experimental data
        sim_data, expt_data = get_sim_and_expt_data(path, run_pydream_filename)

        # loop over observables
        for obs in sim_data['observable'].unique():
            plt.figure(constrained_layout=True, figsize=(6.4 * 0.7, 4.8 * 0.9))
            x_units = None
            y_units = None
            for tc_id in tc_ids:
                # get time units for the x-axis
                time_units = expt_data.loc[
                    (expt_data['expt_id'] == tc_id) & (expt_data['observable'] == obs), 'time_units'].unique()
                if len(time_units) > 1:
                    raise Exception("Time units not consistent for experiment '%s'" % tc_id)
                if x_units is not None and x_units != time_units[0]:
                    raise Exception("Time units not consistent for observable '%s'" % obs)
                x_units = time_units[0]
                x_label = 'Time (%s)' % label_dict.get(x_units, x_units)

                # get observable units for the y-axis
                amount_units = expt_data.loc[
                    (expt_data['expt_id'] == tc_id) & (expt_data['observable'] == obs), 'amount_units'].unique()
                if len(amount_units) > 1:
                    raise Exception("Observable units not consistent for experiment '%s'" % tc_id)
                if y_units is not None and y_units != amount_units[0]:
                    raise Exception("Observable units not consistent for observable '%s'" % obs)
                y_units = amount_units[0]
                y_label = '%s (%s)' % (label_dict.get(obs, obs), label_dict.get(y_units, y_units))

                # plot simulation and experimental data
                add_plot_to_fig(sim_data, expt_data, tc_id, label_dict=label_dict, xlabel=x_label, ylabel=y_label)

            # merge line and fill_between legend handles
            merge_legend()

            # save figure to file
            filename = 'fig_PyDREAM_tc_%s' % obs
            plt.savefig(os.path.join(str(path), filename))

    if show_plot:
        plt.show()
