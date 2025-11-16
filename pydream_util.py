import numpy as np
import matplotlib.pyplot as plt
import os
from math import isclose
import importlib.util
import uuid
import pandas as pd
import math
import warnings
import glob
import pysb
from collections.abc import Iterable
import re


def is_num_pair(obj):
    try:
        # Convert to list so we can check length and reuse elements
        items = list(obj)
    except TypeError:
        # Not iterable
        return False

    return (
            isinstance(obj, Iterable) and
            len(items) == 2 and
            all(isinstance(x, (int, float)) for x in items)
    )


# convert None -> nan
def ensure_float(val):
    return np.array([val], dtype=float)[0]


# convert nan -> None
def nan_to_none(val):
    return None if pd.isnull(val) else val


# Helper function for getting the optimal number of columns for a multi-plot figure
def get_fig_ncols(ndims):
    if not isinstance(ndims, int) or ndims <= 0:
        raise ValueError("'ndims' must be a positive integer")
    r1 = round(math.sqrt(ndims))  # round() returns an int
    r2 = math.ceil(math.sqrt(ndims))  # math.ceil() also returns an int
    while r1 * r2 >= ndims:
        r1 -= 1
        r2 += 1
    return min(r1 + 1, r2 - 1)  # the smaller of the two integers is the # of columns


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


def plot_exp_data(filepath, separate_plots=True, show_plots=True, save_plots=False, **kwargs):

    # check if an array has a single value, return that value if so, throw an Exception if not
    def check_unique(arr):
        if len(arr) == 1:
            return arr[0]
        else:
            raise Exception('More than one value detected for quantity that should be unique:', arr)

    # process kwargs
    fontsizes = kwargs.get('fontsizes', {})
    labels_fs = fontsizes.get('labels', 12)
    ticks_fs = fontsizes.get('ticks', 12)
    legend_fs = fontsizes.get('legend', 10)
    legend_loc = kwargs.get('legend_loc', 'best')
    use_alt_expt_ids = kwargs.get('use_alt_expt_ids', False)

    # read in data
    data = pd.read_csv(filepath)
    observables = data['observable'].unique()
    expt_ids = data['expt_id'].unique()
    '''print(data.columns)
    print(observables)
    print(expt_ids)'''

    figures = []
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
    markers = ['o', '^', 's', 'v', '<', '>', 'D', 'p', 'H', '*']
    # loop over observables and experiments and plot data
    for obs in observables:
        data_obs = data[data['observable'] == obs]
        if not separate_plots:
            fig = plt.figure(obs, constrained_layout=True)
            if save_plots is not False and fig not in figures:
                figures.append(fig)
        for i, expt_id in enumerate(expt_ids):
            data_obs_expt = data_obs[data_obs['expt_id'] == expt_id]
            if separate_plots:
                fig = plt.figure('%s_%s' % (obs, expt_id), constrained_layout=True)
                if save_plots is not False and fig not in figures:
                    figures.append(fig)
            time = data_obs_expt['time']
            average = data_obs_expt['average']
            stderr = data_obs_expt['stderr']
            label = 'expt %s' % expt_id if not use_alt_expt_ids else check_unique(data_obs_expt['alt_expt_id'].unique())
            plt.errorbar(time, average, yerr=stderr, ls='', marker=markers[i % len(markers)],
                         color=colors[i % len(colors)], ms=10, capsize=6, label=label)
            time_units = check_unique(data_obs_expt['time_units'].unique())
            amount_units = check_unique(data_obs_expt['amount_units'].unique())
            plt.xlabel('Time (%s)' % time_units, fontsize=labels_fs)
            plt.ylabel('%s (%s)' % (obs, amount_units), fontsize=labels_fs)
            plt.tick_params(axis='both', which='major', labelsize=ticks_fs)
            plt.legend(loc=legend_loc, fontsize=legend_fs)

    if save_plots is not False:
        outpath = '.' if save_plots is True else save_plots
        prefix = os.path.splitext(os.path.basename(filepath))[0]
        for fig in figures:
            outfile = os.path.join(outpath, prefix + '_' + fig.get_label() + '.pdf')
            fig.savefig(outfile, format='pdf')

    if show_plots:
        plt.show()


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
    if not run_pydream_filename.endswith(".py"):
        run_pydream_filename += ".py"
    run_pydream_file = os.path.normpath(os.path.join(path, run_pydream_filename))
    module_name = f"_dynamic_{uuid.uuid4().hex}"  # unique name (arbitrary)
    spec = importlib.util.spec_from_file_location(module_name, run_pydream_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # get the path to the experimental data file referenced in the run_<...>_pydream.py file that's in the path
    expt_data_file = os.path.normpath(module.expt_data_file) if os.path.isabs(module.expt_data_file) else \
        os.path.normpath(os.path.join(path, module.expt_data_file))
    expt_data = pd.read_csv(expt_data_file) if os.path.exists(expt_data_file) else None
    print('expt_data:', list(expt_data.columns) if expt_data is not None else None)

    sim_data_file = os.path.join(path, 'SIM_DATA.csv')
    sim_data = pd.read_csv(sim_data_file) if os.path.exists(sim_data_file) else None
    print('sim_data:', list(sim_data.columns) if sim_data is not None else None)

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


def plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict=None, show_plot=True,
                          **kwargs):

    # process kwargs
    figsize = kwargs.get('figsize', (6.4 * 0.7, 4.8 * 0.9))

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
                    plt.figure('%d_%d_%s' % (i, j, obs), constrained_layout=True, figsize=figsize)
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
                    # plt.xticks(conc, conc)  # set the x-ticks to the experimental values TODO

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


def plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, label_dict=None, show_plot=True,
                         **kwargs):

    # process kwargs
    figsize = kwargs.get('figsize', (6.4 * 0.7, 4.8 * 0.9))

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
            plt.figure(constrained_layout=True, figsize=figsize)
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


def plot_from_simdata(basepath, directories, run_pydream_filename, expt_doses=None, tc_ids=None, label_dict=None,
                      show_plot=True, **kwargs):
    if expt_doses is not None:
        plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict, show_plot, **kwargs)
    if tc_ids is not None:
        plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, label_dict, show_plot, **kwargs)

    if expt_doses is None and tc_ids is None:
        warnings.warn("No drug doses or timecourse IDs were passed to `plot_from_simdata`")


def plot_pydream_output(dirpath, calibrator, max_iter=None, **kwargs):
    logps_files = glob.glob(os.path.join(dirpath, 'dreamzs*logps*'))
    samples_files = glob.glob(os.path.join(dirpath, 'dreamzs*params*'))
    if max_iter is not None:
        print('max iterations = ', max_iter)
        logps_files = [f for f in logps_files if int(re.search(r'_(\d+).npy$', f).group(1)) <= max_iter]
        samples_files = [f for f in samples_files if int(re.search(r'_(\d+).npy$', f).group(1)) <= max_iter]
    return calibrator.create_figures(logps_files, samples_files, save_plots=dirpath, **kwargs)


def detect_equilibrium(sim, tspan_linspace, rtol=1e-6, show_plot=False):
    fig1, ax1 = plt.subplots(constrained_layout=True, figsize=(6.4 * 1.5, 4.8))
    fig2, ax2 = plt.subplots(constrained_layout=True, figsize=(6.4 * 1.5, 4.8))

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
    initials = None

    t_start = tspan_linspace[0]
    delta_t = tspan_linspace[-1] - tspan_linspace[0]
    n_time_pts = len(tspan_linspace)

    ADD_LABEL = True
    STOP = False
    while not STOP:
        STOP =True
        tspan = np.linspace(t_start, t_start + delta_t, n_time_pts)
        output = sim.run(tspan=tspan, initials=initials)
        for n, name in enumerate(output.all.dtype.names):
            if name[0] == '_':  # this is a species
                conc = output.all[name]
                ax1.plot(tspan, conc, color=colors[n % len(colors)],
                         label=name[2:] if ADD_LABEL else None)  # conc
                d_conc_dt = [(conc[i] - conc[i - 1]) / conc[i - 1] / (tspan[i] - tspan[i - 1])
                             for i in range(1, len(conc))]
                ax2.plot(tspan[1:], d_conc_dt, color=colors[n % len(colors)],
                         label=name[2:] if ADD_LABEL else None)  # d_conc/dt
                # check for equilibrium
                if d_conc_dt[-1] > rtol:
                    STOP = False
        print('Simulated %g time units:' % (tspan[-1]), end=' ')
        if not STOP:
            print('Equilibration NOT detected, continuing...')
            t_start = tspan[-1]
            initials = output.species[-1]
            ADD_LABEL = False
        else:
            print('Equilibration detected!')
    # add a horizontal dashed line indicating equilibration
    ax2.axhline(rtol, linestyle='--', color='r', lw=2)
    ax2.annotate(text='Equil. threshold', xy=(0.7 * tspan[-1], 2.1 * rtol), xycoords='data', color='r',
                 fontsize=12, fontweight='bold')
    # finish up the plots
    for ax, ylabel in zip([ax1, ax2], ['[Conc]', 'Relative d[Conc]/dt']):
        ax.set_xlabel('Time')
        ax.set_ylabel(ylabel)
        ax.set_yscale('log')
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), ncols=3)

    if show_plot:
        plt.show()

    return tspan[-1]


def remove_unneeded_observables(model, obs_to_keep=None):
    obs_to_keep = [] if obs_to_keep is None else [obs_to_keep] if isinstance(obs_to_keep, str) else list(obs_to_keep)
    # remove any entries in 'obs_to_keep' that aren't observable names (e.g., expression names)
    obs_to_keep = [obs_name for obs_name in obs_to_keep if obs_name in [obs.name for obs in model.observables]]
    # get names of Observables in Expressions
    obs_names = []
    for expr in model.expressions:
        symbols = expr.expand_expr().free_symbols
        obs = [s for s in symbols if isinstance(s, pysb.core.Observable)]
        obs_names += [o.name for o in obs if o.name not in obs_names]
    obs_names = list(set(obs_names + obs_to_keep))  # create a set and then cast back to a list to remove duplicates
    # create the new set of observables
    new_observables = pysb.core.ComponentSet([model.observables[obs_name] for obs_name in obs_names])
    model.observables = new_observables
