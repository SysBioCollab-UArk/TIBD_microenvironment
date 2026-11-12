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
import csv
import multiprocessing


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


def plot_expt_data(filepath, figures=None, label_dict=None, add_titles=False, separate_plots=True, show_plots=True,
                   save_plots=False, **kwargs):


    # check if an array has a single value, return that value if so, throw an Exception if not
    def check_unique(arr):
        if len(arr) == 1:
            return arr[0]
        else:
            raise Exception('More than one value detected for quantity that should be unique:', arr)


    if label_dict is None:
        label_dict = {}

    # process kwargs
    fontsizes = kwargs.get('fontsizes', {})
    labels_fs = fontsizes.get('labels', 12)
    ticks_fs = fontsizes.get('ticks', 12)
    legend_fs = fontsizes.get('legend', 10)
    title_fs = fontsizes.get('title', 12)
    use_alt_expt_ids = kwargs.get('use_alt_expt_ids', False)

    legend_loc = kwargs.get('legend_loc', 'best')
    ms = kwargs.get('ms', 10)
    mec = kwargs.get('mec', None)
    mfc = kwargs.get('mfc', None)
    capsize = kwargs.get('capsize', 6)
    ecolor = kwargs.get('ecolor', None)
    # make these variables list-like
    legend_loc, ms, mec, mfc, capsize, ecolor = \
        [(v if np.ndim(v) != 0 else [v]) for v in (legend_loc, ms, mec, mfc, capsize, ecolor)]

    colors = kwargs.get('colors', None)
    markers = kwargs.get('markers', None)

    # read in data
    data = pd.read_csv(filepath)
    observables = data['observable'].unique()
    expt_ids = data['expt_id'].unique()
    '''
    print(data.columns)
    print(observables)
    print(expt_ids)
    '''

    if figures is None:
        figures = []
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
    if markers is None:
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
            if len(data_obs_expt) > 0:
                if separate_plots:
                    fig = plt.figure('%s_%s' % (obs, expt_id), constrained_layout=True)
                    if save_plots is not False and fig not in figures:
                        figures.append(fig)
                time = data_obs_expt['time']
                average = data_obs_expt['average']
                stderr = data_obs_expt['stderr']
                label = 'expt %s' % expt_id if not use_alt_expt_ids else \
                    check_unique(data_obs_expt['alt_expt_id'].unique())
                plt.errorbar(time, average, yerr=stderr, ls='', marker=markers[i % len(markers)],
                             color=colors[i % len(colors)], ms=ms[i % len(ms)], mec=mec[i % len(mec)],
                             mfc=mfc[i % len(mfc)], capsize=capsize[i % len(capsize)], ecolor=ecolor[i % len(ecolor)],
                             label=label_dict.get(label, label))
                time_units = check_unique(data_obs_expt['time_units'].unique())
                amount_units = check_unique(data_obs_expt['amount_units'].unique())
                plt.xlabel('Time (%s)' % label_dict.get(time_units, time_units), fontsize=labels_fs)
                plt.ylabel('%s (%s)' % (label_dict.get(obs, obs), label_dict.get(amount_units, amount_units)),
                           fontsize=labels_fs)
                plt.tick_params(axis='both', which='major', labelsize=ticks_fs)
                plt.legend(loc=legend_loc[i % len(legend_loc)], fontsize=legend_fs)
                if add_titles:
                    title = '%s' % fig.get_label()
                    title = label_dict.get(title, title)
                    plt.title(title, fontsize=title_fs, fontweight='bold')

    if save_plots is not False:
        outpath = '.' if save_plots is True else save_plots
        prefix = os.path.splitext(os.path.basename(filepath))[0]
        for fig in figures:
            outfile = os.path.join(outpath, prefix + '_' + fig.get_label() + '.pdf')
            fig.savefig(outfile, format='pdf')

    if show_plots:
        plt.show()

    return figures


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


def add_plot_to_fig(sim_data, expt_data, expt_id, type='fill_between', xvals=None, label_dict=None, **kwargs):
    """
    This function currently supports 'fill_between' and 'bar' plot types
    """
    _supported = ['fill_between', 'bar']

    if type not in _supported:
        raise ValueError('plot type must be one of', _supported)

    # process kwargs
    fontsizes = kwargs.get('fontsizes', {})
    title_fs = fontsizes.get('title', 14)
    labels_fs = fontsizes.get('labels', 12)
    ticks_fs = fontsizes.get('ticks', 12)
    leg_kwargs = kwargs.get('leg_kwargs', {})
    leg_kwargs['fontsize'] = fontsizes.get('legend', 10)

    if xvals is None:
        xvals_sim = sim_data.loc[sim_data['sim_id'] == expt_id, 'time'].to_numpy()
        xvals_expt = expt_data.loc[expt_data['expt_id'] == expt_id, 'time'].to_numpy()
    else:
        xvals_sim = xvals_expt = np.array(xvals)

    # simulation data
    legend_label = '%s (sim)' % label_dict.get(expt_id, expt_id)
    yval_min = sim_data.loc[sim_data['sim_id'] == expt_id, 'yval_min'].to_numpy()
    yval_max = sim_data.loc[sim_data['sim_id'] == expt_id, 'yval_max'].to_numpy()
    if type == 'fill_between':
        p = plt.plot(xvals_sim, (yval_min + yval_max) / 2, ls='--', label=legend_label)
        plt.fill_between(xvals_sim, yval_min, yval_max, alpha=0.25, color=p[0].get_color(), label='x')
    elif type == 'bar':
        p = plt.bar(xvals_sim - 0.2, (yval_min + yval_max) / 2, yerr=(yval_max - yval_min) / 2, width=0.4, capsize=8,
                    alpha=0.25, label=legend_label)

    # experimental data
    legend_label = '%s (expt)' % label_dict.get(expt_id, expt_id)
    avg = expt_data.loc[expt_data['expt_id'] == expt_id, 'average'].to_numpy()
    stderr = expt_data.loc[expt_data['expt_id'] == expt_id, 'stderr'].to_numpy()
    if type == 'fill_between':
        plt.errorbar(xvals_expt, avg, yerr=stderr, fmt='o', ms=8, capsize=6, color=p[0].get_color(), label=legend_label)
    elif type == 'bar':
        plt.bar(xvals_expt + 0.2, avg, yerr=stderr, width=0.4, capsize=8, color=p[0].get_facecolor(), alpha=1,
                label=legend_label)

    plt.title(kwargs.get('title', None), fontsize=title_fs, fontweight='bold')
    plt.xlabel(kwargs.get('xlabel'), fontsize=labels_fs)
    plt.ylabel(kwargs.get('ylabel'), fontsize=labels_fs)
    plt.xticks(fontsize=ticks_fs)
    plt.yticks(fontsize=ticks_fs)
    plt.legend(loc='best', **leg_kwargs)


def merge_legend(**leg_kwargs):
    # merge line and fill_between legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    n_sims = len([label for label in labels if 'sim' in label])
    n_expt = len([label for label in labels if 'expt' in label])
    n_merge = int((len(labels) - n_expt) / n_sims)

    handles = [tuple([handles[n + m] for m in range(n_merge)]) for n in range(0, n_sims * n_merge, n_merge)] + \
              list(handles[n_sims * n_merge:])
    labels = [labels[n] for n in range(0, n_sims * n_merge, n_merge)] + list(labels[n_sims * n_merge:])

    # reorder handles and labels
    new_handles = []
    new_labels = []
    for n in range(n_sims):
        new_handles += [handles[:n_sims][n], handles[n_sims:][n]]
        new_labels += [labels[:n_sims][n], labels[n_sims:][n]]
    plt.legend(new_handles, new_labels, loc='best', **leg_kwargs)


def plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict=None, show_plot=True,
                          **kwargs):
    # process kwargs
    figsize = kwargs.get('figsize', (6.4 * 0.7, 4.8 * 0.9))

    if label_dict is None:
        label_dict = {}

    if isinstance(directories, str):
        directories = [directories]

    if isinstance(expt_doses, dict):
        expt_doses = [expt_doses]

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
            experiments = [key for key in expt_doses_dict.keys()]  # if key not in ['xlabel', 'title']]
            for k, expt in enumerate(experiments):  # [key for key in expt_doses_dict.keys() if key not in ['xlabel', 'title']]:
                # get dose-response concentrations
                conc = expt_doses_dict[expt]['doses']
                plot_type = expt_doses_dict[expt].get('plot_type', 'fill_between')
                # loop over observables for this experiment (could be different observables for different experiments)
                these_observables = expt_data.loc[expt_data['expt_id'] == expt, 'observable'].unique()
                for obs in these_observables:
                    # for each experiment, create a figure for each observable
                    plt.figure('%d_%d_%s' % (i, k, obs), constrained_layout=True, figsize=figsize)
                    # get units for the observable and make the y-axis label
                    y_units = expt_data.loc[
                        (expt_data['expt_id'] == expt) & (expt_data['observable'] == obs), 'amount_units'].unique()
                    if len(y_units) > 1:
                        raise Exception("Multiple units for observable '%s' in experiment '%s'" % (obs, expt))
                    y_label = '%s (%s)' % (label_dict.get(obs, obs), label_dict.get(y_units[0], y_units[0]))
                    # plot simulation and experimental data
                    add_plot_to_fig(sim_data, expt_data, expt, type=plot_type, xvals=conc, label_dict=label_dict,
                                    title=expt_doses_dict[expt].get('title', None),
                                    xlabel=expt_doses_dict[expt].get('xlabel'), ylabel=y_label, **kwargs)

                # adjust x-limits
                xlim = expt_doses_dict[expt].get('xlim', None)
                if xlim is None and plot_type == 'fill_between':
                    xmin = round_down_nice(min(conc))
                    xmax = round_up_nice(max(conc))
                    buffer = 0.05 * (xmax - xmin)  # 5% buffer
                    xlim = (xmin - buffer, xmax + buffer)
                plt.xlim(xlim)

                # adjust y-limits
                ylim = expt_doses_dict[expt].get('ylim', None)
                if ylim is None:
                    ylocs = plt.yticks()[0]
                    delta_y = ylocs[-1] - ylocs[-2]
                    ymax = plt.ylim()[1]
                    top = ylocs[-1] + delta_y if ymax > ylocs[-1] else ylocs[-1]
                    ylim = (0, 1.05 * top)  # 5% buffer
                plt.ylim(ylim)

                # set x- and y-tick labels
                xticklabels = expt_doses_dict[expt].get('xticklabels', None)
                if xticklabels is not None:
                    plt.xticks(expt_doses_dict[expt]['doses'], xticklabels)

                yticklabels = expt_doses_dict[expt].get('yticklabels', None)
                if yticklabels is not None:
                    locs, _ = plt.yticks()
                    plt.yticks(locs, yticklabels)

                # merge line and fill_between legend handles
                merge_legend(**kwargs.get('leg_kwargs', {}))

                # save figure to file
                filename = 'fig_PyDREAM_DRC_%s_%s' % (expt, obs)
                plt.savefig(os.path.join(str(path), filename))

    if show_plot:
        plt.show()


def plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, plot_type='fill_between', label_dict=None,
                         show_plot=True, **kwargs):

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

        # check for the 'alt_expt_id' column in 'expt_data'. if it doesn't exist, copy the 'expt_id' column.
        if 'alt_expt_id' not in expt_data.columns:
            expt_data['alt_expt_id'] = expt_data['expt_id']

        # loop over observables
        for obs in sim_data['observable'].unique():
            x_units = None
            y_units = None
            # loop over timecourse IDs
            for tc_id in tc_ids:
                if obs not in expt_data[expt_data['expt_id'] == tc_id]['observable'].unique():
                    continue
                plt.figure(obs, constrained_layout=True, figsize=figsize)

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

                # get the alt expt IDs for this expt ID and observable pair
                alt_expt_ids = expt_data.loc[
                    (expt_data['expt_id'] == tc_id) & (expt_data['observable'] == obs), 'alt_expt_id'].unique()

                # get the expt data "chunks" for each alt expt ID
                expt_data_chunks = [
                    expt_data[
                        (expt_data['expt_id'] == tc_id) &
                        (expt_data['observable'] == obs) &
                        (expt_data['alt_expt_id'] == alt_expt_id)
                    ].copy() for alt_expt_id in alt_expt_ids]

                # compute how many rows each chunk has
                chunk_sizes = np.array([len(chunk) for chunk in expt_data_chunks], dtype=int)
                if chunk_sizes.sum() == 0:
                    raise ValueError("All expt_data chunks are empty; cannot infer how to split \'sim_data\'.")

                # split sim_data into matching chunks
                n_sim = len(sim_data)
                chunk_fracs = chunk_sizes / chunk_sizes.sum()
                sim_chunk_sizes = [int(n_sim * f) for f in chunk_fracs]
                start = 0
                sim_data_chunks = []
                for size in sim_chunk_sizes:
                    # use iloc to preserve original row order
                    sim_data_chunks.append(sim_data.iloc[start:start + size].copy())
                    start += size

                # plot simulation and experimental data
                for alt_expt_id, sim_data_chunk, expt_data_chunk in \
                        zip(alt_expt_ids, sim_data_chunks, expt_data_chunks):
                    label_dict[tc_id] = label_dict.get(alt_expt_id, alt_expt_id)  # map alt_expt_id to expt_id
                    add_plot_to_fig(sim_data_chunk, expt_data_chunk, tc_id, type=plot_type, label_dict=label_dict,
                                    xlabel=x_label, ylabel=y_label, **kwargs)

            # merge line and fill_between legend handles
            merge_legend(**kwargs.get('leg_kwargs', {}))

            # save figure to file
            filename = 'fig_PyDREAM_tc_%s' % obs
            plt.savefig(os.path.join(str(path), filename))

    if show_plot:
        plt.show()


def plot_from_simdata(basepath, directories, run_pydream_filename, expt_doses=None, tc_ids=None, label_dict=None,
                      show_plot=True, **kwargs):
    if expt_doses is not None:
        plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict=label_dict,
                              show_plot=show_plot, **kwargs)
    if tc_ids is not None:
        plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, label_dict=label_dict,
                             show_plot=show_plot, **kwargs)

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


def plot_GR_metrics(directory, threshold=1.2, par_names=None, only_unconverged=False, **kwargs):
    # process kwargs
    cmap = kwargs.pop('cmap', 'turbo')  # or 'hsv', 'viridis', etc.

    files = glob.glob(os.path.join(directory, '*Gelman*'))
    iterations = sorted([int(re.search(r"\d\d+", os.path.split(file)[1]).group()) for file in files])
    print('iterations:', iterations)

    GR = []
    for i in iterations:
        file = os.path.join(directory, 'dreamzs_5chain_GelmanRubin_iteration_%d.txt' % i)
        GR.append(np.genfromtxt(file, dtype=None, encoding="utf_8_sig"))
    GR = np.array(GR).T

    # gather GR values and parameter labels (check if 'only_converged'=True)
    GR_and_label = [(GR_array, 'par%d' % i if par_names is None else par_names[i])
                    for i, GR_array in enumerate(GR) if not only_unconverged or GR_array[-1] >= threshold]

    # only proceed if len(GR_and_label) > 0 (won't be true if run is converged and 'only_unconverged'=True)
    if len(GR_and_label) == 0:
        print("All parameters are converged and 'only_unconverged' is True. If want to plot GR metrics for all "
              "parameters for all iterations, set 'only_unconverged' to False.")
        return True
    else:
        max_GR = max([GR[-1] for GR, i in GR_and_label])
        ymax = math.ceil(max_GR * 10) / 10 + 0.1

        # get number of columns and rows for the multi-plot figure
        ncols = get_fig_ncols(len(GR_and_label))
        nrows = math.ceil(len(GR_and_label) / ncols)

        colormap = plt.get_cmap(cmap)
        colors = [colormap(0)] if len(GR_and_label) == 1 else \
            [colormap(i / (len(GR_and_label) - 1)) for i in range(len(GR_and_label))]

        fig = plt.figure(constrained_layout=True, figsize=(6.4 / 3 * ncols, 4.8 / 3.5 * nrows))
        ax_ref = None
        for i, ((GR_array, label), color) in enumerate(zip(GR_and_label, colors)):
            ax = fig.add_subplot(nrows, ncols, i + 1, sharex=ax_ref, sharey=ax_ref)
            if i == 0:
                ax_ref = ax
            ax.plot(iterations, GR_array, '-', color=color)
            ax.set_title(label, color=color, fontweight='bold')
            ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
            ax.label_outer()
            ax.axhline(threshold, ls='--', color='black', lw=1.5)
            ax.set_ylim(top=ymax)
        fig.supxlabel('Iteration')
        fig.supylabel('Gelman-Rubin metric')

        plt.show()

        return True if max_GR <= threshold else False


def load_pydream_history(history_file, nparams, nchains, history_thin=1, nseedchains=None):
    """
    Load a PyDREAM history file and reshape it to (stored_iterations, parameters, chains).

    Parameters
    ----------
    history_file : str
        Path to the .npy history file.
    nparams : int
        Number of adjustable parameters.
    nchains : int
        Number of MCMC chains.
    history_thin : int, optional
        Thinning interval used when saving history. Default is 1.
    nseedchains : int or None, optional
        Number of seed history records. If None, assumes PyDREAM default: nseedchains = 10 * nparams.

    Returns
    -------
    history_3d : np.ndarray
        Array of shape (stored_iterations, nparams, nchains).
    nstored_iterations : int
        Number of stored iterations per chain in the history file.
    niterations_total : int
        Total number of MCMC iterations per chain, inferred from thinning.
    """
    # Load flattened history
    history = np.load(history_file)

    # Default PyDREAM seed history
    if nseedchains is None:
        nseedchains = 10 * nparams

    # Reshape to rows of parameter vectors
    history = history.reshape(-1, nparams)

    # Remove seeded prior draws
    history = history[nseedchains:, :]

    # Number of stored rows remaining
    nrows = history.shape[0]

    if nrows % nchains != 0:
        raise ValueError(
            f"After removing seed history, number of rows ({nrows}) "
            f"is not divisible by nchains ({nchains})."
        )

    # Stored iterations per chain
    nstored_iterations = nrows // nchains

    # Total iterations per chain before thinning
    niterations_total = nstored_iterations * history_thin

    # Reshape to (nstored_iterations, nparams, nchains)
    history_3d = history.reshape(nstored_iterations, nchains, nparams).transpose(0, 2, 1)

    return history_3d, niterations_total


def weighted_percentile_linear_1d(values, weights, percentiles):
    """
    Mimic np.percentile(repeated_values, q, method='linear')
    without explicitly repeating values.

    Parameters
    ----------
    values : array-like, shape (n,)
        Unique data values.
    weights : array-like, shape (n,)
        Multiplicity of each unique value (counts).
    percentiles : scalar or array-like
        Percentiles in [0, 100].

    Returns
    -------
    result : ndarray
        Weighted percentile value(s), matching linear interpolation
        on the repeated sample.
    """
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=np.int64)
    percentiles = np.atleast_1d(percentiles).astype(float)

    if values.ndim != 1:
        raise ValueError("values must be 1D")
    if weights.ndim != 1:
        raise ValueError("weights must be 1D")
    if len(values) != len(weights):
        raise ValueError("values and weights must have the same length")
    if np.any(weights < 0):
        raise ValueError("weights must be nonnegative")
    if np.sum(weights) == 0:
        raise ValueError("sum of weights must be > 0")

    # Sort by value
    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]

    # Collapse any equal values that may still remain
    uniq_vals = []
    uniq_wts = []
    for v, w in zip(values, weights):
        if uniq_vals and v == uniq_vals[-1]:
            uniq_wts[-1] += w
        else:
            uniq_vals.append(v)
            uniq_wts.append(w)

    values = np.asarray(uniq_vals, dtype=float)
    weights = np.asarray(uniq_wts, dtype=np.int64)

    cum = np.cumsum(weights)
    N = int(cum[-1])

    # NumPy linear percentile uses virtual index h = (N - 1) * q / 100
    h = (N - 1) * (percentiles / 100.0)
    lo = np.floor(h).astype(np.int64)
    hi = np.ceil(h).astype(np.int64)
    frac = h - lo


    # Map expanded-sample index -> value without constructing repeated array
    def value_at_index(idx):
        pos = np.searchsorted(cum, idx + 1, side='left')
        return values[pos]


    v_lo = np.array([value_at_index(i) for i in lo], dtype=float)
    v_hi = np.array([value_at_index(i) for i in hi], dtype=float)

    return v_lo + frac * (v_hi - v_lo)


def weighted_percentile_linear_axis0(yvals, weights, percentiles):
    """
    Apply weighted_percentile_linear_1d along axis 0.

    Parameters
    ----------
    yvals : ndarray, shape (n_samples, ...)
        Simulation outputs for unique parameter sets.
    weights : ndarray, shape (n_samples,)
        Counts for each unique parameter set.
    percentiles : scalar or array-like
        Percentiles in [0, 100].

    Returns
    -------
    result : ndarray
        Shape (n_percentiles, ...) if multiple percentiles are requested,
        else shape (...).
    """
    yvals = np.asarray(yvals, dtype=float)
    weights = np.asarray(weights, dtype=np.int64)
    percentiles = np.atleast_1d(percentiles).astype(float)

    if yvals.shape[0] != len(weights):
        raise ValueError("weights length must match yvals.shape[0]")

    trailing_shape = yvals.shape[1:]
    yvals_2d = yvals.reshape(yvals.shape[0], -1)

    out = np.empty((len(percentiles), yvals_2d.shape[1]), dtype=float)

    for j in range(yvals_2d.shape[1]):
        out[:, j] = weighted_percentile_linear_1d(
            yvals_2d[:, j], weights, percentiles
        )

    out = out.reshape((len(percentiles),) + trailing_shape)

    if len(percentiles) == 1:
        return out[0]
    return out


def load_logps_and_param_samples(dirpath):
    logps_files = glob.glob(os.path.join(dirpath, 'dreamzs*logps*'))
    samples_files = glob.glob(os.path.join(dirpath, 'dreamzs*params*'))

    # error check
    if len(logps_files) != len(samples_files):
        raise ValueError("Number of logps (%d) and params (%d) files not equal" %
                         (len(logps_files), len(samples_files)))

    # get chains and iterations from file names
    chains = []
    chains_temp = []
    iterations = []
    iterations_temp = []
    for logps_file, samples_file in zip(logps_files, samples_files):
        m = re.search(r'chain_(\d+)_(\d+).npy$', logps_file)
        chains.append(int(m.group(1)))
        iterations.append(int(m.group(2)))
        m = re.search(r'chain_(\d+)_(\d+).npy$', samples_file)
        chains_temp.append(int(m.group(1)))
        iterations_temp.append(int(m.group(2)))
    chains = np.unique(chains)
    chains_temp = np.unique(chains_temp)
    iterations = np.unique(iterations)
    iterations_temp = np.unique(iterations_temp)

    # error check
    if not all(chains[i] == chains_temp[i] for i in range(len(chains))):
        raise ValueError("chains for logps and samples files not equal")
    if not all(iterations[i] == iterations_temp[i] for i in range(len(iterations))):
        raise ValueError("iterations for logps and samples files not equal")

    # read in logps and parameter samples from files
    log_ps_chain = []
    samples_chain = []
    for chain in chains:
        # get files and sort them numerically by number of iterations (using key=lambda function)
        files = sorted([file for file in logps_files if re.search(r'chain_%d' % chain, file)],
                       key=lambda f: int(re.search(r'(\d+).npy$', f).group(1)))
        log_ps_chain.append(np.concatenate(tuple(np.load(file) for file in files)).flatten())
        files = sorted([file for file in samples_files if re.search(r'chain_%d' % chain, file)],
                       key=lambda f: int(re.search(r'(\d+).npy$', f).group(1)))
        samples_chain.append(np.concatenate(tuple(np.load(file) for file in files)))

    return np.array(log_ps_chain), np.array(samples_chain)


def calculate_pydream_population_ess(dirpath, fburnin=0.5, cutoff=2.0, batch_size=None, min_batches=20, max_iter=None):
    """Estimate an autocorrelation-adjusted ESS for each PyDREAM parameter.

    The saved PyDREAM chains are treated as an interacting population. After removing burn-in independently from every
    chain, this function calculates the across-chain mean at each saved iteration index and estimates its long-run
    variance by non-overlapping batch means. For parameter ``p``, the resulting first-order population ESS is

        ESS_p = T * s_p**2 / Omega_p,

    where ``T`` is the number of retained iteration indices per chain, ``s_p**2`` is the marginal variance across all
    post-burn-in chain samples, and ``Omega_p`` is the estimated long-run variance of the across-chain mean series.

    If ``cutoff`` is not None, the function then applies the same log-posterior-density threshold used by
    ``ParameterCalibration`` and approximates the ESS of the filtered ensemble by multiplying the post-burn-in ESS by
    the fraction of parameter vectors retained. This last step assumes that the sampling efficiency is approximately
    unchanged by the filter and is intended as a descriptive correction for the self-distance calculation, not as a
    formal MCMC diagnostic.

    Parameters
    ----------
    dirpath : str or path-like
        Directory containing the saved ``sampled_params`` and ``logps`` files.
    fburnin : float, optional
        Fraction removed independently from the beginning of every chain.
    cutoff : float or None, optional
        Retain post-burn-in samples whose saved log posterior is greater than ``mean - cutoff * standard_deviation``.
        Use None for no filtering.
    batch_size : int or None, optional
        Number of saved iteration indices per non-overlapping batch. The default is floor(sqrt(T)).
    min_batches : int, optional
        Warn if fewer than this many complete batches are available.
    max_iter : int or None, optional
        Use samples only through the latest saved checkpoint whose cumulative
        iteration number is less than or equal to ``max_iter``. This matches
        the checkpoint-selection behavior of ``plot_pydream_output``.

    Returns
    -------
    dict
        Scalar run metadata and one-dimensional arrays of length ``P`` for ``marginal_variance``, ``long_run_variance``,
        ``ess_postburn``, ``sampling_efficiency``, and ``ess_retained``.
    """
    if not 0 <= fburnin < 1:
        raise ValueError("fburnin must satisfy 0 <= fburnin < 1")
    if cutoff is not None and cutoff < 0:
        raise ValueError("cutoff must be nonnegative or None")
    if min_batches < 2:
        raise ValueError("min_batches must be at least 2")

    log_ps_chain, samples_chain = load_logps_and_param_samples(dirpath)

    selected_checkpoint = samples_chain.shape[1]
    if max_iter is not None:
        if not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
            raise ValueError("max_iter must be a positive integer or None")
        sample_files = glob.glob(os.path.join(dirpath, 'dreamzs*sampled_params_chain_*_*.npy'))
        checkpoints = []
        for filename in sample_files:
            match = re.search(r'chain_\d+_(\d+)\.npy$', filename)
            if match is not None and int(match.group(1)) <= max_iter:
                checkpoints.append(int(match.group(1)))
        if not checkpoints:
            raise ValueError("No saved PyDREAM checkpoint is less than or equal to max_iter=%d" % max_iter)
        selected_checkpoint = max(checkpoints)
        if selected_checkpoint > samples_chain.shape[1]:
            raise ValueError("Saved checkpoint %d exceeds the %d samples loaded per chain" %
                             (selected_checkpoint, samples_chain.shape[1]))
        samples_chain = samples_chain[:, :selected_checkpoint, :]
        log_ps_chain = log_ps_chain[:, :selected_checkpoint]

    if samples_chain.ndim != 3:
        raise ValueError("Expected sampled parameters with shape (chains, iterations, parameters); found %s" %
                         samples_chain.shape)
    if log_ps_chain.ndim != 2:
        raise ValueError("Expected log posterior values with shape (chains, iterations); found %s" % log_ps_chain.shape)
    if log_ps_chain.shape != samples_chain.shape[:2]:
        raise ValueError("Saved log posterior and sampled-parameter arrays do not have matching chain and iteration "
                         "dimensions")

    n_chains, n_iterations, n_parameters = samples_chain.shape
    n_burnin = int(math.floor(fburnin * n_iterations))
    post_samples = samples_chain[:, n_burnin:, :]
    post_log_ps = log_ps_chain[:, n_burnin:]
    n_iterations_postburn = post_samples.shape[1]

    if n_iterations_postburn < 2:
        raise ValueError("At least two saved iterations must remain after burn-in")

    n_raw = n_chains * n_iterations_postburn
    population_mean = np.mean(post_samples, axis=0, dtype=np.float64)
    marginal_variance = np.var(post_samples, axis=(0, 1), ddof=1, dtype=np.float64)

    if batch_size is None:
        batch_size = int(math.floor(math.sqrt(n_iterations_postburn)))
    if not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")

    n_batches = n_iterations_postburn // batch_size
    if n_batches < 2:
        raise ValueError("batch_size=%d leaves fewer than two complete batches" % batch_size)
    if n_batches < min_batches:
        warnings.warn("Only %d complete batches are available; the population ESS estimate may be unstable" % n_batches,
                      RuntimeWarning)

    n_batched_iterations = n_batches * batch_size
    batch_means = population_mean[:n_batched_iterations].reshape(n_batches, batch_size, n_parameters).mean(axis=1)
    long_run_variance = batch_size * np.var(batch_means, axis=0, ddof=1, dtype=np.float64)

    if np.any(marginal_variance <= 0):
        bad = np.flatnonzero(marginal_variance <= 0)
        raise ValueError("Cannot calculate ESS for zero-variance parameter columns: %s" % bad.tolist())

    with np.errstate(divide='ignore', invalid='ignore'):
        ess_postburn = (n_iterations_postburn * marginal_variance / long_run_variance)

    # A zero batch-means variance implies an estimate above the available raw sample count; cap it at that count, as for
    # all finite estimates.
    ess_postburn = np.where(np.isfinite(ess_postburn), ess_postburn, float(n_raw))
    ess_postburn = np.clip(ess_postburn, 1.0, float(n_raw))

    posterior_cutoff_value = None
    if cutoff is None:
        n_retained = n_raw
    else:
        finite_log_ps = np.isfinite(post_log_ps)
        finite_values = post_log_ps[finite_log_ps]
        if finite_values.size == 0:
            raise ValueError("No finite post-burn-in log posterior values found")
        posterior_cutoff_value = (np.mean(finite_values) - cutoff * np.std(finite_values))
        keep_mask = finite_log_ps & (post_log_ps > posterior_cutoff_value)
        n_retained = int(np.count_nonzero(keep_mask))
        if n_retained == 0:
            raise ValueError("The log-posterior-density filter retained no samples")

    retained_fraction = n_retained / float(n_raw)
    sampling_efficiency = ess_postburn / float(n_raw)
    ess_retained = ess_postburn * retained_fraction
    ess_retained = np.clip(ess_retained, 1.0, float(n_retained))

    return {
        'n_chains': n_chains,
        'n_parameters': n_parameters,
        'n_iterations_per_chain': n_iterations,
        'selected_checkpoint': selected_checkpoint,
        'max_iter_requested': max_iter,
        'n_burnin_per_chain': n_burnin,
        'n_iterations_postburn_per_chain': n_iterations_postburn,
        'n_raw_postburn': n_raw,
        'n_retained': n_retained,
        'retained_fraction': retained_fraction,
        'posterior_cutoff': cutoff,
        'posterior_cutoff_value': posterior_cutoff_value,
        'batch_size': batch_size,
        'n_batches': n_batches,
        'n_batched_iterations': n_batched_iterations,
        'marginal_variance': marginal_variance,
        'long_run_variance': long_run_variance,
        'ess_postburn': ess_postburn,
        'sampling_efficiency': sampling_efficiency,
        'ess_retained': ess_retained,
    }


def save_self_distance_comparison(dirpath, param_labels, hist_labels, raw_sample_counts, effective_sample_sizes,
                                  self_distances_ess, max_iter=None, filename=None):
    """Save raw-count and ESS-adjusted self-distance comparisons to a CSV file."""

    param_labels = list(param_labels)
    hist_labels = tuple(hist_labels)
    raw_sample_counts = np.asarray(raw_sample_counts, dtype=float)
    effective_sample_sizes = np.asarray(effective_sample_sizes, dtype=float)
    self_distances_ess = np.asarray(self_distances_ess, dtype=float)

    n_parameters = len(param_labels)
    expected_shape = (2, n_parameters)

    if len(hist_labels) != 2:
        raise ValueError("'hist_labels' must contain exactly two labels")

    if raw_sample_counts.shape != (2,):
        raise ValueError("'raw_sample_counts' must have shape (2,); found %s" % raw_sample_counts.shape)

    if effective_sample_sizes.shape != expected_shape:
        raise ValueError("'effective_sample_sizes' must have shape %s; found %s" %
                         (expected_shape, effective_sample_sizes.shape))

    if self_distances_ess.shape != expected_shape:
        raise ValueError("'self_distances_ess' must have shape %s; found %s" %
                         (expected_shape, self_distances_ess.shape))

    if np.any(raw_sample_counts <= 0):
        raise ValueError("All raw sample counts must be positive")

    if np.any(effective_sample_sizes <= 0):
        raise ValueError("All effective sample sizes must be positive")

    self_distances_raw = self_distances_ess * np.sqrt(effective_sample_sizes / raw_sample_counts[:, np.newaxis])

    if filename is None:
        iteration_label = 'all' if max_iter is None else str(max_iter)
        first_group = str(hist_labels[0]).replace(' ', '_')
        second_group = str(hist_labels[1]).replace(' ', '_')
        filename = 'PyDREAM_self_distance_comparison_%s_%s_%s.csv' % (first_group, second_group, iteration_label)

    output_filepath = os.path.join(dirpath, filename)

    with open(output_filepath, 'w', encoding='utf-8', newline='') as output_file:
        csv_writer = csv.writer(output_file)
        csv_writer.writerow(['parameter', 'group', 'retained_raw_n', 'adjusted_ess', 'ess_raw_ratio',
                             'self_distance_raw', 'self_distance_ess', 'displayed_criterion'])

        for parameter_idx, parameter_label in enumerate(param_labels):
            displayed_criterion = np.max(self_distances_ess[:, parameter_idx])

            for distribution_idx in range(2):
                raw_count = int(raw_sample_counts[distribution_idx])
                adjusted_ess = effective_sample_sizes[distribution_idx, parameter_idx]
                efficiency = adjusted_ess / raw_count
                csv_writer.writerow([parameter_label, hist_labels[distribution_idx], raw_count, '%.1f' % adjusted_ess,
                                     '%.6f' % efficiency, '%.6f' % self_distances_raw[distribution_idx, parameter_idx],
                                     '%.6f' % self_distances_ess[distribution_idx, parameter_idx],
                                     '%.6f' % displayed_criterion])

    print('Saving...')
    print('    %s' % filename)
    print('...to %s' % os.path.abspath(dirpath))

    return output_filepath


def _parse_self_distance_sensitivity_filename(filepath):
    """Return the bandwidth and numerical-grid settings encoded in a CSV filename."""

    filename = os.path.basename(filepath)
    match = re.search(r'_(-?\d+(?:\.\d+)?)_(\d+)_(\d+)\.csv$', filename)
    if match is None:
        raise ValueError("Could not parse bw_adjust, n_self_dist_pts, and n_hist_dist_pts from '%s'" % filename)

    return float(match.group(1)), int(match.group(2)), int(match.group(3))


def _load_self_distance_sensitivity_files(dirpath, csv_files, param_labels=None, group_order=None):
    """Read and validate self-distance sensitivity CSV files."""

    csv_files = list(csv_files)
    if not csv_files:
        raise ValueError("'csv_files' must contain at least one self-distance comparison file")

    required_columns = {
        'parameter', 'group', 'retained_raw_n', 'adjusted_ess', 'self_distance_raw', 'self_distance_ess'
    }
    frames = []
    settings_seen = set()

    for csv_file in csv_files:
        filepath = os.fspath(csv_file)
        if not os.path.exists(filepath):
            filepath = os.path.join(dirpath, filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError("Self-distance comparison file not found: %s" % filepath)

        bw_value, n_self_value, n_hist_value = _parse_self_distance_sensitivity_filename(filepath)
        setting = (bw_value, n_self_value, n_hist_value)
        if setting in settings_seen:
            raise ValueError("More than one self-distance comparison file was supplied for settings %s" %
                             (setting, ))
        settings_seen.add(setting)

        frame = pd.read_csv(filepath)
        missing_columns = required_columns.difference(frame.columns)
        if missing_columns:
            raise ValueError("%s is missing required columns: %s" %
                             (os.path.basename(filepath), ', '.join(sorted(missing_columns))))

        frame = frame.copy()
        frame['parameter'] = frame['parameter'].astype(str)
        frame['group'] = frame['group'].astype(str)
        for column in ('retained_raw_n', 'adjusted_ess', 'self_distance_raw', 'self_distance_ess'):
            frame[column] = pd.to_numeric(frame[column], errors='raise')

        if frame.duplicated(['parameter', 'group']).any():
            raise ValueError("%s contains duplicate parameter/group rows" % os.path.basename(filepath))

        frame['bw_adjust'] = bw_value
        frame['n_self_dist_pts'] = n_self_value
        frame['n_hist_dist_pts'] = n_hist_value
        frame['source_file'] = os.path.basename(filepath)
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    settings = sorted(settings_seen, key=lambda item: (item[0], item[1], item[2]))

    if param_labels is None:
        param_labels = list(dict.fromkeys(frames[0]['parameter'].tolist()))
    else:
        param_labels = list(param_labels)
    if len(param_labels) != len(set(param_labels)):
        raise ValueError("'param_labels' contains duplicate parameter names")

    first_groups = frames[0]['group'].tolist()
    if group_order is None:
        group_order = tuple(dict.fromkeys(first_groups))
    else:
        group_order = tuple(group_order)
    if len(group_order) != 2 or len(set(group_order)) != 2:
        raise ValueError("'group_order' must contain exactly two distinct group labels")

    expected_parameters = set(param_labels)
    expected_groups = set(group_order)
    for frame in frames:
        source_file = frame['source_file'].iloc[0]
        if set(frame['parameter']) != expected_parameters:
            raise ValueError("Parameter names in %s do not match 'param_labels'" % source_file)
        if set(frame['group']) != expected_groups:
            raise ValueError("Group names in %s do not match 'group_order'" % source_file)
        if len(frame) != len(param_labels) * 2:
            raise ValueError("%s does not contain one row per parameter and group" % source_file)

    for (parameter, group), subset in combined.groupby(['parameter', 'group'], sort=False):
        for column in ('retained_raw_n', 'adjusted_ess'):
            values = subset[column].to_numpy(dtype=float)
            if not np.allclose(values, values[0], rtol=0, atol=5e-7):
                raise ValueError("%s varies across sensitivity settings for %s, %s" %
                                 (column, parameter, group))

    return combined, settings, param_labels, group_order


def _self_distance_values(combined, settings, parameter, group, value_column):
    """Return values for one parameter, group, and metric in sensitivity-setting order."""

    values = []
    for bw_value, n_self_value, n_hist_value in settings:
        mask = (
            (combined['parameter'] == parameter) &
            (combined['group'] == group) &
            np.isclose(combined['bw_adjust'], bw_value) &
            (combined['n_self_dist_pts'] == n_self_value) &
            (combined['n_hist_dist_pts'] == n_hist_value)
        )
        matches = combined.loc[mask, value_column]
        if len(matches) != 1:
            raise ValueError("Expected one %s value for %s, %s, settings %s; found %d" %
                             (value_column, parameter, group,
                              (bw_value, n_self_value, n_hist_value), len(matches)))
        values.append(float(matches.iloc[0]))

    return np.asarray(values, dtype=float)


def _find_self_distance_setting_index(settings, selected_settings):
    """Return the index of a selected sensitivity setting, or None when none is selected."""

    if selected_settings is None:
        return None

    selected_bw, selected_n_self, selected_n_hist = selected_settings
    for setting_idx, (bw_value, n_self_value, n_hist_value) in enumerate(settings):
        if (np.isclose(bw_value, selected_bw) and n_self_value == selected_n_self and
                n_hist_value == selected_n_hist):
            return setting_idx

    raise ValueError("The selected settings %s were not found in the supplied CSV files" %
                     (selected_settings, ))


def _write_self_distance_dataset(output_filepath, combined, settings, param_labels, group_order,
                                 group_display_labels, table_param_labels, selected_settings):
    """Write the self-distance sensitivity results and ESS values to an Excel workbook."""

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
        from openpyxl.utils import get_column_letter
    except ImportError as error:
        raise ImportError("Creating Dataset S2 requires openpyxl") from error

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = 'Self-distance'
    worksheet.sheet_view.showGridLines = False

    thin_gray = Side(style='thin', color='B7B7B7')
    medium_black = Side(style='medium', color='000000')
    header_fill = PatternFill(fill_type='solid', fgColor='D9EAF7')
    selected_font = Font(name='Arial', size=10, color='FF0000')
    body_font = Font(name='Arial', size=10)
    bold_font = Font(name='Arial', size=10, bold=True)
    centered = Alignment(horizontal='center', vertical='center')
    left_aligned = Alignment(horizontal='left', vertical='center')

    metric_definitions = (
        ('Raw-count E[D^self]', 'self_distance_raw', group_order[0]),
        ('Raw-count E[D^self]', 'self_distance_raw', group_order[1]),
        ('ESS-adjusted E[D^self]', 'self_distance_ess', group_order[0]),
        ('ESS-adjusted E[D^self]', 'self_distance_ess', group_order[1]),
    )
    first_data_column = 3
    last_data_column = first_data_column + len(settings) - 1

    selected_setting_idx = _find_self_distance_setting_index(settings, selected_settings)
    selected_column = (None if selected_setting_idx is None else
                       first_data_column + selected_setting_idx)

    for parameter_idx, parameter in enumerate(param_labels):
        heading_row = 1 + 6 * parameter_idx
        worksheet.merge_cells(start_row=heading_row, start_column=first_data_column,
                              end_row=heading_row, end_column=last_data_column)
        heading_cell = worksheet.cell(heading_row, first_data_column, table_param_labels.get(parameter, parameter))
        heading_cell.font = body_font
        heading_cell.alignment = centered
        for column_idx in range(first_data_column, last_data_column + 1):
            worksheet.cell(heading_row, column_idx).border = Border(bottom=medium_black)

        for metric_idx, (metric_label, value_column, group) in enumerate(metric_definitions, start=1):
            row_idx = heading_row + metric_idx
            label = '%s, %s' % (metric_label, group_display_labels.get(group, group))
            worksheet.cell(row_idx, 2, label).font = body_font
            worksheet.cell(row_idx, 2).alignment = left_aligned

            values = _self_distance_values(combined, settings, parameter, group, value_column)
            for setting_idx, value in enumerate(values):
                column_idx = first_data_column + setting_idx
                cell = worksheet.cell(row_idx, column_idx, value)
                cell.number_format = '0.000000'
                cell.font = selected_font if column_idx == selected_column else body_font
                cell.alignment = centered

    key_first_column = last_data_column + 2
    key_last_column = key_first_column + 3
    worksheet.merge_cells(start_row=1, start_column=key_first_column, end_row=1, end_column=key_last_column)
    key_title = worksheet.cell(1, key_first_column, 'Sensitivity-setting column key')
    key_title.fill = header_fill
    key_title.font = Font(name='Arial', size=10, bold=True)
    key_title.alignment = centered

    key_headers = ('Column', 'bw_adjust', 'n_self_dist_pts', 'n_hist_dist_pts')
    for offset, header in enumerate(key_headers):
        cell = worksheet.cell(2, key_first_column + offset, header)
        cell.fill = header_fill
        cell.font = bold_font
        cell.alignment = centered
        cell.border = Border(bottom=thin_gray)

    for setting_idx, (bw_value, n_self_value, n_hist_value) in enumerate(settings):
        row_idx = 3 + setting_idx
        data_column = first_data_column + setting_idx
        column_label = get_column_letter(data_column)
        is_selected = data_column == selected_column
        if is_selected:
            column_label = '*' + column_label
        values = (column_label, int(bw_value) if float(bw_value).is_integer() else bw_value,
                  n_self_value, n_hist_value)
        for offset, value in enumerate(values):
            cell = worksheet.cell(row_idx, key_first_column + offset, value)
            cell.font = selected_font if is_selected else body_font
            cell.alignment = centered

    note_row = 3 + len(settings)
    if selected_column is not None:
        note_cell = worksheet.cell(note_row, key_first_column, '*Settings used in Figures 4 and S6')
        note_cell.font = Font(name='Arial', size=10, color='FF0000')
        note_cell.alignment = left_aligned

    ess_title_row = note_row + 2
    ess_last_column = key_first_column + 2
    worksheet.merge_cells(start_row=ess_title_row, start_column=key_first_column,
                          end_row=ess_title_row, end_column=ess_last_column)
    ess_title = worksheet.cell(ess_title_row, key_first_column, 'Adjusted effective sample size (ESS)*')
    ess_title.font = bold_font
    ess_title.alignment = centered
    for column_idx in range(key_first_column, ess_last_column + 1):
        worksheet.cell(ess_title_row, column_idx).border = Border(bottom=medium_black)

    ess_header_row = ess_title_row + 1
    ess_headers = ('Parameter', group_display_labels.get(group_order[0], group_order[0]),
                   group_display_labels.get(group_order[1], group_order[1]))
    for offset, header in enumerate(ess_headers):
        cell = worksheet.cell(ess_header_row, key_first_column + offset, header)
        cell.font = bold_font
        cell.alignment = centered

    raw_counts = []
    for parameter_idx, parameter in enumerate(param_labels):
        row_idx = ess_header_row + 1 + parameter_idx
        parameter_cell = worksheet.cell(row_idx, key_first_column, table_param_labels.get(parameter, parameter))
        parameter_cell.font = body_font
        parameter_cell.alignment = centered

        for group_idx, group in enumerate(group_order, start=1):
            subset = combined[(combined['parameter'] == parameter) & (combined['group'] == group)]
            adjusted_ess = float(subset['adjusted_ess'].iloc[0])
            raw_count = int(round(float(subset['retained_raw_n'].iloc[0])))
            raw_counts.append(raw_count)
            cell = worksheet.cell(row_idx, key_first_column + group_idx, adjusted_ess)
            cell.number_format = '#,##0'
            cell.font = body_font
            cell.alignment = centered

    if len(set(raw_counts)) != 1:
        raise ValueError("The retained raw sample count is not identical across all parameters and groups")
    raw_count_note_row = ess_header_row + 1 + len(param_labels)
    raw_count_note = worksheet.cell(raw_count_note_row, key_first_column,
                                    '*Raw retained parameter set samples = %s' % format(raw_counts[0], ','))
    raw_count_note.font = body_font
    raw_count_note.alignment = left_aligned

    worksheet.column_dimensions['A'].width = 2
    worksheet.column_dimensions['B'].width = 38
    for column_idx in range(first_data_column, last_data_column + 1):
        worksheet.column_dimensions[get_column_letter(column_idx)].width = 13
    worksheet.column_dimensions[get_column_letter(last_data_column + 1)].width = 3
    worksheet.column_dimensions[get_column_letter(key_first_column)].width = 18
    for column_idx in range(key_first_column + 1, key_last_column + 1):
        worksheet.column_dimensions[get_column_letter(column_idx)].width = 19

    for row in worksheet.iter_rows(min_row=1, max_row=worksheet.max_row,
                                   min_col=1, max_col=worksheet.max_column):
        for cell in row:
            if cell.value is not None and cell.font == Font():
                cell.font = body_font

    workbook.save(output_filepath)


def _plot_self_distance_sensitivity(output_filepath, combined, settings, param_labels, group_order,
                                    group_display_labels, figure_param_labels, show_plot,
                                    selected_settings=None, figure_options=None):
    """Create a compact, responsive self-distance sensitivity summary figure."""

    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    defaults = {
        'column_width': 1.50,
        'inter_group_width': 0.30,
        'row_height': 0.90,
        'footer_height': 0.90,
        'minimum_width': 7.40,
        'label_padding': 8,
        'label_fontsize': 15,
        'tick_fontsize': 10,
        'group_fontsize': 12,
        'footer_fontsize': 13,
        'range_height': 0.46,
        'point_size': 23,
        'selected_point_size': 35,
        'selected_point_edgecolor': 'black',
        'selected_point_linewidth': 1.0,
        'maximum_ticks': 2,
        'tick_significant_digits': 4,
        'jitter_span': 0.20,
        'dpi': 350,
        'save_padding': 0.04,
    }
    figure_options = {} if figure_options is None else dict(figure_options)
    unknown_options = set(figure_options).difference(defaults)
    if unknown_options:
        raise ValueError("Unknown figure_options: %s" % ', '.join(sorted(unknown_options)))
    options = defaults.copy()
    options.update(figure_options)

    if not isinstance(options['tick_significant_digits'], (int, np.integer)) or \
            options['tick_significant_digits'] < 1:
        raise ValueError("'tick_significant_digits' must be a positive integer")

    if not param_labels:
        raise ValueError("'param_labels' must contain at least one parameter")
    if not settings:
        raise ValueError("At least one sensitivity setting is required")

    selected_setting_idx = _find_self_distance_setting_index(settings, selected_settings)

    metrics = (
        ('self_distance_raw', group_order[0]),
        ('self_distance_raw', group_order[1]),
        ('self_distance_ess', group_order[0]),
        ('self_distance_ess', group_order[1]),
    )
    values = {
        parameter: [
            _self_distance_values(combined, settings, parameter, group, value_column)
            for value_column, group in metrics
        ]
        for parameter in param_labels
    }

    raw_values = np.concatenate([values[parameter][metric_idx]
                                 for parameter in param_labels for metric_idx in (0, 1)])
    ess_values = np.concatenate([values[parameter][metric_idx]
                                 for parameter in param_labels for metric_idx in (2, 3)])

    def padded_limits(all_values):
        if not np.all(np.isfinite(all_values)):
            raise ValueError('Self-distance values must all be finite')
        data_min = float(np.min(all_values))
        data_max = float(np.max(all_values))
        value_range = data_max - data_min
        scale = max(abs(data_min), abs(data_max), 1e-12)
        padding = max(0.08 * value_range, 0.05 * scale, 1e-12)
        lower = data_min - padding
        if data_min >= 0:
            lower = max(0.0, lower)
        return lower, data_max + padding

    raw_limits = padded_limits(raw_values)
    ess_limits = padded_limits(ess_values)
    if len(settings) == 1:
        jitter = np.zeros(1)
    else:
        jitter = np.linspace(-options['jitter_span'] / 2, options['jitter_span'] / 2, len(settings))

    n_parameters = len(param_labels)
    display_parameter_labels = [figure_param_labels.get(parameter, parameter) for parameter in param_labels]
    display_group_labels = [group_display_labels.get(group, group) for group in group_order]
    footer_labels = (r'Raw-count $\mathrm{E}[D^{\mathrm{self}}]$',
                     r'ESS-adjusted $\mathrm{E}[D^{\mathrm{self}}]$')

    def maximum_text_extent(labels, fontsize, fontweight='normal'):
        measurement_figure = plt.figure(figsize=(1, 1))
        text_artists = [measurement_figure.text(0, 0, str(label), fontsize=fontsize, fontweight=fontweight)
                        for label in labels]
        measurement_figure.canvas.draw()
        renderer = measurement_figure.canvas.get_renderer()
        widths = [text_artist.get_window_extent(renderer).width for text_artist in text_artists]
        heights = [text_artist.get_window_extent(renderer).height for text_artist in text_artists]
        extent = max(widths) / measurement_figure.dpi, max(heights) / measurement_figure.dpi
        plt.close(measurement_figure)
        return extent

    parameter_label_width, parameter_label_height = maximum_text_extent(
        display_parameter_labels, options['label_fontsize'])
    group_label_width, group_label_height = maximum_text_extent(
        display_group_labels, options['group_fontsize'])
    footer_label_width, footer_label_height = maximum_text_extent(
        footer_labels, options['footer_fontsize'], fontweight='bold')

    effective_column_width = max(options['column_width'], group_label_width + 0.16,
                                 footer_label_width / 2 + 0.12)
    effective_row_height = max(options['row_height'], parameter_label_height + 0.35)
    effective_footer_height = max(options['footer_height'], group_label_height + footer_label_height + 0.42)
    plot_width = 4 * effective_column_width + options['inter_group_width']
    figure_width = max(options['minimum_width'], plot_width + parameter_label_width + 0.35)
    figure_height = n_parameters * effective_row_height + effective_footer_height
    figure = plt.figure(figsize=(figure_width, figure_height), layout='constrained', facecolor='white')
    grid = figure.add_gridspec(
        n_parameters + 1,
        5,
        width_ratios=[effective_column_width, effective_column_width, options['inter_group_width'],
                      effective_column_width, effective_column_width],
        height_ratios=[effective_row_height] * n_parameters + [effective_footer_height],
    )
    layout_engine = figure.get_layout_engine()
    if layout_engine is not None:
        layout_engine.set(h_pad=0.02, w_pad=0.02, hspace=0.03, wspace=0.05)

    axis_columns = (0, 1, 3, 4)
    axes = np.empty((n_parameters, 4), dtype=object)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for row_idx, parameter in enumerate(param_labels):
        color = colors[row_idx % len(colors)]
        for column_idx, grid_column in enumerate(axis_columns):
            axis = figure.add_subplot(grid[row_idx, grid_column])
            axes[row_idx, column_idx] = axis
            data = values[parameter][column_idx]
            data_min = np.min(data)
            data_max = np.max(data)
            data_median = np.median(data)
            range_height = options['range_height']

            axis.add_patch(Rectangle((data_min, 1 - range_height / 2), data_max - data_min, range_height,
                                     facecolor=color, edgecolor=color, alpha=0.32, linewidth=1.5, zorder=1))
            axis.plot([data_median, data_median], [1 - range_height / 2, 1 + range_height / 2],
                      color='black', linewidth=1.9, zorder=2)
            axis.scatter(data, 1 + jitter, s=options['point_size'], facecolor='white', edgecolor=color,
                         linewidth=1.0, zorder=3)
            if selected_setting_idx is not None:
                axis.scatter(data[selected_setting_idx], 1 + jitter[selected_setting_idx],
                             s=options['selected_point_size'], facecolor=color,
                             edgecolor=options['selected_point_edgecolor'],
                             linewidth=options['selected_point_linewidth'], zorder=4)

            axis.set_xlim(*(raw_limits if column_idx < 2 else ess_limits))
            axis.set_ylim(0.66, 1.34)
            axis.set_yticks([])
            axis.set_axisbelow(True)
            axis.grid(axis='x', color='#D9D9D9', linewidth=0.7, alpha=0.75)
            axis.spines['top'].set_visible(False)
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_visible(False)
            axis.spines['bottom'].set_color('#8A8A8A')
            axis.xaxis.set_major_locator(MaxNLocator(nbins=options['maximum_ticks'], min_n_ticks=3))
            format_specification = '.%dg' % options['tick_significant_digits']
            formatter = FuncFormatter(lambda value, position, spec=format_specification: format(value, spec))
            axis.xaxis.set_major_formatter(formatter)
            axis.tick_params(axis='x', labelsize=options['tick_fontsize'], width=0.8, length=3,
                             labelbottom=row_idx == n_parameters - 1)
        axes[row_idx, 0].set_ylabel(
            display_parameter_labels[row_idx],
            rotation=0,
            labelpad=options['label_padding'],
            fontsize=options['label_fontsize'],
            ha='right',
            va='center',
        )

    def add_footer(axis, label):
        axis.set_axis_off()
        axis.set_xlim(0, 2)
        axis.set_ylim(0, 1)
        axis.text(0.5, 0.78, display_group_labels[0], ha='center', va='center',
                  fontsize=options['group_fontsize'])
        axis.text(1.5, 0.78, display_group_labels[1], ha='center', va='center',
                  fontsize=options['group_fontsize'])
        axis.plot([0.02, 1.98], [0.43, 0.43], color='black', linewidth=1.1, clip_on=False)
        axis.plot([0.02, 0.02], [0.43, 0.35], color='black', linewidth=1.1, clip_on=False)
        axis.plot([1.98, 1.98], [0.43, 0.35], color='black', linewidth=1.1, clip_on=False)
        axis.text(1.0, 0.12, label, ha='center', va='center',
                  fontsize=options['footer_fontsize'], fontweight='bold')

    raw_footer = figure.add_subplot(grid[-1, 0:2])
    ess_footer = figure.add_subplot(grid[-1, 3:5])
    add_footer(raw_footer, footer_labels[0])
    add_footer(ess_footer, footer_labels[1])

    figure.savefig(output_filepath, dpi=options['dpi'], bbox_inches='tight',
                   pad_inches=options['save_padding'], facecolor='white')
    if show_plot:
        plt.show()
    plt.close(figure)


def save_self_distance_sensitivity_outputs(dirpath, csv_files, param_labels=None, group_order=None,
                                           group_display_labels=None, table_param_labels=None,
                                           figure_param_labels=None, selected_settings=None,
                                           dataset_filename='Dataset_S2.xlsx',
                                           figure_filename='Figure_S7_self_distance_sensitivity.png',
                                           show_plot=False, figure_options=None):
    """Create Dataset S2 and Figure S7 from self-distance sensitivity CSV files.

    When ``param_labels`` or ``group_order`` is None, its order is inferred from first appearance in the first CSV
    file. Supplying either argument remains useful when a different display order is desired. Optional
    ``figure_options`` override the responsive figure-layout and styling defaults. When ``selected_settings`` is
    supplied, its spreadsheet values are highlighted and its figure points are shown as filled markers. A value of
    None leaves all sensitivity settings visually equivalent.
    """

    combined, settings, param_labels, group_order = _load_self_distance_sensitivity_files(
        dirpath=dirpath, csv_files=csv_files, param_labels=param_labels, group_order=group_order)
    group_display_labels = {} if group_display_labels is None else dict(group_display_labels)
    table_param_labels = {} if table_param_labels is None else dict(table_param_labels)
    figure_param_labels = {} if figure_param_labels is None else dict(figure_param_labels)

    os.makedirs(dirpath, exist_ok=True)
    dataset_filepath = os.path.join(dirpath, dataset_filename)
    figure_filepath = os.path.join(dirpath, figure_filename)

    _write_self_distance_dataset(
        output_filepath=dataset_filepath, combined=combined, settings=settings, param_labels=param_labels,
        group_order=group_order, group_display_labels=group_display_labels,
        table_param_labels=table_param_labels, selected_settings=selected_settings)
    _plot_self_distance_sensitivity(
        output_filepath=figure_filepath, combined=combined, settings=settings, param_labels=param_labels,
        group_order=group_order, group_display_labels=group_display_labels,
        figure_param_labels=figure_param_labels, show_plot=show_plot,
        selected_settings=selected_settings, figure_options=figure_options)

    print('Saving...')
    print('    %s' % dataset_filename)
    print('    %s' % figure_filename)
    print('...to %s' % os.path.abspath(dirpath))

    return dataset_filepath, figure_filepath


def _get_pydream_sample_files(dirpath, nchains=None):
    "Return complete PyDREAM sampled-parameter files grouped by chain."
    if nchains is None:
        sample_pattern = 'dreamzs*sampled_params_chain_*_*.npy'
    else:
        sample_pattern = (
            'dreamzs_%dchain_sampled_params_chain_*_*.npy' % nchains
        )

    samples_files = glob.glob(os.path.join(dirpath, sample_pattern))
    if not samples_files:
        raise ValueError(
            "No PyDREAM sampled parameter files found in '%s'" % dirpath
        )

    files_by_chain = {}
    for filename in samples_files:
        match = re.search(r'chain_(\d+)_(\d+)\.npy$', filename)
        if match is None:
            continue

        chain = int(match.group(1))
        iteration = int(match.group(2))
        chain_files = files_by_chain.setdefault(chain, {})

        if iteration in chain_files:
            raise ValueError(
                "Multiple sampled parameter files found for chain %d "
                "at iteration %d" % (chain, iteration)
            )

        chain_files[iteration] = filename

    chains = sorted(files_by_chain)
    if not chains:
        raise ValueError(
            "No valid PyDREAM sampled parameter files found in '%s'" % dirpath
        )

    if nchains is not None:
        expected_chains = list(range(nchains))
        if chains != expected_chains:
            raise ValueError(
                "Expected sampled parameter files for chains %s, "
                "found chains %s" % (expected_chains, chains)
            )

    common_iterations = set(files_by_chain[chains[0]])
    for chain in chains[1:]:
        common_iterations.intersection_update(files_by_chain[chain])

    if not common_iterations:
        raise ValueError(
            "No common sampled-parameter iteration checkpoint found "
            "across chains"
        )

    latest_common_iteration = max(common_iterations)

    # Require every chain to contain the same complete checkpoint sequence
    # through the latest checkpoint common to all chains. Later files from
    # only some chains are ignored because they can be left by an interrupted
    # run, but a missing intermediate checkpoint is treated as an error.
    iterations = sorted(
        iteration for iteration in files_by_chain[chains[0]]
        if iteration <= latest_common_iteration
    )

    for chain in chains[1:]:
        chain_iterations = sorted(
            iteration for iteration in files_by_chain[chain]
            if iteration <= latest_common_iteration
        )
        if chain_iterations != iterations:
            raise ValueError(
                "Sampled parameter checkpoints are not consistent across "
                "chains through iteration %d" % latest_common_iteration
            )

    complete_files = {
        chain: [files_by_chain[chain][iteration] for iteration in iterations]
        for chain in chains
    }

    return chains, iterations, complete_files


def load_pydream_restart_state(dirpath, nchains=None):
    "Load only the iteration count and final point needed to restart PyDREAM."
    chains, iterations, files_by_chain = _get_pydream_sample_files(
        dirpath, nchains=nchains
    )
    latest_iteration = iterations[-1]

    start = []
    nparams = None

    for chain in chains:
        samples = np.load(files_by_chain[chain][-1], mmap_mode='r')

        if samples.ndim != 2 or samples.shape[0] == 0:
            raise ValueError(
                "Invalid sampled parameter array in '%s'" %
                files_by_chain[chain][-1]
            )

        if nparams is None:
            nparams = samples.shape[1]
        elif samples.shape[1] != nparams:
            raise ValueError(
                "Sampled parameter dimensions are not consistent across chains"
            )

        start.append(np.array(samples[-1, :], copy=True))
        del samples

    return latest_iteration, start


def _combine_population_stats(n_a, mean_a, m2_a, n_b, mean_b, var_b):
    "Combine population mean/variance statistics from two sample blocks."
    if n_a == 0:
        return n_b, mean_b.copy(), var_b * n_b

    n_total = n_a + n_b
    delta = mean_b - mean_a
    mean_total = mean_a + delta * (float(n_b) / n_total)
    m2_total = (
        m2_a
        + (var_b * n_b)
        + (delta ** 2) * n_a * n_b / n_total
    )

    return n_total, mean_total, m2_total


def Gelman_Rubin_from_files(dirpath, fburnin=0.5, nchains=None):
    "Calculate PyDREAM Gelman-Rubin statistics from saved sample files."
    if not 0 <= fburnin < 1:
        raise ValueError("fburnin must satisfy 0 <= fburnin < 1")

    chains, _, files_by_chain = _get_pydream_sample_files(
        dirpath, nchains=nchains
    )

    nsamples_by_chain = []
    nparams = None

    for chain in chains:
        nsamples = 0

        for filename in files_by_chain[chain]:
            samples = np.load(filename)

            if samples.ndim != 2:
                raise ValueError(
                    "Invalid sampled parameter array in '%s'" % filename
                )

            if nparams is None:
                nparams = samples.shape[1]
            elif samples.shape[1] != nparams:
                raise ValueError(
                    "Sampled parameter dimensions are not consistent "
                    "across files"
                )

            nsamples += samples.shape[0]
            del samples

        nsamples_by_chain.append(nsamples)

    if len(set(nsamples_by_chain)) != 1:
        raise ValueError(
            "Number of sampled parameter values is not equal across chains"
        )

    nsamples = nsamples_by_chain[0]
    nburnin = int(math.floor(nsamples * fburnin))
    nkeep = nsamples - nburnin

    if nkeep <= 0:
        raise ValueError("No samples remain after burn-in")

    chain_means = []
    chain_vars = []

    for chain in chains:
        n_seen = 0
        n_kept = 0
        mean = np.zeros(nparams, dtype=float)
        m2 = np.zeros(nparams, dtype=float)

        for filename in files_by_chain[chain]:
            samples = np.load(filename)
            file_start = n_seen
            file_end = n_seen + samples.shape[0]

            if file_end > nburnin:
                local_start = max(0, nburnin - file_start)
                block = samples[local_start:, :]
                block_n = block.shape[0]

                if block_n:
                    block_mean = np.mean(
                        block, axis=0, dtype=np.float64
                    )
                    block_var = np.var(
                        block, axis=0, dtype=np.float64
                    )

                    n_kept, mean, m2 = _combine_population_stats(
                        n_kept,
                        mean,
                        m2,
                        block_n,
                        block_mean,
                        block_var,
                    )

                    del block_mean, block_var, block

            n_seen = file_end
            del samples

        if n_kept != nkeep:
            raise ValueError(
                "Unexpected number of post-burn-in samples for chain %d: "
                "expected %d, found %d" % (chain, nkeep, n_kept)
            )

        chain_means.append(mean)
        chain_vars.append(m2 / n_kept)

    W = np.mean(chain_vars, axis=0)
    B = np.var(chain_means, axis=0)
    var_est = W * (1 - (1. / nkeep)) + B

    return np.sqrt(np.divide(var_est, W))


def iter_simulation_outputs_resilient(sample_args, run_fn, run_with_timeout_fn, timeout_seconds=None, n_procs=1):

    batch_size = n_procs if n_procs > 1 else 1
    n_total = len(sample_args)
    i = 0

    while i < n_total:
        batch_end = min(i + batch_size, n_total)
        batch = sample_args[i:batch_end]

        if n_procs == 1:
            for local_j, args in enumerate(batch):
                global_idx = i + local_j
                if timeout_seconds is None:
                    yield global_idx, run_fn(args), False
                else:
                    yield global_idx, *run_with_timeout_fn(args, timeout_seconds)
            i = batch_end
            continue

        import multiprocessing
        pool = multiprocessing.Pool(processes=n_procs)
        async_results = []

        try:
            for args in batch:
                async_results.append(pool.apply_async(run_fn, (args,)))

            batch_outputs = [None] * len(batch)
            completed = [False] * len(batch)
            timeout_triggered = False

            for local_j, async_result in enumerate(async_results):
                try:
                    if timeout_seconds is None:
                        batch_outputs[local_j] = async_result.get()
                    else:
                        batch_outputs[local_j] = async_result.get(timeout=timeout_seconds)
                    completed[local_j] = True
                except multiprocessing.TimeoutError:
                    timeout_triggered = True
                    break

            if not timeout_triggered:
                pool.close()
                pool.join()
                for local_j, output in enumerate(batch_outputs):
                    yield i + local_j, output, False
            else:
                pool.terminate()
                pool.join()

                for local_j, output in enumerate(batch_outputs):
                    if completed[local_j]:
                        yield i + local_j, output, False

                for local_j, args in enumerate(batch):
                    if not completed[local_j]:
                        yield i + local_j, *run_with_timeout_fn(args, timeout_seconds)

        finally:
            try: pool.close()
            except: pass
            try: pool.join()
            except: pass

        i = batch_end


def save_posterior_summary(dirpath, samples, groups, param_labels,
                           group_labels, max_iter=None, filename=None):
    """Save posterior medians and central 90% credible intervals.

    Parameters
    ----------
    dirpath : str or path-like
        Directory in which to save the CSV file.
    samples : array-like
        Retained parameter set ensemble with shape (n_samples, n_parameters). Values must be in log10 space.
    groups : sequence of sequences
        Parameter-column indices for each parameter group.
    param_labels : sequence of sequences
        Parameter names corresponding to the indices in each group.
    group_labels : sequence of str
        Labels describing the parameter groups, such as Bone-Adapted, Parental, and Shared.
    max_iter : int or None, optional
        Maximum PyDREAM iteration used. This is included in the default filename.
    filename : str or None, optional
        Output filename. A default name is generated when None.

    Returns
    -------
    str
        Path to the saved CSV file.
    """
    samples = np.asarray(samples, dtype=float)

    if samples.ndim != 2:
        raise ValueError("'samples' must have shape (n_samples, n_parameters); found %s" % samples.shape)

    if samples.shape[0] == 0:
        raise ValueError("'samples' contains no retained parameter sets")

    if not np.all(np.isfinite(samples)):
        raise ValueError("'samples' contains nonfinite values")

    if not len(groups) == len(param_labels) == len(group_labels):
        raise ValueError("'groups', 'param_labels', and 'group_labels' must have the same length")

    parameter_indices = []

    for group, labels in zip(groups, param_labels):
        if len(group) != len(labels):
            raise ValueError("Each parameter group must have the same length as its corresponding label list")

        for parameter_idx in group:
            if not isinstance(parameter_idx, (int, np.integer)):
                raise TypeError("Parameter indices must be integers; found %r" % parameter_idx)

            if not 0 <= parameter_idx < samples.shape[1]:
                raise IndexError("Parameter index %d is outside the available range 0-%d" %
                                 (parameter_idx, samples.shape[1] - 1))

            parameter_indices.append(int(parameter_idx))

    expected_indices = list(range(samples.shape[1]))

    if sorted(parameter_indices) != expected_indices:
        raise ValueError("'groups' must include every sampled-parameter column exactly once")

    if filename is None:
        iteration_label = 'all' if max_iter is None else str(max_iter)
        filename = 'PyDREAM_posterior_summary_%s.csv' % iteration_label

    output_filepath = os.path.join(dirpath, filename)

    with open(output_filepath, 'w', encoding='utf-8', newline='') as output_file:
        csv_writer = csv.writer(output_file)
        csv_writer.writerow(['parameter', 'group', 'posterior_median', 'credible_interval_5th_percentile',
                             'credible_interval_95th_percentile'])
        for group_label, group, labels in zip(group_labels, groups, param_labels):
            for parameter_idx, parameter_label in zip(group, labels):
                quantiles_log10 = np.percentile(samples[:, parameter_idx], [5.0, 50.0, 95.0])
                lower, median, upper = np.power(10.0, quantiles_log10)
                csv_writer.writerow([parameter_label, group_label, '%.12g' % median, '%.12g' % lower, '%.12g' % upper])

    print('Saving...')
    print('    %s' % filename)
    print('...to %s' % os.path.abspath(dirpath))

    return output_filepath
