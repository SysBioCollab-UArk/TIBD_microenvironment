import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import os
from collections.abc import Iterable
from math import isclose
import importlib


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


def is_list_like(obj):
    """Check if an object is iterable but NOT a string or bytes."""
    return isinstance(obj, Iterable) and not isinstance(obj, (str, bytes))


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


def plot_drc(basepath, directory, run_pydream_filename, expts_doses, label_dict=None, show_plot=True):

    dirs = [directory] if isinstance(directory, str) else directory
    for dir, expt_doses_list in zip(dirs, expts_doses):
        print('Directory:', dir)

        if isinstance(expt_doses_list, dict):
            expt_doses_list = [expt_doses_list]

        # set the 'path' variable to the directory where the SIM_DATA.csv, run_<...>_pydream.py, and expt data file are
        path = os.path.join(basepath, dir)

        # import everything from run_<...>_pydream.py file that's in the path
        run_pydream_file = os.path.join(path, run_pydream_filename)
        import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
        module = importlib.import_module(import_string)  # import the module

        # get the path to the experimental data file referenced in the run_<...>_pydream.py file that's in the path
        exp_data_file = os.path.normpath(module.exp_data_file) if os.path.isabs(module.exp_data_file) else \
            os.path.normpath(os.path.join(path, module.exp_data_file))
        expt_data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
        print('expt_data:', expt_data.dtype.names)

        sim_data_file = os.path.join(path, 'SIM_DATA.csv')
        sim_data = np.genfromtxt(sim_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
        print('sim_data:', sim_data.dtype.names)

        for expt_doses in expt_doses_list:
            plt.figure(constrained_layout=True, figsize=(6.4 * 0.7, 4.8 * 0.9))
            suffix = []
            min_conc = np.inf
            max_conc = 0
            for expt in [key for key in expt_doses.keys() if key not in ['xlabel', 'title']]:
                label = label_dict.get(expt, expt)
                conc = expt_doses[expt]
                if min(conc) < min_conc:
                    min_conc = min(conc)
                if max(conc) > max_conc:
                    max_conc = max(conc)
                # sim data
                yval_min = np.array([d['yval_min'] for d in sim_data if d['sim_id'] == expt])
                yval_max = np.array([d['yval_max'] for d in sim_data if d['sim_id'] == expt])
                legend_label = 'simulation' if label is None else '%s (sim)' % label
                p = plt.plot(conc, (yval_min + yval_max)/2, ls='--', label=legend_label)
                plt.fill_between(conc, yval_min, yval_max, alpha=0.25, color=p[0].get_color(), label='x')
                # expt data
                avg = np.array([d['average'] for d in expt_data if d['expt_id'] == expt])
                stderr = np.array([d['stderr'] for d in expt_data if d['expt_id'] == expt])
                legend_label = 'experiment' if label is None else '%s (expt)' % label
                plt.errorbar(conc, avg, yerr=stderr, fmt='o', ms=8, capsize=6, color=p[0].get_color(), label=legend_label)

                plt.title(expt_doses.get('title', None), fontsize=16, fontweight='bold')
                plt.xlabel(expt_doses.get('xlabel'), fontsize=16)
                yobs = np.unique([d['observable'] for d in expt_data if d['expt_id'] == expt])
                yunits = np.unique([d['amount_units'] for d in expt_data if d['expt_id'] == expt])
                if len(yobs) > 1 or len(yunits) > 1:
                    raise Exception("More than one observable or unit type for experiment '%s'" % expt)
                suffix.append(('%s_%s' % (expt, yobs[0])).split('_'))  # for the filename
                plt.ylabel('%s (%s)' % (label_dict.get(yobs[0], yobs[0]), yunits[0]), fontsize=16)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.ylim(bottom=0)
                plt.legend(loc='best')

            # adjust x-limits
            xmin = round_down_nice(min_conc)
            xmax = round_up_nice(max_conc)
            buffer = 0.05 * (xmax - xmin)  # 5% buffer
            plt.xlim(left=xmin - buffer, right=xmax + buffer)

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

            filename = 'fig_PyDREAM_DRC'
            for i in range(np.array(suffix).shape[1]):
                filename += '_' + '_'.join(np.unique([x[i] for x in suffix]))
            plt.savefig(os.path.join(path, filename))

    if show_plot:
        plt.show()
