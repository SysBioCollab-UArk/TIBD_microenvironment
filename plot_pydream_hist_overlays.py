from param_calibration import *
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.colors import to_rgb
import os
import importlib
import math
from pydream_util import get_fig_ncols


def calc_self_distance(kde, n_samples, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    E_Dself = 0.5 * np.sqrt(2 / n_samples / np.pi) * np.sum(np.sqrt(kde(x_vals)[1:] * dx))
    return E_Dself


def calc_hist_distance(kde, kde_ref, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    D = 0.5 * np.sum(np.abs(kde(x_vals)[1:] - kde_ref(x_vals)[1:]) * dx)
    return D


def get_word_width(fig, word, fontsize, fontweight):
    # Font properties
    fontprops = FontProperties(size=fontsize, weight=fontweight)
    # Measure full word width (includes kerning)
    tp_word = TextPath((0, 0), word, prop=fontprops)
    word_width = tp_word.get_extents().width / 72 / fig.get_size_inches()[0]
    return word_width


def write_multicolor_word(fig, x, y, word, colors, fontsize=10, fontweight='normal', additional_text=None):

    if len(colors) != len(word):
        raise Exception("Length of 'colors' (%d) does not match length of 'word' (%d)" % (len(colors), len(word)))

    word_width = get_word_width(fig, word, fontsize, fontweight)

    # Measure individual letter widths
    letter_widths = []
    for letter in word:
        letter_widths.append(get_word_width(fig, letter, fontsize, fontweight))

    # Estimate average inter-letter spacing from kerning
    letter_spacing = 0 if len(word) <= 1 else (word_width - sum(letter_widths)) / (len(word) - 1)

    # Draw each character
    for i, (letter, color) in enumerate(zip(word, colors)):
        fig.text(x, y, letter, color=color, ha='left', fontsize=fontsize, fontweight=fontweight)
        # Move x-position forward by width of character + additional spacing
        x += letter_widths[i] + letter_spacing

    # Add additional text, if exists
    if isinstance(additional_text, str):
        fig.text(x, y, additional_text, ha='left', fontsize=fontsize, fontweight=fontweight)
    elif isinstance(additional_text, dict):
        text = additional_text.pop('text')
        fs = additional_text.pop('fontsize', fontsize)
        fw = additional_text.pop('fontweight', fontweight)
        fig.text(x, y, text, ha='left', fontsize=fs, fontweight=fw, **additional_text)


def plot_hist_overlays(two_samples, param_labels, hist_labels, E_Dself=None, show_plots=False, save_plots=True,
                       **kwargs):

    # error check
    if not (len(two_samples) == len(hist_labels) == 2):
        raise Exception("The number of sets of parameter samples and histogram labels must equal 2")

    if not (two_samples[0].shape[1] == two_samples[1].shape[1] == len(param_labels)):
        raise Exception("The number of parameter labels must equal the number of columns in 'two_samples'")

    # process kwargs
    fontsizes = kwargs.get('fontsizes', {})
    labels_fs = fontsizes.get('labels', None)
    ticks_fs = fontsizes.get('ticks', None)
    title_fs = fontsizes.get('title', None)
    legend_fs = fontsizes.get('legend', None)
    bw_adjust = kwargs.get('bw_adjust', (1, 1))  # histogram smoothing parameter
    sharex = kwargs.get('sharex', False)
    table_props = kwargs.get('table_props', {})
    table_fs = table_props.get('fontsize', None)
    table_ncols = table_props.get('ncols', 2)
    table_scale = table_props.get('scale', None)  # stretch/squeeze embedded table
    table_nudge = table_props.get('nudge', None)  # move table horizontally and/or vertically

    ### Plot histogram distances with self distances ###
    n_params = len(param_labels)
    ncols = get_fig_ncols(n_params)
    nrows = math.ceil(n_params / ncols)
    labelsize = 10 * max(1, (2 / 5 * np.ceil(nrows / 2)))
    fontsize = 10 * max(1, (3 / 5 * np.ceil(nrows / 2)))
    colors = sns.color_palette(n_colors=n_params)
    fig_width = 0.65 * ncols * 6.4
    fig_overlay = plt.figure(constrained_layout=True, figsize=(fig_width, 0.6 * nrows * 4.8))
    axes = []
    reference_ax = None
    D = []  # histogram distances
    if E_Dself is None:
        E_Dself = [[None] * n_params for _ in range(2)]
    print('Plotting %d histogram overlays' % n_params)
    for n in range(n_params):
        print(n, end=' ')
        # share x-axis with first subplot
        ax = fig_overlay.add_subplot(nrows, ncols, n + 1, sharex=reference_ax)
        if n == 0 and sharex:
            reference_ax = ax
        axes.append(ax)
        for i, (samples, color) in enumerate(zip(two_samples, ['k', colors[n % n_params]])):
            # reference histograms are black, sample histograms are multicolor
            sns.kdeplot(samples[:, n], color=color, fill=True, common_norm=False, ax=ax, bw_adjust=bw_adjust[i])
            x_vals = sorted(ax.collections[i].get_paths()[0].vertices[:, 0])  # get x-axis vals from kdeplot
            # get kernel density estimate (KDE) for calculating histogram distance and self distance
            kde = stats.gaussian_kde(samples[:, n])
            kde.set_bandwidth(kde.factor * bw_adjust[i])
            if i == 0:  # reference histogram
                # calculate self distance (expected value)
                kde_ref = kde
                x_min_ref = x_vals[0]
                x_max_ref = x_vals[-1]
                if E_Dself[0][n] is None:
                    E_Dself[0][n] = calc_self_distance(kde_ref, len(samples[:, n]), x_min_ref, x_max_ref, 1000)
            else:
                # calculate histogram distance relative to the reference (use 2x the points, just to be safe)
                D.append(calc_hist_distance(kde, kde_ref, min(x_vals[0], x_min_ref), max(x_vals[-1], x_max_ref),
                                            2000))
        empty_handle = Line2D([], [], linestyle="none")
        legend = ax.legend([empty_handle, empty_handle],
                           [r'$D$: %.3f' % D[-1],
                            r'$\left\langle D^\mathrm{self} \right\rangle$: %.3f' % E_Dself[0][n]],
                           fontsize=0.9 * labelsize if legend_fs is None else legend_fs,
                           loc='best', handlelength=0, handletextpad=0, labelspacing=0.3, borderpad=0)
        legend.set_frame_on(False)
        ax.set_yticklabels([])
        ax.set_ylabel(None)
        ax.tick_params(axis='x', labelsize=labelsize if ticks_fs is None else ticks_fs)
        # ax.label_outer()
        ax.set_title(param_labels[n], fontsize=labelsize if title_fs is None else title_fs)
    print()
    fig_overlay.supxlabel(r'log$_{10}$ value' + '\n\n', fontsize=fontsize if labels_fs is None else labels_fs)
    fig_overlay.supylabel('Density', fontsize=fontsize if labels_fs is None else labels_fs)

    # Create a common figure legend for the overlay figure
    fontsize = fig_overlay.get_axes()[0].xaxis.label.get_fontsize()
    fig_legend_fs = 1.1 * (fontsize if labels_fs is None else labels_fs)  # make legend text a bit bigger
    space_height = 2 * fig_legend_fs / 72 / fig_overlay.get_size_inches()[0]
    # reference histograms are black
    additional_text = {'text': ": %s " % hist_labels[0], 'fontweight': 'normal'}
    text_width = get_word_width(fig_overlay, 'Black', fig_legend_fs, 'bold')
    text_width += get_word_width(fig_overlay, additional_text['text'], fig_legend_fs, additional_text['fontweight'])
    write_multicolor_word(fig_overlay, 0.5 - text_width - 0.01, 0.45 * space_height, "Black",
                          ['k'] * len("Black"), fontsize=fig_legend_fs, fontweight='bold',
                          additional_text=additional_text)
    # sample histograms are multicolor
    additional_text = {'text': ": %s" % hist_labels[1], 'fontweight': 'normal'}
    write_multicolor_word(fig_overlay, 0.5 + 0.01, 0.45 * space_height, "_Multicolor",
                          [to_rgb('white')] + sns.color_palette(n_colors=len("Multicolor")), fontsize=fig_legend_fs,
                          fontweight='bold', additional_text=additional_text)

    ### Make barplot figure rank-ordered by histogram distance ###
    # Make a table with parameter names rank-ordered by histogram distance
    sorted_idxs = np.argsort(D)[::-1]  # sort from largest to smallest
    sorted_labels = [param_labels[i % n_params] for i in sorted_idxs]
    table_data = []
    for i, col in enumerate(range(table_ncols)):
        start = col * n_params // table_ncols
        end = (col + 1) * n_params // table_ncols
        '''print('start: %d, end: %d, len(sorted_labels): %d' % (start, end, len(sorted_labels)))'''
        for row in range(end - start):
            '''print('row: %d, len(table_data): %d' % (row, len(table_data)))'''
            if row == len(table_data):
                table_data.append([])
            if len(table_data[row]) < col:
                table_data[row] += [""]
            table_data[row] += ["%d. %s" % (start + row, sorted_labels[start + row])]
    '''for row in table_data:
        print(row)'''

    # Create the barplot
    fig_barplot = plt.figure(constrained_layout=True, figsize=(fig_width, fig_width / 6.4 * 4.8 * 0.75))
    plt.bar(np.arange(len(sorted_idxs)), [D[i] for i in sorted_idxs])
    for i, e in enumerate([E_Dself[0][j] for j in sorted_idxs]):
        plt.plot([i - 0.4, i + 0.4], [e, e], color='r', lw=2)
    plt.xlabel('Index', fontsize=labels_fs)
    plt.ylabel('Histogram Distance', fontsize=labels_fs)
    plt.tick_params(axis='both', which='major', labelsize=ticks_fs)

    # Add table to the barplot
    def add_table(loc=None, bbox=None):
        this_table = plt.table(
            cellText=table_data,
            colLabels=None,
            loc=loc,
            bbox=bbox,
            cellLoc='left'
        )
        this_table.set_fontsize(table_fs)  # Set font size
        # Get column widths
        col_widths = [0 for _ in range(table_ncols)]
        for row in range(len(table_data)):
            for col in range(len(table_data[row])):
                word_width = get_word_width(fig_barplot, table_data[row][col], fontsize=table_fs,
                                            fontweight='normal') / fig_barplot.get_axes()[0].get_position().width
                if word_width > col_widths[col]:
                    col_widths[col] = word_width
        # Remove all borders and set column widths
        for (_, col), cell in this_table.get_celld().items():
            cell.set_linewidth(0)
            cell.set_width(col_widths[col])

        return this_table

    table = add_table(loc='upper right')

    # If the user requests the table to be scaled or shifted, get the table location and dimensions, delete the current
    # table, and create a new one with `add_table()` using `bbox` instead of `loc`
    if any([arg is not None] for arg in (table_scale, table_nudge)):
        fig_barplot.canvas.draw()  # ensure renderer is ready
        # Get the bbox to get the current position and size of the table
        bbox = table.get_window_extent(renderer=fig_barplot.canvas.get_renderer())
        bbox = bbox.transformed(fig_barplot.transFigure.inverted())
        table.remove()  # delete the old table
        # Make a new table
        scale_x, scale_y = (1, 1) if table_scale is None else (table_scale[0], table_scale[1])
        dx, dy = (0, 0) if table_nudge is None else (table_nudge[0], table_nudge[1])
        add_table(bbox=(bbox.x0 + dx - bbox.width * (scale_x - 1), bbox.y0 + dy - bbox.height * (scale_y - 1),
                        bbox.width * scale_x, bbox.height * scale_y))

    if save_plots is not False:
        outdir = '.' if save_plots is True else save_plots
        print('Saving...')
        for fig, fig_type in zip([fig_overlay, fig_barplot], ['hist_overlay', 'hist_distances']):
            outfile = 'fig_PyDREAM_%s_%s' % (fig_type, str.join('_', hist_labels))
            print('    %s' % outfile)
            fig.savefig(os.path.join(outdir, outfile))
        print('...to %s' % os.path.abspath(outdir))

    if show_plots:
        plt.show()

    # return figures and E_Dself
    return fig_overlay, fig_barplot, E_Dself


def plot_hist_overlays_from_dirs(dirpath, directories, run_pydream_filename, show_plots=False,
                                             save_plots=True, **kwargs):

    if save_plots is not False:
        outdir = '.' if save_plots is True else save_plots

    directories = np.array(directories)

    # process kwargs
    bw_adjust = kwargs.get('bw_adjust', [1] * len(directories))  # default smoothing parameter

    samples_ALL = []
    for i, directory in enumerate(directories):

        # Set the 'path' variable to the directory where the SIM_DATA.csv, run_<...>_pydream.py, and expt data files are
        path = os.path.join(dirpath, directory)  # os.getcwd()

        # import everything from run_<...>_pydream.py file that's in the path
        run_pydream_file = os.path.join(path, run_pydream_filename)
        import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
        module = importlib.import_module(import_string)  # Import the module

        # get the path to the experimental data file referenced in the run_<...>_pydream.py file that's in the path
        if not os.path.isabs(module.exp_data_file):
            exp_data_file = os.path.normpath(os.path.join(path, module.exp_data_file))
        else:
            exp_data_file = os.path.normpath(module.exp_data_file)

        logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))
        samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))

        calibrator = ParameterCalibration(module.model,
                                          exp_data_file,
                                          module.sim_protocols,
                                          priors=module.custom_priors,
                                          no_sample=module.no_sample,
                                          param_expts_map=module.param_expts_map)

        _, (samples, groups, group_labels), _ = \
            calibrator.create_figures(logps_files, samples_files, obs_labels=None, show_plots=True,
                                      plot_ll_args={'cutoff': 2, 'file_suffix': directory},
                                      plot_pd_args={'sharex': 'all', 'bw_adjust': bw_adjust[i]},
                                      which_plots=2)

        samples_ALL.append(samples)

    # Create figures with parameter histograms overlaid
    n_params_tot = len(calibrator.parameter_idxs)
    # store E_Dself so don't need to recalculate
    E_Dself = np.array([[None] * n_params_tot for _ in range(len(samples_ALL))])
    samples_ALL_idxs = np.arange(len(samples_ALL))

    # Loop over all (ref, sample) pairs
    sample_pairs = set(tuple(sorted((a, b))) for a in samples_ALL_idxs for b in samples_ALL_idxs if a != b)
    for sample_pair in sample_pairs:
        sample_pair = np.array(sample_pair)
        print("Comparing histograms for '%s' and '%s'" % (directories[sample_pair[0]], directories[sample_pair[1]]))

        # Loop over parameter groups
        for g, (group, param_labels) in enumerate(zip(groups, group_labels)):
            two_samples = [samples_ALL[idx][:, group] for idx in sample_pair]
            hist_labels = [directories[idx] for idx in sample_pair]
            E_Dself_g = E_Dself[sample_pair[:, None], np.array([group for _ in range(2)])]
            # Create figures by calling general overlay plotting function
            fig_ov, fig_bp, E_Dself_g = plot_hist_overlays(two_samples, param_labels, hist_labels, E_Dself_g,
                                                           show_plots=False, save_plots=False, **kwargs)
            E_Dself[sample_pair[:, None], np.array([group for _ in range(2)])] = E_Dself_g

            # save overlay figure
            if save_plots is not False:
                outfile = 'fig_PyDREAM_hist_overlay_%s' % str.join('_', directories[sample_pair])
                if len(groups) > 1:
                    outfile += '_group_%d' % g
                fig_ov.savefig(os.path.join(outdir, outfile))

            # save barplot figure
            if save_plots is not False:
                outfile = 'fig_PyDREAM_hist_distances_%s' % str.join('_', directories[sample_pair])
                if len(groups) > 1:
                    outfile += '_group_%d' % g
                fig_bp.savefig(os.path.join(outdir, outfile))

    if show_plots:
        plt.show()


if __name__ == '__main__':

    dirpath = 'SAVED'
    directories = ['Flav_FAD', 'Flav_Fumarate', 'Flav_Fumarate_FAD', 'Flav_Fumarate_FAD_Time']
    run_pydream_filename = 'run_complex_II_pydream.py'

    kwargs = {
        'fontsizes': {'labels': 22, 'ticks': 18, 'title': 18, 'legend': 14},
        'bw_adjust': [3.0, 2.0, 2.0, 2.0],  # histogram smoothing parameters (default = 1, > 1 = smoother)
        'sharex': False
    }

    plot_hist_overlays_from_dirs(dirpath, directories, run_pydream_filename, show_plots=True, **kwargs)
