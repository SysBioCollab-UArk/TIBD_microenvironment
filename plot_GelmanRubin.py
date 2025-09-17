import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import math
from pydream_util import get_fig_ncols


def plot_GR_metrics(directory, threshold=1.2, par_names=None, only_unconverged=False, **kwargs):

    # process kwargs
    cmap = kwargs.pop('cmap', 'turbo')  # or 'hsv', 'viridis', etc.

    files = glob.glob(os.path.join(directory, '*Gelman*'))
    iterations = sorted([int(re.search("\d\d+", os.path.split(file)[1]).group()) for file in files])
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
    else:
        max_GR = max([GR[-1] for GR, i in GR_and_label])
        ymax = math.ceil(max_GR * 10) / 10 + 0.1

        # get number of columns and rows for the multi-plot figure
        ncols = get_fig_ncols(len(GR_and_label))
        nrows = math.ceil(len(GR_and_label) / ncols)

        colormap = plt.get_cmap(cmap)
        colors = [colormap(i / (len(GR_and_label) - 1)) for i in range(len(GR_and_label))]

        fig = plt.figure(constrained_layout=True, figsize=(6.4 / 3 * ncols, 4.8 / 3.5 * nrows))
        ax_ref = None
        for i, ((GR_array, label), color) in enumerate(zip(GR_and_label, colors)):
            ax = fig.add_subplot(nrows, ncols, i+1, sharex=ax_ref, sharey=ax_ref)
            if i == 0:
                ax_ref = ax
            ax.plot(iterations, GR_array, '-', color=color)  #, label=label)
            ax.set_title(label, color=color, fontweight='bold')
            ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
            ax.label_outer()
            ax.axhline(threshold, ls='--', color='black', lw=1.5)
            ax.set_ylim(top=ymax)
            # leg = ax.legend(loc='best', handlelength=0, handletextpad=0, facecolor='white', edgecolor='white', framealpha=1,
            #                 frameon=True, prop={'weight': 'bold'})
            # for text in leg.get_texts():
            #     text.set_color(color)
        fig.supxlabel('Iteration')
        fig.supylabel('Gelman-Rubin metric')

        plt.show()


if __name__ == '__main__':

    # dirpath = os.path.join('MODELS', 'SAVE', 'Leonard')
    # dirname = 'Bennett2024_Johnson2011_stderr_2abs'

    dirpath = os.path.join('..', 'AorticCalcification', 'SAVE')
    dirname = 'Messika_stderr10'

    directory = os.path.join(dirpath, dirname)

    plot_GR_metrics(directory, only_unconverged=True)
