import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import os


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
