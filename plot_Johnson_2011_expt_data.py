import numpy as np
import matplotlib.pyplot as plt
from pydream_util import plot_expt_data
import os

path_orig = os.path.join('DATA', 'TIBD_Johnson2011_orig.csv')
path_stderr = os.path.join('DATA', 'TIBD_Johnson2011_stderr5.csv')
path_tumor = os.path.join('DATA', 'TIBD_Johnson2011_tumor_RFU.csv')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
kwargs = {
    'fontsizes': {'labels': 20, 'ticks': 18, 'legend': 18, 'title': 20},
    'ms': 10,
    'capsize': 10,
    'colors': colors
}

label_dict = {
    'Bone_tot': 'Bone Density',
    'Tumor_tot': 'Tumor Burden',
    'relative % bone vol': 'relative %',
    'expt B': 'Original',
    'expt C': 'Original',
    'expt D': 'Original',
    'Bone_tot_B': 'Johnson (2011): No Drug-Tumor',
    'Bone_tot_C': 'Johnson (2011): ZA-Tumor',
    'Bone_tot_D': 'Johnson (2011): ZA-No Tumor',
    'Tumor_tot_B': 'Johnson (2011): No Drug-Tumor',
    'Tumor_tot_C': 'Johnson (2011): ZA-Tumor'
}

figures = plot_expt_data(path_orig, label_dict=label_dict, add_titles=True, show_plots=False, **kwargs)
for expt_id in ['B', 'C', 'D']:
    label_dict['expt %s' % expt_id] = 'Modified'
plot_expt_data(path_stderr, figures=figures, label_dict=label_dict, add_titles=True, mec='k', mfc='none', ecolor='k',
               **kwargs)

# Exponential fits to tumor fluorescence data
kwargs['colors'] = [colors[1], colors[2], colors[2]]
kwargs['markers'] = ['o', '^', '^']
kwargs['mfc'] = [None, 'none', None]
label_dict = {
    'Tumor_tot': 'Fluorescence',
    'expt B': 'No Drug-Tumor',
    'expt D': 'ZA-Tumor'
}
plot_expt_data(path_tumor, separate_plots=False, label_dict=label_dict, show_plots=False, **kwargs)
time = np.linspace(0, 28, 281)  # days
rfu = 131.25 * np.exp(0.917 * time / 7)
plt.plot(time, rfu, '--', lw=2, color=kwargs['colors'][0])
plt.annotate(r'RFU = 131.25 $\bf e^{%.3g \cdot \mathbf{Day}}$' % (0.917 / 7), (0.23, 0.45), xycoords='axes fraction',
             color=kwargs['colors'][0], fontsize=16, fontweight='bold')
rfu = 131.25 * np.exp(0.680 * time / 7)
plt.plot(time, rfu, '--', lw=2, color=kwargs['colors'][1])
plt.annotate(r'RFU = 131.25 $\bf e^{%.2g \cdot \mathbf{Day}}$' % (0.68 / 7), (0.43, 0.01), xycoords='axes fraction',
             color=kwargs['colors'][1], fontsize=16, fontweight='bold')
plt.xticks([0, 7, 14, 21, 28])

plt.axhline(750, color='k', linestyle='-.')
plt.annotate('Hypothesized Limit of Detection', (0.02, 0.16), xycoords='axes fraction')

handles, labels = plt.gca().get_legend_handles_labels()
plt.legend([handles[0], handles[2]], [labels[0], labels[2]], fontsize=kwargs['fontsizes']['legend'])

plt.show()
