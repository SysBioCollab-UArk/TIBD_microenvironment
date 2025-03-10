import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import ScalarFormatter

path = os.path.join('MODELS', 'SAVE', 'RUN_4')

# Experimental data
expt_data_file = os.path.join(path, 'TIBD_PopD_Data.csv')
expt_data_raw = np.genfromtxt(expt_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(expt_data_raw.dtype.names)
expt_ids = list(dict.fromkeys(expt_data_raw['expt_id']))
print(expt_ids)
# gather experimental data for all alt_expt_ids
expt_data = {}
n_alt_expt_ids = {}
for expt_id in expt_ids:
    alt_expt_ids = list(dict.fromkeys([d['alt_expt_id'] for d in expt_data_raw if d['expt_id'] == expt_id]))
    n_alt_expt_ids[expt_id] = len(alt_expt_ids)
    for alt_expt_id in alt_expt_ids:
        expt_data[alt_expt_id] = \
            np.array([d for d in expt_data_raw if d['expt_id'] == expt_id and d['alt_expt_id'] == alt_expt_id])

print('n_alt_expt_ids:', n_alt_expt_ids)
for key in expt_data.keys():
    print(key)
    # print('    ', expt_data[key].dtype.names)
    # for d in expt_data[key]:
    #     print('    ', d)
# quit()

print('----------')

# Observables
obs_names = np.unique(expt_data_raw['observable'])
print(obs_names)

print('----------')

# Simulation data
sim_data_file = os.path.join(path, 'SIM_DATA.csv')
sim_data_raw = np.genfromtxt(sim_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(sim_data_raw.dtype.names)
sim_ids = list(dict.fromkeys(sim_data_raw['sim_id']))
print(sim_ids)
# gather simulation data for all alt_expt_ids
sim_data = {}
for sim_id in sim_ids:
    if sim_id in expt_ids:
        # get all data points
        data = np.array([d for d in sim_data_raw if d['sim_id'] == sim_id])
        # number of experiments and observables
        n_expts = n_alt_expt_ids[sim_id]
        n_obs = len(obs_names)
        # calculate the number of time points
        n_timepoints = len(data) // (n_expts * n_obs)
        # extract indices for each experiment
        for exp in range(n_expts):
            indices = np.hstack([
                np.arange(start_idx, start_idx + n_timepoints)  # Select contiguous blocks
                for start_idx in range(exp * n_timepoints, len(data), n_expts * n_timepoints)
            ])
            alt_expt_id = list(expt_data.keys())[len(sim_data)]
            sim_data[alt_expt_id] = data[indices]
    else:
        sim_data[sim_id] = np.array([d for d in sim_data_raw if d['sim_id'] == sim_id])

for key in sim_data.keys():
    print(key)
    # print('    ', sim_data[key].dtype.names)
    # for d in sim_data[key][key2]:
    #     print('    ', d)
# quit()

print('----------')

obs_labels = \
    {'Bone_tot': 'Bone Density',
     'C_obs': 'Osteoclasts',
     'OB_tot': 'Osteoblasts',
     'Tumor_tot': 'Tumor Cells'}

legend_labels = \
    {'Bennett2024': 'This work: No Drug-Tumor',
     'Johnson2011 (PBS-Tumor)': 'Johnson (2011): No Drug-Tumor',
     'Johnson2011 (ZA-Tumor)': 'Johnson (2011): ZA-Tumor',
     'Johnson2011 (ZA-NoTumor)': 'Johnson (2011): ZA-No Tumor',
     '3': 'This work: ZA-Tumor'}

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
colors = {'Bennett2024': cycle[0],
          'Johnson2011 (PBS-Tumor)': cycle[1],
          'Johnson2011 (ZA-Tumor)': cycle[2],
          'Johnson2011 (ZA-NoTumor)': cycle[3],
          '3': cycle[4]}

def make_fig(sim_expt_ids, fig_name, expt_data_only=True):
    nrows = 2
    ncols = 2
    fig, axs = plt.subplots(nrows, ncols, sharex='all', figsize=[1*6.4, 1.2*4.8], constrained_layout=True)
    axs_tumor2 = axs[1, 1].twinx()  # 2nd y-axis for Johnson et al. (2011) tumor data
    print(sim_expt_ids)
    xticks = [0, 7, 14, 21, 28] \
        if any('Johnson' in sim_expt_id for sim_expt_id in sim_expt_ids) \
        else [0, 7, 14, 21]
    color = 'black' \
        if any(j_tumor in sim_expt_ids for j_tumor in ['Johnson2011 (PBS-Tumor)', 'Johnson2011 (ZA-Tumor)']) \
        else 'white'
    row = col = 0
    for obs_name in obs_names:
        axs[row, col].set_title(obs_labels[obs_name], fontweight='bold', fontsize=14)
        for sim_expt_id in sim_expt_ids:
            # plot Johnson et al. (2011) tumor data on 2nd y-axis
            ax = axs_tumor2 \
                if obs_name == 'Tumor_tot' and \
                   (sim_expt_id == 'Johnson2011 (PBS-Tumor)' or sim_expt_id == 'Johnson2011 (ZA-Tumor)') \
                else axs[row, col]
            # check if experimental data exists
            expt_time = [] if sim_expt_id not in expt_data.keys() \
                else [d['time'] for d in expt_data[sim_expt_id] if d['observable'] == obs_name]
            # plot simulation data
            if not expt_data_only or expt_data_only and len(expt_time) > 0:
                sim_time = [d['time'] for d in sim_data[sim_expt_id] if d['observable'] == obs_name]
                yval_min = [d['yval_min'] for d in sim_data[sim_expt_id] if d['observable'] == obs_name]
                yval_max = [d['yval_max'] for d in sim_data[sim_expt_id] if d['observable'] == obs_name]
                ax.fill_between(sim_time, yval_min, yval_max, alpha=0.25, label=legend_labels[sim_expt_id],
                                           color=colors[sim_expt_id])
            # plot experimental data
            if len(expt_time) > 0:
                average = [d['average'] for d in expt_data[sim_expt_id] if d['observable'] == obs_name]
                stderr = [d['stderr'] for d in expt_data[sim_expt_id] if d['observable'] == obs_name]
                ax.errorbar(expt_time, average, yerr=stderr, capsize=6, fmt='o', ms=8, color=colors[sim_expt_id])
        if row == nrows - 1:
            axs[row, col].set_xlabel('Time (day)', fontsize=14)
        ylabel = 'Relative BV/TV (%)' if obs_name == 'Bone_tot' else 'Concentration (fM)'
        axs[row, col].set_ylabel(ylabel, fontsize=14)
        axs[row, col].set_xticks(xticks, xticks)
        axs[row, col].set_xlim(right=xticks[-1] + 2)
        axs[row, col].tick_params(axis='both', which='major', labelsize=14)
        # update row and col
        col += 1
        if col == ncols:
            col = 0
            row += 1
    # adjust axis formats
    axs[1, 1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axs[1, 1].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axs_tumor2.set_ylim(bottom=-0.1, top=axs_tumor2.get_ylim()[1])
    axs_tumor2.set_ylabel('Rel. Concentration', fontsize=14, rotation=-90, labelpad=16, color=color)
    axs_tumor2.tick_params(axis='y', which='major', labelsize=14, color=color)
    axs_tumor2.set_yticks([0, 1, 2, 3, 4], [0, 1, 2, 3, 4], color=color)
    # put the legend at the bottom of the figure
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncols=2, bbox_to_anchor=(0.5, 0.06), loc='center', fontsize=12, columnspacing=1)
    left = 0
    bottom = 0.12
    fig.get_layout_engine().set(rect=(left, bottom, 1 - left, 1 - bottom))  # (left, bottom, width, height)
    plt.savefig(fig_name, format='pdf')


if __name__ == '__main__':

    make_fig(['Bennett2024', 'Johnson2011 (PBS-Tumor)', 'Johnson2011 (ZA-Tumor)', 'Johnson2011 (ZA-NoTumor)'],
             'Figure_3.pdf', expt_data_only=True)
    make_fig(['Bennett2024', 'Johnson2011 (PBS-Tumor)'], 'Figure_4.pdf', expt_data_only=False)
    make_fig(['Bennett2024', '3'], 'Figure_5.pdf', expt_data_only=False)
    make_fig(['Johnson2011 (ZA-Tumor)', '3'], 'Figure_ZA.pdf', expt_data_only=False)
    make_fig(['Johnson2011 (PBS-Tumor)', 'Johnson2011 (ZA-Tumor)'],
             'Figure_B_C.pdf', expt_data_only=False)

    plt.show()
