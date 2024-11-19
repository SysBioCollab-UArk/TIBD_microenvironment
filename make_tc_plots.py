import numpy as np
import matplotlib.pyplot as plt
import os

expt_data_file = os.path.join('DATA', 'TIBD_PopD_Data_stderr5.csv')
expt_data = np.genfromtxt(expt_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(expt_data.dtype.names)

sim_data_file = os.path.join('.', 'SIM_DATA.csv')
sim_data = np.genfromtxt(sim_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(sim_data.dtype.names)

obs_names = np.unique(sim_data['observable'])
sim_ids = np.unique(sim_data['sim_id'])
print(obs_names)
print(sim_ids)

obs_labels = \
    {'Bone_tot': 'Bone Density',
     'C_obs': 'Osteoclasts',
     'OB_tot': 'Osteoblasts',
     'Tumor_tot': 'Tumor Cells'}

legend_labels = \
    {'A': 'This work: No Drug-Tumor',
     'B': 'Johnson (2011): No Drug-Tumor',
     'C': 'Johnson (2011): ZA-Tumor',
     'D': 'Johnson (2011): ZA-No Tumor',
     '4': 'This work: ZA-Tumor'}

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
colors = {'A': cycle[0], 'B': cycle[1], 'C': cycle[2], 'D': cycle[3], '4': cycle[4]}

# FIGURE 3: plots for just experimental data
nrows = 2
ncols = 2
fig, axs = plt.subplots(nrows, ncols, sharex='all', figsize=[1*6.4, 1.2*4.8], constrained_layout=True)
row = col = 0
for obs_name in obs_names:
    axs[row, col].set_title(obs_labels[obs_name], fontweight='bold', fontsize=14)
    for sim_id in sim_ids:
        # check if experimental data exists
        expt_time = [d['time'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
        if len(expt_time) > 0:
            # plot simulation data
            sim_time = [d['time'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
            yvals_min = [d['yvals_min'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
            yvals_max = [d['yvals_max'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
            axs[row, col].fill_between(sim_time, yvals_min, yvals_max, alpha=0.25, label=legend_labels[sim_id],
                                       color=colors[sim_id])
            # plot experimental data
            avg = [d['average'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            stderr = [d['stderr'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            axs[row, col].errorbar(expt_time, avg, yerr=stderr, capsize=6, fmt='o', ms=8, color=colors[sim_id])
    if row == nrows - 1:
        axs[row, col].set_xlabel('Time (day)', fontsize=14)
    ylabel = 'Relative BV/TV (%)' if obs_name == 'Bone_tot' else 'Concentration (fM)'
    axs[row, col].set_ylabel(ylabel, fontsize=14)
    axs[row, col].set_xlim(right=30)
    axs[row, col].set_xticks([0, 7, 14, 21, 28], [0, 7, 14, 21, 28])
    axs[row, col].tick_params(axis='both', which='major', labelsize=14)
    # update row and col
    col += 1
    if col == ncols:
        col = 0
        row += 1
# put the legend at the bottom of the figure
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, ncols=2, bbox_to_anchor=(0.5, 0.06), loc='center', fontsize=12, columnspacing=1)
left = 0
bottom = 0.12
fig.get_layout_engine().set(rect=(left, bottom, 1 - left, 1 - bottom))  # (left, bottom, width, height)
plt.savefig('Figure_3.pdf', format='pdf')

# FIGURE 4: comparing PBS-Tumor from this work and Johnson (2011)
nrows = 2
ncols = 2
fig, axs = plt.subplots(nrows, ncols, sharex='all', figsize=[1*6.4, 1.2*4.8], constrained_layout=True)
row = col = 0
for obs_name in obs_names:
    axs[row, col].set_title(obs_labels[obs_name], fontweight='bold', fontsize=14)
    for sim_id in ['A', 'B']:
        # plot simulation data
        sim_time = [d['time'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        yvals_min = [d['yvals_min'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        yvals_max = [d['yvals_max'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        axs[row, col].fill_between(sim_time, yvals_min, yvals_max, alpha=0.25, label=legend_labels[sim_id],
                                   color=colors[sim_id])
        # plot experimental data, if it exists
        expt_time = [d['time'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
        if len(expt_time) > 0:
            avg = [d['average'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            stderr = [d['stderr'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            axs[row, col].errorbar(expt_time, avg, yerr=stderr, capsize=6, fmt='o', ms=8, color=colors[sim_id])
    if row == nrows - 1:
        axs[row, col].set_xlabel('Time (day)', fontsize=14)
    ylabel = 'Relative BV/TV (%)' if obs_name == 'Bone_tot' else 'Concentration (fM)'
    axs[row, col].set_ylabel(ylabel, fontsize=14)
    axs[row, col].set_xlim(right=30)
    axs[row, col].set_xticks([0, 7, 14, 21, 28], [0, 7, 14, 21, 28])
    axs[row, col].tick_params(axis='both', which='major', labelsize=14)
    # update row and col
    col += 1
    if col == ncols:
        col = 0
        row += 1
# put the legend at the bottom of the figure
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, ncols=2, bbox_to_anchor=(0.5, 0.04), loc='center', fontsize=12, columnspacing=1)
left = 0
bottom = 0.08
fig.get_layout_engine().set(rect=(left, bottom, 1 - left, 1 - bottom))  # (left, bottom, width, height)
plt.savefig('Figure_4.pdf', format='pdf')

# FIGURE 5: comparing PBS-Tumor and ZA-Tumor from this work
nrows = 2
ncols = 2
fig, axs = plt.subplots(nrows, ncols, sharex='all', figsize=[1*6.4, 1.2*4.8], constrained_layout=True)
row = col = 0
for obs_name in obs_names:
    axs[row, col].set_title(obs_labels[obs_name], fontweight='bold', fontsize=14)
    for sim_id in ['A', '4']:
        # plot simulation data
        sim_time = [d['time'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        yvals_min = [d['yvals_min'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        yvals_max = [d['yvals_max'] for d in sim_data if d['observable'] == obs_name and d['sim_id'] == sim_id]
        axs[row, col].fill_between(sim_time, yvals_min, yvals_max, alpha=0.25, label=legend_labels[sim_id],
                                   color=colors[sim_id])
        # plot experimental data, if it exists
        expt_time = [d['time'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
        if len(expt_time) > 0:
            avg = [d['average'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            stderr = [d['stderr'] for d in expt_data if d['observable'] == obs_name and d['expt_id'] == sim_id]
            axs[row, col].errorbar(expt_time, avg, yerr=stderr, capsize=6, fmt='o', ms=8, color=colors[sim_id])
    if row == nrows - 1:
        axs[row, col].set_xlabel('Time (day)', fontsize=14)
    ylabel = 'Relative BV/TV (%)' if obs_name == 'Bone_tot' else 'Concentration (fM)'
    axs[row, col].set_ylabel(ylabel, fontsize=14)
    axs[row, col].set_xlim(right=22)
    axs[row, col].set_xticks([0, 7, 14, 21], [0, 7, 14, 21])
    axs[row, col].tick_params(axis='both', which='major', labelsize=14)
    # update row and col
    col += 1
    if col == ncols:
        col = 0
        row += 1
# put the legend at the bottom of the figure
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, ncols=2, bbox_to_anchor=(0.5, 0.04), loc='center', fontsize=12, columnspacing=1)
left = 0
bottom = 0.08
fig.get_layout_engine().set(rect=(left, bottom, 1 - left, 1 - bottom))  # (left, bottom, width, height)
plt.savefig('Figure_5.pdf', format='pdf')

plt.show()