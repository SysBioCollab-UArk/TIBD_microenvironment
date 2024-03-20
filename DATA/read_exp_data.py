import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container

data = np.genfromtxt('Mouse_Data_March2024.csv', dtype=None, delimiter=',', names=True, encoding="utf_8_sig")

print(data.dtype.names)
# print(data)

cond_cellTypes = np.unique([(d['Condition'], d['CellType']) for d in data], axis=0)
# print(cond_cellTypes)

for cc in cond_cellTypes:
    print(cc)
    plt.figure(cc[1], constrained_layout=True)

    # plot all data points
    time = [d['Week'] for d in data if d['Condition'] == cc[0] and d['CellType'] == cc[1]]
    conc = [d['Concentration'] for d in data if d['Condition'] == cc[0] and d['CellType'] == cc[1]]
    weeks = np.unique(time)
    p = plt.plot(time, conc, 'o', mfc='none')  # , label='%s, %s' % (cc[1], cc[0]))

    # plot averages and standard errors
    avg_conc = [np.mean([conc[i] for i in range(len(conc)) if time[i] == w]) for w in weeks]
    se_conc = [[conc[i] for i in range(len(conc)) if time[i] == w] for w in weeks]
    se_conc = [np.std(x, ddof=1) / np.sqrt(len(x)) for x in se_conc]
    plt.errorbar(weeks, avg_conc, se_conc, marker='o', ms=10, capsize=10, lw=3, color=p[0].get_color(),
                 label='%s, %s' % (cc[1], cc[0]))

    plt.xlabel('Week')
    plt.ylabel('Concentration (pM)')
    plt.xticks(ticks=weeks, labels=weeks)
    # Remove error bars from legend
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    plt.legend(handles, labels, loc=0)

plt.show()
