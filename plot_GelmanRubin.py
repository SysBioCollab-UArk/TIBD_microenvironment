import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

dir = 'Bennett2024_Johnson2011_stderr_2abs'
path = os.path.join('MODELS', 'SAVE', 'Leonard', dir)

# dir = 'Messika'
# path = os.path.join('..', 'AorticCalcification', 'SAVE', dir)

files = glob.glob(os.path.join(path, '*Gelman*'))
iterations = sorted([int(re.search("\d\d+", os.path.split(file)[1]).group()) for file in files])
print('iterations:', iterations)

GR = []
for i in iterations:
    file = os.path.join(path, 'dreamzs_5chain_GelmanRubin_iteration_%d.txt' % i)
    GR.append(np.genfromtxt(file, dtype=None, encoding="utf_8_sig"))
GR = np.array(GR).T

cmap = plt.get_cmap('turbo')  # or 'hsv', 'viridis', etc.
colors = [cmap(i / (GR.shape[0] - 1)) for i in range(GR.shape[0])]

labels = ['p%d' % i for i in range(GR.shape[0])]

plt.figure(constrained_layout=True, figsize=(6.4 * 2, 4.8))
for GR_array, color, label in zip(GR, colors, labels):
    plt.plot(iterations, GR_array, '-', color=color, label=label)
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.axhline(1.2, ls='--', color='black', lw=2)
plt.xlabel('Iterations')
plt.ylabel('Gelman-Rubin metric')
plt.legend(loc='best', ncol=4, bbox_to_anchor=(1.02, 1))

plt.show()
