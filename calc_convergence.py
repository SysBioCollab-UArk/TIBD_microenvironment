from pydream.convergence import Gelman_Rubin
import numpy as np
import re
import glob
import os
import itertools

path = '../AorticCalcification/SAVE/Messika'
# '/Users/leonardharris/PycharmProjects/Mitochondrial_Complex_II/SAVED/TEMP'
# 2024_06_06_1expt'

logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))
sample_files = glob.glob(os.path.join(path, 'dreamzs*params*'))

chains = []
iterations = []
for file in sample_files:
    m = re.search(r'chain_(\d+)_(\d+).npy$', file)
    chains.append(int(m.group(1)))
    iterations.append(int(m.group(2)))
chains = np.unique(chains)
all_iterations = np.unique(iterations)

for n in range(1, 3):
    iterations = all_iterations[:n]
    print('iterations:', iterations)
    for fburnin in [0.5]:  # [0.2, 0.4, 0.5, 0.6, 0.8]:
        # print('fburnin:', fburnin)
        combos = itertools.combinations(chains, len(chains))  # - 2)
        for combo in combos:
            print(combo)
            sampled_params = []
            for chain in combo:
                files = []
                for iter in iterations:
                    files += [file for file in sample_files if re.search(r'_%d_%d.npy' % (chain, iter), file)]
                sampled_params.append(np.concatenate(tuple(np.load(file) for file in files)))
            GR = Gelman_Rubin(sampled_params, fburnin)
            print('Gelman-Rubin metric:')
            print(GR)
            print(np.all(GR < 1.2))
            print()
