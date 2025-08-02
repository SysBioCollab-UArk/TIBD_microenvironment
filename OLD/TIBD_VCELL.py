from MODELS.TIBD_PopD_v1 import model

from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os

from util import get_exp_data
# experimental data
exp_data = get_exp_data(os.path.join('OLD', 'Mouse_Data_May2024.csv'))
Sim = ScipyOdeSimulator(model, verbose=True)

#equilibration simulation
output = Sim.run(tspan=[0,500])

# run regular simulation
tspan =np.linspace(0,28, 281)
output = Sim.run(tspan=tspan, initials=output.species[-1])

obs_to_plot = ['C_obs', 'OB_tot', 'Tumor_tot', 'Bone_tot']
names = ['Osteoclasts', 'Osteoblasts', 'Tumor', 'Bone']
for obs, name in zip(obs_to_plot,names):
    plt.plot(tspan, output.observables[obs], lw=2, label=name)
#for obs in model.observables:
#plt.figure(obs.name)
#fig.legend((l1, l2), ('Exp', 'Data'), 'upper left')
    plt.legend(loc=0)
    plt.xlabel('Time (day)', fontsize=20)
    plt.ylabel('Amount', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    #plt.savefig(outfilename, format='pdf')
    plt.savefig(output.observables[obs_to_plot], format='pdf')
#plt.savefig(os.path.join('os.getcwd()', '%s.pdf' % obs.name), format='pdf')
# print(obs, name)
plt.show()




