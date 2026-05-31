from MODELS.TIBD_PopD_v1 import model
from MODULES.perturbations import *
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *
import os

this_dir = os.path.dirname(__file__)
expt_data_file = os.path.join(this_dir, 'DATA', 'Canon_OPGFc_2008.csv')


if __name__ == '__main__':
    import pandas as pd
    import matplotlib.pyplot as plt

    expt_data = pd.read_csv(expt_data_file)

    expt_ids = expt_data['expt_id'].unique()
    alt_expt_ids = [expt_data[expt_data['expt_id'] == expt_id]['alt_expt_id'].unique()[0] for expt_id in expt_ids]

    labels_dict = {
        'OC_score': 'Osteoclasts score',
        'Tumor_BLI': 'Bioluminescence'
    }

    for expt_id, alt_expt_id in zip(expt_ids, alt_expt_ids):
        data_expt_id = expt_data[expt_data['expt_id'] == expt_id]
        observables = data_expt_id['observable'].unique()
        for obs in observables:
            plt.figure(obs, constrained_layout=True)
            time = data_expt_id[data_expt_id['observable'] == obs]['time']
            average = data_expt_id[data_expt_id['observable'] == obs]['average']
            stderr = data_expt_id[data_expt_id['observable'] == obs]['stderr']
            plt.errorbar(time, average, yerr=stderr, fmt='o', capsize=6, label=alt_expt_id)
            time_units = data_expt_id[data_expt_id['observable'] == obs]['time_units'].unique()[0]
            amount_units = data_expt_id[data_expt_id['observable'] == obs]['amount_units'].unique()[0]
            plt.xlabel('Time (%s)' % time_units, fontsize=14)
            plt.ylabel('%s (%s)' % (labels_dict.get(obs, obs), amount_units), fontsize=14)
            plt.tick_params(labelsize=14)
            plt.xlim(left=-1)
            plt.legend(loc='best', fontsize=12)

    plt.show()
