from param_calibration import SimulationProtocol
import numpy as np


class SequentialInjections(SimulationProtocol):
    def __init__(self, solver, t_equil=None, perturb_time_value=None):
        super().__init__(solver, t_equil)
        self.perturb_time_value = perturb_time_value

    def run(self, tspan, param_values):
        if self.perturb_time_value is None:  # just do the default simulation protocol
            output = super().run(tspan, param_values)
        elif isinstance(self.perturb_time_value, dict):
            # Get the sorted perturbation times here, in case any are < tspan[0]
            sorted_ptv_by_time = sorted(list(self.perturb_time_value.items()), key=lambda x: x[1][0])
            # equilibration
            if self.t_equil is not None:
                min_tsim = min(sorted_ptv_by_time[0][1][0], tspan[0])  # min of perturb time and tspan[0]
                out = self.solver.run(tspan=np.linspace(-self.t_equil + min_tsim, min_tsim, 2),
                                      param_values=param_values)
                initials = out.species[-1]
                # if there are NaNs in the initials, just return the current output
                if np.any(np.isnan(initials)):
                    return out.all
            else:
                # set initials for next iteration
                initials = self.solver.initials[0]
            # sort drug treatments by time of application and loop over them
            output = None
            pert_time_last = np.inf
            for i, ptv in enumerate(sorted_ptv_by_time):
                perturb = ptv[0]
                pert_time = ptv[1][0]
                pert_value = ptv[1][1]
                # if time <= tspan[0], don't run a simulation
                if pert_time > tspan[0] or tspan[0] > pert_time > pert_time_last:  # run a simulation
                    if i == 0:  # run a simulation with no perturbation
                        tspan_i = [t for t in tspan if t < pert_time] + [pert_time]
                    else:
                        tspan_i = [pert_time_last] + [t for t in tspan if pert_time_last < t < pert_time] + [pert_time]
                    sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
                    # save output
                    if output is None:
                        # only keep output points for t >= tspan[0]
                        idx_keep = [idx for idx in range(len(tspan_i)) if tspan_i[idx] >= tspan[0]]
                        output = sim_output.all[idx_keep]
                    else:
                        # remove first output point, since time of perturbation is included in last tspan, if it's there
                        output = np.append(output, sim_output.all[1:])
                    # initials for next iteration
                    initials = sim_output.species[-1]
                    # if there are NaNs in the initials, just return the current output
                    if np.any(np.isnan(initials)):
                        return output
                # save the perturbation time for the next iteration
                pert_time_last = pert_time
                # add perturbation to initials or param_values for next iteration
                sp_names = [str(sp) for sp in self.solver.model.species]
                par_names = [p.name for p in self.solver.model.parameters]
                if perturb in sp_names:
                    idx_perturb = sp_names.index(perturb)  # TODO: speed up by saving this index
                    initials[idx_perturb] = pert_value
                elif perturb in par_names:
                    idx_perturb = par_names.index(perturb)  # TODO: speed up by saving this index
                    param_values[idx_perturb] = pert_value
                else:
                    raise Exception("Perturbation '%s' not found in either model.species or model.parameters.")
            # final perturbation
            tspan_i = [pert_time_last] + [t for t in tspan if t > pert_time_last]
            sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
            output = sim_output.all if output is None else np.append(output, sim_output.all[1:])
        else:
            raise Exception("'perturb_time_value' must be either a dict or None.")

        return output

class ParallelExperiments(SimulationProtocol):
    # noinspection PyMissingConstructor
    def __init__(self, solver, t_equil=None, perturb_time_value=None, scale_by_idx=None):
        # if only one dict is passed, make it a list
        if len(np.array(perturb_time_value).shape) == 0:
            perturb_time_value = [perturb_time_value]
        # create SequentialInjections objects to run simulations
        self.sim_protocols = []
        for ptv in perturb_time_value:
            self.sim_protocols.append(SequentialInjections(solver, t_equil, ptv))
        # dict with observables as keys and indices of data points to scale by as values
        self.scale_by_idx = {} if scale_by_idx is None else scale_by_idx

    def run(self, tspan, param_values):
        # run simulations
        output = []
        for protocol in self.sim_protocols:
            output.append(protocol.run(tspan=tspan, param_values=param_values))
        # concatenate output arrays
        output = np.concatenate(output)
        # scale output arrays
        for obs in self.scale_by_idx.keys():
            output[obs] /= output[obs][self.scale_by_idx[obs]]

        return output


if __name__ == "__main__":
    from MODELS.TIBD_PopD_v1 import model
    from MODULES.perturbations import add_bisphosphonate_components
    from pysb.simulator import ScipyOdeSimulator
    import os
    import matplotlib.pyplot as plt

    # read in experimental data
    datafile = os.path.join('..', 'DATA', 'TIBD_PopD_Data.csv')
    data = np.genfromtxt(datafile, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    print(data.dtype.names)

    observables = np.unique([d['observable'] for d in data if d['expt_id'] == 'B'])
    for obs in observables:
        plt.figure(obs, constrained_layout=True)
        expts = np.unique([d['alt_expt_id'] for d in data if d['observable'] == obs and d['expt_id'] == 'B'])
        for expt in expts:
            xvals = [d['time'] for d in data if d['observable'] == obs and d['alt_expt_id'] == expt]
            yvals = [d['average'] for d in data if d['observable'] == obs and d['alt_expt_id'] == expt]
            stderr = [d['stderr'] for d in data if d['observable'] == obs and d['alt_expt_id'] == expt]
            plt.errorbar(xvals, yvals, yerr=stderr, fmt='o', ms=8, capsize=6, label=expt)
        time_units = np.unique([d['time_units'] for d in data if d['observable'] == obs])
        amount_units = np.unique([d['amount_units'] for d in data if d['observable'] == obs])
        plt.xlabel('time (%s)' % time_units[0])
        plt.ylabel('%s (%s)' % (obs, amount_units[0]))
        plt.legend(loc='best')

    # plt.show()
    # quit()

    # run simulations
    add_bisphosphonate_components()

    sim = ScipyOdeSimulator(model)

    # Experiment A
    # tumor_injection = SequentialInjections(sim, t_equil=500, perturb_time_amount={'Tumor()': (0, 1)})

    # Experiment B
    perturb_time_value = [{'Tumor()': (0, 1)}, {'Tumor()': (0, 1), 'Bisphos()': (6, 1)}]
    scale_by_idx = {'Tumor_tot': 3}  # [0, 6, 7, 14, 21, 28, 0, 6, 7, 14, 21, 28], t=14 in expt 1 is idx 3
    multi_exp_injection = ParallelExperiments(sim, t_equil=500, perturb_time_value=perturb_time_value,
                                              scale_by_idx=scale_by_idx)

    result = multi_exp_injection.run(tspan=[0, 6, 7, 14, 21, 28], param_values=sim.param_values[0])

    # Experiment C
    # bisphos_injection = SequentialInjections(sim, t_equil=500, perturb_time_amount={'Bisphos()': (6, 1)})

    tspan_mask = {
        'Bone_tot': [[True, False, True, True, True, True], [True, True, True, True, True, True]],
        'Tumor_tot': [[False, False, False, True, True, True], [False, False, False, True, True, True]],
    }
    tspan = np.array([0, 6, 7, 14, 21, 28])
    for obs in observables:
        plt.figure(obs)
        colors = [line.get_color() for i, line in enumerate(plt.gca().get_lines()) if i % 3 == 0]
        for i, yvals in enumerate(result[obs].reshape(len(tspan_mask[obs]), -1)):
            plt.plot(tspan[tspan_mask[obs][i]], yvals[tspan_mask[obs][i]], lw=2, color=colors[i])

    plt.show()
