from param_calibration import SimulationProtocol
import numpy as np
from collections.abc import Iterable
from math import isclose


def find_exact_index(tspan, target, rel_tol=1e-9):
    """Find the index of an exact floating-point match in a list.

    Raises an Exception if no exact match is found.
    """
    for i, t in enumerate(tspan):
        if isclose(t, target, rel_tol=rel_tol):
            return i
    raise ValueError(f"No exact match found for target {target} in tspan.")


def find_closest_index(tspan, target):
    best_idx = 0
    min_diff = abs(tspan[0] - target)

    for i in range(1, len(tspan)):  # Start from index 1 to avoid redundant checks
        diff = abs(tspan[i] - target)
        if diff < min_diff:
            best_idx, min_diff = i, diff
        else:
            break  # Exit early if difference starts increasing

    return best_idx


def is_list_like(obj):
    """Check if an object is iterable but NOT a string or bytes."""
    return isinstance(obj, Iterable) and not isinstance(obj, (str, bytes))


def validate_time_perturb_value(time_perturb_value):
    """
    Validates and **modifies in place** the given dictionary to ensure it follows the expected format:
    {time: [(perturb_name, value), ...]} or {time: (perturb_name, value)}

    - Keys must be numeric (int or float).
    - Values must be either:
      - A single list-like object of exactly 2 elements
      - A list-like object of list-like objects, where each inner object has 2 elements
    - The first element (perturb_name) must be a string.
    - The second element (value) must be numeric (int or float).

    Returns the modified dictionary with a standardized format.
    """
    if not isinstance(time_perturb_value, dict):
        raise ValueError("Input must be a dictionary.")

    for time in list(time_perturb_value.keys()):  # Iterate over keys safely
        # Check that time keys are numeric
        if not isinstance(time, (int, float)):
            raise ValueError(f"Invalid key {time}: Time keys must be integers or floats.")

        perturbations = time_perturb_value[time]

        # Ensure perturbations is list-like (but not string)
        if not is_list_like(perturbations):
            raise ValueError(f"Invalid value at time {time}: Must be list-like.")

        # If perturbations has length 2 and neither element is list-like, convert to a list of one tuple
        if len(perturbations) == 2 and not any(is_list_like(perturb) for perturb in perturbations):
            time_perturb_value[time] = [perturbations]  # Modify in place

        # Ensure perturbations is now a list-like object of list-like objects
        perturbations = time_perturb_value[time]  # Reassign after modification
        if not all(is_list_like(perturb) and len(perturb) == 2 for perturb in perturbations):
            raise ValueError(
                f"Invalid perturbation format {perturbations} at time {time}: Must be list-like objects with 2 " +
                "elements each."
            )

        for perturb_name, value in perturbations:
            # Check that the first element is a string (parameter/species name)
            if not isinstance(perturb_name, str):
                raise ValueError(f"Invalid parameter name '{perturb_name}' at time {time}: Must be a string.")

            # Check that the second element (value) is numeric (int or float)
            if not isinstance(value, (int, float)):
                raise ValueError(f"Invalid value {value} for '{perturb_name}' at time {time}: Must be a number " +
                                 "(int or float).")

    return time_perturb_value  # ✅ Return the modified dictionary


class SequentialInjections(SimulationProtocol):
    def __init__(self, solver, t_equil=None, time_perturb_value=None):
        super().__init__(solver, t_equil)
        self.time_perturb_value = {} if time_perturb_value is None else validate_time_perturb_value(time_perturb_value)
        # use a set ({}) to get unique perturb names
        perturbs = {perturb[0] for perturbations in time_perturb_value.values() for perturb in perturbations}
        # create a dictionary to store the indices of all perturbations (initials or param_values)
        self.perturb_idx = dict(zip(perturbs, [None for p in perturbs]))
        # store the species and parameter names, so don't have to keep generating these lists over and over
        self.sp_names = [str(sp) for sp in self.solver.model.species]
        self.par_names = [p.name for p in self.solver.model.parameters]

    def run(self, tspan, param_values):

        def save_output():
            if output is None:
                # only keep output points for t >= tspan[0]
                idx_keep = [idx for idx in range(len(tspan_i)) if tspan_i[idx] >= tspan[0]]
                return sim_output.all[idx_keep]
            else:
                # remove first output point, since time of perturbation is included in last tspan, if it's there
                return np.append(output, sim_output.all[1:])

        if len(self.time_perturb_value.keys()) == 0:  # just do the default simulation protocol
            return super().run(tspan, param_values)

        # loop over sorted perturbation times, in case any are < tspan[0]
        sorted_perturb_times = sorted(self.time_perturb_value.keys())
        # equilibration
        if self.t_equil is not None:
            min_tsim = min(sorted_perturb_times[0], tspan[0])  # min of perturb time and tspan[0]
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
        for i, pert_time in enumerate(sorted_perturb_times):
            # if time <= tspan[0], don't run a simulation
            if pert_time > tspan[0] or tspan[0] > pert_time > pert_time_last:  # run a simulation
                if i == 0:  # run a simulation with no perturbation
                    tspan_i = [t for t in tspan if t < pert_time] + [pert_time]
                else:
                    tspan_i = [pert_time_last] + [t for t in tspan if pert_time_last < t < pert_time] + [pert_time]
                sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
                # save output
                output = save_output()
                # initials for next iteration
                initials = sim_output.species[-1]
                # if there are NaNs in the initials, just return the current output
                if np.any(np.isnan(initials)):
                    return output
            # save the perturbation time for the next iteration
            pert_time_last = pert_time
            # add perturbations to initials or param_values for next iteration
            for perturb, pert_value in self.time_perturb_value[pert_time]:
                if perturb in self.sp_names:
                    if self.perturb_idx[perturb] is None:
                        self.perturb_idx[perturb] = self.sp_names.index(perturb)
                    initials[self.perturb_idx[perturb]] = pert_value
                elif perturb in self.par_names:
                    if self.perturb_idx[perturb] is None:
                        self.perturb_idx[perturb] = self.par_names.index(perturb)
                    param_values[self.perturb_idx[perturb]] = pert_value
                else:
                    raise Exception("Perturbation '%s' not found in either model.species or model.parameters." % perturb)
        # final perturbation
        tspan_i = [pert_time_last] + [t for t in tspan if t > pert_time_last]
        sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)

        return save_output()


class ParallelExperiments(SimulationProtocol):
    # noinspection PyMissingConstructor
    def __init__(self, solver, t_equil=None, time_perturb_value=None, scale_by_eidx_time=None):
        # if a dict is passed, wrap it in a list
        if isinstance(time_perturb_value, dict):
            time_perturb_value = [time_perturb_value]
        # create SequentialInjections objects to run simulations
        self.sim_protocols = []
        for tpv in time_perturb_value:
            self.sim_protocols.append(SequentialInjections(solver, t_equil, tpv))
        # dict with observables as keys and tuples with experiment index and times of data points to scale by
        # Example: scale_by_eidx_time = {'pSmad23': (0, 14)}
        self.scale_by_eidx_time = {} if scale_by_eidx_time is None else scale_by_eidx_time

    def run(self, tspan, param_values):
        # run simulations
        output = []
        for protocol in self.sim_protocols:
            output.append(protocol.run(tspan=tspan, param_values=param_values.copy()))
        # scale output arrays
        for obs in self.scale_by_eidx_time.keys():
            e_idx = self.scale_by_eidx_time[obs][0]  # expt index
            t_idx = find_closest_index(tspan, self.scale_by_eidx_time[obs][1])  # time pt index
            scale_val = output[e_idx][obs][t_idx]  # value to scale by
            for out in output:
                out[obs] /= scale_val
        # concatenate output arrays
        output = np.concatenate(output)

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
    # tumor_injection = SequentialInjections(sim, t_equil=500, time_perturb_value={0: ('Tumor()', 1)})

    # Experiment B
    time_perturb_value = [{0: ('Tumor()', 1)},
                          {0: ('Tumor()', 1), 6: ('Bisphos()', 1)}]
    scale_by_eidx_time = {'Tumor_tot': (0, 14)}  # scale by output at t=14 in expt 0
    multi_exp_injection = ParallelExperiments(sim, t_equil=500, time_perturb_value=time_perturb_value,
                                              scale_by_eidx_time=scale_by_eidx_time)

    tspan = np.array([0, 6, 7, 14, 21, 28])
    result = multi_exp_injection.run(tspan=tspan, param_values=sim.param_values[0])

    # Experiment C
    # bisphos_injection = SequentialInjections(sim, t_equil=500, time_perturb_value={6: ('Bisphos()', 1)})

    tspan_mask = {
        'Bone_tot': [[True, False, True, True, True, True], [True, True, True, True, True, True]],
        'Tumor_tot': [[False, False, False, True, True, True], [False, False, False, True, True, True]],
    }

    for obs in observables:
        plt.figure(obs)
        colors = [line.get_color() for i, line in enumerate(plt.gca().get_lines()) if i % 3 == 0]
        for i, yvals in enumerate(result[obs].reshape(len(tspan_mask[obs]), -1)):
            plt.plot(tspan[tspan_mask[obs][i]], yvals[tspan_mask[obs][i]], lw=2, color=colors[i])

    plt.show()
