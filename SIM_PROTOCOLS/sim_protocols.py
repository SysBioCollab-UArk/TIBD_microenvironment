from param_calibration import SimulationProtocol
import numpy as np


class SequentialInjections(SimulationProtocol):
    def __init__(self, solver, t_equil=None, perturb_time_amount=None):
        super().__init__(solver, t_equil)
        self.perturb_time_amount = perturb_time_amount

    def run(self, tspan, param_values):
        if self.perturb_time_amount is None:  # just do the default simulation protocol
            output = super().run(tspan, param_values)
        elif isinstance(self.perturb_time_amount, dict):
            # Get the sorted perturbation times here, in case any are < tspan[0]
            sorted_pta_by_time = sorted(list(self.perturb_time_amount.items()), key=lambda x: x[1][0])
            # equilibration
            if self.t_equil is not None:
                min_tsim = min(sorted_pta_by_time[0][1][0], tspan[0])  # min of perturb time and tspan[0]
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
            for i, pta in enumerate(sorted_pta_by_time):
                perturb = pta[0]
                time = pta[1][0]
                amount = pta[1][1]
                # if time <= tspan[0], don't run a simulation
                if time > tspan[0]:  # run a simulation
                    if i == 0:  # run a simulation with no perturbation
                        tspan_i = [t for t in tspan if t <= time]
                    else:
                        tspan_i += [t for t in tspan if tspan_i[0] < t <= time]
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
                # add perturbation to initials for next iteration
                tspan_i = [time]  # TODO: need to handle the case where multiple perturbations occur before tspan[0]
                # get index of perturbation
                idx_perturb = [str(sp) for sp in self.solver.model.species].index(perturb)  # TODO: speed up by saving this index
                initials[idx_perturb] = amount
            # final perturbation
            tspan_i += [t for t in tspan if t > tspan_i[0]]
            sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
            output = sim_output.all if output is None else np.append(output, sim_output.all[1:])
        else:
            raise Exception("'perturb_time_amount' must be either a dict or None.")

        return output


if __name__ == "__main__":
    from pysb.simulator import ScipyOdeSimulator
    from MODELS.TIBD_PopD_v1 import model

    sim = ScipyOdeSimulator(model)
    protocol = SequentialInjections(sim, t_equil=500, perturb_time_amount={'Tumor()': (0, 1), 'Bisphos()': (6, 1)})
    result = protocol.run(tspan=[0, 6, 7, 14, 21, 28], param_values=sim.param_values[0])

    print(result)
    print(result.shape)
    print(result.dtype.names)
