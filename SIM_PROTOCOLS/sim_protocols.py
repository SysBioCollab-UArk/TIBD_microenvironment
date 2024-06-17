from param_calibration import SimulationProtocol
import numpy as np


class SequentialInjections(SimulationProtocol):
    def __init__(self, solver, equil=None, perturb_day_amount=None):
        super().__init__(solver, equil)
        self.perturb_day_amount = perturb_day_amount

    def run(self, tspan, param_values):
        if self.perturb_day_amount is None:  # just do the default simulation protocol
            output = super().run(tspan, param_values)
        elif isinstance(self.perturb_day_amount, dict):
            # equilibration
            if self.equil is not None:
                out = self.solver.run(tspan=np.linspace(-self.equil, 0, 2), param_values=param_values)
                initials = out.species[-1]
            else:
                initials = None
            # sort drug treatments by time of application and loop over them
            output = None
            for i, pda in enumerate(sorted(list(self.perturb_day_amount.items()), key=lambda x: x[1][0])):
                perturb = pda[0]
                day = pda[1][0]
                amount = pda[1][1]
                if day == 0:  # don't run a simulation
                    # initials for next iteration
                    initials = self.solver.initials[0]
                else:
                    if i == 0 and day > 0:  # run a simulation with no perturbation
                        tspan_i = [t for t in tspan if t <= day]
                    else:
                        tspan_i += [t for t in tspan if tspan_i[0] < t <= day]
                    sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
                    # save output (remove first point from sim_output, since time at which perturbation is added is
                    # included in last tspan, if it's there)
                    output = sim_output.all if output is None else np.append(output, sim_output.all[1:])
                    # initials for next iteration
                    initials = sim_output.species[-1]
                # add perturbation to initials for next iteration
                tspan_i = [day]
                idx_perturb = [str(sp) for sp in self.solver.model.species].index(perturb)  # get index of perturbation
                initials[idx_perturb] = amount
            # final perturbation
            tspan_i += [t for t in tspan if t > tspan_i[0]]
            sim_output = self.solver.run(tspan=tspan_i, param_values=param_values, initials=initials)
            output = sim_output.all if output is None else np.append(output, sim_output.all[1:])
        else:
            raise Exception("'perturb_day_amount' must be either a dict or None.")

        return output
