import numpy as np


def tumor_injection(solver, tspan, param_values):
    # equilibration
    equil = solver.run(tspan=np.linspace(-500, 0, 2), param_values=param_values)
    # add tumor cells
    initials = equil.species[-1]
    idx_tumor = [str(sp) for sp in solver._model.species].index('Tumor()')  # get index of Tumor species
    initials[idx_tumor] = 1  # fM
    output = solver.run(tspan=tspan, param_values=param_values, initials=initials).all

    return output
