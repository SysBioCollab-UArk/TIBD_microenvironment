import numpy as np


def tumor_injection(solver, tspan, param_values):
    # equilibration
    equil = solver.run(tspan=np.linspace(-500, 0, 2), param_values=param_values)
    # add tumor cells
    initials = equil.species[-1]
    idx_tumor = [str(sp) for sp in solver.model.species].index('Tumor()')  # get index of Tumor species
    initials[idx_tumor] = 1  # fM
    output = solver.run(tspan=tspan, param_values=param_values, initials=initials).all

    return output


def tumor_bisphos_injection(solver, tspan, param_values, day=6):
    # equilibration
    equil = solver.run(tspan=np.linspace(-500, 0, 2), param_values=param_values)
    # break tspan up into two parts
    tspan1 = [t for t in tspan if t <= day]
    tspan2 = [day] + [t for t in tspan if t > day]  # add the time drug is added to the beginning of the array
    # add tumor cells
    initials = equil.species[-1]
    idx_tumor = [str(sp) for sp in solver.model.species].index('Tumor()')  # get index of Tumor species
    initials[idx_tumor] = 1  # fM
    output1 = solver.run(tspan=tspan1, param_values=param_values, initials=initials)
    # add bisphosphonate (Zoledronic acid)
    initials = output1.species[-1]
    idx_bisphos = [str(sp) for sp in solver.model.species].index('Bisphos()')  # get index of Bisphos species
    initials[idx_bisphos] = 1  # fM
    output2 = solver.run(tspan=tspan2, param_values=param_values, initials=initials)

    # remove the first point from output2, since the time at which drug is added is included in tspan1, if it's there
    return np.append(output1.all, output2.all[1:])
