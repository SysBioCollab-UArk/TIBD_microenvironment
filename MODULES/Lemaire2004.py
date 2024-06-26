from pysb import *
from pysb.util import alias_model_components


def create_model_elements():

    # This is a reproduction of the computational model presented in Lemaire et al., J. Theor. Biol. 229, 293-309 (2004)
    # https://doi.org/10.1016/j.jtbi.2004.03.023

    # Note that all units have been changed to fM (1 pM = 1000 fM)

    # Vol = 8.7e-9  # liters

    # cell types
    Monomer('R')  # responding osteoblasts
    Monomer('B')  # active osteoblasts
    Monomer('C')  # active osteoclasts

    # initial concentrations
    Parameter('R_0', 0.7734e-3 * 1000)  # * 1e-12 * N_A * Vol)  # fM
    Parameter('B_0', 0.7282e-3 * 1000)  # * 1e-12 * N_A * Vol)  # fM
    Parameter('C_0', 0.9127e-3 * 1000)  # * 1e-12 * N_A * Vol)  # fM

    # rate parameters
    Parameter('Cs', 5e-3 * 1000)  # * 1e-12 * N_A * Vol  # fM
    Parameter('DA', 0.7)  # 0.7) # /day
    Parameter('dB', 0.7)  # /day
    Parameter('DC', 2.1e-3 * 1000)  # * 1e-12 * N_A * Vol  # fM/day
    Parameter('DR', 7e-4 * 1000)  # * 1e-12 * N_A * Vol  # fM/day
    Parameter('f0', 0.05)  # unitless
    Parameter('IL', 0)  # * 1e-12 * N_A * Vol)  # fM/day (range: 0-1e9)
    Parameter('IO', 0)  # * 1e-12 * N_A * Vol)  # fM/day (range: 0-1e9)
    Parameter('IP_const', 0)  # fM/day (range: 0-1e9)
    Parameter('K', 10 * 1000)  # * 1e-12 * N_A * Vol  # fM
    Parameter('k1', 1e-2 / 1000)  # / (1e-12 * N_A * Vol)  # /fM-day
    Parameter('k2', 10)  # /day
    Parameter('k3', 5.8e-4 / 1000)  # / (1e-12 * N_A * Vol)  # /fM-day
    Parameter('k4', 1.7e-2)  # /day
    Parameter('k5', 0.02 / 1000)  # / (1e-12 * N_A * Vol) # /fM-day
    Parameter('k6', 3)  # /day
    Parameter('kB', 0.189)  # 0.189) # /day
    Parameter('KLP', 3e6)  # unitless
    Parameter('kO', 0.35)  # /day
    Parameter('KOP', 2e5)  # /day
    Parameter('kP', 86)  # /day
    Parameter('rL', 1e3 * 1000)  # * 1e-12 * N_A * Vol  # fM/day
    Parameter('SP', 250 * 1000)  # * 1e-12 * N_A * Vol  # fM/day

    alias_model_components()

    # make injected PTH/PTHrP an expression, rather than a constant
    Expression('IP', IP_const)  # * 1e-12 * N_A * Vol)  # fM/day

    # initial conditions
    Initial(R(), R_0)
    Initial(B(), B_0)
    Initial(C(), C_0)

    # observables
    Observable('R_obs', R())
    Observable('B_obs', B())
    Observable('C_obs', C())
    Observable('OB_tot', R() + B())

    # compound constants
    # Pbar = IP/kP
    # P0 = SP/kP  # fM
    # Ps = k6/k5  # fM

    alias_model_components()

    # ratios
    # Expression('pi_P', (IP/kP + P0) / (IP/kP + Ps))  # unitless
    # Expression('pi_P', (IP/kP + SP/kP) / (IP/kP + k6/k5))  # unitless
    Expression('pi_C', (C_obs + f0*Cs) / (C_obs + Cs))  # unitless
    # Expression('pi_L', k3/k4 * KLP*pi_P*B_obs * (1 + IL/rL) / (1 + k3*K/k4 + k1/k2/kO*(KOP/pi_P*R_obs + IO)))  # unitless

    alias_model_components()

    # cell dynamics rules
    Expression('DR_pi_C', DR * pi_C)  # fM/day
    Expression('DB_over_pi_C', f0 * dB / pi_C)  # /day
    # Expression('DC_pi_L', DC * pi_L)  # fM/day
    Expression('DA_pi_C', DA * pi_C)  # /day

    # Expression('DR_pi_C', DR * (C_obs + f0*Cs) / (C_obs + Cs))  # fM/day
    # Expression('DB_over_pi_C', f0 * dB / ((C_obs + f0*Cs) / (C_obs + Cs)))  # /day
    Expression('DC_pi_L', DC * k3/k4 * KLP*((IP/kP + SP/kP) / (IP/kP + k6/k5))*B_obs * (1 + IL/rL) / (1 + k3*K/k4 + k1/k2/kO*(KOP/((IP/kP + SP/kP) / (IP/kP + k6/k5))*R_obs + IO)))  # fM/day
    # Expression('DA_pi_C', DA * (C_obs + f0*Cs) / (C_obs + Cs))  # /day

    alias_model_components()

    Rule('ROB_creation', None >> R(), DR_pi_C)
    Rule('ROB_to_AOB', R() >> B(), DB_over_pi_C)
    Rule('AOB_death', B() >> None, kB)
    Rule('AOC_creation_death', None | C(), DC_pi_L, DA_pi_C)


if __name__ == '__main__':
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()
    create_model_elements()

    # cell injections
    Parameter('kf_AOB', 0)  # fM/day
    Parameter('kf_AOC', 0)  # fM/day
    Parameter('kf_ROB', 0)  # fM/day
    Parameter('kr_AOB', 0)  # /day
    Parameter('kr_AOC', 0)  # /day
    Parameter('kr_ROB', 0)  # /day
    Rule('AOB_injection', None | B(), kf_AOB, kr_AOB)
    Rule('AOC_injection', None | C(), kf_AOC, kr_AOC)
    Rule('ROB_injection', None | R(), kf_ROB, kr_ROB)

    # simulations
    sim = ScipyOdeSimulator(model, verbose=True)

    tspan1 = np.linspace(0, 20, 201)
    tspan2 = np.linspace(20, 80, 601)
    tspan3 = np.linspace(80, 140, 601)

    # perturbations
    perturb = [
        [{'kf_AOB': 0}, {'kf_AOB': 1e-4 * 1000}, {'kf_AOB': 0}],
        [{'kf_AOC': 0}, {'kf_AOC': 1e-4 * 1000}, {'kf_AOC': 0}],
        [{'kf_ROB': 0}, {'kf_ROB': 1e-4 * 1000}, {'kf_ROB': 0}],
        [{'kr_AOB': 0}, {'kr_AOB': 8.3e-2}, {'kr_AOB': 0}],  # value different from paper
        [{'kr_AOC': 0}, {'kr_AOC': 2.9e-1}, {'kr_AOC': 0}],  # value different from paper
        [{'kr_ROB': 0}, {'kr_ROB': 1.2e-1}, {'kr_ROB': 0}],  # value different from paper
        [{'IP_const': 0}, {'IP_const': 1e3 * 1000}, {'IP_const': 0}],
        [{'IO': 0}, {'IO': 2e5 * 1000}, {'IO': 0}],
        [{'IL': 0, 'IO': 0}, {'IL': 1e4 * 1000, 'IO': 0},
         {'IL': 1e4 * 1000, 'IO': 9e4 * 1000}]  # value of IL different from paper
    ]

    colors = ['b', 'r', 'g']
    labels = ['ROB', 'AOB', 'AOC']
    fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, figsize=(10, 8), constrained_layout=True)
    row = 0
    col = 0
    for p in perturb:
        print(p)
        output1 = sim.run(param_values=p[0], tspan=tspan1)
        output2 = sim.run(param_values=p[1], tspan=tspan2, initials=output1.species[-1])
        output3 = sim.run(param_values=p[2], tspan=tspan3, initials=output2.species[-1])
        for obs, color, label in zip(model.observables, colors, labels):
            # convert concentrations back to pM for plotting
            axs[row, col].plot(tspan1, output1.observables[obs.name] / 1000, lw=2, color=color, label=label)
            axs[row, col].plot(tspan2, output2.observables[obs.name] / 1000, lw=2, color=color)
            axs[row, col].plot(tspan3, output3.observables[obs.name] / 1000, lw=2, color=color)
        if row == 2:
            axs[row, col].set_xlabel('time (d)')
        if col == 0:
            axs[row, col].set_ylabel('concentration (pM)')
        axs[row, col].set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
        axs[row, col].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axs[row, col].legend(loc=0)
        if (row+1) % 3 == 0:
            col += 1
            row = 0
        else:
            row += 1

    plt.show()
