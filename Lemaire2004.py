from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
# Parameters
# Initials
# Rules
# Observables
# simulation commands

# Cell types
Monomer('R')  # responding osteoblasts
Monomer('B')  # active osteoblasts
Monomer('C')  # active osteoclasts

# Molecules
Monomer('P', ['pr'])
Monomer('Pr', ['p'])
Monomer('L', ['k'])
Monomer('O', ['l'])
Monomer('K', ['l'])

# Observables
Observable('ROB', R())
Observable('AOB', B())
Observable('AOC', C())
Observable('P_free', P(pr=None))
Observable('Pr_free', Pr(p=None))
Observable('Pr_P_bound', Pr(p=1) % P(pr=1))
Observable('O_free', O(l=None))
Observable('L_free', L(k=None))
Observable('K_free', K(l=None))
Observable('O_L_bound', O(l=1) % L(k=1))
Observable('K_L_bound', K(l=1) % L(k=1))
Observable('L_total', L())
Observable('Pr_total', Pr())

# Expression('piP_actual', Pr_P_bound / (Pr_free + Pr_P_bound))
# Expression('piP_theory', P_free / (P_free + k6/k5))

# Initial cell concentrations
Parameter('R_0', 0.0007734)  # 0.0007734 pM (steady state), Initial conc of ROBs
Parameter('B_0', 0.0007282)  # 0.0007282 pM (steady state), Initial conc of AOBs
Parameter('C_0', 0.0009127)  # 0.0009127 pM (steady state), Initial conc of AOCs
Initial(R(), R_0)
Initial(B(), B_0)
Initial(C(), C_0)

# Parameters for P and Pr
Parameter('Sp', 250)  # pM/day
Parameter('Ip', 0)  # 1e4 # pM/day, range from 0 to 10^6
Parameter('kp', 86)  # 1/day
Parameter('k5', 0.02)  # 1/pM-day (rate of PTH binding with receptor)
Parameter('k6', 3)  # 1/day (rate of PTH unbinding from receptor)

##############################
Rtp = 1e6  # pmol/pmol cells, Number of PTH receptors per osteoblast (ROB + AOB). TODO: need to find this number
P_equil = (Sp.value + Ip.value) / kp.value
Pr_P_equil = P_equil * Rtp * (R_0.value + B_0.value) / (P_equil + k6.value / k5.value)
Pr_equil = Rtp * (R_0.value + B_0.value) - Pr_P_equil
##############################

Parameter('Pr_P_0', Pr_P_equil)
Parameter('Pr_0', Pr_equil)
Parameter('P_0', P_equil)
Initial(P(pr=None), P_0)
Initial(Pr(p=1) % P(pr=1), Pr_P_0)
Initial(Pr(p=None), Pr_0)

# Parameters for O
Expression('pi_P', Pr_P_bound / Pr_total)
Kop = 2e5  # pmol/day per pmol cells, Minimal rate of production of OPG per cell
Expression('So', Kop * ROB / pi_P)
Parameter('Io', 0)  # pM/day range from 0 to 10^6
Parameter('ko', 0.35)  # 1/day (Rate of elimination of OPG)
Parameter('k1', 1e-2)  # 1/pM-day (Rate of OPG-RANKL binding)
Parameter('k2', 10)  # 1/day (Rate of OPG-RANKL unbinding)

##############################
pi_P_equil = Pr_P_equil / (Pr_equil + Pr_P_equil)
O_equil = 1/ko.value * (Kop / pi_P_equil * R_0.value + Io.value)
So_equil = Kop * R_0.value / pi_P_equil
##############################

# Parameters for L
rL = 1e3  # pM/day, Rate of RANKL production and elimination
Klp = 3e6  # pmol/pmol cells, maximum number of RANKL attached to each AOC surface
Parameter('kl', 1e2)  # rate of degradation of RANKL TODO: need to think about this number
Expression('Sl', rL + kl * L_free)
Parameter('Il', 0)  # pM/day range from 0 to 10^6

# Parameters for K
Parameter('k3', 5.8e-4)  # 1/pM-day (Rate of RANK-RANKL binding)
Parameter('k4', 1.7e-2)  # 1/day (Rate of RANK-RANKL unbinding)

##############################
K_equil = 10  # pM, fixed conc of RANK
L_equil = Klp * pi_P_equil * B_0.value / (1 + k3.value / k4.value * K_equil + k1.value / k2.value * O_equil) * \
          (1 + Il.value / rL)
O_L_equil = k1.value / k2.value * O_equil * L_equil
K_L_equil = k3.value / k4.value * L_equil * K_equil
pi_L_equil = K_L_equil / K_equil
Sl_equil = rL + kl.value * L_equil
##############################

# Initial concentrations
Parameter('O_0', O_equil)
Parameter('O_L_0', O_L_equil)
Parameter('L_0', L_equil)
Parameter('K_0', K_equil)
Parameter('K_L_0', K_L_equil)
Initial(O(l=None), O_0)
Initial(O(l=1) % L(k=1), O_L_0)
Initial(L(k=None), L_0)
Initial(K(l=None), K_0)
Initial(K(l=1) % L(k=1), K_L_0)

# Parameters for cell dynamics
Expression('pi_L', K_L_bound / K_free)
f0 = 0.05  # unitless
Cs = 5e-3  # pM
Expression('pi_C', (AOC + f0*Cs) / (AOC + Cs))
DA = 0.7  # /day
dB = 0.7  # /day
DC = 2.1e-3  # pM/day
DR = 7e-4  # pM/day
Expression('DR_pi_C', DR*pi_C)
Expression('DB_over_pi_C', f0*dB/pi_C)
Expression('DC_pi_L', DC*pi_L)
Expression('DA_pi_C', DA*pi_C)
Parameter('kB', 0.189)  # /day

##############################
pi_C_equil = (C_0.value + f0 * Cs) / (C_0.value + Cs)
##############################

# Rules for P + Pr
Rule('P_creation_basal', None >> P(pr=None), Sp)
Rule('P_creation_external', None >> P(pr=None), Ip)
Rule('P_degradation', P(pr=None) >> None, kp)
Rule('P_binds_Pr', P(pr=None) + Pr(p=None) | P(pr=1) % Pr(p=1), k5, k6)

# Rules for O
Rule('O_creation_basal', None >> O(l=None), So)
Rule('O_creation_external', None >> O(l=None), Io)
Rule('O_degradation', O(l=None) >> None, ko)
Rule('O_bind_L', O(l=None) + L(k=None) | O(l=1) % L(k=1), k1, k2)

# Rules for L
Rule('L_creation_basal', None >> L(k=None), Sl)
Rule('L_creation_external', None >> L(k=None), Il)
Rule('L_degradation', L(k=None) >> None, kl)
Expression('k_L_cc', (Sl/L_free - kl) * L_total / (Klp * pi_P * AOB * L_free))
Rule('L_carrying_capacity', L(k=None) + L(k=None) >> L(k=None), k_L_cc)

##############################
k_L_cc_equil = (Sl_equil/L_equil - kl.value) * (L_equil + O_L_equil + K_L_equil) / \
               (Klp * pi_P_equil * B_0.value * L_equil)
##############################

#Rules for K
Rule('K_binds_to_L', K(l=None) + L(k=None) | K(l=1) % L(k=1), k3, k4)

# Cell dynamics rules
Rule('ROB_creation', None >> R(), DR_pi_C)
Rule('ROB_to_AOB', R() >> B(), DB_over_pi_C)
Rule('AOB_death', B() >> None, kB)
Rule('AOC_creation_death', None | C(), DC_pi_L, DA_pi_C)

# cell injections
Parameter('kf_AOB', 0)  # pM/day
Parameter('kf_AOC', 0)  # pM/day
Parameter('kf_ROB', 0)  # pM/day
Parameter('kr_AOB', 0)  # pM/day
Parameter('kr_AOC', 0)  # pM/day
Parameter('kr_ROB', 0)  # pM/day
Rule('AOB_injection', None | B(), kf_AOB, kr_AOB)
Rule('AOC_injection', None | C(), kf_AOC, kr_AOC)
Rule('ROB_injection', None | R(), kf_ROB, kr_ROB)

##############################
# equilibrium relations (ratios should equal 1)
expr_equil = [
    Expression('pP_div_dP', (Sp + Ip) / (kp * P_free)),
    Expression('pO_div_dO', (So + Io) / (ko * O_free)),
    Expression('pL_div_dL', (Sl + Il) / (kl * L_free + k_L_cc * L_free * L_free)),
    Expression('P_Pr_bind_unbind', k5 * P_free * Pr_free / (k6 * Pr_P_bound)),
    Expression('O_L_bind_unbind', k1 * O_free * L_free / (k2 * O_L_bound)),
    Expression('K_L_bind_unbind', k3 * K_free * L_free / (k4 * K_L_bound))
]
##############################

# simulations
tspan =np.linspace(0, 140 ,141)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

# cell populations
obs_cells = {'ROB': R_0, 'AOB': B_0, 'AOC': C_0}
for obs in [o for o in model.observables if o.name in obs_cells.keys()]:
    p = plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.plot(np.linspace(tspan[0], tspan[-1], 10), [obs_cells[obs.name].value]*10, 'o', mfc='none', color=p[0].get_color())
plt.xlabel('time (d)')
plt.ylabel('concentration (pM)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='best')
plt.tight_layout()

# molecule concentrations
fig, axs = plt.subplots(nrows=4, ncols=2, figsize=[6.4, 9.6])  # [6.4, 4.8]
row = 0
col = 0
molec_dict = {P_free: P_equil, Pr_free: Pr_equil, Pr_P_bound: Pr_P_equil, O_free: O_equil, L_free: L_equil,
              K_free: K_equil, O_L_bound: O_L_equil, K_L_bound: K_L_equil}
for obs in molec_dict.keys():
    p = axs[row, col].plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    axs[row, col].plot(np.linspace(tspan[0], tspan[-1], 10), [molec_dict[obs]]*10, 'o', mfc='none',
                       color=p[0].get_color())
    axs[row, col].legend(loc='best')
    if col == 0:
        axs[row, col].set_ylabel('concentration (pM)')
    if row == 3:
        axs[row, col].set_xlabel('time (d)')
    if (col + 1) % 2 == 0:
        row += 1
        col = 0
    else:
        col += 1
plt.tight_layout()

# expressions
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=[6.4, 7.2])  # [6.4, 4.8]
row = 0
col = 0
expr_dict = {'pi_P': pi_P_equil, 'pi_L': pi_L_equil, 'pi_C': pi_C_equil, 'So': So_equil, 'Sl': Sl_equil,
            'k_L_cc': k_L_cc_equil}
for expr in [e for e in model.expressions if e.name in ['pi_P', 'pi_L', 'pi_C', 'So', 'Sl', 'k_L_cc']]:
    p = axs[row, col].plot(tspan, output.expressions[expr.name], lw=2, label=expr.name)
    axs[row, col].plot(np.linspace(tspan[0], tspan[-1], 10), [expr_dict[expr.name]]*10, 'o', mfc='none',
                       color=p[0].get_color())
    axs[row, col].legend(loc='best')
    if col == 0:
        axs[row, col].set_ylabel('value')
    if row == 2:
        axs[row, col].set_xlabel('time (d)')
    if (col + 1) % 2 == 0:
        row += 1
        col = 0
    else:
        col += 1
plt.tight_layout()

# equilibrium ratios (should equal 1)
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=[6.4, 7.2])  # [6.4, 4.8]
row = 0
col = 0
for expr in expr_equil:
    p = axs[row, col].plot(tspan, output.expressions[expr.name], lw=2, label=expr.name)
    axs[row, col].plot(np.linspace(tspan[0], tspan[-1], 10), [1]*10, 'o', mfc='none', color = p[0].get_color())
    axs[row, col].legend(loc='best')
    axs[row, col].set_ylim(bottom=0.99, top=1.02)
    if col == 0:
        axs[row, col].set_ylabel('ratio')
    if row == 2:
        axs[row, col].set_xlabel('time (d)')
    if (col + 1) % 2 == 0:
        row += 1
        col = 0
    else:
        col += 1
plt.tight_layout()

# perturbations
tspan1 = np.linspace(0, 20, 201)
tspan2 = np.linspace(20, 80, 601)
tspan3 = np.linspace(80, 140, 601)

perturb = [
    [{'kf_AOB': 0}, {'kf_AOB': 1e-4}, {'kf_AOB': 0}],
    [{'kf_AOC': 0}, {'kf_AOC': 1e-4}, {'kf_AOC': 0}],
    [{'kf_ROB': 0}, {'kf_ROB': 1e-4}, {'kf_ROB': 0}],
    [{'kr_AOB': 0}, {'kr_AOB': 8.3e-2}, {'kr_AOB': 0}],  # something wrong here
    [{'kr_AOC': 0}, {'kr_AOC': 2.9e-1}, {'kr_AOC': 0}],  # something wrong here
    [{'kr_ROB': 0}, {'kr_ROB': 1.2e-1}, {'kr_ROB': 0}],  # something wrong here
    [{'Ip': 0}, {'Ip': 1e3}, {'Ip': 0}],
    [{'Io': 0}, {'Io': 2e5}, {'Io': 0}],
    [{'Il': 0, 'Io': 0}, {'Il': 1e4, 'Io': 0}, {'Il': 1e4, 'Io': 9e4}]  # something wrong with IL here
]

fig, axs = plt.subplots(nrows=3, ncols=3, sharex='all', figsize=(10,8))
colors = ['b', 'r', 'g']
row = 0
col = 0
for p in perturb:
    print(p)
    output1 = sim.run(param_values=p[0], tspan=tspan1)
    output2 = sim.run(param_values=p[1], tspan=tspan2, initials=output1.species[-1])
    output3 = sim.run(param_values=p[2], tspan=tspan3, initials=output2.species[-1])
    for i, obs in enumerate(['ROB', 'AOB', 'AOC']):
        axs[row, col].plot(tspan1, output1.observables[obs], lw=2, color=colors[i], label=obs)
        axs[row, col].plot(tspan2, output2.observables[obs], lw=2, color=colors[i])
        axs[row, col].plot(tspan3, output3.observables[obs], lw=2, color=colors[i])
    if row == 2:
        axs[row, col].set_xlabel('time (d)')
    if col == 0:
        axs[row, col].set_ylabel('concentration (pM)')
    axs[row, col].set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
    axs[row, col].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axs[row, col].legend(loc=0)
    if (row+1) % 3 == 0:
        col += 1
        row = 0
    else:
        row += 1
plt.tight_layout()

plt.show()
