from pysb import *
from pysb.util import alias_model_components
from sympy import Piecewise


def add_bisphosphonate_components():
    Monomer('Bisphos')
    Parameter('Bisphos_0', 0)  # unitless (ratio of dose to baseline dose = 0.1 mg/mouse)
    Parameter('k_bisphos_AOC', 1)  # 1/day (ZA is unitless, not in fM)
    alias_model_components()
    Initial(Bisphos(), Bisphos_0)
    Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)


def add_RANKLi_components():
    Monomer('RANKLi')
    Parameter('RANKLi_0', 0)  # fM
    # binding and unbinding rate constants for RANKLi + RANKL <-> RANKLi % RANKL
    Parameter('kact_IL', 1)  # kf_IL/kr_IL (1/fM)
    alias_model_components()
    Initial(RANKLi(), RANKLi_0)
    Observable('RANKLi_tot', RANKLi())
    alias_model_components()
    # DC_pi_L.expr = DC * k3/k4 * KLP*((IP/kP + SP/kP) / (IP/kP + k6/k5))*B_obs * (1 + IL/rL) / (1 + k3*K/k4 + k1/k2/kO*(KOP/((IP/kP + SP/kP) / (IP/kP + k6/k5))*R_obs + IO) + kf_IL/kr_IL*RANKLi_tot)  # fM/day
    # Ltot_over_L.expr = 1 + k3*K/k4 + k1/k2/kO*(KOP/((IP/kP + SP/kP) / (IP/kP + k6/k5))*R_obs + IO) + kf_IL/kr_IL*RANKLi_tot
    L_denom.expr += kact_IL * RANKLi_tot


def add_tumor_RANK_components():

    Monomer('TumorK', ['state'], {'state': ['u', 'L']})
    Parameter('TumorK_0', 0)
    Parameter('k_T_adapt', 0.05)  # adaptation rate constant for Tumor -> TumorK cells
    # Rate constants for RANKL binding/unbinding to TumorK cells
    Parameter('kf_TK_L', 0.1)
    Parameter('kr_TK_L', 1)
    alias_model_components()

    Initial(TumorK(state='u'), TumorK_0)
    Observable('TumorK_tot', TumorK())
    Observable('TumorK_u', TumorK(state='u'))
    Observable('TumorK_L', TumorK(state='L'))
    alias_model_components()

    # Redefine 'Tumor_tot' observable
    Tumor_tot.reaction_pattern += TumorK()
    # Adaptation of tumor cells to a RANK expressesing variant
    Rule('Tumor_to_TumorK', Tumor() >> TumorK(state='u'), k_T_adapt)
    # Reversible binding of RANKL to TumorK
    L_numer.expr -= TumorK_L
    Expression('kf_TK_L_times_L', kf_TK_L * L_numer / L_denom)
    alias_model_components()
    Rule('TumorK_binds_RANKL', TumorK(state='u') | TumorK(state='L'), kf_TK_L_times_L, kr_TK_L)

    # 1. Tumor cell division
    Parameter('k_tumorKL_div_basal', k_tumor_div_basal.value)  # initialize with value from unadapted tumor cells
    Parameter('k_tumorKL_div_TGFb', k_tumor_div_TGFb.value)  # initialize with value from unadapted tumor cells
    alias_model_components()
    Expression('k_tumorKL_div', k_tumorKL_div_basal + k_tumorKL_div_TGFb * pi_C)
    alias_model_components()
    Rule('TumorKu_div', TumorK(state='u') >> TumorK(state='u') + TumorK(state='u'), k_tumor_div)
    Rule('TumorKL_div', TumorK(state='L') >> TumorK(state='L') + TumorK(state='L'), k_tumorKL_div)

    # 2. Tumor cell death
    Parameter('k_tumorKL_dth', k_tumor_dth.value)  # initialize with value from unadapted tumor cells
    alias_model_components()
    Rule('TumorKu_death', TumorK(state='u') >> None, k_tumor_dth)
    Rule('TumorKL_death', TumorK(state='L') >> None, k_tumorKL_dth)

    # 3. Density-dependent cell death (carrying capacity)
    Expression('k_tumorKu_cc', Piecewise((0, TumorK_u <= 0),
                                         (CC_ON * (k_tumor_div - k_tumor_dth) / N * Tumor_tot / TumorK_u, True)))
    Expression('k_tumorKL_cc', Piecewise((0, TumorK_L <= 0),
                                         (CC_ON * (k_tumorKL_div - k_tumorKL_dth) / N * Tumor_tot / TumorK_L, True)))
    alias_model_components()
    Rule('TumorKu_cc', TumorK(state='u') + TumorK(state='u') >> TumorK(state='u'), k_tumorKu_cc)
    Rule('TumorKL_cc', TumorK(state='L') + TumorK(state='L') >> TumorK(state='L'), k_tumorKL_cc)

    # 4. PTHrP production by tumor cells
    Parameter('k_tumorKL_PTHrP', k_tumor_PTHrP.value)  # initialize with value from unadapted tumor cells
    alias_model_components()
    IP.expr += (k_tumor_PTHrP * TumorK_u + k_tumorKL_PTHrP * TumorK_L) * pi_C

    # 5. Tumor cells promote OB death
    Parameter('k_tumorKL_OB', k_tumor_OB.value)  # initialize with value from unadapted tumor cells
    alias_model_components()
    Rule('AOB_death_TumorKu', TumorK(state='u') + B() >> TumorK(state='u'), k_tumor_OB)
    Rule('AOB_death_TumorKL', TumorK(state='L') + B() >> TumorK(state='L'), k_tumorKL_OB)

    # 6. Allee effect
    Expression('k_tumorKL_allee', Piecewise(
        (0, Tumor_tot <= 0), (ALLEE_ON * (k_tumorKL_div - k_tumorKL_dth) / (A * Tumor_tot), True)))
    alias_model_components()
    Rule('TumorKu_allee', TumorK(state='u') >> None, k_tumor_allee)
    Rule('TumorKL_allee', TumorK(state='L') >> None, k_tumorKL_allee)
