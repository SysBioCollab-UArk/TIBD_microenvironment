from pysb import *
from sympy import Piecewise, sympify
from pysb.util import alias_model_components


def create_model_elements(OB_OC_model=1):

    # ########## Bone ##########
    Monomer('Bone')
    Parameter('Bone_0', 100)  # percentage

    # Make OB bone production (OC bone consumption) rate increase (decrease) with decreasing (increasing) concentration
    mean_OB = 2.06E-02 * 1000  # 1.06E-02 * 1000  # fM, from data
    mean_OC = 2.19E-03 * 1000  # 9.16E-04 * 1000  # fM, from data

    if OB_OC_model == 1:

        # === MODEL #1 ===
        Parameter('k1_B_builds_bone', 0)  # /day 1.25e8
        Parameter('k2_B_builds_bone', 100)  # fM
        Parameter('nB', 1)  # 1.1
        Parameter('k1_C_consumes_bone', 1e6)  # /day
        Parameter('k2_C_consumes_bone', 1)  # fM
        Parameter('nC', 1)  # 0.8
        alias_model_components()
        Expression('k_B_builds_bone',
                   Piecewise((0, B_obs <= 0),
                             (k1_B_builds_bone * (k2_B_builds_bone ** nB + B_obs ** nB) / B_obs ** nB, True)))
        Expression('k_C_consumes_bone',
                   Piecewise((0, C_obs <= 0),
                             (k1_C_consumes_bone * (k2_C_consumes_bone ** nC + C_obs ** nC) / C_obs ** nC, True)))

        # estimate rate constant for OB -> OB + Bone to get bone density ~ 100
        k1C = k1_C_consumes_bone.value
        k2C = k2_C_consumes_bone.value
        k2B = k2_B_builds_bone.value
        k1_B_builds_bone.value = k1C * (k2C ** nC.value + mean_OC ** nC.value) / mean_OC ** (nC.value-1) * \
                                 Bone_0.value * mean_OB ** (nB.value - 1) / (k2B ** nB.value + mean_OB ** nB.value) * \
                                 1.27

    elif OB_OC_model == 2:

        # === MODEL #2 ===
        Parameter('k1_B_builds_bone', 2.3e6)  # /day
        Parameter('k2_B_builds_bone', 0.4)  # fM 0.8
        Parameter('nB', 1.7)  # 2
        Parameter('kB_star', 1.4e9)  # /day
        Parameter('k1_C_consumes_bone', 1e6)  # /day
        Parameter('k2_C_consumes_bone', 1)  # fM
        Parameter('nC', 1)  # 0.8
        Parameter('kC_star', 2e6)  # /day
        alias_model_components()
        Expression('k_B_builds_bone',
                   k1_B_builds_bone + kB_star * k2_B_builds_bone ** nB / (k2_B_builds_bone ** nB + B_obs ** nB))
        Expression('k_C_consumes_bone',
                   k1_C_consumes_bone + kC_star * k2_C_consumes_bone ** nC / (k2_C_consumes_bone ** nC + C_obs ** nC))

        # estimate rate constant for OB -> OB + Bone to get bone density ~ 100
        k1C = k1_C_consumes_bone.value
        k2C = k2_C_consumes_bone.value
        k2B = k2_B_builds_bone.value
        k1_B_builds_bone.value = ((k1C + kC_star.value * k2C ** nC.value / (k2C ** nC.value + mean_OC ** nC.value)) * \
                                  mean_OC * Bone_0.value / mean_OB - kB_star.value * k2B ** nB.value / \
                                  (k2B ** nB.value + mean_OB ** nB.value)) * \
                                  1.93

    alias_model_components()
    Initial(Bone(), Bone_0)
    Observable('Bone_tot', Bone())
    Rule('B_builds_bone', B() >> B() + Bone(), k_B_builds_bone)
    Rule('C_consumes_bone', C() + Bone() >> C(), k_C_consumes_bone)

    # ######### Tumor #########
    Monomer('Tumor')
    Parameter('Tumor_0', 0)  # fM (1000 fM = 1 pM)

    # standard exponential growth
    Parameter('k_tumor_div_basal', 1)  # 1/day
    Parameter('k_tumor_dth', 0.6)  # 1/day

    # carrying capacity (logistic growth)
    Parameter('CC_ON', 1)
    Parameter('N', 900)  # fM 600

    # allee effect
    Parameter('ALLEE_ON', 0)
    Parameter('A', 1.5)  # 1/fM

    # growth enhancement due to TGF-beta
    Parameter('k_tumor_div_TGFb', 0.1)  # /day
    alias_model_components()
    Expression('k_tumor_div', k_tumor_div_basal + k_tumor_div_TGFb * pi_C)

    Observable('Tumor_tot', Tumor())
    alias_model_components()

    Expression('k_tumor_cc', CC_ON * (k_tumor_div - k_tumor_dth) / N)
    Expression('k_tumor_allee',
               Piecewise((0, Tumor_tot <= 0),
                         (ALLEE_ON * (k_tumor_div - k_tumor_dth) / (A * Tumor_tot), True)))
    alias_model_components()

    Initial(Tumor(), Tumor_0)
    Rule('tumor_division', Tumor() >> Tumor() + Tumor(), k_tumor_div)
    Rule('tumor_death', Tumor() >> None, k_tumor_dth)
    Rule('tumor_cc', Tumor() + Tumor() >> Tumor(), k_tumor_cc)
    Rule('tumor_allee', Tumor() >> None, k_tumor_allee)

    # Make PTH (i.e., PTHrP) synthesis rate a function of the number of tumor cells
    Parameter('k_tumor_PTHrP', 5e2)  # 1/day
    alias_model_components()
    IP.expr = sympify(IP_const + k_tumor_PTHrP * pi_C * Tumor_tot)  # PTHrP expression depends on TGF-beta occupancy pi_C

    # Hypothesize that tumor cells secrete a factor that promotes OB death
    Parameter('k_tumor_OB', 0.002)  # 1/fM-day  # 0.01
    alias_model_components()
    Rule('AOB_death_tumor', Tumor() + B() >> Tumor(), k_tumor_OB)

    # Hypothesize that tumor cells secrete a factor that promotes OC production
    Parameter('k_tumor_OC', 0.004)  # 1/fM-day  # 0.01
    alias_model_components()
    Rule('AOC_creation_tumor', Tumor() >> Tumor() + C(), k_tumor_OC)

    # ######### BISPHOSPHONATES #########

    Monomer('Bisphos')
    Parameter('Bisphos_0', 0)  # fM
    Parameter('k_bisphos_AOC', 1)  # 1/fM-day
    alias_model_components()
    Initial(Bisphos(), Bisphos_0)
    Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)
