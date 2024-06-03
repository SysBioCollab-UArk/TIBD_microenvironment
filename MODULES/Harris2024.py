from pysb import *
from pysb.util import alias_model_components
from sympy import Piecewise, sympify


def create_model_elements(OB_OC_BONE_MODEL=1):

    # ########## BONE ##########

    Monomer('Bone')
    Parameter('Bone_0', 100)  # percentage

    # Make OB bone production and OC bone consumption rates increase with decreasing cell concentrations
    # B >> B + Bone, rate_B
    # C + Bone >> C, rate_C
    # At equilibrium: rate_B * [B] = rate_C * [C] * [Bone]

    # parameters
    Parameter('k1_B_builds_bone')  # /day
    Parameter('k2_B_builds_bone')  # fM
    Parameter('nB')  # unitless
    Parameter('k1_C_consumes_bone')  # /day
    Parameter('k2_C_consumes_bone')  # fM
    Parameter('nC')  # unitless

    alias_model_components()

    if OB_OC_BONE_MODEL == 1:

        # === MODEL #1 ===

        Expression('k_B_builds_bone',
                   Piecewise((0, B_obs <= 0),
                             (k1_B_builds_bone * (k2_B_builds_bone ** nB + B_obs ** nB) / B_obs ** nB, True)))
        Expression('k_C_consumes_bone',
                   Piecewise((0, C_obs <= 0),
                             (k1_C_consumes_bone * (k2_C_consumes_bone ** nC + C_obs ** nC) / C_obs ** nC, True)))

        # some preliminary values so the model will run
        k1_B_builds_bone.value = 3359287  # calculated to get [Bone]_equil = 100
        k2_B_builds_bone.value = 100
        nB.value = 1
        k1_C_consumes_bone.value = 1e6
        k2_C_consumes_bone.value = 1
        nC.value = 1

    elif OB_OC_BONE_MODEL == 2:

        # === MODEL #2 ===

        Parameter('kB_star', 5e9)  # /day 1.4e9
        Parameter('kC_star', 2e6)  # /day 2e6

        alias_model_components()

        Expression('k_B_builds_bone',
                   k1_B_builds_bone + kB_star * k2_B_builds_bone ** nB / (k2_B_builds_bone ** nB + B_obs ** nB))
        Expression('k_C_consumes_bone',
                   k1_C_consumes_bone + kC_star * k2_C_consumes_bone ** nC / (k2_C_consumes_bone ** nC + C_obs ** nC))

        # some preliminary values so the model will run
        k1_B_builds_bone.value = 37207883  # calculated to get [Bone]_equil = 100
        k2_B_builds_bone.value = 0.5
        nB.value = 1.8
        k1_C_consumes_bone.value = 2.2e6
        k2_C_consumes_bone.value = 0.1
        nC.value = 0.1

    alias_model_components()

    Initial(Bone(), Bone_0)
    Observable('Bone_tot', Bone())

    Rule('B_builds_bone', B() >> B() + Bone(), k_B_builds_bone)
    Rule('C_consumes_bone', C() + Bone() >> C(), k_C_consumes_bone)

    # ######### TUMOR #########

    Monomer('Tumor')
    Parameter('Tumor_0', 0)  # fM (1000 fM = 1 pM)
    alias_model_components()
    Initial(Tumor(), Tumor_0)

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

    # rate expressions
    Expression('k_tumor_cc', CC_ON * (k_tumor_div - k_tumor_dth) / N)
    Expression('k_tumor_allee',
               Piecewise((0, Tumor_tot <= 0),
                         (ALLEE_ON * (k_tumor_div - k_tumor_dth) / (A * Tumor_tot), True)))

    alias_model_components()

    Rule('tumor_division', Tumor() >> Tumor() + Tumor(), k_tumor_div)
    Rule('tumor_death', Tumor() >> None, k_tumor_dth)
    Rule('tumor_cc', Tumor() + Tumor() >> Tumor(), k_tumor_cc)
    Rule('tumor_allee', Tumor() >> None, k_tumor_allee)

    # Make PTH (i.e., PTHrP) synthesis rate a function of the number of tumor cells
    Parameter('k_tumor_PTHrP', 5e2)  # 1/day
    alias_model_components()
    # PTHrP expression in tumor cells depends on TGF-beta occupancy pi_C
    IP.expr = sympify(IP_const + k_tumor_PTHrP * pi_C * Tumor_tot)

    # Hypothesize that tumor cells secrete a factor that promotes OB death
    Parameter('k_tumor_OB', 0.002)  # 1/fM-day  # 0.01
    alias_model_components()
    Rule('AOB_death_tumor', Tumor() + B() >> Tumor(), k_tumor_OB)

    # Hypothesize that tumor cells secrete a factor that promotes OC production
    Parameter('k_tumor_OC', 0.004)  # 1/fM-day  # 0.01
    alias_model_components()
    Rule('AOC_creation_tumor', Tumor() >> Tumor() + C(), k_tumor_OC)
