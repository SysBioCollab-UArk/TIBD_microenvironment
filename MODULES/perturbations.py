from pysb import *
from pysb.util import alias_model_components


def add_bisphosphonate_components():
    Monomer('Bisphos')
    Parameter('Bisphos_0', 0)  # unitless (ratio of dose to baseline dose = 0.1 mg/mouse)
    Parameter('k_bisphos_AOC', 1)  # 1/day (ZA is unitless, not in fM)
    alias_model_components()
    Initial(Bisphos(), Bisphos_0)
    Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)

def add_denosumab_components():
    Monomer('Denosumab')
    Parameter('Denosumab_0', 100)  # fM
    Parameter('kf_DL', 1)  # 1/fM-s
    Parameter('kr_DL', 1)  # 1/s
    alias_model_components()
    Initial(Denosumab(), Denosumab_0)
    Observable('Denosumab_tot', Denosumab())
    alias_model_components()
    DC_pi_L.expr = DC * k3/k4 * KLP*((IP/kP + SP/kP) / (IP/kP + k6/k5))*B_obs * (1 + IL/rL) / (1 + k3*K/k4 + k1/k2/kO*(KOP/((IP/kP + SP/kP) / (IP/kP + k6/k5))*R_obs + IO) + kf_DL/kr_DL*Denosumab_tot)  # fM/day

