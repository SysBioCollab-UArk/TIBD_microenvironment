from pysb import *
from pysb.util import alias_model_components


def add_bisphosphonate_components():
    Monomer('Bisphos')
    Parameter('Bisphos_0', 0)  # unitless (ratio of dose to baseline dose = 0.1 mg/mouse)
    Parameter('k_bisphos_AOC', 1)  # 1/day (ZA is unitless, not in fM)
    alias_model_components()
    Initial(Bisphos(), Bisphos_0)
    Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)
