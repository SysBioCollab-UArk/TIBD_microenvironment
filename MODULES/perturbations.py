from pysb import *
from pysb.util import alias_model_components


def add_bisphosphonate_components():
    Monomer('Bisphos')
    Parameter('Bisphos_0', 0)  # fM
    Parameter('k_bisphos_AOC', 1)  # 1/fM-day
    alias_model_components()
    Initial(Bisphos(), Bisphos_0)
    Rule('AOC_death_bisphos', Bisphos() + C() >> Bisphos(), k_bisphos_AOC)
