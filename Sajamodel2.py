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

Monomer('UP')  # uncommitted Progenitors
Monomer('OPG') # Osteoprotegerin
Monomer('PTH', ['rob'])  # (PTH) PTH monomer
Monomer('ROB', ['pth'])  # (ROB) Responding osteoblast monomer
Monomer('AOB', ['state'], {'state': ['X','L']})  # (AOB) Active osteoblast monomer
Monomer('TGFB')  # (TGF-B) TGF-Beta monomer
Monomer('AOC') # Active Osteoclast monomer
Monomer('OCP') # Osteoclast Precursors



Parameter('UP_init', 0)
Initial(UP(), UP_init)

Parameter('OPG_init', 0)
Initial(OPG(), OPG_init)

Parameter('PTH_init', 0)
Initial(PTH(rob=None), PTH_init)

Parameter('ROB_init', 0)
Initial(ROB(pth=None), ROB_init)

Parameter('AOB_init', 0)
Initial(AOB(state='X'), AOB_init)

Parameter('TGFB_init', 0)
Initial(TGFB(), TGFB_init)

Parameter('AOC_init', 0)
Initial(AOC(), AOC_init)

Parameter('krobp_ob', 1)
Rule('ROB_production', UP() + TGFB() >> ROB() + TGFB(), krobp_ob)


Parameter('kf_ROB_PTH',1)
Parameter('kr_ROB_PTH',1)
Rule('ROB_BIND_PTH', ROB(pth=None) + PTH(rob=None) | ROB(pth=1) % PTH(rob=1) ,kf_ROB_PTH, kr_ROB_PTH )


Parameter('kaob_pth', 1)
Rule('AOB_Bind_PTH', AOB() + PTH() >> AOB()+ PTH(), kaob_pth)
