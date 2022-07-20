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

Parameter('kf_UP_ROB', 1)
Rule('ROB_production', UP() + TGFB() >> ROB() + TGFB(), kf_UP_ROB)


Parameter('kf_ROB_PTH',1)
Parameter('kr_ROB_PTH',1)
Rule('ROB_BIND_PTH', ROB(pth=None) + PTH(rob=None) | ROB(pth=1) % PTH(rob=1) ,kf_ROB_PTH, kr_ROB_PTH )


Parameter('kf_AOB_PTH', 1)
Rule('AOB_Bind_PTH', AOB() + PTH() | AOB() + PTH(), kf_AOB_PTH)

Parameter('kf_ROB_AOB',1)
Rule('AOB_Production', ROB() + TGFB() >> AOB() + TGFB(), kf_ROB_AOB)

Parameter('kf_AOB_L',1)
Rule('AOB_EXPRESS_L', AOB(state='X') + PTH() | AOB(state='L') + PTH(), kf_AOB_L)

Parameter('k_AOB_death',1)
Rule('AOB_death', AOB() >> None, k_AOB_death)

Parameter('k_PTH_deg',1)
Rule('PTH_deg', PTH() >> None, k_PTH_deg)

Parameter('k_OPG_deg',1)
Rule('OPG_deg', OPG() >> None, k_OPG_deg)

Parameter('kf_AOB_L_OPG',1)
Rule('AOB_EXPRESS_L', AOB(state='X') + OPG() | AOB(state='L') + OPG(), kf_AOB_L_OPG)

Parameter('kf_AOB_L_PTH_OPG',1)
Rule('AOB_EXPRESS_L', AOB(state='X') + OPG()+ PTH() | AOB(state='L') + OPG() +PTH(), kf_AOB_L_PTH_OPG)


