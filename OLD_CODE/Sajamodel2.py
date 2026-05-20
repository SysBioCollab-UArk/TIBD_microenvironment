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
Monomer('PTH', ['ob'])  # (PTH) PTH monomer
Monomer('ROB', ['pth','tgfb'])  # (ROB) Responding osteoblast monomer
Monomer('AOB', ['pth','state'], {'state': ['X','L']})  # (AOB) Active osteoblast monomer
Monomer('TGFB',['ob'])  # (TGF-B) TGF-Beta monomer
Monomer('AOC') # Active Osteoclast monomer
Monomer('OCP') # Osteoclast Precursors



Parameter('UP_init', 0)
Initial(UP(), UP_init)

Parameter('OPG_init', 0)
Initial(OPG(), OPG_init)

Parameter('PTH_init', 0)
Initial(PTH(ob=None), PTH_init)

Parameter('ROB_init', 0)
Initial(ROB(pth=None,tgfb=None), ROB_init)

Parameter('AOB_init', 0)
Initial(AOB(pth=None,state='X'), AOB_init)

Parameter('TGFB_init', 0)
Initial(TGFB(ob=None), TGFB_init)

Parameter('AOC_init', 0)
Initial(AOC(), AOC_init)

Parameter('kf_UP_ROB', 1)
Rule('ROB_production', UP() + TGFB() >> ROB() + TGFB(), kf_UP_ROB)


Parameter('kf_ROB_PTH',1)
Parameter('kr_ROB_PTH',1)
Rule('ROB_BIND_PTH', ROB(pth=None) + PTH(ob=None) | ROB(pth=1) % PTH(ob=1) ,kf_ROB_PTH, kr_ROB_PTH )

Parameter('kf_AOB_PTH', 1)
Parameter('kr_AOB_PTH', 1)
Rule('AOB_Bind_PTH', AOB(pth=None) + PTH(ob=None) | AOB(pth=1) % PTH(ob=1), kf_AOB_PTH, kr_AOB_PTH)

Parameter('k_AOB_L',1)
Rule('AOB_EXPRESS_L', AOB(pth=1,state='X') % PTH(ob=1) >> AOB(pth=1, state='L') % PTH(ob=1), k_AOB_L)

Parameter('k_AOB_L_X',1)
Rule('AOB_L_TO_X', AOB(pth=None,state='L') >> AOB(pth=None,state='X'), k_AOB_L_X)

Parameter('k_ROB_AOB',1)
Rule('AOB_Production', ROB(pth=None,tgfb=None) >> AOB(pth=None,state='X'), k_ROB_AOB)

Parameter('k_ROB_AOB_pth',1)
Rule('AOB_Production_pth', ROB(pth=ANY,tgfb=None) >> AOB(pth=1,state='X') % PTH(ob=1), k_ROB_AOB_pth)

Parameter('kf_ROB_TGFB',1)
Parameter('kr_ROB_TGFB',1)
Rule('ROB_BIND_TGFB', ROB(tgfb=None) + TGFB(ob=None) | ROB(tgfb=1) % TGFB(ob=1) ,kf_ROB_TGFB, kr_ROB_TGFB )

Parameter('k_AOB_death',1)
Rule('AOB_death', AOB() >> None, k_AOB_death)

Parameter('k_PTH_deg',1)
Rule('PTH_deg', PTH() >> None, k_PTH_deg)

Parameter('k_OPG_deg',1)
Rule('OPG_deg', OPG() >> None, k_OPG_deg)

Parameter('kf_AOB_L_OPG',1)
Parameter('kr_AOB_L_OPG',1)
Rule('AOB_EXPRESS_L_OPG', AOB(state='X') + OPG() | AOB(state='L') + OPG(), kf_AOB_L_OPG, kr_AOB_L_OPG)

Parameter('kf_AOB_L_PTH_OPG',1)
Parameter('kr_AOB_L_PTH_OPG',1)
Rule('AOB_EXPRESS_L_OPG_PTH', AOB(state='X') + OPG()+ PTH() | AOB(state='L') + OPG() +PTH(), kf_AOB_L_PTH_OPG, kr_AOB_L_PTH_OPG)

#Parameter('',1)
#Rule(''
