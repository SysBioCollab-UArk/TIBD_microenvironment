from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot

Model()

# Monomers
# Parameters
# Initials
# Rules
# Observables
# simulation commands


# (T) tumor cell monomer
# (B) osteoblast monomer
# (C) osteoclast monomer
# (P) PTHrP monomer
# (Beta) TGF-Beta monomer


Monomer('T')
Monomer('P')
Monomer('B')
Monomer('C')
Monomer('Beta')

Parameter('T_init',100)
initial(T(T >> T + P) , T_init)


parameter('ktp_tumor', 1)
Rule('tumor' , T >> T + P , ktp_tumor )



#Rule('P_binds_B',