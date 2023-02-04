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

# Monomers
Monomer('P'['pr'],)
Monomer('Pr',['p'])
Monomer('O',['l'])
Monomer('L',['o'],['k'])
Monomer('K',['l'])

Parameters

#Rules
Rule('P_creation',P(None)>>P(1),dp)
Rule('P_destruction',P(1) >> P(None))
Rule('P_binds_Pp',P([p=None]) + Pr([p=None]) | P([pr=1]), Pr([p=1]),k5,k6)

# Rules for O + L
Rule('O_creation',O(None) >> O(1),po)
Rule('O_destruction', o(1) >> O(None),do)
Rule('O_bind_L',O(l=None) + L(o=None) | O(l=1) % L(o=1),k1,k2)
Rule('L_creation',L(None) >> l(1),pl)
Rule('L_destruction',L(1) >>L(None),dl)

#Rules for L + K
# where does K comes from?
Rule('K_binds to L',K(l=None) +L(k=None) | K(l=1) +L(k=1),k3,k4)

