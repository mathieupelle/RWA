# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:39:42 2021

@author: Mathieu Pellé
"""

from LL_Utilities import VortexGeometry, LiftingLine
from BEMT_Utilities import BEMT_execute
from LL_plotting import plot_radial, performance_coefs

#%% Define common parameters between BEM and LL
N_radial = 25
spacing = 'cos'
TSR = 8
U_inf = 10

#%% Call for BEM geometry and results

rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#%% Lifting line


# Lists of: rotor geometries, No. of Radial Elements, Spacing Method, No. of Rotations, Rotor positions
geometry = VortexGeometry([rotor_opt], [2], [0], result_BEM)

#List of rotors and combined geometry
results_LL = LiftingLine([rotor_opt],geometry,result_BEM)


#%% PLotting and results

plot_radial(results_LL, result_BEM, rotor_opt)
performance_coefs(results_LL, result_BEM, rotor_opt)

#%%





