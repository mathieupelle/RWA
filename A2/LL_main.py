# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:39:42 2021

@author: Mathieu Pellé
"""

from LL_Utilities import VortexGeometry, LiftingLine
from BEMT_Utilities import BEMT_execute
from LL_plotting import plot_radial, performance_coefs, plot_radial_2R, performance_coefs_2R


#%% Single rotor case

#Define common parameters between BEM and LL
N_radial = 25
spacing = 'cos'
TSR = 8
U_inf = 10

#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Call for Lifting Line geometry and results

# Lists of: rotor geometries, No. of Radial Elements, No. of Rotations, Rotor positions
geometry = VortexGeometry([rotor_opt], [2], [0], result_BEM)

#List of rotors, combined geometry and BEM results
results_LL = LiftingLine([rotor_opt],geometry,result_BEM)

# PLotting and results
plot_radial(results_LL, result_BEM, rotor_opt)
performance_coefs(results_LL, result_BEM, rotor_opt)

#%% Double rotor case

#Define common parameters for both rotors and for both BEM and LL
N_radial = 25
spacing = 'cos'
TSR = 8
U_inf = 10

#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

# Lists of: rotor geometries, No. of Radial Elements, No. of Rotations, Rotor positions
geometry = VortexGeometry([rotor_opt, rotor_opt], [2, 2], [0, 50*2], result_BEM)

#List of rotors, combined geometry and BEM results
results_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)

# PLotting and results
plot_radial_2R(results_LL_2R, results_LL, rotor_opt, [[1,2,3],[1,2,3]])
performance_coefs_2R(results_LL_2R, results_LL, result_BEM, rotor_opt)
















