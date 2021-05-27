# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:39:42 2021

@author: Mathieu Pellé
"""

from LL_Utilities import VortexGeometry, LiftingLine
from BEMT_Utilities import BEMT_execute
from LL_plotting import plot_radial, performance_coefs, plot_radial_2R, performance_coefs_2R, plot_difference_2R
import numpy as np

#%% Single rotor case

#Define common parameters between BEM and LL
N_radial = 35
spacing = 'cos'
U_inf = 10
TSR = [8]
results_LL = []
results_BEM = []
rotors = []
for i in range(len(TSR)):
    #Call for BEM geometry and results
    rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR[i])
    rotors.append(rotor_opt)
    results_BEM.append(result_BEM)

    #Call for Lifting Line geometry and results

    # Lists of: rotor geometries, No. of Radial Elements, No. of Rotations, Rotor positions
    geometry = VortexGeometry([rotor_opt], [5], [0], result_BEM)

    # List of rotors, combined geometry and BEM results
    result_LL = LiftingLine([rotor_opt],geometry,result_BEM)
    results_LL.append(result_LL)

    # Coefficients
    print('=====TSR '+str(TSR[i])+'=====')
    performance_coefs(result_LL, result_BEM, rotor_opt)


# PLots
plot_radial(results_LL, results_BEM, rotors, TSR)


#%% Double rotor case

#Define common parameters for both rotors and for both BEM and LL
N_radial = 35
spacing = 'cos'
TSR = 8
U_inf = 10

#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Single rotor
geometry = VortexGeometry([rotor_opt], [5], [0], result_BEM)
result_LL = LiftingLine([rotor_opt],geometry,result_BEM)

#Double rotor
geometry = VortexGeometry([rotor_opt, rotor_opt], [2, 2], [0, 50*2], result_BEM)
results_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)

# PLotting and results
plot_radial_2R(results_LL_2R, result_LL, rotor_opt, [[1,2,3],[]])
performance_coefs_2R(results_LL_2R, result_LL, result_BEM, rotor_opt)



#%% Double rotor effect of phase difference


N_radial = 30
spacing = 'cos'
TSR = 8
U_inf = 10

#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Single rotor
geometry = VortexGeometry([rotor_opt], [3], [0], result_BEM)
result_LL = LiftingLine([rotor_opt],geometry,result_BEM)

phase = np.linspace(0,120,11)
CT_lst = [[], [], []]
CP_lst = [[], [], []]
results_LL_2R = []
for i in range(len(phase)):
    print('Computing for phase difference of ' + str(phase[i]))

    geometry = VortexGeometry([rotor_opt, rotor_opt], [3, 3], [0, 50*2], result_BEM, phase[i])
    result_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)
    results_LL_2R.append(result_LL_2R)

    [CT, CP] = performance_coefs_2R(result_LL_2R, result_LL, result_BEM, rotor_opt)
    CT_lst[0].append(CT[0])
    CT_lst[1].append(CT[1])
    CT_lst[2].append(CT[2])
    CP_lst[0].append(CP[0])
    CP_lst[1].append(CP[1])
    CP_lst[2].append(CP[2])

plot_difference_2R(phase, CT_lst, CP_lst, 'phase')

#%% Double rotor effect of rotor distance

N_radial = 35
spacing = 'cos'
TSR = 8
U_inf = 10
D = 2*50
#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Single rotor
geometry = VortexGeometry([rotor_opt], [3], [0], result_BEM)
result_LL = LiftingLine([rotor_opt],geometry,result_BEM)

D_lst = [D, 2*D, 5*D, 10*D]
CT_lst = [[], [], []]
CP_lst = [[], [], []]
results_LL_2R = []
for i in range(len(D_lst)):
    print('Computing for distance of ' + str(D_lst[i]))

    geometry = VortexGeometry([rotor_opt, rotor_opt], [3, 3], [0, D_lst[i]], result_BEM,)
    result_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)
    results_LL_2R.append(result_LL_2R)

    [CT, CP] = performance_coefs_2R(result_LL_2R, result_LL, result_BEM, rotor_opt)
    CT_lst[0].append(CT[0])
    CT_lst[1].append(CT[1])
    CT_lst[2].append(CT[2])
    CP_lst[0].append(CP[0])
    CP_lst[1].append(CP[1])
    CP_lst[2].append(CP[2])

plot_difference_2R(D_lst, CT_lst, CP_lst, 'distance')

#%%