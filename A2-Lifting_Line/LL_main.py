# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:39:42 2021

@author: Mathieu Pellé
"""

from LL_Utilities import VortexGeometry, LiftingLine
from BEMT_Utilities import BEMT_execute
from LL_plotting import plot_radial, performance_coefs, plot_radial_2R, performance_coefs_2R, plot_difference_2R, plot_radial_2R_phase
import numpy as np

#%% Single rotor case - Radial plots, TSR effect, Performance coefficients

#Define common parameters between BEM and LL
N_radial = 35
spacing = 'cos'
U_inf = 10
TSR = [8]
results_LL = []
results_BEM = []
rotors = []
for i in range(len(TSR)):
    # Call for BEM geometry and results
    rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR[i])
    rotors.append(rotor_opt)
    results_BEM.append(result_BEM)

    # Call for Lifting Line geometry and results
    #Lists of: rotor geometries, No. of Radial Elements, No. of Rotations, Rotor positions
    geometry = VortexGeometry([rotor_opt], [5], [0], result_BEM)

    #List of rotors, combined geometry and BEM results
    result_LL = LiftingLine([rotor_opt],geometry,result_BEM)
    results_LL.append(result_LL)

    # Coefficients
    print('=====TSR '+str(TSR[i])+'=====')
    performance_coefs(result_LL, result_BEM, rotor_opt)


# Radial plots
plot_radial(results_LL, results_BEM, rotors, TSR)


#%% Double rotor case - Radial plots and performance coefficients for both rotors and all blades

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
geometry = VortexGeometry([rotor_opt, rotor_opt], [5, 5], [0, 50*2], result_BEM)
results_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)

# Radial plots compared to single rotor case
plot_radial_2R(results_LL_2R, result_LL, rotor_opt, [[1,2,3],[1,2,3]])

# Performance coefficients compared to single rotor case and BEM
performance_coefs_2R(results_LL_2R, result_LL, result_BEM, rotor_opt)



#%% Double rotor effect of phase difference - CT/CP vs phase, specific radial plots per blade, phase, rotor
#COSTLY COMPUTATION!! change phase list or N_radial

#Common parameters
N_radial = 35
spacing = 'cos'
TSR = 8
U_inf = 10

#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Single rotor
geometry = VortexGeometry([rotor_opt], [5], [0], result_BEM)
result_LL = LiftingLine([rotor_opt],geometry,result_BEM)

#Varying over phase angle from 0deg to 120deg
phase = np.linspace(0,120,21)
CT_lst = [[], [], []]
CP_lst = [[], [], []]
results_LL_2R = []
for i in range(len(phase)):
    print('Computing for phase difference of ' + str(phase[i]))

    geometry = VortexGeometry([rotor_opt, rotor_opt], [5, 5], [0, 50*2], result_BEM, phase[i])
    result_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)
    results_LL_2R.append(result_LL_2R)

    [CT, CP] = performance_coefs_2R(result_LL_2R, result_LL, result_BEM, rotor_opt)
    CT_lst[0].append(CT[0])
    CT_lst[1].append(CT[1])
    CT_lst[2].append(CT[2])
    CP_lst[0].append(CP[0])
    CP_lst[1].append(CP[1])
    CP_lst[2].append(CP[2])

#CT/CP vs phase
plot_difference_2R(phase, CT_lst, CP_lst, 'phase')

#Radial plot for different phases on 1 blade
blades_all = [[0],[]]
phase_lst = [0,7,10,14,18]
shift1 = 7
shift2 = -3

plot_radial_2R_phase(results_LL_2R, result_LL, rotor_opt, blades_all, phase_lst, phase, [shift1, shift2])

#Radial plot for multiple blades for 1 phase
blades_all = [[0,1,2],[0,1,2]]
phase_lst = [10]
colors = ['black', 'orangered', 'royalblue']
line = ['o', '-']

plot_radial_2R_phase(results_LL_2R, result_LL, rotor_opt, blades_all, phase_lst, phase, [shift1, shift2], line, colors)


#%% Double rotor effect of rotor distance - CT/CP vs distance
#COSTLY COMPUTATION!! change distance list or N_radial

N_radial = 35
spacing = 'cos'
TSR = 8
U_inf = 10
D = 2*50
#Call for BEM geometry and results
rotor_opt, result_BEM = BEMT_execute(N_radial,spacing,U_inf,TSR)

#Single rotor
geometry = VortexGeometry([rotor_opt], [5], [0], result_BEM)
result_LL = LiftingLine([rotor_opt],geometry,result_BEM)

#Varying for different distances
D_lst = D*np.hstack((np.linspace(1,4,7),np.array([5,7,9,11])))
CT_lst = [[], [], []]
CP_lst = [[], [], []]
results_LL_2R = []
for i in range(len(D_lst)):
    print('Computing for distance of ' + str(D_lst[i]))

    geometry = VortexGeometry([rotor_opt, rotor_opt], [5, 5], [0, D_lst[i]], result_BEM,)
    result_LL_2R = LiftingLine([rotor_opt, rotor_opt],geometry,result_BEM)
    results_LL_2R.append(result_LL_2R)

    [CT, CP] = performance_coefs_2R(result_LL_2R, result_LL, result_BEM, rotor_opt)
    CT_lst[0].append(CT[0])
    CT_lst[1].append(CT[1])
    CT_lst[2].append(CT[2])
    CP_lst[0].append(CP[0])
    CP_lst[1].append(CP[1])
    CP_lst[2].append(CP[2])

# CT/CP vs distance
plot_difference_2R(D_lst, CT_lst, CP_lst, 'distance')

#%%




