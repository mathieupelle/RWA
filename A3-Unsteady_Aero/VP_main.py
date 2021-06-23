# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""


import numpy as np
from VP_plotting import polars, steady_contours, scatter
from VP_utilities import vortex_panel


#%% Steady case - Polar, pressure and velocity contours

alpha = np.arange(-15,15,1)
results = []
for i in range(len(alpha)):
    result = vortex_panel([0], 10, [alpha[i]], [0], U_inf=1)
    results.append(result)


polars(alpha, results, mode='steady')

steady_contours(alpha, results[28], quiver=True)


#%%

k=0.02
omega=k*2*10/0.2
T = 2*np.pi/omega
t=np.linspace(0,T, 50)

theta = 15+10*np.sin(omega*t)
theta_dot = omega*np.cos(omega*t)

result = vortex_panel(t, 20, theta, theta_dot, c=0.2, U_inf=10)

scatter([result])

#STEP
t=np.arange(0,2,0.1)
theta=np.ones(len(t))*5
theta[:2] = 0
theta_dot = np.zeros(len(t))
result = vortex_panel(t, 20, theta, theta_dot, c=0.2, U_inf=10)

#%%

import matplotlib.pyplot as plt
def polars(alpha, result, rho=1.225, mode='steady'):

    if mode=='unsteady':
        C=np.zeros(len(alpha))
        C_theory=np.zeros(len(alpha))
        N_panels = result['N_panels']
        for i in range(len(alpha)):
            U_inf = np.linalg.norm(result['velocity'])
            dL = rho*U_inf*-(result['gamma'][i][:N_panels])
            L = sum(dL)

            C[i] = L/(0.5*rho*U_inf**2*result['chord'])
            C_theory[i]=2*np.pi*np.sin(np.deg2rad(alpha[i]))
            label = '$C_l$ [-]'


        plt.figure()
        plt.plot(alpha, C, 'x' ,label='Vortex panel code')
        plt.plot(alpha, C_theory, '--k', label='2D Airfoil theory')
        plt.grid()
        plt.xlabel(r'$\alpha$ [$^{\circ}$]')
        plt.ylabel(label)
        plt.legend()

        plt.figure()
        plt.plot(result['time'], C, 'x')
        plt.grid()
        plt.xlabel('time [s]')
        plt.ylabel(label)

polars(theta, result, mode='unsteady')


#%%
from VP_utilities import VOR2D

rho = 1.225
N_panels = result['N_panels']
Cl_lst  =[]
time = result['time']
dt = time[1]-time[0]
Cl = np.zeros(len(time))
Cl_theory = np.zeros(len(time))
for t in range(len(time)):
    tan_vec = result['TE'][t]-result['LE'][t]
    tan_vec = tan_vec/np.linalg.norm(tan_vec)
    p = np.zeros(N_panels)
    for j in range(N_panels):
        v=result['velocity']
        for i in range(len(result['vortices'][t])-N_panels):
            i = i+N_panels #TODO CHECK this for other functions
            v = v + VOR2D(result['panels_history'][t][j], result['vortices'][t][i], result['gamma'][t][i,0])

        Q = np.dot(tan_vec.T, v)
        if t==0:
            dgamma_dt = 0 # ??? check
        else:
            dgamma_dt = (sum(result['gamma'][t][:N_panels])-sum(result['gamma'][t-1][:N_panels]))/dt
        p[j] = (Q*result['gamma'][t][j,0]/result['L_panels'] + dgamma_dt)*rho


    L = -sum(p*result['L_panels'])*np.cos(np.deg2rad(theta[t]))
    Cl[t] = L/(0.5*rho*(np.linalg.norm(result['velocity']))**2*result['chord'])
    Cl_theory[t]=2*np.pi*np.sin(np.deg2rad(theta[t]))

plt.scatter(theta, Cl, marker='.')
#plt.plot(theta, Cl_theory, '-k', alpha=0.8)
plt.grid()
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('$C_l$ [-]')
