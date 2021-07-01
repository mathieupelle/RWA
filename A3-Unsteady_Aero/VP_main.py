# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""


import numpy as np
from VP_plotting import steady_polars, contours, scatter, unsteady_polars, flap_analysis, step_response
from VP_utilities import vortex_panel,Sensitivity_NPanels, Sensitivity_DeltaT


#%% Steady case - Polars

alpha = np.arange(-10,22,2)
results = []
for i in range(len(alpha)):
    result = vortex_panel([0], 10, [alpha[i]], [0])
    results.append(result)

steady_polars(alpha, results)

#%% Steady case - Contours
alpha=20
result = vortex_panel([0], 30, [alpha], [0], U_inf_vec=[np.array([[1],[0]])], c=1)
#scatter([result])
bound = [[0.2,1.5],[-0.5,0.8],[1,-2]]
contours(result, streamlines=True, condition=alpha, bound=bound)

#%% Unsteady case - oscillations

k=[0.01, 0.05, 0.1]
U_inf = 1
c = 1
results = []
t_frame=5

for i in range(len(k)):

    omega=k[i]*2*U_inf/c
    T = 2*np.pi/omega
    N = 191
    dt = 2*T/N
    t = np.linspace(0,2*T, N)
    idx_frames = ((np.abs(t - t_frame)).argmin())
    U_inf_vec = [U_inf*np.array([[1],[0]])]*len(t)
    theta = 15+10*np.sin(omega*t)
    theta_dot = omega*np.cos(omega*t)

    result = vortex_panel(t, 2, theta, theta_dot, c=c, U_inf_vec=U_inf_vec)
    results.append(result)

    #scatter([result])
    unsteady_polars(theta, result, rho=1.225, condition=str(k[i]))
    contours(result, rho=1.225, streamlines=True, frames=[idx_frames], condition=str(k[i]))


#%% Unsteady case - gust
U0 = np.array([[1],[0]])
U1 = np.array([[1],[0.1]])

t=np.arange(0,2,0.05)
t_step = 1
U_inf_vec = [U0]*int(t_step)+[U1]*int(len(t)-t_step)


theta=np.zeros(len(t))
theta_dot = np.zeros(len(t))
result = vortex_panel(t, 5, theta, theta_dot, c=0.1, U_inf_vec=U_inf_vec)

theta = []
for i in range(len(U_inf_vec)):
    U = U_inf_vec[i]
    theta.append(np.rad2deg(np.arctan(U[1,0]/U[0,0])))

step_response(theta, result, 'gust')
contours(result, rho=1.225, streamlines=True, frames=[6], condition='gust')

#%% Unsteady case - step in angle of attack

t=np.arange(0,2,0.02)
theta=np.ones(len(t))*5
theta_dot = np.zeros(len(t))
U_inf_vec= [np.array([[1],[0]])]*len(t)

result = vortex_panel(t, 5, theta, theta_dot, c=0.1, U_inf_vec=U_inf_vec)

step_response(theta, result, step='angle', idx=0)

contours(result, rho=1.225, streamlines=True, frames=[10], condition='step')

#%% Steady case - flap polars

flap_angle = np.linspace(0,24,5)
alpha = np.arange(-5,20,2)
flap_results = {}
flaps = []
for j in range(len(flap_angle)):
    flap = {'length':0.2, 'angle':flap_angle[j], 'N_panels':5}
    flaps.append(flap)
    results = []
    for i in range(len(alpha)):
        result = vortex_panel([0], 10, [alpha[i]], [0], c=0.8, flap=flap)
        results.append(result)
    flap_results[j] = results

flap_analysis(flap_results, flaps, alpha, 'flap_angle', theory=True)

flap_length = np.linspace(0,1,5)
alpha = np.arange(-5,20,2)
flap_results = {}
flaps = []
for j in range(len(flap_length)):
    flap = {'length':flap_length[j], 'angle':10, 'N_panels':5}
    flaps.append(flap)
    results = []
    for i in range(len(alpha)):
        result = vortex_panel([0], 10, [alpha[i]], [0], c=1, flap=flap)
        results.append(result)
    flap_results[j] = results

flap_analysis(flap_results, flaps, alpha, 'flap_length',theory=True)



#%% Steady case - flap contours
alpha = 20

flap = {'length':0.25, 'angle':18, 'N_panels':3}
result = vortex_panel([0], 10, [alpha], [0], c=1, U_inf_vec=[np.array([[1],[0]])], flap=flap)

contours(result, streamlines=True, flap=flap, condition='flap'+str(alpha))


#%% Unsteady case - flap contours

t=np.arange(0,2,0.1)
theta=np.ones(len(t))*5
theta_dot = np.zeros(len(t))
U_inf_vec= [np.array([[1],[0]])]*len(t)

flap = {'length':0.2, 'angle':45, 'N_panels':3}

result = vortex_panel(t, 5, theta, theta_dot, c=1, U_inf_vec=U_inf_vec, flap=flap, shed_loc_ratio=0.3)

scatter([result])
contours(result, rho=1.225, streamlines=True, flap=flap, frames=[19])

#%% Sensitivity study

#Number of panels to be tested
N_panels = np.logspace(0, 2, 15) 
Sensitivity_NPanels(N_panels)

#%% Effect of dt
dt = np.log10(np.logspace((0.02),(0.5),6))
Sensitivity_DeltaT(dt)