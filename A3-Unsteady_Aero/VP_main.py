# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""


import numpy as np
from VP_plotting import steady_polars, contours, scatter, unsteady_polars, flap_analysis, step_response
from VP_utilities import vortex_panel



#%% Steady case - Polars

alpha = np.arange(-20,20,2)
results = []
for i in range(len(alpha)):
    result = vortex_panel([0], 10, [alpha[i]], [0])
    results.append(result)

steady_polars(alpha, results)

#%% Steady case - Contours

result = vortex_panel([0], 30, [15], [0], U_inf_vec=[np.array([[1],[0]])], c=1)
scatter([result])
contours(result, streamlines=True)

#%% Unsteady case - oscillations

#TODO Check orientation of arrows and make plots clearer

k=[0.01, 0.05, 0.1]
U_inf = 1
c = 0.1

for i in range(len(k)):

    omega=k[i]*2*U_inf/c
    T = 2*np.pi/omega
    N = 81
    dt = 2*T/N
    t = np.linspace(0,2*T, N)
    U_inf_vec = [U_inf*np.array([[1],[0]])]*len(t)
    theta = 10+5*np.sin(omega*t)
    theta_dot = omega*np.cos(omega*t)

    result = vortex_panel(t, 10, theta, theta_dot, c=c, U_inf_vec=U_inf_vec)

    #scatter([result])
    unsteady_polars(theta, result, quasi=True, inputs=True, arrows=True)
    contours(result, rho=1.225, streamlines=True, frames=[5])

#%% Unsteady case - gust
#TODO Check if gust and angle of attack input make sense

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
contours(result, rho=1.225, streamlines=True, frames=[10])

#%% Unsteady case - step in angle of attack

t=np.arange(0,2,0.02)
theta=np.ones(len(t))*5
theta_dot = np.zeros(len(t))
U_inf_vec= [np.array([[1],[0]])]*len(t)

result = vortex_panel(t, 5, theta, theta_dot, c=0.1, U_inf_vec=U_inf_vec)

step_response(theta, result, step='angle', idx=0)
contours(result, rho=1.225, streamlines=True, frames=[10])

#%% Steady case - flap polars

flap_angle = np.linspace(0,25,5)
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

flap = {'length':0.4, 'angle':10, 'N_panels':3}
result = vortex_panel([0], 10, [15], [0], c=1, U_inf_vec=[np.array([[2],[0]])], flap=flap)

contours(result, streamlines=True, flap=flap)


