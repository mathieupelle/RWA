# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""


import numpy as np
from VP_plotting import steady_polars, steady_contours, scatter, unsteady_polars
from VP_utilities import vortex_panel



#%% Steady case - Polars

alpha = np.arange(-20,20,2)
results = []
for i in range(len(alpha)):
    result = vortex_panel([0], 10, [alpha[i]], [0])
    results.append(result)

steady_polars(alpha, results)

#%% Steady case - Contours

result = vortex_panel([0], 30, [15], [0], U_inf=1, c=1)
scatter([result])
steady_contours(alpha, result, streamlines=True)

#%% Unsteady case - oscillations

k=0.01
U_inf = 5
c = 0.05
omega=k*2*U_inf/c
T = 2*np.pi/omega
N = 81
dt = 2*T/N
t = np.linspace(0,2*T, N)

theta = 10+5*np.sin(omega*t)
theta_dot = omega*np.cos(omega*t)

result = vortex_panel(t, 10, theta, theta_dot, c=c, U_inf=U_inf)

#scatter([result])

unsteady_polars(theta, result, basic_method=True, inputs=True)

#%% Unsteady case - step

# U_inf = 5
# c = 0.1
# t=np.arange(0,2,0.1)
# theta=np.ones(len(t))*5
# theta[:5] = 0
# theta_dot = np.zeros(len(t))
# result = vortex_panel(t, 20, theta, theta_dot,  c=c, U_inf=U_inf)

# unsteady_polars(theta, result, inputs=True, arrows=False)

#%% Steady case - flap

flap_angle = np.linspace(0,45,10)
alpha = np.arange(0,15,1)
flap_results = {}

for j in range(len(flap_angle)):
    flap = {'length':0.2, 'angle':flap_angle[j], 'N_panels':5}
    results = []
    for i in range(len(alpha)):
        result = vortex_panel([0], 10, [alpha[i]], [0], c=0.8, U_inf=2, flap=flap)
        results.append(result)
    flap_results['flap_angle'+str(flap_angle[j])] = results


flap_length = np.linspace(0.1,1,10)
alpha = np.arange(0,15,1)
flap_results = {}

for j in range(len(flap_length)):
    flap = {'length':flap_length[j], 'angle':10, 'N_panels':5}
    results = []
    for i in range(len(alpha)):
        result = vortex_panel([0], 10, [alpha[i]], [0], c=1, U_inf=2, flap=flap)
        results.append(result)
    flap_results['flap_length'+str(flap_length[j])] = results

#%%
import matplotlib.pyplot as plt


def flap_analysis(flap_results, flap, alpha, parameter, parameter_lst, contour=False, rho=1.225):
    Cl = np.zeros((len(parameter_lst),len(alpha)))
    for j in range(len(parameter_lst)):
        results = flap_results[parameter+str(parameter_lst[j])]
        for i in range(len(alpha)):
            result = results[i]
            c = result['chord']+flap['length']
            U_inf = np.linalg.norm(result['velocity'])
            dL = rho*U_inf*(result['gamma'][0])
            L = sum(dL)
            Cl[j,i] = L/(0.5*rho*U_inf**2*c)

        if j%2==0:
            plt.plot(alpha, Cl[j,:], label=r'$\beta$='+str(round(parameter_lst[j]))+'$^\circ$')
        plt.legend()
        plt.xlabel(r'$\alpha$ [$^\circ$]')
        plt.ylabel('$C_l$ [-]')
        plt.grid()


    if contour:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        X = alpha
        Y = parameter_lst
        X, Y = np.meshgrid(X, Y)

        ax.plot_surface(X, Y, Cl, cmap='viridis')
        ax.set_ylabel(r'$\beta$ [$^\circ$]')
        ax.set_xlabel(r'$\alpha$ [$^\circ$]')
        ax.set_zlabel('$C_l$ [-]')


#flap_analysis(flap_results, flap, alpha, 'flap_angle', flap_angle, rho=1.225)
flap_analysis(flap_results, flap, alpha, 'flap_length', flap_length, rho=1.225)