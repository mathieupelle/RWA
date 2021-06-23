# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""


import numpy as np
from VP_plotting import steady_polars, steady_contours, scatter, unsteady_polars
from VP_utilities import vortex_panel


#%% Steady case - Polar, pressure and velocity contours

#Polars
alpha = np.arange(-15,15,1)
results = []
for i in range(len(alpha)):
    result = vortex_panel([0], 10, [alpha[i]], [0], U_inf=2)
    results.append(result)

steady_polars(alpha, results)

#Contours
result = vortex_panel([0], 10, [10], [0], U_inf=1)
steady_contours(alpha, result, quiver=False, streamlines=True)

#%% Unsteady case - oscillations

k=0.02
U_inf = 5
c = 0.1
omega=k*2*U_inf/c
T = 2*np.pi/omega
t = np.hstack((np.linspace(0, T, 40),np.linspace(0, T, 40)))

theta = 10+10*np.sin(omega*t)
theta_dot = omega*np.cos(omega*t)

result = vortex_panel(t, 10, theta, theta_dot, c=c, U_inf=U_inf)

#scatter([result])

unsteady_polars(theta, result, basic_method=False)

#%% Unsteady case - step

# U_inf = 5
# c = 0.1
# t=np.arange(0,2,0.1)
# theta=np.ones(len(t))*5
# theta[:5] = 0
# theta_dot = np.zeros(len(t))
# result = vortex_panel(t, 20, theta, theta_dot,  c=c, U_inf=U_inf)

# unsteady_polars(theta, result, inputs=True, arrows=False)

#%%
