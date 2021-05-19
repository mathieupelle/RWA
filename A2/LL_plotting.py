# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:24:49 2021

@author: Mathieu Pellé
"""

import matplotlib.pyplot as plt
import math as m
import numpy as np


x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
plt.rc('font', family='serif')

save=False # Save or not


def plot_radial(LL, BEM, rotor):
    LL = LL[0] #CHANGE
    models = ['LL', 'BEMT']
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for j in range(2):
            if j==0:
                dic = LL
                idx1 = 0
                idx2 = int(len(LL.mu)/rotor.n_blades)
            else:
                dic = BEM
                idx1 = 0
                idx2 = -1
            if var[i]=='f_tan' or var[i]=='f_nor':
                Z=getattr(dic, str(var[i]))/(0.5*rotor.wind_speed**2*rotor.radius)
            elif var[i]=='circulation':
                Z=getattr(dic, str(var[i]))/((m.pi*rotor.wind_speed**2/(rotor.n_blades*rotor.omega)))
            else:
                Z=getattr(dic, str(var[i]))

            plt.plot(dic.mu[idx1:idx2], Z[idx1:idx2], label = models[j])

            plt.legend()
            if save==True:
                plt.savefig('figures/TSR_'+str(var[i])+'.pdf')




def performance_coefs(LL, BEM, rotor):

    LL = LL[0]
    N = len(rotor.mu)
    dr = np.zeros(N-1)
    for i in range(N-1):
        dr[i] = (rotor.mu[i+1]-rotor.mu[i])*rotor.radius
    CT = 0
    CP = 0
    for j in range(rotor.n_blades):
        idx1 = j*(N-1)
        idx2 = idx1 + N -1
        f_nor = LL.f_nor[idx1:idx2]
        f_tan = LL.f_tan[idx1:idx2]
        mu = LL.mu[idx1:idx2]
        CT += np.dot(f_nor, dr)/(0.5*rotor.wind_speed**2*m.pi*rotor.radius**2)
        CP += np.dot(dr, f_tan*mu*rotor.radius)*rotor.omega/(0.5*rotor.wind_speed**3*m.pi*rotor.radius**2)
    print('BEMT CT = ' + str(BEM.CT), 'LL CT = ' + str(CT))
    print('BEMT CP = ' + str(BEM.CP), 'LL CT = ' + str(CP))