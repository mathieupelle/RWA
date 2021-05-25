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
   # var=['circulation']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for j in range(2):
            if j==0:
                dic = LL
                blade = 0
                idx1 = blade*(len(rotor.mu)-1)
                idx2 = idx1 + len(rotor.mu) -1
            else:
                dic = BEM
                idx1 = 0
                idx2 = len(dic.mu)
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
    print('BEMT CP = ' + str(BEM.CP), 'LL CP = ' + str(CP))



def plot_radial_2R(LL_2R, LL, rotor, blades_all):

    var=['alpha','phi','a','ap','f_tan','f_nor','circulation']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])

        for j in range(len(LL_2R)+1):
            if j==0:
                dic = LL[0]
                blade = 0
                idx1 = blade*(len(rotor.mu)-1)
                idx2 = idx1 + len(rotor.mu) -1
                if var[i]=='f_tan' or var[i]=='f_nor':
                    Z=getattr(dic, str(var[i]))/(0.5*rotor.wind_speed**2*rotor.radius)
                elif var[i]=='circulation':
                    Z=getattr(dic, str(var[i]))/((m.pi*rotor.wind_speed**2/(rotor.n_blades*rotor.omega)))
                else:
                    Z=getattr(dic, str(var[i]))

                plt.plot(dic.mu[idx1:idx2], Z[idx1:idx2], '--k', label = '1R')

            else:
                dic = LL_2R[j-1]
                blades = blades_all[j-1]
                for k in range(len(blades)):
                    b =  blades[k]
                    idx1 = b*(len(rotor.mu)-1)
                    idx2 = idx1 + len(rotor.mu) -1
                    if var[i]=='f_tan' or var[i]=='f_nor':
                        Z=getattr(dic, str(var[i]))/(0.5*rotor.wind_speed**2*rotor.radius)
                    elif var[i]=='circulation':
                        Z=getattr(dic, str(var[i]))/((m.pi*rotor.wind_speed**2/(rotor.n_blades*rotor.omega)))
                    else:
                        Z=getattr(dic, str(var[i]))

                    plt.plot(dic.mu[idx1:idx2], Z[idx1:idx2], label = '2R'+str(j)+' B'+str(b))

            plt.legend()
            if save==True:
                plt.savefig('figures/TSR_'+str(var[i])+'.pdf')


def performance_coefs_2R(LL_2R, LL, BEM, rotor):

    for i in range(len(LL_2R)+1):
        if i==0:
            LL = LL[0]
        else:
            LL = LL_2R[i-1]
        N = len(rotor.mu)
        dr = np.zeros(N-1)
        for n in range(N-1):
            dr[n] = (rotor.mu[n+1]-rotor.mu[n])*rotor.radius
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
        if i==0:
            print('LL 1R CT = ' + str(CT), 'LL 1R CP = ' + str(CP))
        else:
            print('LL 2R'+str(i)+' CT = ' + str(CT), 'LL 2R'+str(i)+' CP = ' + str(CP))

    print('BEMT CT = ' + str(BEM.CT), 'BEMT CP = ' + str(BEM.CP))