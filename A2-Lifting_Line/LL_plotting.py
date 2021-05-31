# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:24:49 2021

@author: Mathieu Pellé
"""

import matplotlib.pyplot as plt
import math as m
import numpy as np
# import logging
# logging.getLogger().setLevel(logging.CRITICAL)

x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
plt.rc('font', family='serif')

save=True # Save or not

#Radial plots for single rotor and multiple TSRs
def plot_radial(result_LL, result_BEM, rotors, TSR):
    #models = ['BEMT', 'LL']
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]']
    colours = ['deepskyblue', 'firebrick', 'mediumpurple']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for l in range(len(result_LL)):
            LL = result_LL[l][0]
            BEM = result_BEM[l]
            rotor = rotors[l]
            for j in range(2):
                if j==0:
                    dic = BEM
                    idx1 = 0
                    idx2 = len(dic.mu)
                else:
                    dic = LL
                    blade = 0
                    idx1 = blade*(len(rotor.mu)-1)
                    idx2 = idx1 + len(rotor.mu) -1
                if var[i]=='f_tan' or var[i]=='f_nor':
                    Z=getattr(dic, str(var[i]))/(0.5*rotor.wind_speed**2*rotor.radius)
                elif var[i]=='circulation':
                    Z=getattr(dic, str(var[i]))/((m.pi*rotor.wind_speed**2/(rotor.n_blades*rotor.omega)))

                else:
                    Z=getattr(dic, str(var[i]))

                if j==0:
                    plt.plot(dic.mu[idx1:idx2], Z[idx1:idx2],'--', color=colours[l])

                else:
                    plt.plot(dic.mu[idx1:idx2], Z[idx1:idx2], label = '$\lambda$='+str(TSR[l]), color=colours[l])
                plt.legend()
                if save==True:
                    plt.savefig('figures/TSR_'+str(var[i])+'.pdf')



#Performance coefficients (CT/CP) for single rotor
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


#Radial plots for double rotor and comparing to single rotor
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
                plt.savefig('figures/2R_'+str(var[i])+'.pdf')

#Performance coefficients (CT/CP) for double rotor and comparing to single rotor
def performance_coefs_2R(LL_2R, LL, BEM, rotor):
    CT_lst = []
    CP_lst = []
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
        CT_lst.append(CT)
        CP_lst.append(CP)
    print('BEMT CT = ' + str(BEM.CT), 'BEMT CP = ' + str(BEM.CP))

    return CT_lst, CP_lst

#CP/CT vs phase or distance
def plot_difference_2R(x, CT_lst, CP_lst, var):
    if var=='phase':
        xlab = 'Phase difference [deg]'

    else:
        xlab = 'Distance [m]'

    for i in range(2):
        if i==0:
            lst = CT_lst
            ylab = '$\Delta C_T$ [%]'
            name = 'CT'
        else:
            lst = CP_lst
            ylab = '$\Delta C_P$ [%]'
            name = 'CP'

        ref = lst[0][0]
        y1 = (lst[1]-ref)/ref*100
        y2 = (lst[2]-ref)/ref*100
        plt.figure()
        plt.plot(x, y1,'-o', label='Rotor 1')
        plt.plot(x, y2, label='Rotor 2')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.grid()
        plt.legend()

        if save==True:
            plt.savefig('figures/'+str(var)+str(name)+'.pdf')

#Radial plots for any blade or rotor for any phase
def plot_radial_2R_phase(LL_2Rs, LL, rotor, blades_all, phase_lst, phase, shift, line=None, colors=None):
    var = ['circulation','alpha']
    labels = ['$\Delta \Gamma$ [%]', '$\Delta$'+r'$\alpha$ [%]']
    shift1 = shift[0]
    shift2 = shift[1]
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])

        for j in range(3):
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
                ref = Z
                #plt.plot(dic.mu[idx1+shift1:idx2+shift2], Z[idx1+shift1:idx2+shift2], '--k', label = '1R')

            else:
                for p in range(len(phase_lst)):
                    LL_2R = LL_2Rs[phase_lst[p]]
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
                        Z = (Z-ref)/ref*100

                        if line:
                            name = 'phase_'+str(phase[phase_lst[p]])
                            plt.plot(dic.mu[idx1+shift1:idx2+shift2], Z[idx1+shift1:idx2+shift2],line[j-1], label='R'+str(j)+' B'+str(b+1), color=colors[b])
                        else:
                            name = 'phase_all'
                            plt.plot(dic.mu[idx1+shift1:idx2+shift2], Z[idx1+shift1:idx2+shift2], label='$\Delta \Phi$ = '+str(phase[phase_lst[p]])+'$\degree$')


        plt.legend()
        if save==True:
            plt.savefig('figures/2R_'+name+'_'+str(var[i])+'.pdf')

def plot_convergence_TSR(results_LL, TSR):
    colours = ['deepskyblue', 'firebrick', 'mediumpurple']
    #fig = plt.figure()
    ax = plt.gca()
    ax.grid()
    plt.xlabel('Number of iterations [-]')
    plt.ylabel('Error')
    ax.set_yscale('log')
    ax.set_xscale('log')
    for i in range(len(TSR)):
        LL = results_LL[i][0]
        error = getattr(LL, 'error')
        ax.plot(np.linspace(0,len(error),len(error)), error, label='$\lambda = $'+str(TSR[i]), color=colours[i])

    plt.legend()
    if save==True:
        plt.savefig('figures/convergence_TSR.pdf')