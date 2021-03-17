# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:03:39 2021

@author: MathieuPelle
"""

import matplotlib.pyplot as plt
import numpy as np


x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# input Res_org and TSR_lst=[6,8,10]
def plot_TSR(results,TSR_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor',] #add Cq and circulation
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for j in range(len(TSR_lst)):
            dic=results['TSR'+str(TSR_lst[j])+'_yaw'+str(0)]
            plt.plot(dic.mu, getattr(dic, str(var[i])),label= '$\lambda$=' +str(TSR_lst[j]))

        plt.legend()
        #plt.savefig('figures/TSR_'+str(var[j])+'.pdf')

# input Res_org and yaw_lst=[0,15,30]
def plot_yaw(results,yaw_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor',] #add Cq
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]']
    for i in range(len(yaw_lst)):
        for j in range(len(var)):
            dic=results['TSR'+str(8)+'_yaw'+str(yaw_lst[i])]
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            #ax.set_theta_zero_location("N")
            if yaw_lst[i]==0:
                dic2=results['TSR'+str(8)+'_yaw'+str(yaw_lst[1])]
                r = np.hstack((dic2.mu,dic2.mu[:,0].reshape(len(dic2.mu[:,0]),1)))
                psi = np.hstack((dic2.azimuth,2*np.pi*np.ones((len(r[:,0]),1))))
                pp=ax.contourf(psi,r , np.tile(getattr(dic, str(var[j])),len(psi[0])))
            else:
                r = np.hstack((dic.mu,dic.mu[:,0].reshape(len(dic.mu[:,0]),1)))
                psi = np.hstack((dic.azimuth,2*np.pi*np.ones((len(r[:,0]),1))))
                Z = getattr(dic, str(var[j]))
                Z = np.hstack((Z,Z[:,0].reshape(len(Z[:,0]),1)))
                pp=ax.contourf(psi,r,Z)
            cbar = plt.colorbar(pp, orientation='vertical', pad=0.1)
            cbar.ax.set_ylabel(labels[j])
            plt.set_cmap('viridis')
            #plt.savefig('figures/Yaw_'+str(var[j])+'.pdf')

# input Rotor_org
def plot_polars(dic):
    plt.figure()
    plt.grid()
    plt.xlabel(r'$\alpha [deg]$')
    plt.ylabel('$C_l [-]$')
    plt.xlim([-20,30])
    plt.plot(dic.polars.alpha,dic.polars.Cl)
    plt.figure()
    plt.grid()
    plt.xlabel('$C_d [-]$')
    plt.ylabel('$C_l [-]$')
    plt.xlim([0,.1])
    plt.plot(dic.polars.Cd,dic.polars.Cl)



#Input Res_org, TSR and Yaw case (Res_org,8,0)
def plot_correction(results, TSR, Yaw):
    dic=results['TSR'+str(TSR)+'_yaw'+str(Yaw)]
    plt.figure()
    plt.grid()
    plt.plot(dic.mu, dic.f,color='k', label='Prandtl',zorder=1)
    plt.scatter(dic.mu, dic.f_tip,marker='.',s=80, label='Prandtl tip', zorder=2)
    plt.scatter(dic.mu, dic.f_root,marker='.',s=80, label='Prandtl root', zorder=3)
    plt.xlabel('r/R')
    plt.ylabel('f [-]')
    plt.legend()
