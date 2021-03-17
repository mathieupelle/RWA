# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:03:39 2021

@author: MathieuPelle
"""

import matplotlib.pyplot as plt
import numpy as np

'''
Put this in main, change save or not save option and uncomment latex default text
↓↓↓↓↓↓↓↓↓↓↓↓↓↓
'''
#plot_TSR(Res_org,Rotor_org,[6,8,10])
#plot_yaw(Res_org,Rotor_org,[15,30])
# plot_polars(Rotor_org)
# plot_correction(Res_org, 8, 0)
#plot_enthalpy(Res_org,Rotor_org,8,0)


x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
save=False

def plot_TSR(results,param,TSR_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        for j in range(len(TSR_lst)):
            dic=results['TSR'+str(TSR_lst[j])+'_yaw'+str(0)]
            if var[i]=='f_tan' or var[i]=='f_nor':
                Z=getattr(dic, str(var[i]))/(0.5*param.rho*param.wind_speed**2*param.radius)
            elif var[i]=='circulation':
                Z=getattr(dic, str(var[i]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
            else:
                Z=getattr(dic, str(var[i]))
            plt.plot(dic.mu,Z,label='$\lambda$=' +str(TSR_lst[j]))

        plt.legend()
        if save==True:
            plt.savefig('figures/TSR_'+str(var[j])+'.pdf')

def plot_yaw(results,param,yaw_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$']
    for i in range(len(yaw_lst)):
        for j in range(len(var)):
            dic=results['TSR'+str(8)+'_yaw'+str(yaw_lst[i])]
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            Z = getattr(dic, str(var[j]))
            if var[j]=='f_tan' or var[j]=='f_nor':
                Z=getattr(dic, str(var[j]))/(0.5*param.rho*param.wind_speed**2*param.radius)
            elif var[i]=='circulation':
                Z=getattr(dic, str(var[j]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
            else:
                Z=getattr(dic, str(var[j]))
            if yaw_lst[i]==0:
                dic2=results['TSR'+str(8)+'_yaw'+str(yaw_lst[1])]
                r = np.hstack((dic2.mu,dic2.mu[:,0].reshape(len(dic2.mu[:,0]),1)))
                psi=dic2.azimuth
                psi=psi-psi[0]*np.ones(len(psi))
                psi=np.append(psi,2*np.pi)
                psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
                pp=ax.contourf(psi,r , np.tile(Z,len(psi[0])),20)
            else:
                ax.set_theta_zero_location('N') ## CHECK!!!!!!!
                r = np.hstack((dic.mu,dic.mu[:,0].reshape(len(dic.mu[:,0]),1)))
                psi=dic.azimuth
                psi=psi-psi[0]*np.ones(len(psi))
                psi=np.append(psi,2*np.pi)
                psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
                Z = np.hstack((Z,Z[:,0].reshape(len(Z[:,0]),1)))
                pp=ax.contourf(psi,r,Z,20)
            ax.set_rlabel_position(135) ## Change depending on visualisation
            rlabels = ax.get_ymajorticklabels()
            for label in rlabels:
                label.set_color('white')
            cbar = plt.colorbar(pp, orientation='vertical', pad=0.1)
            cbar.ax.set_ylabel(labels[j])
            plt.set_cmap('viridis')
            if save==True:
                plt.savefig('figures/Yaw'+str(yaw_lst[i])+'_'+str(var[j])+'.pdf')

def plot_polars(dic):
    plt.figure()
    plt.grid()
    plt.xlabel(r'$\alpha [deg]$')
    plt.ylabel('$C_l [-]$')
    plt.xlim([-20,30])
    plt.plot(dic.polars.alpha,dic.polars.Cl)
    if save==True:
        plt.savefig('figures/Cl.pdf')
    plt.figure()
    plt.grid()
    plt.xlabel('$C_d [-]$')
    plt.ylabel('$C_l [-]$')
    plt.xlim([0,.1])
    plt.plot(dic.polars.Cd,dic.polars.Cl)
    if save==True:
        plt.savefig('figures/Cd.pdf')

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
    if save==True:
        plt.savefig('figures/tip_corrections_'+str(TSR)+'_'+str(Yaw)+'.pdf')

def plot_enthalpy(results, param, TSR, Yaw):
    dic=results['TSR'+str(TSR)+'_yaw'+str(Yaw)]
    dic2=param
    enthalpy_1=np.ones((len(dic.mu),1))*(101325/dic2.rho+dic2.wind_speed**2/2)
    enthalpy_2=enthalpy_1
    enthalpy_3=dic.enthalpy_3+np.ones((len(dic.enthalpy_3),1))*101325/dic2.rho
    enthalpy_4=np.ones((len(dic.mu),1))*((dic2.wind_speed*(1-2*dic.a_global))**2/2+101325/dic2.rho)
    plt.figure()
    plt.grid()
    plt.plot(dic.mu, enthalpy_1, 'k', linewidth=2.5, label='Infinity upstream',zorder=1)
    plt.plot(dic.mu, enthalpy_2, 'm--',label='Rotor upwind',zorder=2)
    plt.plot(dic.mu, enthalpy_3, label='Rotor downwind',zorder=3)
    plt.plot(dic.mu, enthalpy_4, label='Infinity downwind',zorder=4)
    plt.xlabel('r/R')
    plt.ylabel(r'$h_0$ [-]')
    plt.legend(bbox_to_anchor = [0.58, 0.5])
    if save==True:
        plt.savefig('figures/enthalpy_'+str(TSR)+'_'+str(Yaw)+'.pdf')


# Need data for radial spacing on axial induction

# def plot_spacing(results,TSR,Yaw):
#     dic=results['TSR'+str(TSR)+'_yaw'+str(Yaw)]
#     plt.figure()
#     plt.grid()
#     plt.plot(dic.mu, dic.a,'b-x', label='Constant spacing N='+str(len(dic.mu)),zorder=2)
#     plt.plot(dic.mu, dic.a,'g-x', label='Constant spacing N='+str(len(dic.mu)),zorder=1)
#     plt.plot(dic.mu, dic.a,'r--x', label='Cosine spacing',zorder=3)
#     plt.xlabel('r/R')
#     plt.ylabel('a [-]')
#     plt.legend()
#     if save==True:
#         plt.savefig('figures/tip_corrections_'+str(TSR)+'_'+str(Yaw)+'.pdf')


