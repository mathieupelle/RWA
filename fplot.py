# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:03:39 2021

@author: MathieuPelle
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm


'''
1) Uncomment latex defaul text
2) Create folder in directory called 'Figures'
3) Put this in main
↓↓↓↓↓↓↓↓↓↓↓↓↓↓
'''

#from fplot import*

#plot_TSR(Res_org,Rotor_org,[6,8,10])
#plot_yaw(Res_org,Rotor_org,[0,15,30])
# plot_polars(Rotor_org)
# plot_correction(Res_org, 8, 0)
# plot_enthalpy(Res_org,Rotor_org,8,0)
#enthalpy_tube(Res_org,Rotor_org,8,0)


x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')


save=True # Save or not

def plot_TSR(results,param,TSR_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ','local_CT']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$', '$C_T [-]$']
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
            plt.savefig('figures/TSR_'+str(var[i])+'.pdf')

def plot_yaw(results,param,yaw_lst):
    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ', 'local_CT']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$','$C_T [-]$']

    for j in range(len(var)):
        fig, axs = plt.subplots(1, len(yaw_lst), figsize=(12,5),subplot_kw=dict(projection='polar'))
        fig.tight_layout()
        cmax=-1e12
        cmin=1e12
        for i in range(len(yaw_lst)):
            dic=results['TSR'+str(8)+'_yaw'+str(yaw_lst[i])]
            Z = getattr(dic, str(var[j]))
            if var[j]=='f_tan' or var[j]=='f_nor':
                Z=getattr(dic, str(var[j]))/(0.5*param.rho*param.wind_speed**2*param.radius)
            elif var[i]=='circulation':
                Z=getattr(dic, str(var[j]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
            else:
                Z=getattr(dic, str(var[j]))

            if cmax<Z.max():
                cmax=Z.max()
            if cmin>Z.min():
                cmin=Z.min()
            if cmin<0:
                cmin=0
            if var[j]=='a':
                cmax=0.5
            if var[j]=='ap':
                cmax=0.05

        for i in range(len(yaw_lst)):
            dic=results['TSR'+str(8)+'_yaw'+str(yaw_lst[i])]
            axs[i].set_theta_zero_location('E')
            axs[i].set_title('Yaw angle: '+str(yaw_lst[i])+'$^\circ$')
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
                psi=np.append(psi,2*np.pi+psi[0])
                psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
                axs[i].contourf(psi,r , np.tile(Z,len(psi[0])),20,vmin=cmin, vmax=cmax)
            else:

                r = np.hstack((dic.mu,dic.mu[:,0].reshape(len(dic.mu[:,0]),1)))
                psi=dic.azimuth
                psi=np.append(psi,2*np.pi+psi[0])
                psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
                Z = np.hstack((Z,Z[:,0].reshape(len(Z[:,0]),1)))
                axs[i].contourf(psi,r,Z,20,vmin=cmin, vmax=cmax)

            axs[i].set_rlabel_position(225)
            axs[i].tick_params(axis='x', which='major', labelsize=12)
            rlabels = axs[i].get_ymajorticklabels()
            for label in rlabels:
                label.set_color('white')

        m = plt.cm.ScalarMappable(cmap=cm.viridis)
        cbar=plt.colorbar(m, ax=[axs[0:len(yaw_lst)]],orientation='horizontal', boundaries=np.linspace(cmin, cmax, 50))
        cbar.ax.set_xlabel(labels[j], fontsize=16)
        cbar.ax.tick_params(labelsize=14)
        if save==True:
            plt.savefig('figures/polar_'+str(var[j])+'.pdf',bbox_inches = "tight")


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
        plt.figure()
    plt.grid()
    plt.xlabel(r'$\alpha [deg]$')
    plt.ylabel('$C_d [-]$')
    plt.xlim([-20,30])
    plt.plot(dic.polars.alpha,dic.polars.Cd)
    if save==True:
        plt.savefig('figures/Cdalpha.pdf')


def plot_correction(results, param, TSR, Yaw,Res_prandtl,Res_no_prandtl):
    dic=results['TSR'+str(TSR)+'_yaw'+str(Yaw)]
    plt.figure()
    plt.grid()
    plt.plot(dic.mu, dic.f,color='k', label='Prandtl total',zorder=1)
    plt.scatter(dic.mu, dic.f_tip,marker='.',s=80, label='Prandtl tip', zorder=2)
    plt.scatter(dic.mu, dic.f_root,marker='.',s=80, label='Prandtl root', zorder=3)
    plt.xlabel('r/R [-]')
    plt.ylabel('f [-]')
    plt.legend()
    if save==True:
        plt.savefig('figures/prandtl_correction/tip_corrections_'+str(TSR)+'_'+str(Yaw)+'.pdf')

    var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ','local_CT']
    labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$','$C_T [-]$']
    for i in range(len(var)):
        plt.figure()
        plt.grid()
        plt.xlabel(r'Radius $\frac{r}{R}$ [-]')
        plt.ylabel(labels[i])
        if var[i]=='f_tan' or var[i]=='f_nor':
            Z = getattr(Res_prandtl, str(var[i]))/(0.5*param.rho*param.wind_speed**2*param.radius)
            Z_np = getattr(Res_no_prandtl, str(var[i]))/(0.5*param.rho*param.wind_speed**2*param.radius)

        elif var[i]=='circulation':
            Z=getattr(Res_prandtl, str(var[i]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
            Z_np=getattr(Res_no_prandtl, str(var[i]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))

        else:
            Z=getattr(Res_prandtl, str(var[i]))
            Z_np=getattr(Res_no_prandtl, str(var[i]))

        plt.plot(dic.mu,Z,label='Tip and loss correction applied')
        plt.plot(dic.mu,Z_np,label='No tip and loss correction')


        plt.legend()
        if save==True:
            plt.savefig('figures/prandtl_correction/'+str(var[i])+'.pdf')




def plot_enthalpy_tube(results, param, TSR, Yaw):


    dic=results['TSR'+str(TSR)+'_yaw'+str(Yaw)]
    dic2=param

    p_inf=0


    mu1=np.zeros((len(dic.mu),1))
    mu2=dic.mu
    mu3=mu2
    mu4=np.zeros((len(dic.mu),1))
    h1=np.zeros((len(dic.mu),1))
    h3=np.zeros((len(dic.mu),1))
    for i in range(len(mu2)):
        mu1[i]=mu2[i]*np.sqrt((1-dic.a[i]*dic.f[i]))
        mu4[i]=mu2[i]*np.sqrt((1-dic.a[i]*dic.f[i])/(1-2*dic.a[i]*dic.f[i]))
        h1[i]=p_inf/param.rho+0.5*param.wind_speed**2
        h3[i]=p_inf/param.rho+0.5*param.wind_speed**2*(1-2*dic.a[i]*dic.f[i])**2/h1[i]

    h1=h1/h1[0]
    h2=h1
    h4=h3
    Z=np.ones((len(dic.mu),1))


    u = np.linspace(0,  2*np.pi, 100)

    pos=[0,400,700,1200]

    x1 = np.outer(mu1, np.cos(u))
    y1 = np.outer(mu1, np.sin(u))
    z1 = np.zeros((len(mu1),len(u)))

    x2 = np.outer(mu2, np.cos(u))
    y2 = np.outer(mu2, np.sin(u))
    z2 = pos[1]*np.ones((len(mu1),len(u)))

    x3 = np.outer(mu3, np.cos(u))
    y3 = np.outer(mu3, np.sin(u))
    z3 = pos[2]*np.ones((len(mu1),len(u)))

    x4 = np.outer(mu4, np.cos(u))
    y4 = np.outer(mu4, np.sin(u))
    z4 = pos[3]*np.ones((len(mu1),len(u)))

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x1,z1,y1, color='r')
    ax.plot_surface(x2,z2,y2,  color='r')
    scamap = plt.cm.ScalarMappable(cmap='viridis')
    fcolors = scamap.to_rgba(np.tile(h3,100))
    ax.plot_surface(x3,z3,y3, facecolors=fcolors, cmap='viridis',vmin=0, vmax=50)
    ax.plot_surface(x4,z4,y4, facecolors=fcolors, cmap='viridis',vmin=0, vmax=50)

    cb = fig.colorbar(scamap)
    cb.set_label(r'$h_0$ [-]', rotation=0,fontsize=18,labelpad=25)
    cb.ax.tick_params(labelsize=12)
    fake2Dline = mpl.lines.Line2D([0],[0], linestyle="none", c='r', marker = 'o')
    ax.legend([fake2Dline], [r'$h_0$ = 1'], numpoints = 1, fontsize=18)

    ax.set_yticks(pos)
    ax.set_yticklabels(['Infinity Upwind','Rotor Upwind','Rotor Downwind','Infinity Downwind'],rotation=-45,)
    ax.set_xlabel('x', fontsize=16)
    ax.set_zlabel('y', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='y', which='major', labelsize=12)

    ax.view_init(11,-18)
    if save==True:
        plt.savefig('figures/enthalpy1_'+str(TSR)+'_'+str(Yaw)+'.pdf')

    #Ugly plot 1
    plt.figure()
    plt.grid()
    plt.plot(mu1,h1,label='Infinity Upwind')
    plt.plot(mu2,h2,label='Rotor Upwind')
    plt.plot(mu3,h3,label='Rotor Downwind')
    plt.plot(mu4,h4,label='Infinity Downwind')
    plt.xlabel('r/R [-]')
    plt.ylabel('h_0 [-]')
    plt.legend()

    if save==True:
        plt.savefig('figures/enthalpy2_'+str(TSR)+'_'+str(Yaw)+'.pdf')

    #Ugly plot 2
    # plt.grid()
    # plt.plot(Z*0,mu1)
    # plt.plot(Z,mu2)
    # plt.plot(2*Z,mu3)
    # plt.plot(3*Z,mu4)
    # plt.plot(Z*0,-mu1)
    # plt.plot(Z,-mu2)
    # plt.plot(2*Z,-mu3)
    # plt.plot(3*Z,-mu4)


# = = = = = = = = = =   OLD  = = = = = = = = = = = = = #


# def plot_yaw(results,param,yaw_lst):
#     var=['alpha','phi','a','ap','f_tan','f_nor','circulation','local_CQ']
#     labels=[r'$\alpha$ [deg]','$\phi$ [deg]', 'a [-]','$a^,[-]$', '$C_t$ [-]', '$C_n$ [-]','$\Gamma$ [-]','$C_q [-]$']
#     for i in range(len(yaw_lst)):
#         for j in range(len(var)):
#             dic=results['TSR'+str(8)+'_yaw'+str(yaw_lst[i])]
#             fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
#             Z = getattr(dic, str(var[j]))
#             if var[j]=='f_tan' or var[j]=='f_nor':
#                 Z=getattr(dic, str(var[j]))/(0.5*param.rho*param.wind_speed**2*param.radius)
#             elif var[i]=='circulation':
#                 Z=getattr(dic, str(var[j]))/((np.pi*param.wind_speed**2/(param.n_blades*param.omega)))
#             else:
#                 Z=getattr(dic, str(var[j]))
#             if yaw_lst[i]==0:
#                 dic2=results['TSR'+str(8)+'_yaw'+str(yaw_lst[1])]
#                 r = np.hstack((dic2.mu,dic2.mu[:,0].reshape(len(dic2.mu[:,0]),1)))
#                 psi=dic2.azimuth
#                 psi=np.append(psi,2*np.pi+psi[0])
#                 psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
#                 pp=ax.contourf(psi,r , np.tile(Z,len(psi[0])),20)
#                 ax.set_rlabel_position(225)
#             else:
#                 ax.set_theta_zero_location('N') ## CHECK!!!!!!!
#                 r = np.hstack((dic.mu,dic.mu[:,0].reshape(len(dic.mu[:,0]),1)))
#                 psi=dic.azimuth
#                 psi=np.append(psi,2*np.pi+psi[0])
#                 psi=np.tile(psi.transpose(),(len(r[:,0]), 1))
#                 Z = np.hstack((Z,Z[:,0].reshape(len(Z[:,0]),1)))
#                 pp=ax.contourf(psi,r,Z,20)
#                 ax.set_rlabel_position(135)
#             rlabels = ax.get_ymajorticklabels()
#             for label in rlabels:
#                 label.set_color('white')
#             cbar = plt.colorbar(pp, pad=0.1)
#             cbar.ax.set_ylabel(labels[j], rotation=270, labelpad=15)
#             plt.set_cmap('viridis')
#             if save==True:
#                 plt.savefig('figures/Yaw'+str(yaw_lst[i])+'_'+str(var[j])+'.pdf')


# NOT NEEDED now
# Spacing plot, uncomment if you have the put the spacing_data file in the same folder

# def plot_spacing(TSR,Yaw):
#     key='TSR'+str(TSR)+'_yaw'+str(Yaw)
#     plt.figure()
#     plt.grid()
#     points_N=[20,40]
#     marker=['-o','-x']
#     for i in range(len(points_N)):
#         results1 = np.load('spacing_data\dicN'+str(points_N[i])+'.npy',allow_pickle='TRUE').item()
#         dic1=results1[key]
#         plt.plot(dic1.mu, dic1.a, marker[i],markersize=7, label='Constant spacing N='+str(len(dic1.mu)))
#     plt.xlabel('r/R')
#     plt.ylabel('a [-]')
#     plt.legend()
#     if save==True:
#         plt.savefig('figures/spacing1_'+str(TSR)+'_'+str(Yaw)+'.pdf')

#     points_cos=[40,40]
#     plt.figure()
#     plt.grid()
#     results1 = np.load('spacing_data\dicN'+str(points_cos[i])+'.npy',allow_pickle='TRUE').item()
#     results2 = np.load('spacing_data\diccos'+str(points_cos[i])+'.npy',allow_pickle='TRUE').item()
#     dic1=results1[key]
#     dic2=results2[key]
#     plt.plot(dic2.mu, dic2.a,'-o',markersize=5, label='Cosine spacing N='+str(len(dic2.mu)))
#     plt.plot(dic1.mu, dic1.a,'-x',markersize=8, label='Constant spacing N='+str(len(dic2.mu)))
#     plt.xlabel('r/R')
#     plt.ylabel('a [-]')
#     plt.legend()
#     if save==True:
#         plt.savefig('figures/spacing2_'+str(TSR)+'_'+str(Yaw)+'.pdf')


# plot_spacing(8,0)
