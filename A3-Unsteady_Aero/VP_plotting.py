# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 23:53:54 2021

@author: Mathieu Pellé
"""

import matplotlib.pyplot as plt
import numpy as np
from VP_utilities import VOR2D

def steady_polars(alpha, results, rho=1.225, moment=True, flap=False):
    if moment:
        N=2
    else:
        N=1
    for j in range(N):
        C=np.zeros(len(alpha))
        C_theory=np.zeros(len(alpha))
        for i in range(len(alpha)):
            result = results[i]
            c = result['chord']
            if flap:
                c = c+flap['length']
            U_inf = np.linalg.norm(result['velocity'][0])
            dL = rho*U_inf*(result['gamma'][0])
            L = sum(dL)

            if j==0:
                C[i] = L/(0.5*rho*U_inf**2*c)
                if flap:
                    hinge = result['chord']/(result['chord']+flap['length'])
                    theta_k=np.arccos(1-2*hinge)
                    C_theory[i] = 2*np.pi*np.sin(np.deg2rad(alpha[i]+flap['angle']*((1-theta_k/np.pi)+1/np.pi*np.sin(theta_k))))

                else:
                    C_theory[i]=2*np.pi*np.sin(np.deg2rad(alpha[i]))
                label = '$C_l$ [-]'
            else:
                x_lst = []
                for p in range(result['N_panels']):
                    x = result['panels'][p][0,0]-0.5*result['L_panels']
                    x_lst.append(x)
                M = -np.dot(x_lst, dL)
                C[i] = M/(0.5*rho*U_inf**2*result['chord']**2)
                C_theory[i]=-np.pi/2*np.sin(np.deg2rad(alpha[i]))
                label = '$C_m$ [-]'


        plt.figure()
        plt.plot(alpha, C_theory, 'o', label='2D Airfoil theory')
        plt.plot(alpha, C,'-k', label='Vortex panel code')
        plt.grid()
        plt.xlabel(r'$\alpha$ [$^{\circ}$]')
        plt.ylabel(label)
        plt.legend()


def contours(result, rho=1.225, streamlines=True, flap=False, frames=[0]):
    c = result['chord']

    for t in range(len(frames)):
        idx = frames[t]
        TE = result['TE'][idx]
        LE = result['LE'][idx]

        # TE = TE-LE
        # LE = np.array([[0],[0]])

        if flap:
            TE_flap = result['flap_TE']
        u_inf = np.linalg.norm(result['velocity'][idx])

        #Finer mesh - not really working
        # x_pre = np.linspace(-1.5*c,LE[0,0],20)
        # x_aft = np.linspace(TE[0,0],1.5*c,10)
        # x_finer = np.linspace(LE[0,0],TE[0,0],5)
        # x = np.hstack((x_pre, x_finer, x_aft))
        # z_pre = np.linspace(-1.5*c,TE[1,0],20)
        # z_aft = np.linspace(LE[1,0],1.5*c,10)
        # z_finer = np.linspace(TE[1,0],LE[1,0],5)
        # z = np.hstack((z_pre, z_finer, z_aft))

        z = np.linspace(-1*c+LE[0,0],2*c,32)
        x = np.linspace(-1.5*c,1.5*c,32)
        V_mag = np.zeros((len(z), len(x)))
        p = np.zeros((len(z), len(x)))
        U = np.zeros((len(z), len(x)))
        V = np.zeros((len(z), len(x)))
        for i in range(len(z)):
            for j in range(len(x)):
                v = np.array([[0], [0]])
                for k in range(len(result['vortices'][idx])):
                    v = v + VOR2D(np.array([[z[i]], [x[j]]]), result['vortices'][idx][k], result['gamma'][idx][k,0])
                vel = v+result['velocity'][idx]
                V_mag[j, i] = np.sqrt(vel[0]**2+vel[1]**2)
                U[j, i] = vel[0,0]
                V[j, i] = vel[1,0]
                p[j, i] = 0.5*rho*(-V_mag[j, i]**2+u_inf**2)

        plt.figure()
        xx, zz = np.meshgrid(z,x)
        cp = plt.contourf(xx, zz, U/u_inf, 200, cmap='jet')
        cbar = plt.colorbar(cp)
        cbar.ax.set_ylabel('$U/U_{\infty}$ [-]', rotation=270, labelpad=15)
        plt.ylabel('z/c [-]')
        plt.xlabel('x/c [-]')
        if flap:
            plt.plot([TE[0,0], TE_flap[0,0]], [TE[1,0], TE_flap[1,0]], 'k',  linewidth=3.2)
        plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3.2)
        #plt.scatter(zz,xx, marker='+', color='white', alpha=0.5)

        plt.figure()
        cp = plt.contourf(xx, zz, V/u_inf, 200, cmap='jet')
        cbar = plt.colorbar(cp)
        cbar.ax.set_ylabel('$V/U_{\infty}$ [-]', rotation=270, labelpad=15)
        plt.ylabel('z/c [-]')
        plt.xlabel('x/c [-]')
        plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3.2)
        if flap:
            plt.plot([TE[0,0], TE_flap[0,0]], [TE[1,0], TE_flap[1,0]], 'k',  linewidth=3.2)


        if streamlines:
            plt.figure()
            #seed_points = np.linspace((-1.5*c,-1.5*c),(-1.5*c,1.5*c),20)
            plt.streamplot(z,x,U,V)#, start_points=seed_points)
            plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3)
            if flap:
                plt.plot([TE[0,0], TE_flap[0,0]], [TE[1,0], TE_flap[1,0]], 'k',  linewidth=3.2)

        plt.figure()
        cp = plt.contourf(xx, zz, p/(0.5*rho*u_inf**2), 200, cmap='jet')
        cbar = plt.colorbar(cp)
        cbar.ax.set_ylabel('$\Delta q/q_{\infty}$', rotation=270, labelpad=15)
        plt.ylabel('z/c [-]')
        plt.xlabel('x/c [-]')
        plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3.2)
        if flap:
            plt.plot([TE[0,0], TE_flap[0,0]], [TE[1,0], TE_flap[1,0]], 'k',  linewidth=3.2)



def scatter(results, TE=True, LE=True, Vortex=True):
    markers = ['.', 'x', '+']
    plt.figure()
    plt.grid()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    if Vortex:
        for i in range(len(results)):
            result = results[i]
            N = len(result['time'])
            coords = result['vortices'][N-1]
            x = []
            y =[]
            for j in range(len(coords)):
                x.append(coords[j][0,0])
                y.append(coords[j][1,0])

            plt.scatter(x,y,marker=markers[i])

    if TE:
        result = results[0]
        coords = result['TE']
        x = []
        y =[]
        for i in range(len(coords)):
            x.append(coords[i][0,0])
            y.append(coords[i][1,0])
        plt.scatter(x,y, marker='+', color='black')

    if LE:
        result = results[0]
        coords = result['LE']
        x = []
        y =[]
        for i in range(len(coords)):
            x.append(coords[i][0,0])
            y.append(coords[i][1,0])
        plt.scatter(x,y, marker='+', color='purple')

def get_CL(result, rho, theta):
    N_panels = result['N_panels']
    time = result['time']
    dt = time[1]-time[0]
    Cl = np.zeros(len(time))
    Cl_NC = np.zeros(len(time))
    Cl_theory = np.zeros(len(time))
    for t in range(len(time)):
        tan_vec = result['TE'][t]-result['LE'][t]
        tan_vec = tan_vec/np.linalg.norm(tan_vec)
        p = np.zeros(N_panels)
        p_NC = np.zeros(N_panels)
        for j in range(N_panels):
            v=result['velocity'][t]
            for i in range(len(result['vortices'][t])-N_panels):
                i = i+N_panels
                v = v + VOR2D(result['panels_history'][t][j], result['vortices'][t][i], result['gamma'][t][i,0])

            Q = np.dot(tan_vec.T, v)
            if t==0:
                dgamma_dt = 0
            else:
                dgamma_dt = (sum(result['gamma'][t][:N_panels])-sum(result['gamma'][t-1][:N_panels]))/dt
            p[j] = (Q*result['gamma'][t][j,0]/result['L_panels'] + dgamma_dt)*rho
            p_NC[j] = dgamma_dt*rho


        L = sum(p*result['L_panels'])*np.cos(np.deg2rad(theta[t]))
        L_NC = sum(p_NC*result['L_panels'])*np.cos(np.deg2rad(theta[t]))
        Cl[t] = L/(0.5*rho*(np.linalg.norm(result['velocity'][t]))**2*result['chord'])
        Cl_NC[t] = L_NC/(0.5*rho*(np.linalg.norm(result['velocity'][t]))**2*result['chord'])
        Cl_theory[t]=2*np.pi*np.sin(np.deg2rad(theta[t]))
    return Cl, Cl_theory, Cl_NC

def unsteady_polars(theta, result, rho=1.225, quasi=False, inputs=False, arrows=True):

    Cl, Cl_theory, Cl_NC = get_CL(result, rho, theta)

    time = result['time']
    N_panels = result['N_panels']

    idx = int(len(time)/2)
    idx2 = int(3*len(time)/4)

    plt.figure()
    plt.plot(theta[idx:], Cl[idx:], '--b', label='$C_{l_{us}}$')
    plt.plot(theta[idx:], Cl_NC[idx:],'-k', label='$C_{l_{NC}}$')
    plt.plot(theta[idx:], Cl[idx:]-Cl_NC[idx:],'--r', label='$C_{l_{C}}$')
    plt.plot(theta[idx:], Cl_theory[idx:],'k', label='$C_{s}$', alpha=0.7)


    if arrows:
        p1 = (theta[idx+2], Cl[idx+2]+0.2)
        p2 = (theta[idx], Cl[idx]+0.2)
        d1=dict(color='black', shrink=0.05, width=.5, headwidth=3, headlength=4, linewidth=.2)
        plt.annotate('', xy=p1, xytext=p2, arrowprops=d1)
        p1 = (theta[idx2+2], Cl[idx2+2]-0.2)
        p2 = (theta[idx2], Cl[idx2]-0.2)
        d1=dict(color='black', shrink=0.05, width=.5, headwidth=3, headlength=4, linewidth=.2)
        plt.annotate('', xy=p1, xytext=p2, arrowprops=d1)



    if quasi:
        C=np.zeros(len(theta))
        for i in range(len(theta)):
            U_inf = np.linalg.norm(result['velocity'][i])
            dL = rho*U_inf*(result['gamma'][i][:N_panels])
            L = sum(dL)
            C[i] = L/(0.5*rho*U_inf**2*result['chord'])

        plt.plot(theta[idx:], C[idx:], '--g', label='$C_{l_{qs}}$')
    plt.grid()
    plt.tight_layout()
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_l$ [-]')
    plt.legend()

    if inputs:
        plt.figure()
        s = []
        for i in range(len(time)):
            semichord = 2*time[i]*np.linalg.norm(result['velocity'][i])/result['chord']
            s.append(semichord)

        s = s[idx:] - s[idx]

        plt.plot(s, Cl[idx:], '--b', label='$C_l$')
        plt.plot(s, Cl[idx:]-Cl_NC[idx:], '-.r', label='$C_{l_{C}}$')
        plt.plot(s, Cl_NC[idx:], '-k', label='$C_{l_{NC}}$')
        plt.xlim([0,max(s)])
        plt.grid()
        plt.legend()
        plt.xlabel('s [-]')
        plt.ylabel('$C_l$')

def flap_analysis(flap_results, flaps, alpha, parameter, rho=1.225, theory=False):
    Cl = np.zeros((len(flap_results),len(alpha)))
    Cl_theory = np.zeros((len(flap_results),len(alpha)))
    plt.figure()
    for j in range(len(flap_results)):
        results = flap_results[j]
        for i in range(len(alpha)):
            result = results[i]
            c = result['chord']+flaps[j]['length']
            U_inf = np.linalg.norm(result['velocity'])
            dL = rho*U_inf*(result['gamma'][0])
            L = sum(dL)
            Cl[j,i] = L/(0.5*rho*U_inf**2*c)
            hinge = result['chord']/(result['chord']+flaps[j]['length'])
            theta_k=np.arccos(1-2*hinge)
            Cl_theory[j,i] = 2*np.pi*np.sin(np.deg2rad(alpha[i]+flaps[j]['angle']*((1-theta_k/np.pi)+1/np.pi*np.sin(theta_k))))

        if j%1==0:
            if parameter == 'flap_length':
                lab = '$c_{f}/c$='+str(round(flaps[j]['length']/result['chord'],2))
            else:
                lab = r'$\beta$='+str(round(flaps[j]['angle']))+'$^\circ$'

            plt.plot(alpha, Cl[j,:], label=lab)
            if theory:
                if j==len(flap_results)-1:
                    plt.plot(alpha, Cl_theory[j,:],'.k', label='Theoretical')
                else:
                    plt.plot(alpha, Cl_theory[j,:],'.k')
    plt.legend()
    plt.xlabel(r'$\alpha$ [$^\circ$]')
    plt.ylabel('$C_l$ [-]')
    plt.grid()


def step_response(theta, result, step, rho=1.225, idx=0):

    Cl, Cl_theory, Cl_NC = get_CL(result, rho, theta)
    time = result['time']
    s = []
    for i in range(len(time)):
        semichord = 2*time[i]*np.linalg.norm(result['velocity'][i])/result['chord']
        s.append(semichord)

    s = s[idx:] - s[idx]
    if step=='gust':
        Cl_function = (1-0.5*(np.exp(-0.13*s)+np.exp(-s)))*Cl_theory[idx:]
        lab = 'Kussner function'
    else:
        Cl_function = (1-0.165*np.exp(-0.045*s)-0.335*np.exp(-0.3*s))*Cl_theory[idx:]
        lab = 'Wagner function'
    plt.figure()
    plt.plot(s,(Cl[idx:]-Cl_NC[idx:])/max(Cl_theory), label='Vortex Panel Code')
    plt.plot(s,Cl_function/max(Cl_theory), '--k', label=lab)

    plt.grid()
    plt.legend()
    plt.xlabel('s [-]')
    plt.ylabel('$C_l/C_{l_{qs}}$')