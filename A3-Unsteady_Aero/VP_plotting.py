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
            U_inf = np.linalg.norm(result['velocity'])
            dL = rho*U_inf*(result['gamma'][0]) #CHECK!!!!
            L = sum(dL)

            if j==0:
                C[i] = L/(0.5*rho*U_inf**2*c)
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
        plt.plot(alpha, C, label='Vortex panel code')
        plt.plot(alpha, C_theory, 'x', label='2D Airfoil theory')
        plt.grid()
        plt.xlabel(r'$\alpha$ [$^{\circ}$]')
        plt.ylabel(label)
        plt.legend()


def steady_contours(alpha, result, rho=1.225, quiver=True, streamlines=True):
    c = result['chord']
    TE = result['TE'][0]
    LE = result['LE'][0]

    #For contours
    x = np.linspace(-1.5*c,1.5*c,20)
    xx, zz = np.meshgrid(x,x)
    V_mag = np.zeros((len(x), len(x)))
    p = np.zeros((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            v = np.array([[0], [0]])
            for k in range(len(result['vortices'][0])):
                v = v + VOR2D(np.array([[x[i]], [x[j]]]), result['vortices'][0][k], result['gamma'][0][k,0])
            vel = v+result['velocity']
            V_mag[i, j] = np.sqrt(vel[0]**2+vel[1]**2)
            p[i, j] = 0.5*rho*(-V_mag[i, j]**2+np.linalg.norm(result['velocity'])**2)

    plt.figure()
    cp = plt.contourf(xx, zz, V_mag, 150, cmap='jet')

    #For quiver
    x = np.linspace(-1.5*c, 1.5*c, 20)
    Vx = np.zeros((len(x), len(x)))
    Vz = np.zeros((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            v = np.array([[0], [0]])
            for k in range(len(result['vortices'][0])):
                v = v + VOR2D(np.array([[x[i]], [x[j]]]), result['vortices'][0][k], result['gamma'][0][k,0])
            vel = v+result['velocity']
            Vx[i, j] = vel[0,0]
            Vz[i, j] = vel[1,0]
    if quiver:
        plt.quiver(x,x,Vx,Vz, color='black')
    cbar = plt.colorbar(cp)
    cbar.ax.set_ylabel('Velocity magnitude [m/s]', rotation=270, labelpad=15)
    plt.ylabel('z [m]')
    plt.xlabel('x [m]')
    plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3, alpha=0.6)
    #plt.scatter(zz,xx, marker='+')

    if streamlines:
        plt.figure()
        #seed_points = np.linspace((-1.5*c,-1.5*c),(-1.5*c,1.5*c),20)
        plt.streamplot(x,x,Vx,Vz)#, start_points=seed_points)
        plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3, alpha=0.6)

    plt.figure()
    cp = plt.contourf(xx, zz, p, 150, cmap='jet')
    cbar = plt.colorbar(cp)
    cbar.ax.set_ylabel('Dynamic pressure [Pa]', rotation=270, labelpad=15)
    plt.ylabel('z [m]')
    plt.xlabel('x [m]')
    plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3, alpha=0.6)


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
    Cl_theory = np.zeros(len(time))
    for t in range(len(time)):
        tan_vec = result['TE'][t]-result['LE'][t]
        tan_vec = tan_vec/np.linalg.norm(tan_vec)
        p = np.zeros(N_panels)
        for j in range(N_panels):
            v=result['velocity']
            for i in range(len(result['vortices'][t])-N_panels):
                i = i+N_panels
                v = v + VOR2D(result['panels_history'][t][j], result['vortices'][t][i], result['gamma'][t][i,0])

            Q = np.dot(tan_vec.T, v)
            if t==0:
                dgamma_dt = 0 # ??? check
            else:
                dgamma_dt = (sum(result['gamma'][t][:N_panels])-sum(result['gamma'][t-1][:N_panels]))/dt
            p[j] = (Q*result['gamma'][t][j,0]/result['L_panels'] + dgamma_dt)*rho


        L = sum(p*result['L_panels'])*np.cos(np.deg2rad(theta[t]))
        Cl[t] = L/(0.5*rho*(np.linalg.norm(result['velocity']))**2*result['chord'])
        Cl_theory[t]=2*np.pi*np.sin(np.deg2rad(theta[t]))
    return Cl, Cl_theory

def unsteady_polars(theta, result, rho=1.225, basic_method=False, inputs=False, arrows=True):

    Cl, Cl_theory = get_CL(result, rho, theta)

    time = result['time']
    N_panels = result['N_panels']

    idx = int(len(time)/2)
    idx2 = int(3*len(time)/4)

    plt.figure()
    plt.plot(theta[idx:], Cl[idx:], marker='.')
    #plt.plot(theta[idx:], Cl_theory[idx:], '--k', label='US', alpha=0.6)
    if arrows:
        p1 = (theta[idx+2], Cl[idx+2]+0.2)
        p2 = (theta[idx], Cl[idx]+0.2)
        d1=dict(color='black', shrink=0.05, width=.5, headwidth=3, headlength=4, linewidth=.2)
        plt.annotate('', xy=p1, xytext=p2, arrowprops=d1)
        p1 = (theta[idx2+2], Cl[idx2+2]-0.2)
        p2 = (theta[idx2], Cl[idx2]-0.2)
        d1=dict(color='black', shrink=0.05, width=.5, headwidth=3, headlength=4, linewidth=.2)
        plt.annotate('', xy=p1, xytext=p2, arrowprops=d1)

    plt.grid()
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_l$ [-]')

    if basic_method:
        C=np.zeros(len(theta))
        for i in range(len(theta)):
            U_inf = np.linalg.norm(result['velocity'])
            dL = rho*U_inf*(result['gamma'][i][:N_panels])
            L = sum(dL)
            C[i] = L/(0.5*rho*U_inf**2*result['chord'])

        plt.scatter(theta[idx:], C[idx:], marker='.', color='black', label='S')
        plt.legend()

    if inputs:
        plt.figure()
        theta_norm = theta/max(theta)
        cl_norm = Cl/max(Cl)
        plt.plot(result['time'][idx:], cl_norm[idx:], label='$C_l$')
        plt.plot(result['time'][idx:], theta_norm[idx:], label=r'$\alpha$')
        plt.grid()
        plt.legend()
        plt.xlabel('time [s]')
        plt.ylabel(r'$C_l/C_{l_{max}}$, $\alpha$/$\alpha_{max}$ [-]')

