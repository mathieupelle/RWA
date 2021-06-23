# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 23:53:54 2021

@author: Mathieu Pellé
"""

import matplotlib.pyplot as plt
import numpy as np
from VP_utilities import VOR2D

def polars(alpha, results, rho=1.225, mode='steady'):

    if mode=='steady':
        for j in range(2):
            C=np.zeros(len(alpha))
            C_theory=np.zeros(len(alpha))
            for i in range(len(alpha)):
                result = results[i]
                U_inf = np.linalg.norm(result['velocity'])
                dL = rho*U_inf*-(result['gamma'][0]) #CHECK!!!!
                L = sum(dL)

                if j==0:
                    C[i] = L/(0.5*rho*U_inf**2*result['chord'])
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


def steady_contours(alpha, result, rho=1.225, quiver=True):
    c = result['chord']
    TE = result['TE'][0]
    LE = result['LE'][0]

    x1 = np.linspace(-0.4*c, 0.4*c, 20)
    x2 = np.linspace(-1.5*c, -0.4*c, 15)
    x3 = np.linspace(0.4*c, 1.5*c, 15)
    x = np.hstack((x2,x1,x3))
    x = np.linspace(-1.5*c,1.5*c,30)
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
            p[i, j] = 0.5*rho*V_mag[i, j]**2

    plt.figure()
    cp = plt.contourf(xx, zz, V_mag, cmap='jet')

    if quiver:
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

        plt.quiver(x,x,Vx,Vz, color='black')
    cbar = plt.colorbar(cp)
    cbar.ax.set_ylabel('Velocity magnitude [m/s]', rotation=270, labelpad=15)
    plt.ylabel('z [m]')
    plt.xlabel('x [m]')
    plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k',  linewidth=3, alpha=0.6)
    #plt.scatter(zz,xx, marker='+')

    plt.figure()
    cp = plt.contourf(xx, zz, p, cmap='jet')
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