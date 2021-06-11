# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 12:28:41 2021

@author: Mathieu Pellé
"""

import math as m
import numpy as np

def VOR2D(point, vortex, gamma):
    """
      Computes induced velocities from vortex in 2D

      Parameters
      ----------
      point: [2x1 array] Point coordinates at which to compute velocities
      vortex: [2x1 array] Vortex element coordinates
      gamma: [float/int] Circulation of vortex

      Returns
      -------
      v : [2x1 array] Induced velocities vector

    """
    r = np.linalg.norm(point-vortex)

    if r<1e-4:
        v = 0
    else:
        v = gamma/(2*m.pi*r**2)*np.matrix([[0, 1], [-1, 0]])*(point-vortex)
    return v

def transform_coords(point, angle, speed, dt):
    """
      Transforms coordinates with speed and rotation of airfoil. Local --> Global

      Parameters
      ----------
      point: [2x1 array] Point coordinates to transform
      angle: [float/int] Rotation angle
      speed: [2x1 array] Speed of translation
      dt: [float/int] Time step

      Returns
      -------
      X : [2x1 array] Transormed coordinates

    """
    T = np.matrix([[m.cos(angle), m.sin(angle)], [-m.sin(angle), m.cos(angle)]])
    X = T*point - speed*dt

    return X

def transform_vel(speed, angle, rotational_speed, x_pos):

    T = np.matrix([[m.cos(angle), -m.sin(angle)], [m.sin(angle), m.cos(angle)]])
    V = T*-speed + rotational_speed*np.array([[0], [x_pos]])

    return V

def influence_matrix(colloc_lst, vortex_lst, normal, N_panels):
    """
      Computes influcence matrix using panel vortices and latest wake vortex

      Parameters
      ----------
      colloc_lst: [list] Collocation points
      vortex_lst: [list] All vortices
      normal: [1x2 array] Normal vector in global reference frame
      N_panels: [float/int] Numer of panels

      Returns
      -------
      a : [matrix] Influence matrix

    """
    if len(vortex_lst)>N_panels:
        a = np.zeros((N_panels, N_panels+1))
        vortex_lst_LHS = vortex_lst[:N_panels]
        vortex_lst_LHS.append(vortex_lst[-1])
    else:
        a = np.zeros((N_panels, N_panels))
        vortex_lst_LHS = vortex_lst


    for i in range(len(colloc_lst)):
        for j in range(len(vortex_lst_LHS)):
            v = VOR2D(colloc_lst[i], vortex_lst_LHS[j], 1)
            a[i, j] = np.dot(normal, v)

    return a

def RHS_vector(colloc_lst, vortex_lst, gamma, N_panels, velocity_vec, normal_vec):
    """
      Computes RHS vector using all wake vortices but latest one

      Parameters
      ----------
      colloc_lst: [list] Collocation points
      vortex_lst: [list] All vortices
      gamma: [list] Circulation of all vortices
      N_panels: [float/int] Numer of panels
      velocity_vec: [2x1 array] Velocity vecotr in local reference frame ???
      normal_vec: [1x2 array] Normal vector in local reference frame

      Returns
      -------
      RHS: [N_panelsx1 array] RHS vector

    """
    RHS = np.zeros((len(colloc_lst),1))
    wake_vortex_lst = vortex_lst[N_panels:-1]
    gamma_wake_vortex = gamma[N_panels:]
    for i in range(len(colloc_lst)):
        for j in range(len(wake_vortex_lst)):
            v = VOR2D(colloc_lst[i], wake_vortex_lst[j], gamma_wake_vortex[j,0])
            velocity_vec[i] = velocity_vec[i] + v

        RHS[i,0] = -np.dot(normal_vec, velocity_vec[i])
    return RHS

def vortex_wake_rollup(vortex_lst, gamma, dt, N_panels):
    """
      Updates position of wake vortices based on induced velocity at each location

      Parameters
      ----------
      vortex_lst: [list] All vortices
      gamma: [list] Circulation of all vortices
      N_panels: [float/int] Numer of panels

      Returns
      -------
      vortex_lst : [list] Updated positinos of vortices

    """
    wake_vortex_lst = vortex_lst[N_panels:]
    for i in range(len(wake_vortex_lst)):
        v = np.zeros((2,1))
        for j in range(len(vortex_lst)):

            v += VOR2D(wake_vortex_lst[i], vortex_lst[j], gamma[j,0])

        wake_vortex_lst[i] = wake_vortex_lst[i] + v*dt
    vortex_lst[N_panels:] = wake_vortex_lst
    return vortex_lst

def vortex_panel(time, N_panels, theta=0, omega=0):

    if theta:
        print('>> Steady, AOA '+str(theta)+' deg')
    elif omega:
        print('>> Unsteady, Oscillation freq. '+str(omega)+' rad/s')
    else:
         print('>> Invalid inputs')



    U_inf = 1
    U_inf_vec = U_inf*np.array([[1], [0]])

    c = 1 #chord
    TE_loc = np.array([[c], [0]]) #Trailing edge location

    #Creating airfoil geometry
    L = c/N_panels
    vortex_lst = []
    colloc_panels = []
    colloc_lst = []
    for n in range(N_panels):
        vortex = np.array([[1/4], [0]])*L+n*L*np.array([[1], [0]]) #Bound vortex
        colloc = np.array([[3/4], [0]])*L+n*L*np.array([[1], [0]]) #Collocation point
        vortex_lst.append(vortex)
        colloc_lst.append(colloc)
        colloc_panels.append(colloc)

    # Local normal vector for uncambered airfoil
    normal_vec = np.array([[0, 1]])

    # List for change in circulation with time
    gamma_lst = []

    # Time loop
    for t in range(len(time)):

        if t==0:
            dt=0
        else:
            dt=time[t]-time[t-1]

        V_origin = -U_inf_vec

        # Pitching rotation and speed
        if theta:
            theta = np.deg2rad(theta)
            theta_dot = 0
        else:
            theta = m.sin(omega*time[t])
            theta_dot = omega*m.cos(omega*time[t])

        normal_vec_global = transform_coords(normal_vec.T, theta, U_inf_vec, dt)
        normal_vec_global = normal_vec_global.T

        #Updated positions of collocations pts and vortices on airfoil only
        for i in range(N_panels):
            colloc_lst[i] = transform_coords(colloc_lst[i], theta, U_inf_vec, dt)

        for j in range(len(vortex_lst)):
            vortex_lst[j] = transform_coords(vortex_lst[j], theta, U_inf_vec, dt)

        #For initial time, no shed vortex
        if t==0:
            a = influence_matrix(colloc_lst, vortex_lst, normal_vec_global, N_panels)
            velocity_vec = []
            for i in range(N_panels):
                velocity_vec.append(transform_vel(V_origin, theta, theta_dot, colloc_panels[i][0][0]))
            f = RHS_vector(colloc_lst, vortex_lst, np.zeros((2,1)), N_panels, velocity_vec, normal_vec)

        #For other time, shed vortices
        else:
            shed_vortex_loc = transform_coords(1.2*TE_loc, theta, U_inf_vec, dt) #position of newly shed vortex
            vortex_lst.append(shed_vortex_loc) #store position

            # influence matrix for each collocation pts and latest shed vortex
            a = influence_matrix(colloc_lst, vortex_lst, normal_vec_global, N_panels)
            a = np.vstack((a, np.ones((1, N_panels+1)))) #Last row of ones for Kelvin equation

            # RHS matrix with kinematic vel. and induced vel. of all but latest shed vortices
            velocity_vec = []
            for i in range(N_panels):
                velocity_vec.append(transform_vel(V_origin, theta, theta_dot, colloc_panels[i][0][0]))
            f = RHS_vector(colloc_lst, vortex_lst, gamma_lst[t-1], N_panels, velocity_vec, normal_vec)
            gamma_previous = np.sum(gamma_lst[t-1])
            f = np.vstack((f, gamma_previous)) #Total circulation at previous time step

        #Solve system
        gamma = np.linalg.inv(np.asmatrix(a))*f #Gamma for collocation points and wake vortices ordered in time
        if t>=2:
            gamma_all = np.vstack((gamma[:N_panels], gamma_lst[t-1][N_panels:]))
            gamma_all = np.vstack((gamma_all, gamma[-1]))
            gamma = gamma_all
        gamma_lst.append(gamma) #Storing for each time step

        #Updating position of wake vortices using induced velocity from all other vortices
        vortex_lst = vortex_wake_rollup(vortex_lst, gamma, dt, N_panels)

    return gamma_lst, vortex_lst


#%%

import matplotlib.pyplot as plt

rho = 1.225
U_inf = 1
c = 1
alpha = np.arange(0,20,1)
Cl = np.zeros(len(alpha))
for i in range(len(alpha)):
    gamma_lst, vortex_lst = vortex_panel([0], 2, theta=alpha[i])

    dL = rho*U_inf*gamma_lst[0]
    L = sum(dL)
    Cl[i] = L/(0.5*rho*U_inf**2*c)


plt.plot(alpha, Cl, label='Vortex panel code')
plt.plot(alpha, 2*np.pi*np.sin(np.deg2rad(alpha)), 'x', label='2D Airfoil theory')
plt.grid()
plt.xlabel(r'$\alpha$ [$^{\circ}$]')
plt.ylabel('$C_l$ [-]')
plt.legend()

#%%


U_inf_vec = np.array([[1], [0]])
alpha = 5
LE = transform_coords( np.array([[0], [0]]), np.deg2rad(alpha), 0, 0)
TE = transform_coords( np.array([[c], [0]]), np.deg2rad(alpha), 0, 0)

gamma_lst, vortex_lst = vortex_panel([0], 2, theta=alpha)


x = np.linspace(-1.5*c, 1.5*c, 50)
y = np.linspace(-1.5*c, 1.5*c, 50)
V_mag = np.zeros((len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        v = np.array([[0], [0]])
        for k in range(len(vortex_lst)):
            v = v + VOR2D(np.array([[x[i]], [y[j]]]), vortex_lst[k], gamma_lst[0][k,0])
        vel = v+U_inf_vec
        V_mag[i, j] = np.sqrt(vel[0]**2+vel[1]**2)

xx, yy = np.meshgrid(x,y)
cp = plt.contourf(xx, yy, V_mag)
cbar = plt.colorbar(cp)
cbar.ax.set_ylabel('Velocity magnitude [m/s]', rotation=270, labelpad=15)
plt.ylabel('y [m]')
plt.ylabel('x [m]')
plt.plot([LE[0,0], TE[0,0]], [LE[1,0], TE[1,0]], 'k')
plt.grid()
