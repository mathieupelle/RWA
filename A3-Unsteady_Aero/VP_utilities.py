# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 23:54:43 2021

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

    if r<1e-8:
        v = 0
    else:
        v = gamma/(2*m.pi*r**2)*np.matrix([[0, 1], [-1, 0]])*(point-vortex)
    return v

def transform_coords(point, angle, LE):
    """
      Transforms coordinates with speed and rotation of aerofoil. Local --> Global

      Parameters
      ----------
      point: [2x1 array] Point coordinates to transform
      angle: [float/int] Rotation angle
      speed: [2x1 array] Speed of translation
      dt: [float/int] Time step

      Returns
      -------
      X : [2x1 array] Transformed coordinates

    """
    T = np.matrix([[m.cos(angle), m.sin(angle)], [-m.sin(angle), m.cos(angle)]])
    X = T*point + LE

    return X

def transform_vel(speed, angle, rotational_speed, x_pos):

    T = np.matrix([[m.cos(angle), -m.sin(angle)], [m.sin(angle), m.cos(angle)]])
    V = T*-speed + rotational_speed*np.array([[0], [x_pos]])

    return V

def influence_matrix(colloc_lst, vortex_lst, normal, N_panels, flap=False):
    """
      Computes influence matrix using panel vortices and latest wake vortex

      Parameters
      ----------
      colloc_lst: [list] Collocation points
      vortex_lst: [list] All vortices
      normal: [1x2 array] Normal vector in global reference frame
      N_panels: [float/int] Number of panels

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
            if flap:
                if i>N_panels-flap['N_panels']-1:
                    normal = flap['norm_vec_global']
            a[i, j] = np.dot(normal, v)

    return a

def RHS_vector(colloc_lst, vortex_lst, gamma, N_panels, velocity_vec, normal_vec, flap=False):
    """
      Computes RHS vector using all wake vortices but latest one

      Parameters
      ----------
      colloc_lst: [list] Collocation points
      vortex_lst: [list] All vortices
      gamma: [list] Circulation of all vortices
      N_panels: [float/int] Number of panels
      velocity_vec: [2x1 array] Velocity vector in local reference frame ???
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

        if flap:
            if i>N_panels-flap['N_panels']-1:
                normal_vec = flap['norm_vec']

        RHS[i,0] = -np.dot(normal_vec, velocity_vec[i])
    return RHS

def vortex_wake_rollup(vortex_lst, gamma, dt, N_panels):
    """
      Updates position of wake vortices based on induced velocity at each location

      Parameters
      ----------
      vortex_lst: [list] All vortices
      gamma: [list] Circulation of all vortices
      N_panels: [float/int] Number of panels

      Returns
      -------
      vortex_lst : [list] Updated positions of vortices

    """
    wake_vortex_lst = vortex_lst[N_panels:]
    for i in range(len(wake_vortex_lst)):
        v = np.zeros((2,1))
        for j in range(len(vortex_lst)):
            v += VOR2D(wake_vortex_lst[i], vortex_lst[j], gamma[j,0])

        wake_vortex_lst[i] = wake_vortex_lst[i] + v*dt
    lst = vortex_lst[:N_panels]+wake_vortex_lst
    return lst


def vortex_panel(time, N_panels, theta, theta_dot, c=1, U_inf_vec=[np.array([[1],[0]])], flap=False):
    if len(theta)==1:
        mode = 'steady'
    elif len(theta)>1:
        mode = 'unsteady'
    else:
         print('>> Invalid inputs')

    TE_loc = np.array([[c], [0]]) #Trailing edge
    LE_loc = np.array([[0],[0]]) #Leading edge
    theta = np.deg2rad(theta) #angle of attack
    normal_vec = np.array([[0, 1]]) #Normal vector for uncambered aerofoil (local frame)

    #Creating panels
    L = c/N_panels #Panel length

    colloc_panels = [] #Original collocation points coordinates
    vortex_panels = [] #Original vortex coordinates
    colloc_lst = [] #List of collocation points in local reference frame
    vortex_lst = [] #List of vortices in local reference frame
    vortex_history = []
    colloc_history = []

    #Loop over number of panels
    for n in range(N_panels):
        vortex = np.array([[1/4], [0]])*L+n*L*np.array([[1], [0]]) #Bound vortices
        colloc = np.array([[3/4], [0]])*L+n*L*np.array([[1], [0]]) #Collocation points
        vortex_lst.append(vortex)
        vortex_panels.append(vortex)
        colloc_lst.append(colloc)
        colloc_panels.append(colloc)


    if flap:
        flap_angle = np.deg2rad(flap['angle'])
        L = flap['length']/flap['N_panels']
        if flap['length']==0:
            flap['N_panels'] = 0
        N_panels = flap['N_panels']+ N_panels
        for n in range(flap['N_panels']):
            vortex = np.array([[1/4], [0]])*L+n*L*np.array([[1], [0]])
            vortex = transform_coords(vortex, flap_angle, TE_loc)
            colloc = np.array([[3/4], [0]])*L+n*L*np.array([[1], [0]])
            colloc = transform_coords(colloc, flap_angle, TE_loc)

            vortex_lst.append(vortex)
            vortex_panels.append(np.asarray(vortex))
            colloc_lst.append(colloc)
            colloc_panels.append(np.asarray(colloc))


    gamma_lst = []
    LE_loc_lst = [] #Tracking leading edge location
    TE_loc_lst = [] #Tracking trailing edge location

    #Steady code
    if mode=='steady':
        print('>> Steady, AOA '+str(round(np.rad2deg(theta[0]),2))+' deg')

        normal_vec_global = transform_coords(normal_vec.T, theta[0], np.zeros((2,1))) #Transforming normal vector
        normal_vec_global = normal_vec_global.T #Transposing

        V_origin = -U_inf_vec[0]

        if flap:
            flap_normal_vec = transform_coords(normal_vec.T, flap_angle, np.zeros((2,1))) #Transforming normal vector
            flap_normal_vec_global = transform_coords(normal_vec.T, theta[0]+flap_angle, np.zeros((2,1))) #Transforming normal vector
            flap['norm_vec']=flap_normal_vec.T
            flap['norm_vec_global']=flap_normal_vec_global.T


        TE_loc_T = transform_coords(TE_loc, theta[0], LE_loc) #Transforming trailing edge location
        TE_loc_lst.append(TE_loc_T) #Storing trailing edge position
        LE_loc_lst.append(LE_loc) #Storing leading edge position

        if flap:
            flap_TE_loc = np.array([[flap['length']],[0]])
            flap_TE_loc_T = transform_coords(flap_TE_loc, theta[0]+flap_angle, TE_loc_T)

        #Updated positions of collocations pts and vortices on aerofoil only
        for i in range(N_panels):
            colloc_lst[i] = transform_coords(colloc_panels[i], theta[0], LE_loc)
            vortex_lst[i] = transform_coords(vortex_panels[i], theta[0], LE_loc)

        a = influence_matrix(colloc_lst, vortex_lst, normal_vec_global, N_panels, flap=flap) #Influence matrix
        velocity_vec = []
        #Computing velocity vector at each collocation point
        for i in range(N_panels):
            velocity_vec.append(transform_vel(V_origin, theta[0], theta_dot[0], colloc_panels[i][0][0]))

        f = RHS_vector(colloc_lst, vortex_lst, np.zeros((N_panels,1)), N_panels, velocity_vec, normal_vec, flap=flap) #RHS vector
        gamma = np.linalg.inv(np.asmatrix(a))*f #solving system
        gamma_lst.append(gamma) #storing circulation
        vortex_history.append(vortex_lst)
        colloc_history.append(colloc_lst)

    #Unsteady code
    else:
        if theta[-1]==theta[-2]:
            print('>> Unsteady, Step')
        else:
            print('>> Unsteady, Oscillations')

        for t in range(len(time)):
            if t==0:
                dt=0
            else:
                dt=time[2]-time[1]

            normal_vec_global = transform_coords(normal_vec.T, theta[t], np.zeros((2,1))) #Transforming normal vector
            normal_vec_global = normal_vec_global.T #Transposing

            V_origin = -U_inf_vec[t]

            LE_loc = LE_loc - U_inf_vec[t]*dt #Updating leading edge (origin) position. Aerofoil moves left.
            TE_loc_T = transform_coords(TE_loc, theta[t], LE_loc) #Transforming trailing edge location and shifting based on LE.
            TE_loc_lst.append(TE_loc_T) #Storing trailing edge position
            LE_loc_lst.append(LE_loc) #Storing leading edge position

            #shed_loc = TE_loc_T+np.array([[0.2],[0]])*c #shedding/new vortex location
            if t==0:
                shed_loc = TE_loc_T+np.array([[0.25],[0]])*c #shedding/new vortex location
            else:
                shed_loc = TE_loc_T-(TE_loc_T-TE_loc_lst[t-1])*0.25 #shedding/new vortex location
            vortex_lst.append(shed_loc) #Storing shedding/new vortex position

            #Updated positions of collocations pts and vortices on aerofoil only
            for i in range(N_panels):
                colloc_lst[i] = transform_coords(colloc_panels[i], theta[t], LE_loc)
                vortex_lst[i] = transform_coords(vortex_panels[i], theta[t], LE_loc)

            # influence matrix for each collocation pts and latest shed vortex
            a = influence_matrix(colloc_lst, vortex_lst, normal_vec_global, N_panels)
            a = np.vstack((a, np.ones((1, N_panels+1)))) #Last row of ones for Kelvin equation

            # RHS matrix with kinematic vel. and induced vel. of all but latest shed vortices
            velocity_vec = []
            for i in range(N_panels):
                velocity_vec.append(transform_vel(V_origin, theta[t], theta_dot[t], colloc_panels[i][0][0]))

            if t==0:
                gamma_temp = [np.array([[0],[0]])]*2 #Using gamma=0 for initial time step
                f = RHS_vector(colloc_lst, vortex_lst, gamma_temp, N_panels, velocity_vec, normal_vec)
                gamma_previous = 0
            else:
                f = RHS_vector(colloc_lst, vortex_lst, gamma_lst[t-1], N_panels, velocity_vec, normal_vec)
                gamma_previous = -np.sum(gamma_lst[t-1][N_panels:]) #Total circulation at previous time step
            f = np.vstack((f, gamma_previous)) #RHS vector

            gamma = np.linalg.inv(np.asmatrix(a))*f #solving system

            #Re-arranging circulation --> [Bound vortices, shed vortices from earliest to latest]
            if t>=1:
                gamma_all = np.vstack((gamma[:N_panels], gamma_lst[t-1][N_panels:]))
                gamma_all = np.vstack((gamma_all, gamma[-1]))
                gamma = gamma_all

            gamma_lst.append(gamma) #Storing circulation
            vortex_lst = vortex_wake_rollup(vortex_lst, gamma, dt, N_panels) #Shifting wake vortices positions
            vortex_history.append(vortex_lst[:])
            colloc_history.append(colloc_lst[:])

    results = {'N_panels':N_panels, 'L_panels':L, 'chord':c, 'velocity':U_inf_vec, 'panels': colloc_panels,
               'vortices':vortex_history, 'gamma': gamma_lst, 'LE':LE_loc_lst, 'TE':TE_loc_lst,
               'time': time, 'panels_history': colloc_history}

    if flap:
        results['flap_TE']=flap_TE_loc_T
    return results