# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 12:39:17 2021

@author: Mathieu Pellé
"""

import numpy as np
import math as m

def SingleFilament(ctrl_point,point1,point2,gamma):

    [xp,yp,zp] = ctrl_point
    [x1,y1,z1] = point1
    [x2,y2,z2] = point2

    lim = 1e-5

    R1=m.sqrt((xp-x1)**2+(yp-y1)**2+(zp-z1)**2)
    R2=m.sqrt((xp-x2)**2+(yp-y2)**2+(zp-z2)**2)

    R12x=(yp-y1)*(zp-z2)-(zp-z1)*(yp-y2)
    R12y=-(xp-x1)*(zp-z2)+(zp-z1)*(xp-x2)
    R12z=(xp-x1)*(yp-y2)-(yp-y1)*(xp-x2)

    R12sq=R12x**2+R12y**2+R12z**2
    R01=(x2-x1)*(xp-x1)+(y2-y1)*(yp-y1)+(z2-z1)*(zp-z1)
    R02=(x2-x1)*(xp-x2)+(y2-y1)*(yp-y2)+(z2-z1)*(zp-z2)

    if R12sq<lim:
        R12sq=lim**2;

    if R1<lim:
        R1=lim

    if R2<lim:
        R2=lim

    K=gamma/(4*m.pi*R12sq)*(R01/R1-R02/R2)

    U=K*R12x
    V=K*R12y
    W=K*R12z

    return [U,V,W]


class Geometry:
  def __init__(self, s_lst, alpha, L_w):


     self.alpha = alpha
     self.L_W = L_w
     self.control_pts = np.zeros(len(s_lst)-1)
     self.HS_vortex = []
     for i in range(len(s_lst)-1):
         self.control_pts[i] = (s_lst[i+1]-s_lst[i])/2 + s_lst[i] # Control points position

         # Position and circulation of 5 filaments making up horseshoe vortex
         VF1 = {'pos1':[L_w,s_lst[i],L_w*m.sin(alpha)], 'pos2':[1.25,s_lst[i],0], 'Gamma':1}
         VF2 = {'pos1':[1.25,s_lst[i],0], 'pos2':[0,s_lst[i],0], 'Gamma':1}
         VF3 = {'pos1':[0,s_lst[i],0], 'pos2':[0,s_lst[i+1],0], 'Gamma':1}
         VF4 = {'pos1':[0,s_lst[i+1],0], 'pos2':[1.25,s_lst[i+1],0], 'Gamma':1}
         VF5 = {'pos1':[1.25,s_lst[i+1],0], 'pos2':[L_w,s_lst[i+1],L_w*m.sin(alpha)], 'Gamma':1}
         self.HS_vortex.append({'VF1':VF1,'VF2':VF2,'VF3':VF3,'VF4':VF4,'VF5':VF5})


def InducedVelocities(geometry):
    N = len(geometry.control_pts)
    u_ind_mat = np.zeros((N,N))
    v_ind_mat = np.zeros((N,N))
    w_ind_mat = np.zeros((N,N))
    # For every control point
    for i in range(N):
        x_cp = geometry.control_pts[i]

        # For every HS vortex
        for j in range(N):
            hs = geometry.HS_vortex[j]

            # For ever filament of HS vortex
            u_ind = 0
            v_ind = 0
            w_ind = 0
            for k in range(5):
                fil = hs['VF'+str(k+1)]

                # Compute induced velocity
                w = SingleFilament([0,x_cp,0],fil['pos1'],fil['pos2'],fil['Gamma'])
                u_ind+=w[0]
                v_ind+=w[1]
                w_ind+=w[2]

            # Store induced velocities per control point and HS vortex
            u_ind_mat[i,j] = u_ind
            v_ind_mat[i,j] = v_ind
            w_ind_mat[i,j] = w_ind

    return u_ind_mat, v_ind_mat, w_ind_mat

def LiftingLine(u_ind_mat, v_ind_mat, w_ind_mat, geometry):

    alpha = geometry.alpha
    N = len(u_ind_mat)
    gamma_new = np.zeros((N))
    gamma = np.zeros((N))
    CL_lst = np.zeros((N))

    #Iteration controls
    it = 0
    it_max =1000
    error = 1e12
    limit = 1e-4
    UR = 0.3
    while it<it_max and error>limit:
        for i in range(N):
            u = 0
            v = 0
            w = 0
            for j in range(N):
                u = u + u_ind_mat[i,j]*gamma[j]
                v = v + v_ind_mat[i,j]*gamma[j]
                w = w + w_ind_mat[i,j]*gamma[j]

            V = [1 + u , v, m.sin(alpha) + w]
            phi = m.atan(V[2]/V[0])
            CL_lst[i] = 2*m.pi*m.sin(phi)
            Vmag = m.sqrt(np.dot(V,V))
            gamma_new[i] = 0.5*Vmag*CL_lst[i]

        error=max(abs(gamma_new-gamma))
        UR = 0.2

        ## Carlos' (weird) way to control the iteration
        # referror=max(abs(gamma_new))
        # referror=max(referror,0.001)
        # error=max(abs(gamma_new-gamma))
        # diff = error
        # error = error/referror
        # UR = max((1-error)*0.3,0.1)

        gamma = UR*gamma_new + (1-UR)*gamma
        it+=1

    return CL_lst

AR = 5
s_lst = np.linspace(0,np.pi,10)
L_w = 1000*AR
alpha = np.deg2rad(5)

geometry = Geometry(s_lst,alpha,L_w)
[u_ind_mat, v_ind_mat, w_ind_mat] = InducedVelocities(geometry)
CL_lst = LiftingLine(u_ind_mat, v_ind_mat, w_ind_mat, geometry)











