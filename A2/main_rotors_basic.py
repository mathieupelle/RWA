# -*- coding: utf-8 -*-
"""
Created on Tue May  4 10:17:22 2021

@author: Mathieu Pellé
"""

import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt


class Rotor:
  def __init__(self, Radius, TipSpeedRatio, NumberofBlades, GlobalPitch, WindSpeedVector):

     self.R = Radius
     self.TSR = TipSpeedRatio
     self.N_blades = NumberofBlades
     self.theta_p = GlobalPitch
     self.V_inf = WindSpeedVector
     self.omega = self.TSR*self.V_inf[0]/self.R
     self.polars = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])

  def Polars(self,alpha):
      cl = np.interp(alpha,np.array(self.polars['alpha']),np.array(self.polars['Cl']))
      cd = np.interp(alpha,np.array(self.polars['alpha']),np.array(self.polars['Cd']))
      return cl, cd

  def BladeGeometry(self, r):

      c = 3*(1-r/self.R)+1
      beta = -14*(1-r/self.R)
      theta = np.deg2rad(self.theta_p + beta)

      return c, theta

class VortexGeometry:
    def __init__(self, rotors, NumberofRadiusPositions, Spacing, NumberofRotations, RotorLocations):

        N_WT = len(rotors)
        self.D_lst = RotorLocations
        self.WT_lst = {}

        for t in range(N_WT):
            self.WT_lst['WT'+str(t+1)] = {'r_lst':[], 'psi_lst':[],
                                          'N_pts':0, 'N_hs':0, 'N_rot':0, 'N_rad':0,
                                          'blades':{'control_pts':[], 'HS_vortex':[], 'panels':[]}}

            rotor = rotors[t]
            N_blades = rotor.N_blades
            TSR = rotor.TSR
            R = rotor.R
            N_rad = NumberofRadiusPositions[t]
            D = self.D_lst[t]

            mu_start = 0.2
            if Spacing == 'cosine':
                s_lst = np.zeros(N_rad)
                angle=np.linspace(0,np.pi,N_rad)
                for i in range(N_rad):
                    s_lst[i]=(0.4*(1-m.cos(angle[i]))+mu_start)*R
            else:
                s_lst = np.linspace(mu_start*R,R,N_rad)

            N_rot = NumberofRotations[t]
            psi_lst = np.linspace(0,2*m.pi*N_rot,2*N_rot*10+1)

            pts_lst = []
            hs_lst = []
            panel_lst = []
            r_lst = []

            #For each blade
            for b in range(N_blades):
                theta_r = 2*m.pi*b/N_blades #Rotation angle for coordinate transform

                #For each element
                for i in range(N_rad-1):
                    fil = []
                    ## Control points ##
                    r = (s_lst[i+1]+s_lst[i])/2 #Radial position
                    r_lst.append(r)
                    [c, theta] = rotor.BladeGeometry(r) #Chord and pitch
                    vec_n = [m.cos(theta),0,-m.sin(theta)]
                    vec_n = VortexGeometry.Transform(vec_n, theta_r, D) #Transformed normal vector
                    vec_t = [-m.sin(theta),0,-m.cos(theta)]
                    vec_t = VortexGeometry.Transform(vec_t, theta_r, D) #Transformed tangential vector
                    pts_lst.append({'coords': VortexGeometry.Transform([0,r,0],theta_r,D), 'chord': c, 'normal': vec_n, 'tangent': vec_t})

                    ## Horseshoe vortices ##
                    #VF1: filament in span direction at 0.25c
                    VF = {'pos1':[0,s_lst[i],0], 'pos2':[0,s_lst[i+1],0], 'Gamma':1}
                    fil.append(VF)

                    #VF2: filament in upstream direction at 1st element position
                    [c, theta] = rotor.BladeGeometry(s_lst[i])
                    VF = {'pos1':[c*m.sin(-theta),s_lst[i],-c*m.cos(theta)], 'pos2':[0,s_lst[i],0], 'Gamma':1}
                    fil.append(VF)

                    #All filaments downstream of VF2 up to limit defined by number of rotations
                    for j in range(len(psi_lst)-1):
                        x = fil[-1]['pos1'][0]
                        y = fil[-1]['pos1'][1]
                        z = fil[-1]['pos1'][2]
                        dx = (psi_lst[j+1]-psi_lst[j])/TSR*R
                        dy = (m.cos(-psi_lst[j+1])-m.cos(-psi_lst[j]))*s_lst[i]
                        dz = (m.sin(-psi_lst[j+1])-m.sin(-psi_lst[j]))*s_lst[i]
                        VF = {'pos1':[x+dx,y+dy,z+dz], 'pos2':[x,y,z], 'Gamma':1}
                        fil.append(VF)

                    #VF3: filament in downstream direction at 2nd element position
                    [c, theta] = rotor.BladeGeometry(s_lst[i+1])
                    VF = {'pos1':[0,s_lst[i+1],0], 'pos2':[c*m.sin(-theta),s_lst[i+1],-c*m.cos(theta)], 'Gamma':1}
                    fil.append(VF)
                    #All filaments downstream of VF3 up to limit defined by number of rotations
                    for j in range(len(psi_lst)-1):
                        x = fil[-1]['pos2'][0]
                        y = fil[-1]['pos2'][1]
                        z = fil[-1]['pos2'][2]
                        dx = (psi_lst[j+1]-psi_lst[j])/TSR*R
                        dy = (m.cos(-psi_lst[j+1])-m.cos(-psi_lst[j]))*s_lst[i+1]
                        dz = (m.sin(-psi_lst[j+1])-m.sin(-psi_lst[j]))*s_lst[i+1]
                        VF = {'pos1':[x,y,z], 'pos2':[x+dx,y+dy,z+dz], 'Gamma':1}
                        fil.append(VF)

                    #Transforming all filaments coordinates
                    for j in range(len(fil)):
                        fil[j]['pos1'] = VortexGeometry.Transform(fil[j]['pos1'], theta_r, D)
                        fil[j]['pos2'] = VortexGeometry.Transform(fil[j]['pos2'], theta_r, D)

                    ## Panels ##
                    [c1, theta1] = rotor.BladeGeometry(s_lst[i])
                    [c2, theta2] = rotor.BladeGeometry(s_lst[i+1])
                    p1 = VortexGeometry.Transform([-0.25*c1*m.sin(-theta1)*0 , s_lst[i], 0.25*c1*m.cos(theta1)],theta_r,D)
                    p2 = VortexGeometry.Transform([-0.25*c2*m.sin(-theta2)*0 , s_lst[i+1], 0.25*c2*m.cos(theta2)],theta_r,D)
                    p3 = VortexGeometry.Transform([0.75*c2*m.sin(-theta2)*0 , s_lst[i+1], -0.75*c2*m.cos(theta2)],theta_r,D)
                    p4 = VortexGeometry.Transform([0.75*c1*m.sin(-theta1)*0 , s_lst[i], -0.75*c1*m.cos(theta1)],theta_r,D)

                    panel_lst.append({'p1':p1, 'p2':p2, 'p3':p3, 'p4':p4})
                    hs_lst.append(fil)

            #Storing all data
            self.WT_lst['WT'+str(t+1)]['r_lst'] = r_lst
            self.WT_lst['WT'+str(t+1)]['psi_lst'] = psi_lst
            self.WT_lst['WT'+str(t+1)]['N_pts'] = len(pts_lst)
            self.WT_lst['WT'+str(t+1)]['N_hs'] = len(hs_lst)
            self.WT_lst['WT'+str(t+1)]['N_rot'] = N_rot
            self.WT_lst['WT'+str(t+1)]['N_rad'] = N_rad
            self.WT_lst['WT'+str(t+1)]['blades']['control_pts'] = pts_lst
            self.WT_lst['WT'+str(t+1)]['blades']['HS_vortex'] = hs_lst
            self.WT_lst['WT'+str(t+1)]['blades']['panels'] = panel_lst

    def Transform(vec,theta_r,D):
        x = vec[0]
        y = vec[1]*m.cos(theta_r)-vec[2]*m.sin(theta_r)
        z = vec[1]*m.sin(theta_r)+vec[2]*m.cos(theta_r)
        return [x,y,z+D]


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

def InducedVelocities(geometry):

    rotors = geometry.WT_lst
    N_pts = 0
    N_hs = 0
    ctrl_pts = []
    hs_vortex = []
    idx_lst = [0]
    for t in range(len(rotors)):
        N_pts += rotors['WT'+str(t+1)]['N_pts']
        N_hs += rotors['WT'+str(t+1)]['N_hs']
        idx_lst.append(N_hs)
        ctrl_pts = ctrl_pts + rotors['WT'+str(t+1)]['blades']['control_pts']
        hs_vortex =  hs_vortex + rotors['WT'+str(t+1)]['blades']['HS_vortex']

    u_ind_mat = np.zeros((N_pts,N_hs))
    v_ind_mat = np.zeros((N_pts,N_hs))
    w_ind_mat = np.zeros((N_pts,N_hs))

    # For every control point
    for i in range(N_pts):
        coords = ctrl_pts[i]['coords']

        # For every HS vortex
        for j in range(N_hs):
            hs = hs_vortex[j]
            u_ind = 0
            v_ind = 0
            w_ind = 0
            # For ever filament of HS vortex
            for k in range(len(hs)):
                fil = hs[k]

                # Compute induced velocity
                w = SingleFilament(coords,fil['pos1'],fil['pos2'],fil['Gamma'])
                u_ind+=w[0]
                v_ind+=w[1]
                w_ind+=w[2]

            # Store induced velocities per control point and HS vortex
            u_ind_mat[i,j] = u_ind
            v_ind_mat[i,j] = v_ind
            w_ind_mat[i,j] = w_ind

    return u_ind_mat, v_ind_mat, w_ind_mat, idx_lst





def LiftingLine(rotors,geometry):

    #Induced velocity matrices
    [u_ind_mat, v_ind_mat, w_ind_mat,idx_lst] = InducedVelocities(geometry)

    gamma = np.zeros((len(u_ind_mat[:,0])))
    D_lst = geometry.D_lst
    results = {}

    #Iteration controls
    it = 0
    it_max =1000
    error = 1e12
    limit = 1e-6
    UR = 0.2
    while it<it_max and error>limit:
        u_all = (u_ind_mat*gamma).sum(axis=1)
        v_all = (v_ind_mat*gamma).sum(axis=1)
        w_all = (w_ind_mat*gamma).sum(axis=1)
        gamma_new = []

        for t in range(len(rotors)):
            results['WT'+str(t+1)] = {'r':[], 'a':[], 'at':[], 'F_ax':[], 'F_az':[], 'Gamma':[]}
            rotor = rotors[t]
            rotor_geo = geometry.WT_lst['WT'+str(t+1)]
            omega = rotor.omega
            V_inf = rotor.V_inf
            N_pts = rotor_geo['N_pts']

            u = u_all[idx_lst[t]:idx_lst[t+1]]
            v = v_all[idx_lst[t]:idx_lst[t+1]]
            w = w_all[idx_lst[t]:idx_lst[t+1]]

            r_lst = rotor_geo['r_lst']

            a = np.zeros((N_pts))
            at = np.zeros((N_pts))
            F_ax = np.zeros((N_pts))
            F_az = np.zeros((N_pts))
            gamma_new_wt = np.zeros((N_pts))

            for i in range(N_pts):
                ctrl_pt = rotor_geo['blades']['control_pts'][i]['coords']
                coords = [ctrl_pt[0], ctrl_pt[1], ctrl_pt[2]- D_lst[t]]
                r = r_lst[i]
                omega_vec = np.cross([-omega,0,0],coords) #Rotational velocity at point
                # Velocity in x,y,z
                V = [V_inf[0] + u[i] + omega_vec[0], V_inf[1] + v[i] + omega_vec[1], V_inf[2] + w[i] + omega_vec[2]]

                azim_vec = np.cross([-1/r,0,0],coords)
                V_az = np.dot(azim_vec,V)
                axial_vec = [1,0,0]
                V_ax = np.dot(axial_vec,V)
                a[i] = (-u[i] + omega_vec[0])/V_inf[0]
                at[i] = V_az/(omega*r)-1

                V_mag = m.sqrt(V_ax**2+V_az**2)
                phi = m.atan(V_ax/V_az)
                [c, theta] = rotor.BladeGeometry(r)
                alpha = np.rad2deg(theta + phi)
                [Cl,Cd] = rotor.Polars(alpha)

                L = 0.5*V_mag**2*Cl*c
                D = 0.5*V_mag**2*Cd*c

                F_ax[i] = L*m.cos(phi)+D*m.sin(phi)
                F_az[i] = L*m.sin(phi)-D*m.cos(phi)
                gamma_new_wt[i] = 0.5*V_mag*Cl*c
                gamma_new.append(gamma_new_wt[i])

            results['WT'+str(t+1)]['r'] = np.array(r_lst)
            results['WT'+str(t+1)]['a'] = a
            results['WT'+str(t+1)]['at'] = at
            results['WT'+str(t+1)]['F_ax'] = F_ax
            results['WT'+str(t+1)]['F_az'] = F_az
            results['WT'+str(t+1)]['F_ax'] = F_ax
            results['WT'+str(t+1)]['Gamma'] = gamma_new_wt

        gamma_new = np.array(gamma_new)
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

    return results


#%% Single rotor

rotor = Rotor(50, 8, 3, 2, [1,0,0])

geometry = VortexGeometry([rotor], [40], 'ns', [4], [0])

results = LiftingLine([rotor],geometry)

dt = pd.read_csv('results.txt')
results_CSF=dt.to_numpy()


plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:39]/rotor.R,results['WT1']['F_ax'][0:39]/(0.5*rotor.V_inf[0]**2*rotor.R), label=r'Fnorm (LL)')
plt.plot(results['WT1']['r'][0:39]/rotor.R,results['WT1']['F_az'][0:39]/(0.5*rotor.V_inf[0]**2*rotor.R), label=r'Ftan (LL)')

plt.plot(results_CSF[:,2], results_CSF[:,3]/(0.5*rotor.V_inf[0]**2*rotor.R), 'b', label=r'Fnorm (BEM)')
plt.plot(results_CSF[:,2], results_CSF[:,4]/(0.5*rotor.V_inf[0]**2*rotor.R), 'g', label=r'Ftan (BEM)')
plt.legend()
plt.xlabel('r/R')


plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:39]/rotor.R,results['WT1']['a'][0:39], label=r'a (LL)')
plt.plot(results['WT1']['r'][0:39]/rotor.R,results['WT1']['at'][0:39], label=r'at (LL)')

plt.plot(results_CSF[:,2], results_CSF[:,0], 'b', label=r'a (BEM)')
plt.plot(results_CSF[:,2], results_CSF[:,1], 'g', label=r'at (BEM)')
plt.legend()
plt.xlabel('r/R')


plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:39]/rotor.R,results['WT1']['Gamma'][0:39])



#%% Double rotor

# Radius, TipSpeedRatio, NumberofBlades, GlobalPitch, WindSpeedVector
rotor1 = Rotor(50, 8, 3, 2, [1,0,0])
rotor2 = Rotor(50, 8, 3, 2, [1,0,0])

# # Rotor parameters, Number of Radius Positions, Spacing method, Number of Rotations, z-origin of rotors
geometry = VortexGeometry([rotor1,rotor2], [20,20], 'ns', [2,2], [0,100])
results = LiftingLine([rotor1,rotor2],geometry)


dt = pd.read_csv('results.txt')
results_CSF=dt.to_numpy()

plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:19]/rotor1.R,results['WT1']['F_ax'][0:19]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Fnorm (LL) WT1 1')
plt.plot(results['WT2']['r'][0:19]/rotor1.R,results['WT2']['F_ax'][0:19]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Fnorm (LL) WT2 1')
plt.plot(results['WT2']['r'][19:38]/rotor1.R,results['WT2']['F_ax'][19:38]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Fnorm (LL) WT2 2')
plt.plot(results['WT2']['r'][38:57]/rotor1.R,results['WT2']['F_ax'][38:57]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Fnorm (LL) WT2 3')
plt.plot(results_CSF[:,2], results_CSF[:,3]/(0.5*rotor1.V_inf[0]**2*rotor1.R), 'b', label=r'Fnorm (BEM)')
plt.legend()
plt.xlabel('r/R')

plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:19]/rotor1.R,results['WT1']['F_az'][0:19]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Ftan (LL) WT1 1')
plt.plot(results['WT2']['r'][0:19]/rotor1.R,results['WT2']['F_az'][0:19]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Ftan (LL) WT2 1')
plt.plot(results['WT2']['r'][19:38]/rotor1.R,results['WT2']['F_az'][19:38]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Ftan (LL) WT2 2')
plt.plot(results['WT2']['r'][38:57]/rotor1.R,results['WT2']['F_az'][38:57]/(0.5*rotor1.V_inf[0]**2*rotor1.R), label=r'Ftan (LL) WT2 3')
plt.plot(results_CSF[:,2], results_CSF[:,4]/(0.5*rotor1.V_inf[0]**2*rotor1.R), 'g', label=r'Ftan (BEM)')
plt.legend()
plt.xlabel('r/R')


plt.figure()
plt.grid()
plt.plot(results['WT1']['r'][0:19]/rotor1.R,results['WT1']['Gamma'][0:19], label='WT1 1')
plt.plot(results['WT2']['r'][0:19]/rotor1.R,results['WT2']['Gamma'][0:19], label='WT2 1')
plt.plot(results['WT1']['r'][19:38]/rotor1.R,results['WT1']['Gamma'][19:38], label='WT1 2')
plt.plot(results['WT2']['r'][19:38]/rotor1.R,results['WT2']['Gamma'][19:38], label='WT2 2')
plt.plot(results['WT1']['r'][38:57]/rotor1.R,results['WT1']['Gamma'][38:57], label='WT1 3')
plt.plot(results['WT2']['r'][38:57]/rotor1.R,results['WT2']['Gamma'][38:57], label='WT2 3')
plt.legend()
plt.xlabel('r/R')


