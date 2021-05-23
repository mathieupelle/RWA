# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:38:44 2021

@author: Mathieu Pellé
"""


import numpy as np
import math as m


class VortexGeometry:
    def __init__(self, rotors, NumberofRotations, RotorLocations, result_BEM):
        """
          Class that defines the entire vortex system with filaments and control points

          Parameters
          ----------
          rotors : List with different rotor geometries imported from BEM code
          NumberofRotations : Number of rotations for wake discretisation
          RotorLocations : List with locations of rotors in z. First rotor set at 0.
          result_BEM : BEM results used to correct the wake propagation speed locally

        """

        N_WT = len(rotors) #Number of wind turbines
        self.D_lst = RotorLocations
        self.WT_lst = {}
        self.a_BEM = result_BEM.a[:,0]
        self.mu_BEM = result_BEM.mu[:,0]

        # Looping for every wind turbine
        for t in range(N_WT):

            # Empty dictionnay for storing data
            self.WT_lst['WT'+str(t+1)] = {'r_lst':[], 'psi_lst':[],
                                          'N_pts':0, 'N_hs':0, 'N_rot':0, 'N_rad':0,
                                          'blades':{'control_pts':[], 'HS_vortex':[], 'panels':[]}}

            # Specific parameters for one rotor
            rotor = rotors[t] #Rotor class
            N_blades = rotor.n_blades #Number of blades
            TSR = rotor.TSR #Tip speed ratio
            R = rotor.radius #Radius
            N_rad = len(rotor.mu) #Number of radial positions
            D = self.D_lst[t] #Z shift of rotor center
            s_lst = rotor.mu*rotor.radius #Array of radial positions
            N_rot = NumberofRotations[t] # Number of rotations
            psi_lst = np.linspace(0,2*m.pi*N_rot,2*N_rot*30+1) #Arrray of discretised rotation

            #Emptry lists to store geometry parameters
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
                    r = (s_lst[i+1]+s_lst[i])/2 #Radial position in the middle of each element
                    r_lst.append(r)
                    [c, theta] = VortexGeometry.geometry_interp(r, rotor) #Chord and pitch
                    vec_n = [m.cos(theta),0,-m.sin(theta)] #Normal vector
                    vec_n = VortexGeometry.Transform(vec_n, theta_r, D) #Transformed normal vector
                    vec_t = [-m.sin(theta),0,-m.cos(theta)]
                    vec_t = VortexGeometry.Transform(vec_t, theta_r, D) #Transformed tangential vector
                    pts_lst.append({'coords': VortexGeometry.Transform([0,r,0],theta_r,D), 'chord': c, 'normal': vec_n, 'tangent': vec_t})

                    ## Horseshoe vortices ##
                    #VF1: filament in span direction at 0.25c
                    VF = {'pos1':[0,s_lst[i],0], 'pos2':[0,s_lst[i+1],0], 'Gamma':1}
                    fil.append(VF)

                    #VF2: filament in upstream direction at 1st radial element position
                    [c, theta] = VortexGeometry.geometry_interp(s_lst[i], rotor)
                    VF = {'pos1':[c*m.sin(-theta),s_lst[i],-c*m.cos(theta)], 'pos2':[0,s_lst[i],0], 'Gamma':1}
                    fil.append(VF)

                    #All filaments downstream of VF2 up to limit defined by number of rotations
                    for j in range(len(psi_lst)-1):
                        # 1st point (downstream)
                        x = fil[-1]['pos1'][0]
                        y = fil[-1]['pos1'][1]
                        z = fil[-1]['pos1'][2]

                        # 2nd point (upstream)
                        #a_temp = 1/3 # Average wake velocity
                        a_temp = np.interp(s_lst[i]/rotor.radius, self.mu_BEM, self.a_BEM) # Local wake velocity
                        dx = (psi_lst[j+1]-psi_lst[j])/TSR*R*(1-a_temp)
                        dy = (m.cos(-psi_lst[j+1])-m.cos(-psi_lst[j]))*s_lst[i]
                        dz = (m.sin(-psi_lst[j+1])-m.sin(-psi_lst[j]))*s_lst[i]
                        VF = {'pos1':[x+dx,y+dy,z+dz], 'pos2':[x,y,z], 'Gamma':1}
                        fil.append(VF)

                    #VF3: filament in downstream direction at 2nd radial element position
                    [c, theta] = VortexGeometry.geometry_interp(s_lst[i+1], rotor)
                    VF = {'pos1':[0,s_lst[i+1],0], 'pos2':[c*m.sin(-theta),s_lst[i+1],-c*m.cos(theta)], 'Gamma':1}
                    fil.append(VF)
                    #All filaments downstream of VF3 up to limit defined by number of rotations
                    for j in range(len(psi_lst)-1):
                        # 1st point (upstream)
                        x = fil[-1]['pos2'][0]
                        y = fil[-1]['pos2'][1]
                        z = fil[-1]['pos2'][2]

                        # 2nd point (downstream)
                        #a_temp = 1/3
                        a_temp = np.interp(s_lst[i+1]/rotor.radius, self.mu_BEM, self.a_BEM)
                        dx = (psi_lst[j+1]-psi_lst[j])/TSR*R*(1-a_temp)
                        dy = (m.cos(-psi_lst[j+1])-m.cos(-psi_lst[j]))*s_lst[i+1]
                        dz = (m.sin(-psi_lst[j+1])-m.sin(-psi_lst[j]))*s_lst[i+1]
                        VF = {'pos1':[x,y,z], 'pos2':[x+dx,y+dy,z+dz], 'Gamma':1}
                        fil.append(VF)

                    #Transforming all filaments coordinates
                    for j in range(len(fil)):
                        fil[j]['pos1'] = VortexGeometry.Transform(fil[j]['pos1'], theta_r, D)
                        fil[j]['pos2'] = VortexGeometry.Transform(fil[j]['pos2'], theta_r, D)

                    # ## Panels ##
                    # [c1, theta1] = VortexGeometry.geometry_interp(s_lst[i], rotor)
                    # [c2, theta2] = VortexGeometry.geometry_interp(s_lst[i+1], rotor)
                    # p1 = VortexGeometry.Transform([-0.25*c1*m.sin(-theta1)*0 , s_lst[i], 0.25*c1*m.cos(theta1)],theta_r,D)
                    # p2 = VortexGeometry.Transform([-0.25*c2*m.sin(-theta2)*0 , s_lst[i+1], 0.25*c2*m.cos(theta2)],theta_r,D)
                    # p3 = VortexGeometry.Transform([0.75*c2*m.sin(-theta2)*0 , s_lst[i+1], -0.75*c2*m.cos(theta2)],theta_r,D)
                    # p4 = VortexGeometry.Transform([0.75*c1*m.sin(-theta1)*0 , s_lst[i], -0.75*c1*m.cos(theta1)],theta_r,D)

                    # panel_lst.append({'p1':p1, 'p2':p2, 'p3':p3, 'p4':p4})
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
        """
          Coordinate transformation from x,y,z to correct blade position with z-shift of rotor

          Parameters
          ----------
          vec : Vector to transform
          theta_r : Rotation angle corresponding to blade position
          D : z-coordinate of rotor

          Returns
          -------
          [x,y,z+D] : Transformed vector

        """
        x = vec[0]
        y = vec[1]*m.cos(theta_r)-vec[2]*m.sin(theta_r)
        z = vec[1]*m.sin(theta_r)+vec[2]*m.cos(theta_r)
        return [x,y,z+D]

    def geometry_interp(r, rotor):
        c = np.interp(r/rotor.radius, rotor.mu, rotor.chord)
        theta = np.deg2rad(np.interp(r/rotor.radius, rotor.mu, rotor.beta) + rotor.theta)
        return c, theta

    def polar_interp(alpha, rotor):
        cl = np.interp(alpha,np.array(rotor.polars['alpha']),np.array(rotor.polars['Cl']))
        cd = np.interp(alpha,np.array(rotor.polars['alpha']),np.array(rotor.polars['Cd']))
        return cl, cd


def SingleFilament(ctrl_point,point1,point2,gamma):
    """
      Computes induced velocities from vortex filament using Biot-Savart

      Parameters
      ----------
      ctrl_point : Point at which to compute velocities (control point)
      point1 : Vortex filament first point coordinates
      point2 : Vortex filament second point coordinates
      gamma : Circulation of vortex

      Returns
      -------
      [U,V,W] : Induced velocities vector

    """

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
    """
      Computes induced velocity matrix (all control points and all horshoe vortices) for unit circulation

      Parameters
      ----------
      geometry : Class with vortex system geometry (control points and vortices positions)

      Returns
      -------
      u_ind_mat : x induced velocity matrix
      v_ind_mat : y induced velocity matrix
      w_ind_mat : z induced velocity matrix
      idx_lst : list of indices to know when new rotor starts

    """

    # Combining multiple rotor system
    rotors = geometry.WT_lst
    N_pts = 0
    N_hs = 0
    ctrl_pts = []
    hs_vortex = []
    idx_lst = [0]
    # Loop over rotors
    for t in range(len(rotors)):
        N_pts += rotors['WT'+str(t+1)]['N_pts'] #Total number of control points of system
        N_hs += rotors['WT'+str(t+1)]['N_hs'] #Total number of vortices of system
        idx_lst.append(N_hs) #Index of where circulation switches rotors
        ctrl_pts = ctrl_pts + rotors['WT'+str(t+1)]['blades']['control_pts'] #Stacking control points
        hs_vortex =  hs_vortex + rotors['WT'+str(t+1)]['blades']['HS_vortex'] #Stacking vortices

    # Empty induced velocity matrices
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
                # Sum effect of all filaments at control point for 1 vortex
                u_ind+=w[0]
                v_ind+=w[1]
                w_ind+=w[2]

            # Store induced velocities per control point and HS vortex
            u_ind_mat[i,j] = u_ind
            v_ind_mat[i,j] = v_ind
            w_ind_mat[i,j] = w_ind

    return u_ind_mat, v_ind_mat, w_ind_mat, idx_lst





def LiftingLine(rotors,geometry,result_BEM):
    """
      Iterative procedure for finding the circulation of system and rotor performance

      Parameters
      ----------
      rotors : List with rotors containing operating conditions and basic geometry
      geometry : Class with vortex system geometry (control points and vortices positions)
      result_BEM : BEM results used as starting guess for circulation to speed up convegence

      Returns
      -------
      results : Class with all relevant results

    """

    #Induced velocity matrices
    [u_ind_mat, v_ind_mat, w_ind_mat,idx_lst] = InducedVelocities(geometry)

    #Arranging BEM  circulation results to fit required format (N blades, M rotors)
    array = np.append(result_BEM.circulation,result_BEM.circulation)
    gamma = np.append(array,result_BEM.circulation)
    #gamma = np.zeros((len(u_ind_mat[:,0])))

    D_lst = geometry.D_lst #z distance list

    #Initialising list with multiple results classes
    results = []
    for i in range(len(rotors)):
        results.append(WT_Result(i))

    ## Lifting line iterative loop ##
    #Iteration controls
    it = 0
    it_max = 1000 #Max iteration number
    error = 1e12
    limit = 1e-6 #Error convergence criteria
    UR = 0.2 #Under relaxation factor
    while it<it_max and error>limit:
        if it == it_max - 1:
            print('Max. iterations reached')

        # Multiplying induced velocity matrix with circulation
        u_all = (u_ind_mat*gamma).sum(axis=1)
        v_all = (v_ind_mat*gamma).sum(axis=1)
        w_all = (w_ind_mat*gamma).sum(axis=1)

        gamma_new = []
        # For every rotor
        for t in range(len(rotors)):
            rotor = rotors[t] #Specific rotor parameters
            rotor_geo = geometry.WT_lst['WT'+str(t+1)] #Specific vortex system geometry
            omega = rotor.omega #Rotation speed
            V_inf = [rotor.wind_speed, 0, 0] #Wind velocity
            N_pts = rotor_geo['N_pts'] #Number of radial points

            #Unpacking induced velocities for current rotor
            u = u_all[idx_lst[t]:idx_lst[t+1]]
            v = v_all[idx_lst[t]:idx_lst[t+1]]
            w = w_all[idx_lst[t]:idx_lst[t+1]]

            r_lst = rotor_geo['r_lst'] #Radial positions (middle of element)

            # Empty lists for storage
            a = np.zeros((N_pts))
            at = np.zeros((N_pts))
            F_ax = np.zeros((N_pts))
            F_az = np.zeros((N_pts))
            phi = np.zeros((N_pts))
            alpha = np.zeros((N_pts))
            gamma_new_wt = np.zeros((N_pts))

            # For every radial position
            for i in range(N_pts):
                ctrl_pt = rotor_geo['blades']['control_pts'][i]['coords'] #Control point coordinates
                coords = [ctrl_pt[0], ctrl_pt[1], ctrl_pt[2]- D_lst[t]] #Correcting coordinate to local rotor axis
                r = r_lst[i] #Radius at control point
                omega_vec = np.cross([-omega,0,0],coords) #Rotational velocity at point

                # Velocity vector
                V = [V_inf[0] + u[i] + omega_vec[0], V_inf[1] + v[i] + omega_vec[1], V_inf[2] + w[i] + omega_vec[2]]
                azim_vec = np.cross([-1/r,0,0],coords) #Azimuthal vector
                V_az = np.dot(azim_vec,V) #Azimuthal velocity
                axial_vec = [1,0,0] #Axial vector
                V_ax = np.dot(axial_vec,V) #Axial velocity

                # BEM equations
                #a[i] = -(u[i] + omega_vec[0])/V_inf[0]
                a[i] = 1-V_ax/V_inf[0] #Axial induction
                at[i] = V_az/(omega*r)-1 #Azimuthal induction
                V_mag = m.sqrt(V_ax**2+V_az**2) #Velocity magnitude
                phi[i] = m.atan(V_ax/V_az) #Flow angle
                [c, theta] = VortexGeometry.geometry_interp(r, rotor)
                alpha[i] = np.rad2deg(phi[i] - theta) #Angle of attack
                [Cl,Cd] =  VortexGeometry.polar_interp(alpha[i], rotor) #Lift and drag coefficients

                L = 0.5*V_mag**2*Cl*c #Lift
                D = 0.5*V_mag**2*Cd*c #Drag

                F_ax[i] = L*m.cos(phi[i])+D*m.sin(phi[i]) #Axial force
                F_az[i] = L*m.sin(phi[i])-D*m.cos(phi[i]) #Azimuthal force
                gamma_new_wt[i] = 0.5*V_mag*Cl*c #Updated circulation
                gamma_new.append(gamma_new_wt[i]) #Stacking circulation to combine for all rotors

            # Storing all results
            results[t].mu = np.array(r_lst)/rotor.radius
            results[t].a = a
            results[t].ap = at
            results[t].f_nor = F_ax
            results[t].f_tan = F_az
            results[t].circulation = gamma_new_wt
            results[t].alpha = alpha
            results[t].phi = np.rad2deg(phi)

        # Cimputing error between new and old circulation values
        gamma_new = np.array(gamma_new)
        error=np.linalg.norm(gamma_new - gamma)

        ## Carlos' (weird) way to control the iteration
        # referror=max(abs(gamma_new))
        # referror=max(referror,0.001)
        # error=max(abs(gamma_new-gamma))
        # error = error/referror
        # UR = max((1-error)*0.3,0.1)

        gamma = UR*gamma_new + (1-UR)*gamma #Applying under relaxation
        it+=1

    return results


class WT_Result(object):
    """
      Empty class for storing results

    """
    def __init__(self,name):
        self.a = 0
        self.ap = 0
        self.phi = 0
        self.alpha = 0
        self.f_nor = 0
        self.f_tan = 0
        self.mu = 0
        self.circulation = 0







