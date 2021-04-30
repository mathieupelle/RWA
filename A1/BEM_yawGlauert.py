# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:31:00 2021

@author: Mathieu Pellé
"""

import numpy as np
import pandas as pd
import math as m
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy.solvers import solve
from sympy import Symbol

def Glauert(x, sigma, Cn, phi, F):
    """
    param x: Unknown
    param sigma: Solidity [-]
    param Cn: Normal force coefficient [-]
    param phi: Flow angle [-]
    param F: Prandtl tip loss correction factor [-]
    return CT: Thrust force coefficient [-]
    return aa: Axial induction factor [-] """

    return [x[0] - (1 - x[1])**2 * sigma * Cn / m.sin(phi)**2,
            x[0] - 4 * x[1] * (1 - 0.25 * (5 - 3 * x[1]) * x[1]) * F]



class blade:
  def __init__(self, radius, number_of_blades, blade_start_ratio, global_pitch):
    self.R = radius
    self.B = number_of_blades
    self.r0 = blade_start_ratio*self.R
    self.theta0 = np.deg2rad(global_pitch)

  def chord(self,r):
    """
    Computes value of chord for radial position
    param r: radial positino [m]
    return c: chord length [m]
    """
    c=3*(1-r/self.R)+1
    return c

  def twist(self,r):
    """
    Computes value of twist for radial position
    param r: radial positino [m]
    return t: twist length [rad]
    """
    t=np.deg2rad(14*(1-r/self.R))
    return t

  def polar(self,alpha,df):
    """
    Computes interpolated lift and drag coefficients from angle of attack
    param alpha: Angle of attack [deg]
    return cl: Lift coefficient [-]
    return cd: Drag coefficient [-]
    """
    idx = df.iloc[(df['alpha'] - alpha).abs().argsort()[:2]]
    gradcl = (idx.cl.iloc[1] - idx.cl.iloc[0]) / (idx.alpha.iloc[1] - idx.alpha.iloc[0])
    gradcd = (idx.cd.iloc[1] - idx.cd.iloc[0]) / (idx.alpha.iloc[1] - idx.alpha.iloc[0])
    cl = gradcl * (alpha - idx.alpha.iloc[0]) + idx.cl.iloc[0]
    cd = gradcd * (alpha - idx.alpha.iloc[0]) + idx.cd.iloc[0]
    return cl, cd


class conditions:
    def __init__(self, wind_velocity, tip_speed_ratio, density, yaw_angle):
        self.U0 = wind_velocity
        self.TSR = tip_speed_ratio
        self.omega = self.TSR*self.U0/blade.R
        self.rho = density
        self.gamma = np.deg2rad(yaw_angle)

class functions:

    def loads(r,psi,aa,at):
        """
        Computes value of twist for radial position
        param r: Radial position [m]
        param aa: Axial induction factor [-]
        param at: Tangential induction factor [-]
        return F_az: Azimuthal Force [-]
        return F_ax: Axial Force [-]
        return aa: Corrected axial induction factor [-]
        """
        xsi=(0.6*aa+1)*cond.gamma
        K=m.tan(xsi/2)

        #====Basic Glauert Yaw====#
        aa=aa*(1+K*r/blade.R*m.sin(cond.gamma))
        V_a = cond.U0*(m.cos(cond.gamma)-aa)
        V_t = cond.omega*r*(1+at)+cond.U0*m.sin(cond.gamma)*m.sin(psi)

        #====COleman-Glauert Yaw====#
        # F=0.5*(r/blade.R+0.4*(r/blade.R)**3+0.4*(r/blade.R)**5)
        # V_a = cond.U0*(m.cos(cond.gamma)-aa*(1+F*K*m.sin(psi)))+cond.omega*r*at*m.cos(psi)*m.sin(xsi)*(1+m.sin(psi)*m.sin(xsi))
        # V_t = cond.omega*r*(1+at*m.cos(xsi)*(1+m.sin(psi)*m.sin(xsi)))+cond.U0*m.cos(cond.gamma)*(aa*m.tan(xsi/2)*(1+F*K*m.sin(psi))-m.sin(cond.gamma))

        V = m.sqrt(V_a**2+V_t**2)
        phi = m.atan(V_a/V_t)
        alpha = phi - (blade.theta0 + blade.twist(r))
        cl, cd = blade.polar(np.rad2deg(alpha),df)
        L = 0.5*cond.rho*V**2*blade.chord(r)*cl
        D = 0.5*cond.rho*V**2*blade.chord(r)*cd
        F_az = L*m.sin(phi) - D*m.cos(phi) #per unit length
        F_ax = L*m.cos(phi) + D*m.sin(phi) #per unit length
        CT = F_ax*blade.B/(0.5*cond.rho*cond.U0**2*2*m.pi*r)


        #====Basic Glauert Yaw====#
        CT1=1.816
        CT2=2*m.sqrt(CT1)-CT1
        if CT>=CT2:
            aa = 1 + (CT-CT1)/(4*m.sqrt(CT1)-4)
        else:
            aa = 0.5 - m.sqrt(1-CT)/2

        #====Wilson & Walker====#
        # F = 2 / m.pi * m.acos(m.exp(-blade.B / 2 * (blade.R - r) / r / m.sin(abs(phi))))
        # sigma=blade.chord(r)*blade.B/(2*m.pi*r)
        # Cn=F_az/(0.5*cond.rho*V**2*blade.chord(r))
        # if aa > 0.2:
        #     ac = 0.2
        #     x = Symbol('x')
        #     CT = (1 - x)**2 * sigma * Cn / m.sin(phi)**2
        #     WW = CT - 4 * (ac**2 + (1 - 2 * ac) * x) * F
        #     sol = solve(WW, x)
        #     aa = sol[0]
        #====Glauert correction====#
        # if aa > 1 / 3:
        #     CT, aa = fsolve(Glauert, [1, aa], args=(sigma, Cn, phi, F))
        # else:
        #     aa = 1 / (4 * F * m.sin(phi)**2 / sigma / Cn + 1)

        return F_az, F_ax, CT, aa

    def Prandtl(r, aa):
        """
        Computes Prandtl correction factor
        param r: Radial position [m]
        param aa: Axial induction factor [-]
        return f_total: Prandtl correction factor [-]
        """
        mu = r/blade.R
        mu_root = blade.r0/blade.R
        f_tip = 2/m.pi*m.acos(m.exp(-blade.B/2*((1-mu)/mu)*m.sqrt(1+(cond.TSR*mu)**2/(1-aa)**2)))
        f_root = 2/m.pi*m.acos(m.exp(-blade.B/2*((mu-mu_root)/mu)*m.sqrt(1+(cond.TSR*mu)**2/(1-aa)**2)))
        f_total = f_tip*f_root
        if f_total<1e-6:
            f_total=1e-6
        return f_total


    def iteration(r,psi):
        """
        Iterative procedure to compute final induction factors
        param r: Radial position [m]
        return aa: Axial induction factor [-]
        return at: Tangential induction factor [-]
        """
        aa = 0.3
        at = 0.02
        error = 1e-6
        iterations = 0
        iterations_max = 1000
        residual = 1e4
        UR_factor = 0.25
        while iterations<iterations_max and residual>error:
            F_az, F_ax, CT, aa_new = functions.loads(r,psi,aa,at)
            f_total = functions.Prandtl(r,aa_new)
            aa_new = aa_new/f_total
            aa = (1-UR_factor)*aa + UR_factor*aa_new


            at = F_az*blade.B/(2*cond.rho*2*m.pi*r*cond.U0**2*(1-aa)*cond.TSR*r/blade.R)
            at = at/f_total

            residual = abs(aa - aa_new)
            iterations+=1
        return aa, at

    def outputs(r_lst,psi_lst):
        aa_lst = np.zeros((len(r_lst),len(psi_lst)))
        at_lst = np.zeros((len(r_lst),len(psi_lst)))
        Faz_lst = np.ones((len(r_lst),len(psi_lst)))
        Fax_lst = np.ones((len(r_lst),len(psi_lst)))
        Q_lst = np.ones((len(r_lst),len(psi_lst)))

        for i in range(len(r_lst)):
            for j in range(len(psi_lst)):
                aa_lst[i,j], at_lst[i,j] = functions.iteration(r_lst[i],psi_lst[j])
                Faz_lst[i,j], Fax_lst[i,j], CT, aa = functions.loads(r_lst[i],psi_lst[j],aa_lst[i,j], at_lst[i,j])
                Q_lst[i,j] = Faz_lst[i,j]*r_lst[i]*dr*blade.B
            if i in [len(r_lst)/4,len(r_lst)/2,len(r_lst)*3/4,len(r_lst)-1]:
                print('Completed ' +str(25*round(i/len(r_lst)*100/25))+' %')

        Faz_lst=Faz_lst*dpsi/2/m.pi
        Fax_lst=Fax_lst*dpsi/2/m.pi
        return Fax_lst, Faz_lst, aa_lst, at_lst, Q_lst

    def coefficients(Fax_lst,Faz_lst):
        cT=dr*np.sum(Fax_lst)*blade.B/(0.5*cond.rho*cond.U0**2*m.pi*blade.R**2)
        cP=dr*np.dot(np.sum(Faz_lst,axis=1),r_lst)*blade.B*cond.omega/(0.5*cond.rho*cond.U0**3*m.pi*blade.R**2)
        return cT, cP

class plotting:

    def induction(r_lst,psi_lst,aa_lst,at_lst,style):
        plt.figure()
        if style=='lines':
            for i in range(len(aa_lst[0])):
                plt.plot(r_lst/blade.R, aa_lst[:,i])
               # plt.plot(r_lst/blade.R, at_lst[:,i])
            plt.grid()
        elif style=='polar':
            theta, r = np.meshgrid(psi_lst,r_lst)
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            pp=ax.contourf(theta, r, aa_lst)
            cbar = plt.colorbar(pp, orientation='vertical')
            cbar.ax.set_ylabel('Normalised Loads')
        else:
            plt.plot(r_lst/blade.R, np.mean(aa_lst,axis=1), label=r'$a$')
            plt.plot(r_lst/blade.R, np.mean(at_lst,axis=1), label=r'$a^,$')
            plt.legend()
            plt.grid()
            plt.xlabel('r/R')


    def forces(r_lst,psi_lst,Fax_lst,Faz_lst,style):
        plt.figure()
        Fax_lst=Fax_lst/(0.5*cond.U0**2*blade.R*cond.rho)
        Faz_lst=Faz_lst/(0.5*cond.U0**2*blade.R*cond.rho)
        if style=='lines':
            for i in range(len(Fax_lst[0])):
                plt.plot(r_lst/blade.R, Fax_lst[:,i], label=r'Fnorm')
                plt.plot(r_lst/blade.R, Faz_lst[:,i], label=r'Ftan')
        elif style=='polar':
            theta, r = np.meshgrid(psi_lst,r_lst)
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            pp = ax.contourf(theta, r, Fax_lst)
            cbar = plt.colorbar(pp, orientation='vertical')
            cbar.ax.set_ylabel('Normalised Loads')
        else:
            plt.plot(r_lst/blade.R, np.sum(Fax_lst,axis=1), label=r'Fnorm')
            plt.plot(r_lst/blade.R, np.sum(Faz_lst,axis=1), label=r'Ftan')
            plt.legend()
            plt.grid()
            plt.xlabel('r/R')

    def torque(r_lst,psi_lst,Q_lst,style):
        plt.figure()
        Q_lst=Q_lst/(0.5*cond.rho*cond.U0**2*blade.R**2)
        if style=='lines':
            for i in range(len(Fax_lst[0])):
                plt.plot(r_lst/blade.R, Q_lst[:,i])
        elif style=='polar':
            theta, r = np.meshgrid(psi_lst,r_lst)
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            pp = ax.contourf(theta, r, Q_lst)
            cbar = plt.colorbar(pp, orientation='vertical')
            cbar.ax.set_ylabel('Torque coefficient')
        else:
            plt.plot(r_lst/blade.R, np.sum(Q_lst,axis=1), label=r'Cq')
            plt.legend()
            plt.grid()
            plt.xlabel('r/R')



#%% Inputs and Calculations

#Import airfoil file
df = pd.read_excel('polar DU95W180 (3).xlsx',  header=3)
df.columns=['alpha','cl','cd','cm']

#Blade geometry inputs
blade = blade(50,3,0.2,-2)

#Operating condions
cond = conditions(1,8,1.225,30)

#Radial and azimuthal finite differencing
r_lst = np.linspace(blade.r0+blade.R*0.005,blade.R*0.995,30)
dr = r_lst[1] - r_lst[0]
#psi_lst = np.arange(0,2*m.pi,2*m.pi/30)
psi_lst = np.linspace(0,2*m.pi,20)
dpsi = psi_lst[1] - psi_lst[0]

#BEM code calculations
Fax_lst, Faz_lst, aa_lst, at_lst, Q_lst = functions.outputs(r_lst,psi_lst)
cT, cP = functions.coefficients(Fax_lst, Faz_lst)
print('Thrust coefficient: ' +str(round(cT,3)))
print('Power coefficient: ' +str(round(cP,3)))

#Plotting
plotting.induction(r_lst,psi_lst,aa_lst,at_lst,'lines')
plotting.induction(r_lst,psi_lst,aa_lst,at_lst,'polar')
plotting.forces(r_lst,psi_lst,Fax_lst,Faz_lst,'polar')
plotting.torque(r_lst,psi_lst,Q_lst,'polar')



