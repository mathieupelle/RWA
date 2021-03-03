# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:31:00 2021

@author: Mathieu Pellé
"""

import numpy as np
import pandas as pd
import math as m
import matplotlib.pyplot as plt


class blade:
  def __init__(self, radius, number_of_blades, blade_start_ratio, global_pitch, airfoil_file_name):
    self.R = radius
    self.B = number_of_blades
    self.r0 = blade_start_ratio*self.R
    self.theta0 = np.deg2rad(global_pitch)
    self.file = airfoil_file_name

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

  def polar(self,alpha):
    """
    Computes interpolated lift and drag coefficients from angle of attack
    param alpha: Angle of attack [deg]
    return cl: Lift coefficient [-]
    return cd: Drag coefficient [-]
    """
    df = pd.read_excel(self.file,  header=3)
    df.columns=['alpha','cl','cd','cm']
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

    def loads(r,aa,at):
        """
        Computes value of twist for radial position
        param r: Radial position [m]
        param aa: Axial induction factor [-]
        param at: Tangential induction factor [-]
        return F_az: Azimuthal Force [-]
        return F_ax: Axial Force [-]
        return aa: Corrected axial induction factor [-]
        """
        V_a = cond.U0*(m.cos(cond.gamma)-aa)
        V_t = cond.omega*r*(1+at)+cond.U0*m.sin(cond.gamma)
        V = m.sqrt(V_a**2+V_t**2)
        phi = m.atan(V_a/V_t)
        alpha = phi - (blade.theta0 + blade.twist(r))
        cl, cd = blade.polar(np.rad2deg(alpha))
        L = 0.5*cond.rho*V**2*blade.chord(r)*cl
        D = 0.5*cond.rho*V**2*blade.chord(r)*cd
        F_az = L*m.sin(phi) - D*m.cos(phi)
        F_ax = L*m.cos(phi) + D*m.sin(phi)
        CT = F_ax*blade.B/(0.5*cond.rho*cond.U0**2*2*m.pi*r)
        #CT updtade
        CT1=1.816
        CT2=2*m.sqrt(CT1)-CT1
        if CT>=CT2:
            aa = 1 + (CT-CT1)/(4*m.sqrt(CT1)-4)
        else:
            aa = 0.5 - m.sqrt(1-CT)/2
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


    def iteration(r):
        """
        Iterative procedure to compute final induction factors
        param r: Radial position [m]
        return aa: Axial induction factor [-]
        return at: Tangential induction factor [-]
        """
        aa = 0.2
        at = 0.2
        error = 1e-6
        iterations = 0
        iterations_max = 1000
        residual = 1e4
        UR_factor = 0.25
        while iterations<iterations_max and residual>error:
            F_az, F_ax, CT, aa_new = functions.loads(r,aa,at)
            f_total = functions.Prandtl(r,aa_new)
            aa_new = aa_new/f_total
            aa = (1-UR_factor)*aa + UR_factor*aa_new

            at = F_az*blade.B/(2*cond.rho*2*m.pi*r*cond.U0**2*(1-aa)*cond.TSR*r/blade.R)
            at = at/f_total

            residual = abs(aa - aa_new)
            iterations+=1
        return aa, at

#%% Inputs and Calculations

blade = blade(50,3,0.2,-2, 'polar DU95W180 (3).xlsx')
cond = conditions(1,8,1.225,0)

r_lst = np.linspace(blade.r0+blade.R*0.005,blade.R*0.995,80)
dr = r_lst[1] - r_lst[0]


induction = np.zeros((len(r_lst),2))
output=np.zeros((len(r_lst),4))
for i in range(len(r_lst)):
    induction[i,:] = functions.iteration(r_lst[i])
    output[i,:] = functions.loads(r_lst[i],induction[i,0],induction[i,1])
    if i in [len(r_lst)/4,len(r_lst)/2,len(r_lst)*3/4,len(r_lst)-1]:
        print('Completed ' +str(25*round(i/len(r_lst)*100/25))+' %')
cT=dr*np.sum(output[:,1])*blade.B/(0.5*cond.rho*cond.U0**2*m.pi*blade.R**2)
cP=dr*np.dot(output[:,0],r_lst)*blade.B*cond.omega/(0.5*cond.rho*cond.U0**3*m.pi*blade.R**2)
print('Thrust coefficient: ' +str(round(cT,3)))
print('Power coefficient: ' +str(round(cP,3)))

#%% Plotting

dt = pd.read_csv('Validation/results.txt')
results=dt.to_numpy()

fig1 = plt.figure(figsize=(12, 6))
plt.title('Axial and tangential induction')
plt.plot(r_lst/blade.R, induction[:,0], 'r-', label=r'$a$')
plt.plot(r_lst/blade.R, induction[:,1], 'g-', label=r'$a^,$')
plt.plot(results[:,2], results[:,0], 'r--', label=r'$a_{val}$')
plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,_{val}$')
plt.grid()
plt.xlabel('r/R')
plt.legend()

fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Normal and tagential force, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$')
plt.plot(r_lst/blade.R, output[:,1]/(0.5*cond.U0**2*blade.R*cond.rho), 'r-', label=r'Fnorm')
plt.plot(r_lst/blade.R, output[:,0]/(0.5*cond.U0**2*blade.R*cond.rho), 'g-', label=r'Ftan')
plt.plot(results[:,2], results[:,3]/(0.5*cond.U0**2*blade.R), 'r--', label=r'Fnorm (val)')
plt.plot(results[:,2], results[:,4]/(0.5*cond.U0**2*blade.R), 'g--', label=r'Ftan (val)')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()