## BEMT code using classes

import numpy as np
#import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import math as m


class Rotor:
  def __init__(self):
      
    #Geometric data of the blade
    self.radius = 50 #[m]
    self.n_blades = 3
    self.theta = -2 #Pitch angle [deg]
    self.N_radial = 60 #Number of sections
    self.mu = np.linspace(0.2,1,self.N_radial)
    self.beta = 14*(1-self.mu) #Twist angle in degrees
    self.chord = 3*(1-self.mu)+1 #Chord length in meters
    
    self.N_azimuth = 15 #Number of angular sections
    self.azimuth = np.linspace(0,2*np.pi,self.N_azimuth)
    
    #Polar data
    self.polars = pd.read_excel('polar DU95W180 (3).xlsx',header = 3,names=['alpha','Cl','Cd','Cm'])

  def SetOperationalData(self,wind_speed,TSR,yaw,rho=1.225):
    self.wind_speed = wind_speed
    self.TSR = TSR
    self.yaw = yaw*np.pi/180 #Input yaw should be in degrees!
    self.omega = wind_speed*TSR/self.radius
    self.rho = rho
    
    
class Results: #Create the variables to store the results from BEMT
    def __init__(self,N_radial,N_azimuth):
        self.a,self.ap,self.phi,self.alpha,self.cl,self.cd,self.f_nor,self.f_tan,self.f,self.ite,self.mu,self.chord,self.beta =  np.zeros((13,N_radial-1,N_azimuth-1))
       # self.CT, self.CP, self.CQ, self.P, self.T, self.Q = np.zeros((6,1))
        
    def Integrate(self,Rotor):
        #Calculate global CT and CP
        d_r = (Rotor.mu[2]-Rotor.mu[1])*Rotor.radius
        d_azi = 2*np.pi/np.size(self.a,1)
        self.CT = np.sum(self.f_nor*Rotor.n_blades*d_azi/2/np.pi*d_r)/(0.5*Rotor.rho*Rotor.wind_speed**2*np.pi*Rotor.radius**2)
        
        dTorque = self.f_tan*d_r*self.mu*Rotor.radius*d_azi/2/np.pi
        self.CP = np.sum(dTorque*Rotor.n_blades*Rotor.omega)/(0.5*Rotor.rho*Rotor.wind_speed**3*np.pi*Rotor.radius**2)

            
        
 #       CT = np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2))
#CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))
    
class BEMT:
    def __init__(self,Rotor):
        self.Rotor = Rotor
    
    def RelativeVelocities(self,a,ap,mu,azimuth=0):
        u_tan = self.Rotor.omega*self.Rotor.radius*mu*(1+ap) + self.Rotor.wind_speed*np.sin(self.Rotor.yaw)*np.sin(azimuth)
        u_nor = self.Rotor.wind_speed*(np.cos(self.Rotor.yaw)-a)
        
        # Testing Glauert correction
        psi = (0.6*a+1)*self.Rotor.yaw
        K = 2*np.tan(psi/2)
        
      # u_nor = self.Rotor.wind_speed*(np.cos(self.Rotor.yaw)-a*K*mu*np.sin(self.Rotor.yaw))

        u_rel = np.sqrt(u_tan**2 + u_nor**2)
        phi = np.arctan2(u_nor,u_tan)
        
        return u_tan,u_nor,u_rel,phi 
    
    
    def AirfoilCoefficients(self,alpha):
        cl = np.interp(alpha*180/np.pi,np.array(self.Rotor.polars['alpha']),np.array(self.Rotor.polars['Cl']))
        cd = np.interp(alpha*180/np.pi,np.array(self.Rotor.polars['alpha']),np.array(self.Rotor.polars['Cd']))    
   
        return cl,cd
    
    def Forces(self,chord,phi,u_rel,cl,cd):
        lift = 0.5*self.Rotor.rho*u_rel**2*chord*cl
        drag = 0.5*self.Rotor.rho*u_rel**2*chord*cd
        
        f_tan = lift*np.sin(phi) - drag*np.cos(phi)
        f_nor = lift*np.cos(phi) + drag*np.sin(phi)
        
        return lift,drag,f_tan,f_nor
    
    def NewInductionFactor(self,CT,yaw,a):
        if yaw == 0:
            CT_1 = 1.816 #Constants for the empirical calculation
            CT_2 = 2*np.sqrt(CT_1) - CT_1
            if CT < CT_2:
                a_new = 0.5 - np.sqrt(1-CT)/2
            else: 
                a_new = 1 + (CT-CT_1)/(4*np.sqrt(CT_1)-4)               
        else:
                a_new = CT/(4*np.sqrt(1-a*(2*np.cos(yaw)-a)))
                
        return a_new
    
    
    def PrandtlTipCorrection(self,mu,a_new):
        mu_root = self.Rotor.mu[0]
        #Tip correction
        exp = np.exp(-self.Rotor.n_blades/2 * ((1-mu)/mu) * np.sqrt(1+self.Rotor.TSR**2*mu**2/(1-a_new)**2))
        f_tip = 2/np.pi * np.arccos(exp)
        #Root correction
        exp = np.exp(-self.Rotor.n_blades/2 * ((mu-mu_root)/mu) * np.sqrt(1+self.Rotor.TSR**2*mu**2/(1-a_new)**2))
        f_root = 2/np.pi * np.arccos(exp)
        #Combined correction
        f = f_tip*f_root
        
        if f < 0.0001:
            f = 0.0001
        
        return f
    
    
    def Solver(self,N_iter_max = 100,delta=1e-5):
        
        if self.Rotor.yaw == 0:
            N_azimuth = 2
        else:
            N_azimuth = self.Rotor.N_azimuth
            
        self.Results = Results(self.Rotor.N_radial,N_azimuth)

                
        for i in range(self.Rotor.N_radial-1): #Loop of each blade section
            #Calculate blade parameters in this section
            mu = (self.Rotor.mu[i]+self.Rotor.mu[i+1])/2
            chord = np.interp(mu,self.Rotor.mu,self.Rotor.chord)
            beta = np.interp(mu,self.Rotor.mu,self.Rotor.beta)
            
            [self.Results.mu[i],self.Results.chord[i],self.Results.beta[i]] = [mu,chord,beta]
    
                
            for j in range(N_azimuth-1):
                    azimuth = (self.Rotor.azimuth[j]+self.Rotor.azimuth[j+1])/2
                    
                    a,ap = (0.2,0.2) #Initialize induction factors
                    for ite in range(N_iter_max):
                        #Velocities and angles
                        [u_tan,u_nor,u_rel,phi] = self.RelativeVelocities(a,ap,mu,azimuth)
                        alpha = phi - (self.Rotor.beta[i] + self.Rotor.theta)*np.pi/180
                        
                        #Airfoil forces
                        [cl,cd] = self.AirfoilCoefficients(alpha)
                        [lift,drag,f_tan,f_nor] = self.Forces(chord,phi,u_rel,cl,cd)
                        
                        #Thrust coefficient
                        CT = f_nor*self.Rotor.n_blades/(0.5*self.Rotor.rho*self.Rotor.wind_speed**2*2*np.pi*mu*self.Rotor.radius)
                        
                        #Get new value of axial induction factor
                        a_new = self.NewInductionFactor(CT,self.Rotor.yaw,a)
                        
                        #Apply the tip loss correction factor
                        f = self.PrandtlTipCorrection(mu,a_new)
                        a_new = a_new/f
                        
                        #Induction factor for the next iteration
                        a = 0.75*a + 0.25*a_new
                        
                        #Calculate tangential induction
                        ap_new = f_tan*self.Rotor.n_blades/(2*self.Rotor.rho*2*np.pi*mu*self.Rotor.radius*self.Rotor.wind_speed**2*(1-a)*self.Rotor.TSR*mu*f)
                        ap = 0.75*ap + 0.25*ap_new
                        
                        #Check convergency
                        if np.abs(a_new-a) < delta and np.abs(ap_new-ap):
                            break
                        
                    #Store all the results
                    [self.Results.a[i,j],self.Results.ap[i,j],self.Results.phi[i,j],self.Results.alpha[i,j],self.Results.cl[i,j],
                     self.Results.cd[i,j],self.Results.f_nor[i,j],self.Results.f_tan[i,j],self.Results.f[i,j],self.Results.ite[i,j]] = \
                        [a,ap,phi*180/np.pi,alpha*180/np.pi,cl,cd,f_nor,f_tan,f,ite]
              
                    
                        
def Plotting(Rotor,Results,Validation):
    
    mu = np.zeros((Rotor.N_radial-1))
    for i in range(Rotor.N_radial-1): #Loop of each blade section
       mu[i] = (Rotor.mu[i]+Rotor.mu[i+1])/2

    
    fig = plt.figure()
    plt.plot(mu,Results.a)
    plt.plot(Validation['r_R'],Validation['a'])
    plt.legend(['a','a (validation)'])
    
    fig = plt.figure()
    plt.plot(mu,Results.ap)
    plt.plot(Validation['r_R'],Validation['aline'])
    plt.legend(['ap','ap (validation)'])

    fig = plt.figure()
    plt.plot(mu,Results.f_nor)
    plt.plot(Validation['r_R'],Validation['fnorm']*1.225)
    plt.legend(['fnorm','fnorm (validation)'])
                        
    fig = plt.figure()
    plt.plot(mu,Results.f_tan)
    plt.plot(Validation['r_R'],Validation['ftan']*1.225)
    plt.legend(['ftan','ftan (validation)'])                        


class Optimizer:
    def __init__(self, Rotor_original, a):
        self.a = a
        self.R = Rotor_original.radius
        self.TSR = Rotor_original.TSR
        self.B = Rotor_original.n_blades
        self.mu = Rotor_original.mu
        
        #Calculate optimal Cl and E
        Cl = Rotor_original.polars['Cl']
        Cd = Rotor_original.polars['Cd']
        self.E = max(Cl/Cd)
        self.cl = Cl[np.argmax(Cl/Cd)]
        
        
        
    def residuals(self,x):
        
        c,ap = x #Unpack the input
        
        #Flow angle
        phi = m.atan((1-self.a)*self.R/((1+ap)*self.r*self.TSR))
        
        #Tip loss
        f = self.B * (self.R-self.r)/(2*self.r*np.sin(phi))
        F = 2*m.acos(np.exp(-f))/np.pi
        
        #Force coefficients
        Cy = self.cl * np.cos(phi) + self.cl/self.E*np.sin(phi)
        Cx = self.cl * np.sin(phi) - self.cl/self.E*np.cos(phi)
        
        #Solidity
        sigma = c*self.B/(2*np.pi*self.r)
     
        #Get residual c and ap
        res_c = 4*np.pi*self.r*m.sin(phi)**2*F*2*self.a/(Cy*self.B*(1-self.a)) - c
        res_ap = 1/(4*F*np.sin(phi)*np.cos(phi)/(sigma*Cx)-1) - ap
        
        return res_c,res_ap
    
    
    def ChordOpt(self):
        [self.c,self.ap] = np.zeros((2,len(self.mu)))
        for i in range(len(self.mu)):
            self.r = self.mu[i]*self.R #Radial position
            x0 = [3,0.001] #Initial guess
            bounds = ((0.0,0),(7,1)) #Lower and upper bounds
            results = least_squares(self.residuals,x0,bounds=bounds) #Calculate with the least-squares method the chord and a'
            self.c[i],self.ap[i] = results.x
            
            
        
        
    
    
            
Blade = Rotor()

Blade.SetOperationalData(1,8,30)
    
Blade_BEMT = BEMT(Blade)

Blade_BEMT.Solver()

#Integrate to get total CP
Blade_BEMT.Results.Integrate(Blade)
print('Thrust Coefficient=', Blade_BEMT.Results.CT)
print('Power Coefficient=', Blade_BEMT.Results.CP)



Validation = pd.read_csv('Validation/results.txt')

Plotting(Blade,Blade_BEMT.Results,Validation)


## Optimization
# We want to optimize for a CT = 0.75
CT_opt = 0.75
a_opt = 1/2 - np.sqrt(1-CT_opt)/2

Opti = Optimizer(Blade, a_opt)
Opti.ChordOpt()

fig = plt.figure()
plt.plot(Opti.mu*Opti.R,Opti.c)
plt.plot(Blade.mu*Blade.radius,Blade.chord)
plt.xlabel('Span [m]')
plt.ylabel('Chord [m]')
plt.grid()
plt.legend(['Optimal','Original'])
    
    
