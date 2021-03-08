## BEMT code using classes

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


class Rotor:
  def __init__(self):
      
    #Geometric data of the blade
    self.radius = 50 #[m]
    self.n_blades = 3
    self.theta = -2 #Pitch angle [deg]
    self.N_radial = 80 #Number of sections
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
        self.a,self.ap,self.phi,self.alpha,self.cl,self.cd,self.f_nor,self.f_tan,self.f,self.ite =  np.zeros((10,N_radial-1,N_azimuth-1))
        self.CT, self.CP, self.CQ, self.P, self.T, self.Q = np.zeros((6,1))
            
    
class BEMT:
    def __init__(self,Rotor):
        self.Rotor = Rotor
    
    def RelativeVelocities(self,a,ap,mu,azimuth=0):
        u_tan = self.Rotor.omega*self.Rotor.radius*mu*(1+ap) + self.Rotor.wind_speed*np.sin(self.Rotor.yaw)*np.sin(azimuth)
        u_nor = self.Rotor.wind_speed*(np.cos(self.Rotor.yaw)-a)
        
        # Testing Glauert correction
        psi = (0.6*a+1)*self.Rotor.yaw
        K = 2*np.tan(psi/2)
        
        #u_nor = self.wind_speed*(np.cos(self.yaw)-a*K*mu*np.sin(self.yaw))

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

            
            
Blade = Rotor()

Blade.SetOperationalData(1,8,0)
    
Blade_BEMT = BEMT(Blade)

Blade_BEMT.Solver()

Validation = pd.read_csv('Validation/results.txt')

Plotting(Blade,Blade_BEMT.Results,Validation)

    
    
