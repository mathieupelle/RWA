## BEMT code using classes
import numpy as np
#import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import math as m
from BEMT_Utilities import Rotor,BEMT,Results,Optimizer,Plotting
#from Mathieu import NicePlots
    
    
            
    
#Initialize the rotor that we will analyze
Rotor_org = Rotor()

#Create an optimized rotor at CT=0.75
CT_opt = 0.75
a_opt = 1/2 - np.sqrt(1-CT_opt)/2 #Calculate the corresponding value of a

Optimal_geo = Optimizer(Rotor_org, a_opt, TSR = 8)
Rotor_opt = Rotor(Optimal_geo)

#Test different operational conditions
TSR_list =  [6, 8, 10]
yaw_angles_list = [0, 15, 30]
wind_speed = 10

Res_org = {} #Initializize dictionary to store results
Res_opt = {}

idx = 1


for TSR in TSR_list:
    for yaw_angle in yaw_angles_list:
        
        #Set operational data to each rotor
        Rotor_org.SetOperationalData(wind_speed,TSR,yaw_angle)
        Rotor_opt.SetOperationalData(wind_speed,TSR,yaw_angle)
        
        #Compute BEMT
        print('Computing: TSR = ', TSR, 'Yaw = ', yaw_angle, '-->[Case ', idx, '/',len(TSR_list)*len(yaw_angles_list),']')
        var = ["TSR" + str(TSR) + "_yaw" + str(yaw_angle)]
        BEMT_org = BEMT(Rotor_org)
        BEMT_opt = BEMT(Rotor_opt)
        BEMT_org.Solver()
        BEMT_opt.Solver()
        
        Res_org[var[0]] = BEMT_org.Results
        Res_opt[var[0]] = BEMT_opt.Results
        
        idx = idx+1
        


#Generate CP-Pitch-Lambda plots for both turbines
TSR_list =  list(np.linspace(6,12,10))
theta_list = list(np.linspace(-7,0,10))

CpLambda_org = BEMT_org.CpLambda(TSR_list,theta_list)
CpLambda_opt = BEMT_opt.CpLambda(TSR_list,theta_list)

       
        

Plotting(Rotor_org,Rotor_opt,Res_org,Res_opt,CpLambda_org,CpLambda_opt)

## Optimization
# We want to optimize for a CT = 0.75


# fig = plt.figure()
# plt.plot(Opti.mu*Opti.R,Opti.c)
# plt.plot(Blade.mu*Blade.radius,Blade.chord)
# plt.xlabel('Span [m]')
# plt.ylabel('Chord [m]')
# plt.grid()
# plt.legend(['Optimal','Original'])
    
    
