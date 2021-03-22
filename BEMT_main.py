## Blade Element Momentum Code
## TU Delft - Rotor Wake Aerdoynamics - Assignment 1
## Group 1: Mathieu Pellé, Nils Gaukroger, Guillem Vergés i Plaza

#%% Import the necessary libraries and functions

import numpy as np
import matplotlib.pyplot as plt
from BEMT_Utilities import Rotor,BEMT,Optimizer,MeshSensitivity,plot_optimized_geometry,plot_mesh_sensitivity
from fplot import plot_enthalpy_tube,plot_yaw,plot_TSR,plot_polars,plot_correction
    
    
#%% Define the original rotor and calculate optimized geometry

#Initialize the rotor that we will analyze
Rotor_org = Rotor()

#Create an optimized rotor at CT=0.75
CT_opt = 0.75
a_opt = 1/2 - np.sqrt(1-CT_opt)/2 #Calculate the corresponding value of a as an initial guess

for i in range(100):
    #Get the blade geometry
    Optimal_geo = Optimizer(Rotor_org, a_opt, TSR = 8)
    
    #Creat a rotor class with that geometry
    Rotor_opt = Rotor(Optimized_geometry=Optimal_geo)
    
    #Solve the BEMT with this rotor
    BEMT_opt = BEMT(Rotor_opt)
    BEMT_opt.Solver()
    
    #Evaluate the CT obtained
    CT_res = BEMT_opt.Results.CT
    if abs(CT_opt-CT_res)<0.001:
        break
    
    #If it's not the desired one, correct the axial induction factor for the next iteration
    if i>0:  
        x = (CT_res-CT_opt)/(CT_res-CT_old/(a_opt-BEMT_opt.Results.a_global))
        a_opt = a_opt + x
    else:
        a_opt = a_opt**2/BEMT_opt.Results.a_global
       
    CT_old = CT_res
    
    fig = plt.figure('Convergency')
    plt.plot(i,CT_res,'o')

    print(CT_res)
    print('New a_opt',a_opt,'Result a_global',BEMT_opt.Results.a_global)
    

#%% Evalute both rotors under different operational conditions

#Define operational conditions to be analyzed
TSR_list =  [6, 8, 10]
yaw_angles_list = [0, 15, 30]
wind_speed = 10

#Initializize dictionary to store results
Res_org = {} 
Res_opt = {}

#Initialize index for status message during the loop
idx = 1 

#Main loop: solving all operational conditions for both rotors
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
        

#%% Generate CP-Pitch-Lambda plots for both turbines

TSR_list =  list(np.linspace(5,12,10))
theta_list = list(np.linspace(-6,0,10))

CpLambda_org = BEMT_org.CpLambda(TSR_list,theta_list)
CpLambda_opt = BEMT_opt.CpLambda(TSR_list,theta_list)


#%% Mesh Sensitivity Analysis

N_array = np.geomspace(10,400,10,dtype=int) #Number of radial points that we will test (integrer dtype to avoid decimals)

[CT_lin,err_lin,N_chosen_lin,execution_time_lin] = MeshSensitivity(N_array, Spacing = 'lin')
[CT_cos,err_cos,N_chosen_cos,execution_time_cos] = MeshSensitivity(N_array, Spacing = 'cos')

#%% Effect of the Prandt'l tip and root correction

Rotor_org.SetOperationalData(wind_speed = 10,TSR = 8,yaw = 0) #Set the default operational conditions
BEMT_NoPrandtl = BEMT(Rotor_org) #Initialize BEMT instance
BEMT_NoPrandtl.Solver(Prandtl_correction = False) #Solve the rotor
Res_NoPrandtl = BEMT_NoPrandtl.Results #Store the no-prandtl correction results
Res_Prandtl = Res_org['TSR8_yaw0'] #Store the corrected results
      
        
#%% Plotting functions

x = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * x), 33.11 * .5**(.5 * x)]   )
#plt.rc('text', usetex=False)
plt.rc('font', family='serif')
global save
save=True


plot_optimized_geometry(Rotor_org,Rotor_opt,Res_org,Res_opt,CpLambda_org,CpLambda_opt)
plot_mesh_sensitivity(N_array,CT_lin,CT_cos,N_chosen_lin,N_chosen_cos,execution_time_lin,execution_time_cos,err_lin,err_cos)
plot_TSR(Res_org,Rotor_org,[6,8,10])
plot_yaw(Res_org,Rotor_org,[0,15,30])
plot_polars(Rotor_org)
plot_correction(Res_org, Rotor_org, 8, 0,Res_Prandtl,Res_NoPrandtl)
plot_enthalpy_tube(Res_org,Rotor_org,8,0)

