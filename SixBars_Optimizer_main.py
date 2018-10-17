######### 

# Author: Raed Bsili 
# Supervisors: Alberto Parmiggiani & Giorgio Metta     
# Istituto Italiano di Tecnologia (IIT), iCub
# Octobre 2018

#################################################################################
import numpy as np
import matplotlib.pyplot as plt
import time
import sympy as sy
from sympy import *
from sympy import Matrix, sin, cos, tan, pi, symbols, solve, nsolve
init_printing(use_unicode=True)

from matplotlib.mlab import griddata
from numpy import linalg as LA
 

########
from scipy.optimize import root, minimize, differential_evolution
import scipy.optimize as sc
#######

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#####
import sympy as sy
from sympy import *
from sympy import Matrix, sin, cos, tan, pi, symbols, solve, nsolve
init_printing(use_unicode=True)
###
import math
from scipy.spatial import ConvexHull
from array import array 
import pandas as pd
import os
from scipy.optimize._differentialevolution import DifferentialEvolutionSolver as DES 
###
f_Energies = array('d')
best_l1,best_l2,best_l3,best_alpha = array('d'),array('d'),array('d'),array('d')
import SixBars_OptimizerFunctions as OF 
###


#best_parameters.csv log the best candidate per generation, this is useful for ploting the pareto front for example
#all_parameters.csv log all information about the every candidate per generation. 
########################### I N I T I A L I Z I N G  #########################
Amplitude = 65
dirname,filename_best,filename_all = 'Log_dir','/best_parameters','/all_parameters'
N_contour = 50 # This defines the contour ploting resolution, keep this constant! 
N_grid = 25 # This defines the square grid input resolution 


# Preparing the file to log into it optimization information 
iThetaL,iThetaR = OF.Generate_Input(Amplitude,N_grid)
OF.createFolder(dirname)
TargetFile_best = dirname + filename_best
TargetFile_all = dirname + filename_all
###############

########################### U N I F O R M I T Y     I N D E X     D E S I G N  #########################
#Define here your definition to the mechanism uniformity, J is an mxn matrix! 
def unformity_function(J):
    
    m = np.shape(J)[0] # m is the dimension of the JJ.T matrix 
    M = np.power(LA.det(np.dot(J,J.T)),1/m)
    P = (1/m)*np.trace(np.dot(J,J.T))
    index  = M/P 
    
    return index 

########################### L O S S     F U N C T I O N     D E S I G N  #########################

# Design the loss function here, keep in mind that delta_ should be a vector  
def loss_function(delta_,Global_Iso,WS_Number):
     
    loss = 1 - np.min(delta_)*Global_Iso#*WS_Number

    return loss

########################### O P T I M I Z E R     D E S I G N  #########################
minimum_bnd, maximum_bnd = np.radians(9),np.radians(85)
maximum_generations_number = 100 

Bnds = [(minimum_bnd,maximum_bnd),(minimum_bnd,maximum_bnd),
       (minimum_bnd,maximum_bnd),(minimum_bnd,maximum_bnd)]

data = (TargetFile_all,Amplitude,N_contour,np.radians(iThetaL),np.radians(iThetaR),
        loss_function,unformity_function)

# Tune the DE optimizer here 
solver  = DES(OF.f_isotropy_evo,bounds=Bnds,args=data,
                              polish=True,
                              init='latinhypercube',
                              mutation = (1.5,1.9),
                              recombination = 0.1,
                              maxiter = maximum_generations_number, 
                              disp=True)



########################### S T A R T I N G     O P T I M I Z E R   #########################
print("Running Optimizer . . .")


time_start = time.clock()


for i in range(maximum_generations_number):
    P, e = next(solver)
    
    best_l1.append(np.degrees(P[0]))
    best_l2.append(np.degrees(P[1]))
    best_l3.append(np.degrees(P[2]))
    best_alpha.append(np.degrees(P[3]))
    
    f_Energies.append(e)
    print("[Step,Energy] = ", [i,e])
    df = pd.DataFrame(data = {'f_Energy_Pareto':f_Energies,
                              'best_l1':best_l1,'best_l2':best_l2,
                              'best_l3':best_l3,'best_alpha':best_alpha},
                     columns = ['f_Energy_Pareto', 'best_l1', 'best_l2', 'best_l3','best_alpha'])
    df.to_csv(TargetFile_best)
    if i >= maximum_generations_number - 1 :
        print("Final Optimum is: ", np.degrees(P))
        print("With a final energy: ", e)
     


    ######################   
    

print("Computation time [s] = " (time.clock() - time_start) )    