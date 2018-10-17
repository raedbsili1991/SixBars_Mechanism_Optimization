
# coding: utf-8

# In[2]:

#################################################################################
import numpy as np
import matplotlib.pyplot as plt

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


#################################################################################
fx_array,NoG_array,NumberOfGeneration = array('d'),array('d'),1

l1_array,l2_array,l3_array,alpha_array = array('d'),array('d'),array('d'),array('d')

isotropy_array, ws_array = array('d'),array('d')

#################################################################################
#------------------------------------------------------------------------------
# Decalre symbols

l1, l2, l3, alpha = symbols('l1 l2 l3 alpha') # link arc angles

# Inputs angle
theta_L, theta_R = symbols('theta_L theta_R')
#Zeros 
theta_L0, theta_R0 = symbols('theta_L0 theta_R0')
uL0, vL0 = symbols('uL0 vL0')
uR0, vR0 = symbols('uR0 vR0')
# Passive Angles
u_L, u_R, v_L, v_R= symbols('u_L u_R v_L v_R')  #passive angles 
# Euler angles
x, y, z = symbols('x y z') # platform coordinates
theta,psi,phi = symbols('theta psi phi') # platform roll, pitch, yaw 
#------------------------------------------------------------------------------

Rl = symbols('Rl')

#-------------------------------------------------------------------------------

# RotationSym Function
def RotationSym(angle, axis):
    
    if axis == 'Z':
        ROT = Matrix([[cos(angle), -sin(angle), 0, 0], [sin(angle), cos(angle), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        
    elif axis == 'Y':
        ROT = Matrix([[cos(angle), 0, sin(angle), 0], [0, 1, 0, 0], [-sin(angle), 0, cos(angle), 0], [0, 0, 0, 1]])
        
    elif axis == 'X':
        ROT = Matrix([[1, 0, 0, 0], [0, cos(angle), -sin(angle), 0], [0, sin(angle), cos(angle), 0], [0, 0, 0, 1]])
    
    else:
        print('Incorrect axis definition. Use: "X,Y,Z"')
    
    return ROT
#------------------------------------------------------------------------------
# TranslationSym Function
def TranslationSym(position):
    
    x, y, z = position[0], position[1], position[2]
    TRS = Matrix([[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, z], [0, 0, 0, 1]])
    
    return TRS

def VectorSym(position):
    x, y, z = position[0], position[1], position[2]
    return Matrix([[x], [y], [z]])
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Compute Matrix C for Platform
RZc = RotationSym(phi, 'Z')
RYc = RotationSym(theta, 'Y')
RXc = RotationSym(psi, 'X')

ROTc = RYc*RZc
POSc = TranslationSym([x, y, z])

C = POSc*ROTc
#------------------------------------------------------------------------------



theta_pitch, phi_yaw = symbols('theta_pitch phi_yaw')

RO_L = RotationSym(theta_pitch,'Y')*RotationSym(phi_yaw,'Z')*RotationSym(-pi/2,'Z')*RotationSym(pi/2 - l3,'X')*RotationSym(alpha,'Z')
VL = (RO_L[:3,:3])[:,1]
UL = ((RotationSym(theta_L, 'Y')*RotationSym(l1, 'X'))[:3,:3])[:,1]

#--

RO_R = RotationSym(theta_pitch,'Y')*RotationSym(phi_yaw,'Z')*RotationSym(-pi/2,'Z')*RotationSym(pi/2 - l3,'X')*RotationSym(-alpha,'Z')
VR = (RO_R[:3,:3])[:,1]
UR = ((RotationSym(pi,'Z')*RotationSym(theta_R, 'Y')*RotationSym(l1, 'X'))[:3,:3])[:,1]

fL = (Transpose(VL)*UL)[0]-cos(l2)
fR = (Transpose(VR)*UR)[0]-cos(l2)

fL = collect(trigsimp(fL,deep='True'),sin(theta_L-theta_pitch))
fR = collect(trigsimp(fR,deep='True'),cos(theta_R+theta_pitch))
f3 = (Transpose(VL)*VR)[0]-cos(2*alpha)

dataL = (theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

UL_,VL_ = lambdify(dataL,UL),lambdify(dataL,VL)
UR_,VR_ = lambdify(dataL,UR),lambdify(dataL,VR)

FL_ = lambdify(dataL,fL)
FR_ = lambdify(dataL,fR)
F3_ = lambdify(dataL,f3)


gL = (Transpose(UL.cross(VL))*VectorSym([0,1,0]))[0,0]
gL = simplify(gL)
gR = (Transpose(UR.cross(VR))*VectorSym([0,1,0]))[0,0]
gR = simplify(gR)


jacL_constraint = Matrix([gL]).jacobian([theta_pitch,phi_yaw])
jacR_constraint = Matrix([gR]).jacobian([theta_pitch,phi_yaw])

GL_ = lambdify(dataL,gL)
GR_ = lambdify(dataL,gR)
JLcs_ = lambdify(dataL,jacL_constraint)
JRcs_ = lambdify(dataL,jacR_constraint)

F = pow(fL,2) + pow(fR,2)
Jf = Matrix([F]).jacobian([theta_pitch,phi_yaw])
Hf = hessian(F,(theta_pitch,phi_yaw))

F_ = lambdify(dataL,F)
JF_ = lambdify(dataL,Jf)
HF_ = lambdify(dataL,Hf)

d_FL_L,d_FL_R = diff(fL,theta_L),diff(fL,theta_R)
d_FR_L,d_FR_R = diff(fR,theta_L),diff(fR,theta_R)

d_FL_P,d_FL_Y = simplify(diff(fL,theta_pitch)),diff(fL,phi_yaw)
d_FR_P,d_FR_Y = simplify(diff(fR,theta_pitch)),diff(fR,phi_yaw)

Ji = Matrix([ [d_FL_L,d_FL_R],
        [d_FR_L,d_FR_R] ])
Jo = Matrix([ [d_FL_P,d_FL_Y],
        [d_FR_P,d_FR_Y] ])
dataL = (theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

Ji_ = lambdify(dataL,Ji)
Jo_ = lambdify(dataL,Jo)

J = -Jo.inv()*Ji

J_ = lambdify(dataL,J)

#################################################################################


def grid(x, y, z, resX=100, resY=100):

    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

def SolveFK(E,data):
    
    theta_pitch,phi_yaw = E[0],E[1]
    theta_L,theta_R,l1,l2,l3,alpha = data 
    
    cost = [FL_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha),
           FR_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)]
     
    return cost

def SolveFK_min(E,data):
    
    theta_pitch,phi_yaw = E[0],E[1]
    theta_L,theta_R,l1,l2,l3,alpha = data 

    return F_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

def Jac_F(E,data):
    theta_pitch,phi_yaw = E[0],E[1]
    theta_L,theta_R,l1,l2,l3,alpha = data 
    return JF_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

def Hess_F(E,data):
    theta_pitch,phi_yaw = E[0],E[1]
    theta_L,theta_R,l1,l2,l3,alpha = data 
    return F_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)


def f_cons(E,theta_L,theta_R,l1,l2,l3,alpha):
    
    theta_pitch,phi_yaw = E[0],E[1]
    
    ul = UL_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)
    vl = VL_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

    C = np.cross(ul.reshape(1,3), vl.reshape(1,3))[0][1]
     
    return C

def fl_cons(E,theta_L,theta_R,l1,l2,l3,alpha):
    theta_pitch,phi_yaw = E[0],E[1]
    return GL_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

def fr_cons(E,theta_L,theta_R,l1,l2,l3,alpha):
    theta_pitch,phi_yaw = E[0],E[1]   
    return GR_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)

def Jac_fl_cons(E,theta_L,theta_R,l1,l2,l3,alpha):
    theta_pitch,phi_yaw = E[0],E[1]
    return JLcs_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)
def Jac_fr_cons(E,theta_L,theta_R,l1,l2,l3,alpha):
    theta_pitch,phi_yaw = E[0],E[1]
    return JRcs_(theta_pitch,phi_yaw,theta_L,theta_R,l1,l2,l3,alpha)




#################################################################################



def GetEuler(thetaL_corr,thetaR_corr,dataP):

    
    bnd = [ (np.radians(-180),np.radians(180)),
       (np.radians(-180),np.radians(180)) ]
    
    LL1,LL2,LL3,ALPHA = dataP
    
    m,M = 0,np.size(thetaL_corr)
    pitch_,yaw_ = np.zeros(np.size(np.arange(m,M,1))),np.zeros(np.size(np.arange(m,M,1)))
    
    guess0 = np.radians([0,0])
    
    for k in np.arange(m,M,1):

        datak = [thetaL_corr[k],thetaR_corr[k],LL1,LL2,LL3,ALPHA]
   
        cns = ({'type': 'ineq', 'fun':fl_cons,'jac':Jac_fl_cons,'args':datak},
              {'type': 'ineq', 'fun':fr_cons,'jac':Jac_fr_cons,'args':datak})

        sol = minimize(SolveFK_min,guess0,args=datak,jac=Jac_F,
                       bounds=bnd,constraints=cns,tol=1e-10)

        pitch,yaw = sol.x[0],sol.x[1]


        pitch_[k],yaw_[k] = pitch,yaw
        
    return pitch_,yaw_

def SolveThetas0(theta_L0,data):
    
    l1,l2,l3,alpha = data
    fL0 = np.sin(alpha)*np.cos(l1) + np.sin(l1)*np.sin(l3)*np.sin(theta_L0)*np.cos(alpha) + np.sin(l1)*np.cos(alpha)*np.cos(l3)*np.cos(theta_L0) - np.cos(l2)
    
    return fL0

def GetThetas0(data):
    
    LL1,LL2,LL3,ALPHA = data
    sol = sc.root(SolveThetas0,0,args=[LL1,LL2,LL3,ALPHA],method='broyden2')
    THETAL0 = math.fmod(sol.x,math.pi)
    THETAR0 = - THETAL0
    
    return THETAL0,THETAR0

def Get_delta(J): 
    m = 2
    M = np.power(LA.det(np.dot(J,J.T)),1/m)
    P = (1/m)*np.trace(np.dot(J,J.T))
    return M/P

def show_contours(theta1cr,theta5cr,delta,xtitle,ytitle,title_,N):
    
    T1,T5,DELTA = grid(theta1cr,theta5cr,delta, resX=100, resY=100)
    cs = plt.contourf(T1,T5,DELTA,N,cmap=plt.cm.jet)
    plt.title(title_)
    plt.colorbar()
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.show()

    #plt.show()
    return cs

def get_contours(theta1cr,theta5cr,delta,N):
    
    T1,T5,DELTA = grid(theta1cr,theta5cr,delta, resX=100, resY=100)
    cs = plt.contourf(T1,T5,DELTA,N,cmap=plt.cm.jet)

    return cs

def Get_Area(x_axis,y_axis):
    return np.abs(0.5*np.sum(y_axis[:-1]*np.diff(x_axis) - x_axis[:-1]*np.diff(y_axis)))

def GetGlobalIsotropy(cs):
    
    d = cs.levels
    ro, d_areas = np.zeros(np.size(d)), np.zeros(np.size(d))

    for i in range(np.size(d)-1): 
            x= cs.collections[i].get_paths()[0].vertices[:,0]
            y = cs.collections[i].get_paths()[0].vertices[:,1]
            points = np.array([x,y])
            d_areas[i] = ConvexHull(points.T).volume
                
    ro = np.sum(np.multiply(d,d_areas))/np.sum(d_areas)
    
    return ro 

def GetIsotropy(pitch_,yaw_,thetaL_corr,thetaR_corr,data,get_uniformity):
    
    LL1,LL2,LL3,ALPHA = data
    m,M = 0,np.size(pitch_)
    delta_ = np.zeros(np.size(np.arange(m,M,1)))

    for k in np.arange(m,M,1):
        Jac = J_(pitch_[k],yaw_[k],thetaL_corr[k],thetaR_corr[k],LL1,LL2,LL3,ALPHA)
        delta_[k] = get_uniformity(Jac)
        
    return delta_


def Generate_Input(M,n):
    
    tL,tR =  np.radians(np.linspace(-M,M,n)),np.radians(np.linspace(-M,M,n))
    # -------------------
    TETAL, TETAR = np.meshgrid(tL,tR,indexing='ij')
    tL_,tR_  = np.ravel(TETAL),np.ravel(TETAR)
   
    return np.degrees(tL_),np.degrees(tR_)




#################################################################################


def IsoCons_1(P):
    l1,l2,l3,alpha = P[0],P[1],P[2],P[3]
    return np.sin(l1)*np.cos(alpha) - (np.cos(l2)-np.sin(alpha)*np.cos(l1))

def f_isotropy_evo(P,*data):
            
    LL1,LL2,LL3,ALPHA = P[0],P[1],P[2],P[3]
    filename_all,Input,N,iThetaL,iThetaR,f,get_uniformity = data
    dataP = (LL1,LL2,LL3,ALPHA)
    
    global fx_array, NoG_array, NumberOfGeneration
    global l1_array,l2_array,l3_array,alpha_array
    global isotropy_array, ws_array
    
    if IsoCons_1([LL1,LL2,LL3,ALPHA]) >= 0:
        
        # Solve Correction Angle 
        print("Solving Candidate-LL1,LL2,LL3,ALPHA",np.degrees([LL1,LL2,LL3,ALPHA]))
        
        THETAL0,THETAR0 = GetThetas0([LL1,LL2,LL3,ALPHA])
        iThetaL_corr,iThetaR_corr = iThetaL+THETAL0, iThetaR+THETAR0

        #FK
        pitch_,yaw_ = GetEuler(iThetaL_corr,iThetaR_corr,dataP)
        
        if np.isnan(pitch_).any() or np.isnan(yaw_).any() :
            Loss= np.inf 
            DELTA_NUMBER = 0
            OMEGA_NUMBER = 0            
        else: 
            

        
            # Isotropy
            delta_ = GetIsotropy(pitch_,yaw_,iThetaL_corr,iThetaR_corr,dataP,get_uniformity)
                            
            # Workspace Area Evaluation
            points = np.array([np.degrees(pitch_),np.degrees(yaw_)])
            WS_Area_Number = ConvexHull(points.T).volume/(2*Input)**2
            
            if (np.isnan(delta_).any()) or (WS_Area_Number > 1): 
                Loss= np.inf
                DELTA_NUMBER = 0
                OMEGA_NUMBER = 0
            else: 
                
                # Global Isotropy
                csi = get_contours(np.degrees(iThetaL),np.degrees(iThetaR),delta_,N)
                Global_Iso = GetGlobalIsotropy(csi)
        
                #Lossto minimize
                DELTA_NUMBER = np.min(delta_)*Global_Iso*np.max(delta_)

                Loss = f(delta_,Global_Iso,WS_Area_Number) 
                                
                fx_array.append(Loss)
                l1_array.append(np.degrees(LL1))
                l2_array.append(np.degrees(LL2))  
                l3_array.append(np.degrees(LL3))  
                alpha_array.append(np.degrees(ALPHA))
    
                isotropy_array.append(DELTA_NUMBER)
                ws_array.append(WS_Area_Number)
    
    
                print("Min/Global Iso: ",[np.min(delta_),Global_Iso])
                print("WN = ",WS_Area_Number)
                print("\n Current Loss: ",Loss)
  
    else: 
        
        Loss= np.inf 
        DELTA_NUMBER = 0
        OMEGA_NUMBER = 0
    

    
    df = pd.DataFrame({'fx':fx_array,'l1':l1_array,'l2':l2_array,'l3':l3_array,'Alpha':alpha_array,
                       'ISO':isotropy_array,'WS':ws_array},
                     columns = ['fx', 'l1', 'l2', 'l3','Alpha','ISO','WS'])
    
    df.to_csv(filename_all)
            
    
    return Loss

#################################################################################

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)
        
        

