# -*- coding: utf-8 -*-
"""
Created on Sun May 18 12:34:17 2014

@author: redhotkenny
"""

# --------------------------------------------------------------
# -- Load modules
import numpy as np
import time as time
#from matplotlib import pyplot as plt
from time import *
from math import *
#from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
from scipy.sparse.linalg import cgs
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------
# -------------------- Declare functions ------------------------
# ---------------------------------------------------------------


# -- Van Genuchten functions

def theta_curve(h,soil_carac,snum):
    #h is pressure head
    #soil_carac is the list of vectors of soil characteristics. Each vector has as the same size as the number on nodes.
    #snum is the number of the node at the moment   #snum had to be added in rder to have soil layers with different properties
    h = float(h)
    n = soil_carac['n'][snum]
    theta_r = soil_carac['theta_r'][snum]
    theta_s = soil_carac['theta_s'][snum]
    alpha = soil_carac['alpha'][snum]
    m = 1. - 1./n
    if h < 0 : 
        return( theta_r + (theta_s - theta_r) * (1. + abs(alpha*h)**n)**(-m) )
    else :
        return( theta_s )


def C_curve(h,soil_carac,snum):
    #h is pressure head
    #soil_carac is the list of vectors of soil characteristics. Each vector has as the same size as the number on nodes.
    #snum is the number of the node at the moment   #snum had to be added in rder to have soil layers with different properties
    h = float(h)
    n = soil_carac['n'][snum]
    theta_r = soil_carac['theta_r'][snum]
    theta_s = soil_carac['theta_s'][snum]
    alpha = soil_carac['alpha'][snum]
    m = 1. - 1./n
    if h < 0 : 
        return(  m * n * alpha**n * abs(h) **(n-1.) * (theta_s-theta_r) * (1. + abs(alpha*h)**n)**(-m-1.)   )
    else :
        return(0)


def Kr_curve(h,soil_carac,snum):
    #h is pressure head
    #soil_carac is the list of vectors of soil characteristics. Each vector has as the same size as the number on nodes.
    #snum is the number of the node at the moment   #snum had to be added in rder to have soil layers with different properties
    h = float(h)
    n = soil_carac['n'][snum]
    theta_r = soil_carac['theta_r'][snum]
    theta_s = soil_carac['theta_s'][snum]
    alpha = soil_carac['alpha'][snum]
    m = 1. - 1./n
    if h < 0 : 
        return( ( 1. - abs(alpha*h)**(n-1.)*(1.+abs(alpha*h)**n)**(-m) )**2. / 
        (( 1. + abs(alpha*h)**n  )**(0.5*m)) )
    else : 
        return(1.)



# --------------------------------------------------------------
# --------------------------------------------------------------

# Vector K, hydraulic conductivity
def get_K(h,soil_carac):
    K = np.zeros((h.size,))
    for i in range(h.size):
        K[i] =  soil_carac['Ksat'][i]*Kr_curve(h[i],soil_carac,i)
    return(K)

# Vector K, hydraulic conductivity
def get_K(h,soil_carac):
    K = np.zeros((h.size,))
    for i in range(h.size):
        K[i] =  soil_carac['Ksat'][i]*Kr_curve(h[i],soil_carac,i)
    return(K)


# vector Kp, hydraulic cond. of interface
# Following Clement et al. 1994 (J. Hydrol.)
# Equivalent Kp is arithmetic mean
# 0 0->1 1->2 2->3 3
def get_Kp(K):
    Kp = np.zeros((K.size+1,))
    Kp[0]  = K[0]
    Kp[-1] = K[-1]
    for i in range(1,K.size):
        Kp[i] = ( K[i-1] + K[i] ) / 2.
    return(Kp)

# vector theta
def get_theta(h,soil_carac):
    theta = np.zeros((h.size,))
    for i in range(h.size):
        theta[i] =  theta_curve(h[i],soil_carac,i)
    return(theta)

# vector C
def get_C(h, soil_carac):
    C = np.zeros((h.size,))
    for i in range(h.size):
        C[i] =  C_curve(h[i],soil_carac,i)
    return(C)

#Root Distribution
    #defines the weight of the effect of roots at each node depending of the chosen distribution
def get_RD(L, I, LR, z, dz, type):    #type can be homogeneous or follow Hoffman and van Genuchten, 1983 distribution (Trapezoidal) 
    RD=np.zeros(I)
    if type =='homogeneous':
        for i in range(I):
            if z[i] < (L-LR):
                RD[i]= 0.
            else:
                RD[i]=1.*I/(I-1)
    elif type == 'distribution':
        for i in range(I):
            if z[i] > (L-0.2*LR):
                RD[i] = 1.667/LR
            elif z[i] < (L-LR):
                RD[i]=0.
            else:
                RD[i]=2.0833/LR*(1-(L-z[i])/LR)
    RD=RD/np.sum(RD)
    return(RD)


# coefficient m1 (a in Clement et al., 1994)
def get_m1(i,Kp):
    return( Kp[i]/dz**2 )

# coefficient m2 (c in Clement et al., 1994)
def get_m2(i,Kp,C,theta):
    return( - ( get_m1(i,Kp) + get_m3(i,Kp) + get_b1(i,C) + get_b2(i,theta,soil_carac)  ) )

# coefficient m3 (e in Clement et al., 1994)
def get_m3(i,Kp):
    return( Kp[i+1]/dz**2 )

# coefficient b1 (f1 in Clement et al., 1994)
def get_b1(i,C):
    return(C[i] / dt)

# coefficient b2 in (Clement et al., 1994) 
def get_b2(i,theta,soil_carac):
    Ss = soil_carac['Ss'][i]
    eta = soil_carac['eta'][i]
    return( (Ss*theta[i]) / (eta*dt) )

# coefficient b3 (g in Clement et al., 1994)
def get_b3(i,K):
    #return(-(Kp[i+1]-Kp[i])/dz)
    return( -( K[i+1]-K[i-1] ) / (2.*dz) )
    #if i == 0 : 
        #return( -( K[i+1]-K[i] ) / (dz) )
    #return( -( K[i+2]-K[i] ) / (2.*dz) )
    #else :     
        #return( -( K[i+1]-K[i-1] ) / (2.*dz) )

# coefficient b4 (h in Clement et al., 1994)
def get_b4(i,theta0,theta):
    return( (theta[i] - theta0[i] )/dt )

#build Potential Water Uptake Distribution
def w_uptake(theta0, RD, q, soil_carac, dz):  #from Feddes et al., 1978
    #defines the weight of the effect of root uptake depending on stress conditions following the Feddes distribution 
    WU= np.zeros(I)
    for i in range(I):
        #if theta0[i]>=soil_carac['FC2']:
        #    red=(soil_carac['theta_s']-theta0[i])/(soil_carac['theta_s']-soil_carac['FC2'])
        if theta0[i]>=soil_carac['FC1'][i]:
            red=1.
        elif theta0[i]>=soil_carac['theta_r'][i]:
            red=(-soil_carac['theta_r'][i]+theta0[i])/(-soil_carac['theta_r'][i]+soil_carac['FC1'][i])
        elif theta0[i]<soil_carac['theta_r'][i]:
            red=0.
        WU[i]=red*RD[i]*q*I/((I-1)*dz)
    return(WU)

#build Flux             #Follow Darcy Law and mass conservancy equation from Celia, as used in Hydrus
def get_flux(Kp,h,dz,theta,theta0,dts,WU):
    #Defines the flux at each node, it is only used for the mass balance!
    Fl = np.zeros(h.size)
    #Fl[0] = -Kp[0]*((h[1]-h[0])/dz+1)
    Fl[0] = -Kp[1]*((h[1]-h[0])/dz+1)-dz/2*((theta[0]-theta0[0])/dts+WU[0]) 
    Fl[-1] = -Kp[-2]*((h[-1]-h[-2])/dz+1)-dz/2*((theta[-1]-theta0[-1])/dts+WU[-1])
    for i in range(1,(h.size-1)): 
        Fl[i]=(-Kp[i+1]*((h[i+1]-h[i])/dz+1)-Kp[i]*((h[i]-h[i-1])/dz+1))/2
    return (Fl)

#build SUM functions            #Arithmetic mean
def get_SUM(x,dz,I):
    #Is similar to the function Kp but is used with the other variables to get the total sum, maybe should be optimized
    X=np.zeros(x.size-1)
    for i in range(0,I-1):
        X[i]=(x[i+1]+x[i])/2.
    return(X.sum())


##Get the soil characteristics matrix (independant, only used for data entry) 
def get_soil_carac(numlay,L,dz,LC,soil_data,z,pp):   
    #numlay: number of layers: (1) one soil layer. (3) Ksat distribution (any other value) Two soil layers. ###Should be optimized for more than 2 layers
    #L: total depth
    #dz: spatial discretization
    #LC: deep when respective layer starts
    #soil_data: List with soil characteristics
    #z: vector with coordinates of nodes
    #pp: parameter used only if variable K distribution is chosen

    I = int(round(L/dz)+1)
    if numlay==1:
        soil_matrix={'Ksat':np.array([soil_data['Ksat']]*I),'Ss':np.array([soil_data['Ss']]*I),'eta':np.array([soil_data['eta']]*I), 'theta_r':np.array([soil_data['theta_r']]*I), 'theta_s':np.array([soil_data['theta_s']]*I), 'n':np.array([soil_data['n']]*I), 'alpha':np.array([soil_data['alpha']]*I), 'FC1':np.array([soil_data['FC1']]*I), 'FC2':np.array([soil_data['FC2']]*I)}
    elif numlay==10:
        soil_matrix={'Ksat':np.zeros(I),'Ss':np.zeros(I),'eta':np.zeros(I), 'theta_r':np.zeros(I), 'theta_s':np.zeros(I), 'n':np.zeros(I), 'alpha':np.zeros(I), 'FC1':np.zeros(I), 'FC2':np.zeros(I)}
        for i in range(I):
            soil_matrix['Ss'][i]=soil_data['Ss']
            soil_matrix['eta'][i]=soil_data['eta']
            soil_matrix['theta_r'][i]=soil_data['theta_r']
            soil_matrix['theta_s'][i]=soil_data['theta_s']
            soil_matrix['n'][i]=soil_data['n']
            soil_matrix['alpha'][i]=soil_data['alpha']
            soil_matrix['FC1'][i]=soil_data['FC1']
            soil_matrix['FC2'][i]=soil_data['FC2']
            if i!=0:
                soil_matrix['Ksat'][i]=soil_data['Ksat']*(1-exp(-pp*z[i]))
            else:
                soil_matrix['Ksat'][i]=soil_data['Ksat']*(1-exp(-pp*z[i+1]))
    elif numlay==2:
        soil_matrix={'Ksat':np.zeros(I),'Ss':np.zeros(I),'eta':np.zeros(I), 'theta_r':np.zeros(I), 'theta_s':np.zeros(I), 'n':np.zeros(I), 'alpha':np.zeros(I), 'FC1':np.zeros(I), 'FC2':np.zeros(I)}    
        LCC=(L-LC)/dz
        for i in range(I):
            if i<=int(LCC[-1]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-1]
                soil_matrix['Ss'][i]=soil_data['Ss'][-1]
                soil_matrix['eta'][i]=soil_data['eta'][-1]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-1]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-1]
                soil_matrix['n'][i]=soil_data['n'][-1]
                soil_matrix['alpha'][i]=soil_data['alpha'][-1]
                soil_matrix['FC1'][i]=soil_data['FC1'][-1]
                soil_matrix['FC2'][i]=soil_data['FC2'][-1]
            else:
                soil_matrix['Ksat'][i]=soil_data['Ksat'][0]
                soil_matrix['Ss'][i]=soil_data['Ss'][0]
                soil_matrix['eta'][i]=soil_data['eta'][0]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][0]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][0]
                soil_matrix['n'][i]=soil_data['n'][0]
                soil_matrix['alpha'][i]=soil_data['alpha'][0]
                soil_matrix['FC1'][i]=soil_data['FC1'][0]
                soil_matrix['FC2'][i]=soil_data['FC2'][0]
    elif numlay==3:
        soil_matrix={'Ksat':np.zeros(I),'Ss':np.zeros(I),'eta':np.zeros(I), 'theta_r':np.zeros(I), 'theta_s':np.zeros(I), 'n':np.zeros(I), 'alpha':np.zeros(I), 'FC1':np.zeros(I), 'FC2':np.zeros(I)}    
        LCC=(L-LC)/dz
        for i in range(I):
            if i<=int(LCC[-1]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-2]
                soil_matrix['Ss'][i]=soil_data['Ss'][-2]
                soil_matrix['eta'][i]=soil_data['eta'][-2]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-2]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-2]
                soil_matrix['n'][i]=soil_data['n'][-2]
                soil_matrix['alpha'][i]=soil_data['alpha'][-2]
                soil_matrix['FC1'][i]=soil_data['FC1'][-2]
                soil_matrix['FC2'][i]=soil_data['FC2'][-2]
            elif i>=int(LCC[-1]) and i<=int(LCC[-2]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-1]
                soil_matrix['Ss'][i]=soil_data['Ss'][-1]
                soil_matrix['eta'][i]=soil_data['eta'][-1]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-1]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-1]
                soil_matrix['n'][i]=soil_data['n'][-1]
                soil_matrix['alpha'][i]=soil_data['alpha'][-1]
                soil_matrix['FC1'][i]=soil_data['FC1'][-1]
                soil_matrix['FC2'][i]=soil_data['FC2'][-1]                
            else:
                soil_matrix['Ksat'][i]=soil_data['Ksat'][0]
                soil_matrix['Ss'][i]=soil_data['Ss'][0]
                soil_matrix['eta'][i]=soil_data['eta'][0]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][0]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][0]
                soil_matrix['n'][i]=soil_data['n'][0]
                soil_matrix['alpha'][i]=soil_data['alpha'][0]
                soil_matrix['FC1'][i]=soil_data['FC1'][0]
                soil_matrix['FC2'][i]=soil_data['FC2'][0]
    else:
        soil_matrix={'Ksat':np.zeros(I),'Ss':np.zeros(I),'eta':np.zeros(I), 'theta_r':np.zeros(I), 'theta_s':np.zeros(I), 'n':np.zeros(I), 'alpha':np.zeros(I), 'FC1':np.zeros(I), 'FC2':np.zeros(I)}    
        LCC=(L-LC)/dz
        for i in range(I):
            if i<=int(LCC[-1]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-1]
                soil_matrix['Ss'][i]=soil_data['Ss'][-1]
                soil_matrix['eta'][i]=soil_data['eta'][-1]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-1]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-1]
                soil_matrix['n'][i]=soil_data['n'][-1]
                soil_matrix['alpha'][i]=soil_data['alpha'][-1]
                soil_matrix['FC1'][i]=soil_data['FC1'][-1]
                soil_matrix['FC2'][i]=soil_data['FC2'][-1]
            elif i>int(LCC[-1])and i<=int(LCC[-2]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-2]
                soil_matrix['Ss'][i]=soil_data['Ss'][-2]
                soil_matrix['eta'][i]=soil_data['eta'][-2]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-2]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-2]
                soil_matrix['n'][i]=soil_data['n'][-2]
                soil_matrix['alpha'][i]=soil_data['alpha'][-2]
                soil_matrix['FC1'][i]=soil_data['FC1'][-2]
                soil_matrix['FC2'][i]=soil_data['FC2'][-2]                
            elif i>int(LCC[-2])and i<=int(LCC[-3]):
                soil_matrix['Ksat'][i]=soil_data['Ksat'][-3]
                soil_matrix['Ss'][i]=soil_data['Ss'][-3]
                soil_matrix['eta'][i]=soil_data['eta'][-3]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][-3]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][-3]
                soil_matrix['n'][i]=soil_data['n'][-3]
                soil_matrix['alpha'][i]=soil_data['alpha'][-3]
                soil_matrix['FC1'][i]=soil_data['FC1'][-3]
                soil_matrix['FC2'][i]=soil_data['FC2'][-3]                            
            else:
                soil_matrix['Ksat'][i]=soil_data['Ksat'][0]
                soil_matrix['Ss'][i]=soil_data['Ss'][0]
                soil_matrix['eta'][i]=soil_data['eta'][0]
                soil_matrix['theta_r'][i]=soil_data['theta_r'][0]
                soil_matrix['theta_s'][i]=soil_data['theta_s'][0]
                soil_matrix['n'][i]=soil_data['n'][0]
                soil_matrix['alpha'][i]=soil_data['alpha'][0]
                soil_matrix['FC1'][i]=soil_data['FC1'][0]
                soil_matrix['FC2'][i]=soil_data['FC2'][0]
    return(soil_matrix)


# build system of equation
def build_system(h0,h,theta0, bc, BC, q, soil_carac, ww, bcff):
    #h0 previous or initial profile of pressure head depending on the iteration
    #h actual profile of pressure head
    #theta0 previous ir initial profile of water content
    #bc boundary conditions
    #BC top boundary condition, its value depends on the time interval of the current iteration. 
            #Ex. Considering bc['top']['fixed_flow']=([11,4,1]) at a time bc['time']=([0,5,10]). If the current time of the iteration is 1.5
            #then the value of BC is 11. If the current time of the iteration is 9.9 then the value of BC is 4.
    #q is the source term, in our case is the potential transpiration estimated from the interception model as a flux [m/s]
    #soil_carac is the list of vectors of soil characteristics. Each vector has as the same size as the number on nodes.
    #ww is the index that represents the current time interval and it isused at the 'boundary conditions' and 'q' at the current time step.
            #Ex. Considering bc['top']['fixed_flow']=([11,4,1]) at a time bc['time']=([0,5,10]) and that the time step t[i] is variable.
            # At a time iteration t[10]=1.5, the corresponding ww is '0' so bc['top']['fixed_flow'][0]=11.
            # At a time iteration t[48]=5.892, the corresponding ww is '1' so bc['top']['fixed_flow'][1]=4.
            # At a time iteration t[73]=9.786, the corresponding ww is '1' so bc['top']['fixed_flow'][1]=4.
            # At a time iteration t[75]=10 the corresponding ww is '2' so bc['top']['fixed_flow'][1]=1.
    #bcff is a condition. When bcff='True', the 'top fixed flow' boundary condition is not longer considered in the iteration and it is changed by
        #a 'top fixed head' boundary condition, where the value of pressure head is 0 (it is used in the case of runoff and the current inflow is
        # is larger than the infiltration capacity of the soil) 

    # init vectors
    theta = get_theta(h,soil_carac)
    C     = get_C(h,soil_carac)
    K     = get_K(h,soil_carac)
    Kp    = get_Kp(K)
    WU    = w_uptake(theta0, RD, q[ww], soil_carac,dz)
    M     = np.zeros((h.size,h.size))
    B     = np.zeros((h.size,))
    # fill matrix M and forcing vector B
    for i in range(1,I-1):
       M[i,i-1] = get_m1(i,Kp)
       M[i,i]   = get_m2(i,Kp,C,theta)
       M[i,i+1] = get_m3(i,Kp)
       B[i] = - get_b1(i,C)*h[i] - get_b2(i,theta,soil_carac)*h0[i] \
        + get_b3(i,K) + get_b4(i,theta0,theta) + WU[i]
    # set boundary conditions
    # fixed head
    # top
    if bc['top'][0] == 'fixed_head':
       M[-1,-1] = 1. 
       M[-1,-2] = 0.  
       B[-1]   = bc['top'][1][ww]
    # bot
    if bc['bot'][0] == 'fixed_head':
       M[0,0] = 1. 
       M[0,1] = 0.  
       B[0]   = bc['bot'][1][ww]
    # fixed flow
    # top
    if bc['top'][0] == 'fixed_flow':
            if bcff=='True':
                M[-1,-1] = 1. 
                M[-1,-2] = 0.  
                B[-1]   = 0.
            else:
                #M[-1,-2] = get_m1(I,Kp)
                #M[-1,-1] = - get_m1(I, Kp) - get_b1(I-1,C) - get_b2(I-1,theta,soil_carac)
                #B[-1] = - get_b1(I-1,C)*h[-1] - get_b2(I-1,theta,soil_carac)*h0[-1] \
                #+ (K[I-1]-K[I-2])/dz + get_b4(I-1,theta0,theta) + WU[I-1] + get_m3(I-1,Kp)*dz*(bc['top'][1][ww]/Kp[-1] +1)
                M[-1,-2] = get_m1(-2,Kp)
                M[-1,-1] = -get_b1(-1,C)/2. - get_m1(-2,Kp) #- get_b2(-1,theta,soil_carac)
                B[-1] = -get_b1(-1,C)*h[-1]/2. + get_b4(-1,theta0,theta)/2. + Kp[-2]/dz + WU[-1]/2.+BC/dz   #this bc should be re-check eventhough it appers to work fine and it is in accordance to the bc used at Hydrus
                       
                
        #change get_m1(I-1,Kp)*dz*(bc['top'][1][n]/Kp[-1] +1) by get_m3(I-1,Kp)*dz*(bc['top'][1][n]/Kp[-1] +1)
        # bot
        #if bc['bot'][0] == 'fixed_flow':
        #M[0,1] = get_m1(0,Kp)
        #M[0,0] = get_m2(0,Kp,C,theta) + get_m3(0,Kp) 
        #B[0] = - get_b1(0,C)*h[0] - get_b2(I-1,theta,soil_carac)*h0[0] \
        #   + get_b3(I-2,K) + get_b4(I-1,theta0,theta) + q[0,n]*dt/dz + get_m1(I-1,Kp)*dz*(bc['top'][1][n]/Kp[-1] +1)
        #Change the values... get_b3 not sure what will result
    if bc['bot'][0] == 'fixed_flow':
        M[0,1] = get_m3(0,Kp)
        M[0,0] = -get_b1(0,C)/2.-get_m3(0,Kp)
        B[0] = - get_b1(0,C)*h[0]/2. + get_b4(0,theta0,theta)/2. -Kp[1]/dz +WU[0]/2.- (bc['bot'][1][ww])/dz 
        #Change the values... get_b3 not sure what will result

       # free drainage at the bottom
    if bc['bot'][0] == 'free_drainage':
        M[0,1] = get_m3(0,Kp)
        M[0,0] = -get_b1(0,C)/2.-get_m3(0,Kp)
        B[0] = - get_b1(0,C)*h[0]/2. + get_b4(0,theta0,theta)/2. -Kp[1]/dz +WU[0]/2.- (-K[0])/dz 

     
    # return M and B
    return(np.mat(M), np.transpose(np.mat(B)))

# -- Global dt variable
def gl_val(dts):
    global dt
    dt=dts


# model
def run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac, PICmax = 20, CRIT_CONV = 1.e-2, Runoff='True'):
    # L : column length
    # T : simulation duration
    # dz : vertical discretization
    # tstep : initial time step
    # h_init : initial profile of pressure head, must be of size I = int(round(L/dz))
    # bc : dictionary of boundary conditions : e.g. {'top':['fixed_head',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}
    #      allowed bc at the top : 'fixed_head','fixed_flow',
    #      allowed bc at the bottom : 'fixed_head','fixed_flow','free_drainage'
    # q : source term. in our case is the potential transpiration estimated from the interception model as a flux [m/s]
    #soil_carac is the list of vectors of soil characteristics. Each vector has as the same size as the number of nodes.
    # PICmax : maximum number of Picard iteration
    # CRIT_CONV : convergence criteria.
    # Runoff='True' --> Runoff will be calculated instead of a pounding condition
    # -- model initialization
    # init vars
    start = time()  #estimate time of calculation
    ttmin=100.  #minimum time step
    tmin=ttmin   #minimum time step (variable)
    tmax=900.  #maximum time step
    dtss=tmax
    if tstep>T: #conditions to change the initial time step if necessary
        tstep=T
    if tstep>tmax:
        tstep=900.
    if tstep<tmin:
        tstep=tmin

    t=np.array(0)   #create t, just initial value
    tt=t    #create tt, where time value shoul be changed at each time step of the prescribed boundary conditions.
    #t = np.hstack((0,np.cumsum(dts))) # time array
    N = int(round(T/tmin)) # number of time steps
    I = int(round(L/dz)+1) # number of cells
    z = np.linspace(0,L,I) # z coordinates of nodes
    # check h_init
    if len(h_init) != I:
        print('ERROR: check dimension of h_init')
        return(0,0,0)
    # -- check input data
    # check bc : 
    if bc['top'][0] not in ['fixed_head','fixed_flow','free_drainage'] :
        print('ERROR: check top boundary condition')
    if bc['bot'][0] not in ['fixed_head','fixed_flow','free_drainage'] :
        print('ERROR: check bottom boundary condition')
    #create variables as vectors     #It should be notice that several variables are commented because they are only used for the mass balance.
    RO = np.zeros((1))    #create vector RO, where runoff will be saved at each time step of the prescribed boundary conditions.
    INF = np.zeros((1))   #create vector INF, where infiltration will be saved at each time step of the prescribed boundary conditions.
    PER = np.zeros((1))   #create vector PER, where deep percolation will be saved at each time step of the prescribed boundary conditions.
    TA = np.zeros((1))    #create vector TA, where actual transpiration will be saved at each time step of the prescribed boundary conditions.
    VOL = np.zeros((1))   #create vector VOL, where soil water volume will be saved at each time step of the prescribed boundary conditions.
    

    # -- run initization
    # initialize output matrices and variables
    S=np.mat(np.zeros(shape=(I,1)))
    VOL[0] = (get_SUM(get_theta(h_init,soil_carac),dz,I))*dz
    Ta=0.
    Ro=0.
    Infil=0.
    Per=0.

     #include the time 0
    # Store initial condition
    S[:,0]= np.transpose(np.mat(h_init))
    h = h0 = h_init 
    h_list = []

    # create matrix of variables
    Theta = np.mat(np.zeros(shape=(I,1)))
    Theta[:,0]=np.transpose(np.mat(get_theta(h_init,soil_carac)))
    Ka = np.mat(np.zeros(shape=(I,1)))
    K0 = get_K(h,soil_carac)
    Ka[:,0]=np.transpose(np.mat(K0))
    Flow = np.mat(np.zeros(shape=(I,1)))

  

# iterate over time
    dts=tstep   #dts: the time step at the current iteration, should change depending of the conditions
    dt=dts      #dt: maybe no longer necessary, not sure
    t=np.zeros([1]) #create vector t, where time will be saved at each iteration.
    tt=t
    www=1
    iterlast=0
    alert=0.
    alert2=0.
    for n in range(1,N):    #range 
        t=np.hstack((t,(t[-1]+dts)))    #add new time step to the vector t
        for w in range(www,len(bc['time'])):  #condition to estimate ww, so the bc used are in accordance to its time interval
            if t[-1]<=bc['time'][w]:
                ww=w-1
                www=w
                break
            else:
                ww=w
        #if t[-1]==bc['time'][ww] or t[-1]==T:   #with this conditions we make sure that at the current time we are using bc and q of the previous time interval
            #www=ww-1
        gl_val(dts) #make dt=dts, maybe no longer necessary, not sure
        BC=bc['top'][1][ww] #choose the appropriate bc at the time t
        theta0 = get_theta(h0,soil_carac)   #get theta0 before the picard iterations
        # Picard iteration 
        h00=h   #save the initial value of h if runoff occurs so the top bc could be changed        
        for m in range(PICmax):
            Iter=m  #save the number of iterations, used for the time optimization
            M , B = build_system(h0, h, theta0 , bc, BC, q, soil_carac, ww, bcff='False')   #first with prescribed bc
            # solver linear system
            h1 = np.transpose(np.matrix(cgs(csr_matrix(M),B)[0]))
            #h1 = inv(M)*B
            if Iter==(PICmax-1):
                alert=1.
            if np.max(abs(h1-h)) < CRIT_CONV:
                break
            h=h1
        if abs((h00[-1]-h1[-1])/min(abs(h1[-1]),abs(h00[-1])))>0.3:
            alert=1.
        if alert==1.:
            print('Interval:' + str(tt[-1]) + '     timestep:' + str(t[-1]) + 'dt:' + str(dts) + '    BC:' + str(BC) + '      htop:' + str(h00[-1]))
            dts=1.
            gl_val(dts)
            #tmin=1.
            t[-1]=t[-2]+dts
            www=www-1
            for w in range(www,len(bc['time'])):  #condition to estimate ww, so the bc used are in accordance to its time interval
                if t[-1]<=bc['time'][w]:
                    ww=w-1
                    www=w
                    break
                else:
                    ww=w
            h=h00
             for mmm in range((PICmax)):
                Iter=mmm
                M , B = build_system(h0, h, theta0 , bc, BC, q, soil_carac, ww, bcff='False')   #first with prescribed bc
                # solver linear system
                h1 = np.transpose(np.matrix(cgs(csr_matrix(M),B)[0]))
                #h1 = inv(M)*B
                if Iter==(PICmax-1):
                    alert2=1.    
                if np.max(abs(h1-h)) < CRIT_CONV:
                    break
                h=h1     
            if alert2==1.:
                print('Interval:' + str(tt[-1]) + '     timestep:' + str(t[-1]) + 'dt:' + str(dts) + '     BC:' + str(BC) + '      htop:' + str(h00[-1])+ 'Per:' + str(Per))
                print('valio paloma')
            iterlast=0.          
            h=h1
            alert=0.
            alert2=0.
        if Runoff=='True':  #If runoff should be estimated
            if h1[-1]> 0 :  #if the top node is equal or higher than 0, saturation conditions should be expected and runoff is activated
                #Include the Runoff Modification: these are optional conditions to maintain equilibrium on the solution, not always effective                      
                #dts=tmin
                #dt=dts                        
                h=h00   #in such case the last estimation of runoff is forgot
                print('Runoff modification')
                for mm in range (PICmax):   #and a new estimation starts
                    Iter=mm
                    M , B = build_system(h0, h, theta0 , bc, BC, q, soil_carac, ww, bcff='True')    #the changed bc is used
                    # solver linear system
                    h1 = np.transpose(np.matrix(cgs(csr_matrix(M),B)[0])) # see also cgs, csr, spsolve
                    if np.max(abs(h1-h)) < CRIT_CONV:   
                        #print( 'Runoff Modification PIC iteration ='+str(mm))
                        break
                    h=h1
        h=h1   
            
        #variables estimation    
        theta = get_theta(h,soil_carac)
        C     = get_C(h,soil_carac)
        K     = get_K(h,soil_carac)
        Kp    = get_Kp(K)
        WU    = w_uptake(theta0, RD, q[ww], soil_carac, dz)
        Ta    = get_SUM(WU,dz,I)*dts*dz+Ta
        Flux  = get_flux(Kp,h,dz,theta,theta0,dts,WU)
        Per = Flux[0]*dts+Per     #indirect estimation        
        # Per = ((-(soil_carac['Ksat'][1]*Kr_curve(h[1],soil_carac,1)+soil_carac['Ksat'][0]*Kr_curve(h[0],soil_carac,1))/2.)*((h[1]-h[0])/dz+1.)-dz/2.*((theta_curve(h[0],soil_carac,0)-theta0[0])/dts+WU[0]))*dts+Per    #direct estimation
        if h[-1]>=0.:
            Infil = Flux[-1]*dts+Infil
            Ro = Ro+(bc['top'][1][ww]-Flux[-1])*dts
        else:
            Infil = bc['top'][1][ww]*dts+Infil
            Ro=Ro+0

        if t[-1]==bc['time'][www] or t[-1]==T:   #check if time is at a time step of the bc, if yes, all variables estimation are saved for this time step 
            S = np.c_[S,np.zeros(I)]
            S[:,-1] = np.mat(h)
            tt = np.hstack((tt,bc['time'][www]))
            Theta = np.c_[Theta,np.zeros(I)]
            Theta[:,-1] = np.transpose(np.mat(theta))
            Ka = np.c_[Ka,np.zeros(I)]
            Ka[:,-1] = np.transpose(np.mat(K)) 
            Flow = np.c_[Flow,np.zeros(I)]
            Flow[:,-1] = np.transpose(np.mat(Flux))
            TA = np.hstack((TA, Ta))
            RO = np.hstack((RO, Ro))
            INF = np.hstack((INF,Infil))
            PER = np.hstack((PER, Per))
            VOL = np.hstack((VOL, (get_SUM(theta,dz,I))*dz)) 
            #initial conditions to the next iteration
            Ro=0.
            Infil=0.
            Per=0.
            Ta=0. 


       ##Time optimization: 
            #it follows the methodology used in Hydrus
            #change the time step dts from tmin to bc timestep
            #time can be between bc prescribed intervals but it will always have to coincide with the intervals
            #some limitations were made in order to achieve equilibrium
        tmin=min(dts,ttmin)
        if t[-1]==T:    #first check if time is the total time, in such case the calculation is over
            tt[-1] = T
            break
        h0 = h  #new h0 for next iteration
        if t[-1]==bc['time'][www]:   #check if time is equal to a bc time step
            if (BC-bc['top'][1][www])>0.:    #condition for large flux or large changes in bc
                if BC==0. and abs(bc['top'][1][www])>=3./(1000*15*60):
                    dts=min(dts,100.)
                    tmin=dts
                    iterlast=0.
                elif BC!=0. and abs((BC-bc['top'][1][www])/BC)>=10.:
                    dts=min(dts,100.)
                    tmin=dts
                    iterlast=0.
                else:                
                    if iterlast==1.:
                        dts=dtss
                        iterlast=0.
                        if dts>(tmax-tmin):
                            dts=tmax
                    else:
                        dts=dts 
            else:
                if iterlast==1.:
                    dts=dtss
                    iterlast=0.
                    if dts>(tmax-tmin):
                        dts=tmax
                else:
                    dts=dts        
        else:
            if t[-1]>bc['time'][(len(bc['time'])-1)]:   #this is the case of the bc time step before the total time T
                if Iter<=3:
                    dts=max(1.3*dts,tmin)
                    if (T-t[-1])<=(dts+tmin):
                        dtss=min(tmax,dts)
                        iterlast=1.
                elif Iter>3 and Iter<=7:
                    dts=dts
                    if (T-t[-1])<=(dts+tmin):
                        dts=max(tmin,(T-t[-1]))
                elif Iter>7 and Iter<PICmax:
                    dts=max(0.7*dts,tmin)
                    if (T-t[-1])<=(dts+tmin):
                        dts=max(tmin,(T-t[-1]))
                else:
                    dts=max(dts/3,tmin)
                    if (T-t[-1])<=(dts+tmin):
                        dts=max(tmin,(T-t[-1]))
            else:        #this is the case for any other time step
                if Iter<=3:    
                    dts=max(1.3*dts,tmin)

                    if (bc['time'][www]-t[-1])<=(dts+tmin):
                        dtss=min(tmax,dts)

                        if tmin>=ttmin:
                            dts=max(tmin,(bc['time'][www]-t[-1]))
                        else:
                            dts=(bc['time'][www]-t[-1])
                        iterlast=1.

                        if (dts+tmin)>(bc['time'][www]-t[-1]):
                            dts= (bc['time'][www]-t[-1])

                elif Iter>3 and Iter<=7:
                    dts=dts
                    if (bc['time'][www]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][www]-t[-1]))
                        if (dts+tmin)>(bc['time'][www]-t[-1]):
                            dts= (bc['time'][www]-t[-1])

                elif Iter>7 and Iter<PICmax:
                    dts=max(0.7*dts,tmin)
                    if (bc['time'][www]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][www]-t[-1]))
                        if (dts+tmin)>(bc['time'][www]-t[-1]):
                            dts= (bc['time'][www]-t[-1])
                else:
                    dts=max(dts/3,tmin)
                    if (bc['time'][www]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][www]-t[-1]))
                        if (dts+tmin)>(bc['time'][www]-t[-1]):
                            dts= (bc['time'][www]-t[-1])
        #if t[-1]>=25269300 and t[-1]<=25271100:
        #    dts=1.
        #    tmin=1.
        #    iterlast=0   
        #if t[-1]>=25784100 and t[-1]<=25785900:
        #    dts=1.
        #    tmin=1.
        #    iterlast=0    
        if float(www)%96. ==0 and t[-1]==bc['time'][www]:  #uncomment if print daily values is wanted
            print('time:'+str(www/96)+'    iteration ' + str(n) + ' terminated.')
        #print('iteration ' + str(n) + ' terminated.')
    # return simulation results
    return(t,tt,z,S,Theta,Ka,Flow,INF,RO,PER,TA,VOL)
##end of varsat

##Mass Balance


# def mass_balance(TimeLevel,t,INF,PER,TA,RO,VOL):
#     if TimeLevel[0] not in t[:]:
#         print('ERROR: Time levels are not correct. Choose the prescribed values in vector "t" ')
#     if TimeLevel[-1] not in t[:]:
#         print('ERROR: Time levels are not correct. Choose the prescribed values in vector "t" ')
    
#     for i in range (0,t.size):
#         if TimeLevel[0]==t[i]:
#             a=i
#             break
#     for i in range (0,t.size): 
#         if TimeLevel[-1]==t[i]:
#             b=i
#             break
#     INFcom=sum(INF[a+1:b+1])
#     PERcom=sum(PER[a+1:b+1])
#     TAcom=sum(TA[a+1:b+1])
#     ROcom=sum(RO[a+1:b+1])
#     Erabs=VOL[b]-VOL[a]-(-INFcom+PERcom+ROcom-TAcom)
#     Errel=abs(Erabs)*100./abs(VOL[b]-VOL[a])
#     if Errel>0.1:
#         warning='Attention. Error too high.'
#     else:
#         warning=''
#     print('---------Mass Balance-----------')
#     print('')
#     print('')
#     print('Interval t1:' + str(t[a]) + '  &  t2:' +str(t[b]))
#     print('')
#     print('Initial volume:  ' + str(VOL[a]) + '  [L]')
#     print('Final volume:    ' + str(VOL[b]) + '  [L]')
#     print('______')
#     print('')
#     print('Total infiltration:          ' + str(-INFcom) + '  [L]')
#     print('Total percolation:           ' + str(-PERcom) + '  [L]')
#     print('Total actual transpiration:  ' + str(TAcom) + '  [L]')
#     print('Total runoff:                ' + str(-ROcom) + '  [L]')
#     print('______')
#     print('')
#     print('Total input:     ' + str(-INFcom) + '  [L]')
#     print('Total output:    ' + str(-(PERcom+ROcom-TAcom)) + '  [L]')
#     print('______')
#     print('')
#     print('Input')
#     print('Relative Rainfall: ' + str(INFcom*100./INFcom) + '  [%]')
#     print('')
#     print('Output')
#     print('Relative percolation:            ' + str(PERcom*100./(PERcom+ROcom-TAcom)) + '  [%]')
#     print('Relative actual transpiration:   ' + str(-TAcom*100./(PERcom+ROcom-TAcom)) + '  [%]')
#     print('Relative runoff:                 ' + str(ROcom*100./(PERcom+ROcom-TAcom)) + '  [%]')
#     print('______')
#     print('')
#     print('Absolute error: ' + str(Erabs) + '  [L]')
#     print('Relative error: ' + str(Errel) + '  [%]' + '   ' + str(warning))
