# -*- coding: utf-8 -*-
"""
Created on Sun May 18 12:34:17 2014

@author: redhotkenny
"""

# --------------------------------------------------------------
# -- Load modules
import numpy as np
from matplotlib import pyplot as plt
from time import *
from math import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
from scipy.sparse.linalg import cgs
from scipy.sparse import csr_matrix


# --------------------------------------------------------------
# -------------------- Declare functions ------------------------
# --------------------------------------------------------------


# -- Van Genuchten functions

def theta_curve(h,soil_carac):
    h = float(h)
    n = soil_carac['n']
    theta_r = soil_carac['theta_r']
    theta_s = soil_carac['theta_s']
    alpha = soil_carac['alpha']
    m = 1 - 1/n
    if h < 0 : 
        return( theta_r + (theta_s - theta_r) * (1 + abs(alpha*h)**n)**(-m) )
    else :
        return( theta_s )


def C_curve(h,soil_carac):
    h = float(h)
    n = soil_carac['n']
    theta_r = soil_carac['theta_r']
    theta_s = soil_carac['theta_s']
    alpha = soil_carac['alpha']
    m = 1 - 1./n
    if h < 0 : 
        return(  m * n * alpha**n * abs(h) **(n-1) * (theta_s-theta_r) * (1 + abs(alpha*h)**n)**(-m-1)   )
    else :
        return(0)


def Kr_curve(h,soil_carac):
    h = float(h)
    n = soil_carac['n']
    theta_r = soil_carac['theta_r']
    theta_s = soil_carac['theta_s']
    alpha = soil_carac['alpha']
    m = 1 - 1./n
    if h < 0 : 
        return( ( 1 - abs(alpha*h)**(n-1)*(1+abs(alpha*h)**n)**(-m) )**2 / 
        (( 1 + abs(alpha*h)**n  )**(0.5*m)) )
    else : 
        return(1)



## ------Check and visualize soil properties --------------------

##Check and visualize soil properties
#h_vals = np.arange( -5, 1, 0.05 )
#h_vals = np.arange( -1, 0, 0.05 )
#
#theta_vals = [ theta_curve(h_val,soil_carac) for h_val in h_vals ]
#C_vals = [ C_curve(h_val,soil_carac) for h_val in h_vals ]
#Kr_vals = [ Kr_curve(h_val,soil_carac) for h_val in h_vals ]
#
## plot theta-psi and C-psi plot
#plt.ion()
#plt.plot(h_vals,theta_vals,label='theta(h)')
##plt.plot(h_vals,C_vals,label='C(h)')
#plt.xlim(-1.5,0.5)
#plt.show()
#
#
##plot theta_Kr
#plt.plot(h_vals,Kr_vals,label='Kr(h)')
#plt.xlim(-1.5,0.5)
#plt.legend()
#
#
# --------------------------------------------------------------
# --------------------------------------------------------------

# Vector K, hydraulic conductivity
def get_K(h,soil_carac):
    Ksat  = soil_carac['Ksat']
    K = np.zeros((h.size,))
    for i in range(h.size):
        K[i] =  Ksat*Kr_curve(h[i],soil_carac)
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
        theta[i] =  theta_curve(h[i],soil_carac)
    return(theta)

# vector C
def get_C(h, soil_carac):
    C = np.zeros((h.size,))
    for i in range(h.size):
        C[i] =  C_curve(h[i],soil_carac)
    return(C)

#Root Distribution
def get_RD(L, I, LR, z, dz, type):    #type can be homogeneous or follow Hoffman and van Genuchten, 1983 distribution (Trapecio) 
    RD=np.zeros(I)
    if type =='homogeneous':
        for i in range(I):
            if z[i] <= (L-LR):
                RD[i]= 0
            else:
                RD[i]=1
    elif type == 'distribution':
        for i in range(I):
            if z[i] > (L-0.2*LR):
                RD[i] = 1.667/LR
            elif z[i] < (L-LR):
                RD[i]=0
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
    Ss = soil_carac['Ss']
    eta = soil_carac['eta']
    return( (Ss*theta[i]) / (eta*dt) )

# coefficient b3 (g in Clement et al., 1994)
# coefficient b3 (g in Clement et al., 1994)
def get_b3(i,Kp):
    #if i == 0 : 
    return(-(Kp[i+1]-Kp[i])/dz)
    #return( 0 )
        #return( -( K[i+1]-K[i] ) / (2.*dz) )
    #return( -( K[i+2]-K[i] ) / (2.*dz) )
    #else :     
    #return( -( K[i+1]-K[i-1] ) / (2.*dz) )


# coefficient b4 (h in Clement et al., 1994)
def get_b4(i,theta0,theta):
    return( (theta[i] - theta0[i] )/dt )

#build Potential Water Uptake Distribution
def w_uptake(theta0, RD, q, soil_carac):  #from Feddes et al., 1978 
    WU= np.zeros(I)
    for i in range(I):
        if theta0[i]>=soil_carac['FC2']:
            red=(soil_carac['theta_s']-theta0[i])/(soil_carac['theta_s']-soil_carac['FC2'])
        elif theta0[i]>=soil_carac['FC2']:
            red=1
        elif theta0[i]>=soil_carac['theta_r']:
            red=(-soil_carac['theta_r']+theta0[i])/(-soil_carac['theta_r']+soil_carac['FC2'])
        elif theta0[i]<soil_carac['theta_r']:
            red=0
        else:
            WU[i]=red*RD[i]*q[ww]*I/(I-1)
    return(WU)


# build system of equation
def build_system(h0,h,theta0, bc, q, soil_carac, ww, bcff):
    # init vectors
    theta = get_theta(h,soil_carac)
    C     = get_C(h,soil_carac)
    K     = get_K(h,soil_carac)
    Kp    = get_Kp(K)
    WU    = w_uptake(theta0, RD, q, soil_carac)
    M     = np.zeros((h.size,h.size))
    B     = np.zeros((h.size,))
    # fill matrix M and forcing vector B
    for i in range(1,I-1):
       M[i,i-1] = get_m1(i,Kp)
       M[i,i]   = get_m2(i,Kp,C,theta)
       M[i,i+1] = get_m3(i,Kp)
       B[i] = - get_b1(i,C)*h[i] - get_b2(i,theta,soil_carac)*h0[i] \
        + get_b3(i,Kp) + get_b4(i,theta0,theta) + WU[i]
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
                M[-1,-2] = get_m1(I-1,Kp)
                M[-1,-1] = get_m2(I-1,Kp,C,theta) + get_m3(I-1,Kp)
                B[-1] = - get_b1(I-1,C)*h[-1] - get_b2(I-1,theta,soil_carac)*h0[-1] \
                + get_b3(I-2,Kp) + get_b4(I-1,theta0,theta) + WU[I-1] + get_m3(I-1,Kp)*dz*(bc['top'][1][ww]/Kp[-1] +1)
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
        M[0,0] = get_m2(0,Kp,C,theta) + get_m1(0,Kp)
        B[0] = - get_b1(0,C)*h[0] - get_b2(0,theta,soil_carac)*h0[0] \
        + get_b3(0,Kp) + get_b4(0,theta0,theta) + WU[0] - get_m1(0,Kp)*dz*(bc['bot'][1][ww]/Kp[0] +1)
        #Change the values... get_b3 not sure what will result

    # free drainage at the bottom
    #if bc['bot'][0] == 'free_drainage':
    #M[0,1] = get_m1(0,Kp)
    #M[0,0] = get_m2(0,Kp,C,theta) + get_m3(0,Kp)
    #B[0] = - get_b1(0,C)*h[0] - get_b2(I-1,theta,soil_carac)*h0[0] \
        #+ get_b3(I-2,K) + get_b4(I-1,theta0,theta) + q[0,n]*dt/dz 

   # free drainage at the bottom
    if bc['bot'][0] == 'free_drainage':
        M[0,1] = get_m3(0,Kp)
        M[0,0] = get_m2(0,Kp,C,theta) + get_m1(0,Kp)
        B[0] = - get_b1(0,C)*h[0] - get_b2(0,theta,soil_carac)*h0[0] \
        + get_b3(0,Kp) + get_b4(0,theta0,theta) + WU[0]
     
    # return M and B
    return(np.mat(M), np.transpose(np.mat(B)), theta, C, K, Kp, WU )

# -- Global dt variable
def gl_val(dts):
    global dt
    dt=dts

# model
def run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac, PICmax = 10, CRIT_CONV = 1e-2, Runoff='True'):
    # L : column length
    # T : simulation duration
    # dz : vertical discretization
    # dts : time step, which may be constant or variable. In the latter case, we must have : np.cumsum(dts) = T
    # h_init : initial profile of pressure head, must be of size I = int(round(L/dz))
    # bc : dictionary of boundary conditions : e.g. {'top':['fixed_head',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}
    #      allowed bc at the top : 'fixed_head','fixed_flow',
    #      allowed bc at the bottom : 'fixed_head','fixed_flow','free_drainage'
    # q : source term. May be constant homogeneous (scalar), constant (size I array), variable in space and time (size I*N array)
    # soil_carac : e.g. soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35}
    # PICmax : maximum number of Picard iteration
    # CRIT_CONV : convergence criteria.
    # Runoff='True' --> Runoff will be calculated instead of a pounding condition
    # -- model initialization
    # init vars
    tmin=1.
    tmax=3600.
    if tstep>T:
        tstep=T
    if tstep>tmax:
        tstep=3600.
    if tstep<tmin:
        tstep=1.

    #if hasattr(dts,"__len__") == False:
        #dts  = np.array( [dts] * int(round(T/dts)) )
    t=np.array(0)
    #t = np.hstack((0,np.cumsum(dts))) # time array
    N = int(round(T/tmin)) # number of time steps
    I = int(round(L/dz)+1) # number of cells
    z = np.linspace(0,L,I) # z coordinates of nodes
    # check h_init
    if len(h_init) != I:
        print('ERROR: check dimension of h_init')
        return(0,0,0)
    #  check q (not fully implemented)
    #if hasattr(q,"__len__") == False: #  constant and homogeneous q
    #    q  = np.array( [ [q] * N] * I  )
    #elif np.array(q).shape == (I,) :  # transient homogeneous q
    #    q = np.transpose(np.array(  [list(q)]*N ) ) 
    # -- check input data
    # check bc : 
    if bc['top'][0] not in ['fixed_head','fixed_flow','free_drainage'] :
        print('ERROR: check top boundary condition')
    if bc['bot'][0] not in ['fixed_head','fixed_flow','free_drainage'] :
        print('ERROR: check bottom boundary condition')
    #initialize Runoff and Infiltration
    RO = np.zeros((1))
    INF = np.zeros((1))
    
    # -- run initization
    # initialize output matrices
    S=np.mat(np.zeros(shape=(I,1))) 
    #include the time 0
    # Store initial condition
    S[:,0]= np.transpose(np.mat(h_init))
    h = h0 = h_init 
    h_list = []
    # iterate over time
    dts=tstep
    dt=dts
    t=np.zeros([1])
    for n in range(1,N):
        t=np.hstack((t,(t[-1]+dts)))
        for w in range(1,len(bc['time'])):
            if t[-1]<bc['time'][w]:
                ww=w-1
                break
            else:
                ww=w
        gl_val(dts)
        theta0 = get_theta(h0,soil_carac)
        RO=np.hstack((RO,0.))
        INF=np.hstack((INF,bc['top'][1][ww]*dt)) #we must choose if we want to calculate the infiltration even with fixed head, now it only calculates to fixed flow
        # Picard iteration 
        h00=h
        for m in range(PICmax):
            Iter=m
            M , B, = build_system(h0, h, theta0 , bc, q, soil_carac, ww, bcff='False')
            # solver linear system
            h1 = np.transpose(np.matrix(cgs(csr_matrix(M),B)[0]))
            if np.max(h1-h) < CRIT_CONV:
                print('PIC iteration = '+ str(m))
                break
            h=h1
        if Runoff=='True':
            if h1[-1]> 0 :
                #Include the Runoff Modification                        
                h=h00
                for mm in range (PICmax):
                    Iter=mm
                    M , B = build_system(h0, h, theta0 , bc, q, soil_carac, ww, bcff='True')
                    # solver linear system
                    h1 = np.transpose(np.matrix(cg(csr_matrix(M),B)[0])) # see also cgs, csr, spsolve
                    if np.max(h1-h) < CRIT_CONV:   
                        print( 'Runoff Modification PIC iteration ='+str(mm))
                        break
                    h=h1
                h=h1
                Infil = Kp[-1]*(((get_m1(I-1,Kp)*h1[-2]+get_b2(I-1,theta,soil_carac)*h0[-1] \
                    -get_b3(I-2,Kp)-get_b4(I-1,theta0,theta))/(get_m3(I-1,Kp)*dz))-1)
                RO=np.hstack((RO, RO[-1] = (bc['top'][1][ww]-Infil[0,0])*dt))
                INF=np.hstack((INF,Infil[0,0]*dt))      
        h=h1
        h0 = h
        if t[-1]==bc['time'][ww]:
            dts=dts #key feature, if flux input vary variable it would be better to use dts=tmin
        else:
            if t[-1]>bc['time'][(len(bc['time'])-1)]:
                if Iter<=3:
                    dts=max(1.3*dts,tmin)
                    if (T-t[-1])<=(dts+tmin):
                        dts=max(tmin,(T-t[-1]))
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
            else:
                if Iter<=3:
                    dts=max(1.3*dts,tmin)
                    if (bc['time'][ww+1]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][ww+1]-t[-1]))
                elif Iter>3 and Iter<=7:
                    dts=dts
                    if (bc['time'][ww+1]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][ww+1]-t[-1]))
                elif Iter>7 and Iter<PICmax:
                    dts=max(0.7*dts,tmin)
                    if (bc['time'][ww+1]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][ww+1]-t[-1]))
                else:
                    dts=max(dts/3,tmin)
                    if (bc['time'][ww+1]-t[-1])<=(dts+tmin):
                        dts=max(tmin,(bc['time'][ww+1]-t[-1]))
        print('iteration ' + str(n) + ' terminated.')
        if t[-1]==T:
            break
    # return simulation results

    return(t,z,S,INF,RO)
