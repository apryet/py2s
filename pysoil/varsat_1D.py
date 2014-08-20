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
            if z[i] < (L-LR):
                RD[i]= 0
            else:
                RD[i]=1*I/(I-1)
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
    WU= np.zeros(I)
    for i in range(I):
        if theta0[i]>=soil_carac['FC2']:
            red=(soil_carac['theta_s']-theta0[i])/(soil_carac['theta_s']-soil_carac['FC2'])
        elif theta0[i]>=soil_carac['FC1']:
            red=1
        elif theta0[i]>=soil_carac['theta_r']:
            red=(-soil_carac['theta_r']+theta0[i])/(-soil_carac['theta_r']+soil_carac['FC1'])
        elif theta0[i]<soil_carac['theta_r']:
            red=0
        WU[i]=red*RD[i]*q*I/((I-1)*dz)
    return(WU)

#build Flux             #Follow Darcy Law and mass conservancy equation from Celia, as used in Hydrus p 112 del manual
def get_flux(Kp,h,dz,theta,theta0,dts,WU):
    Fl = np.zeros(h.size)
    #Fl[0] = -Kp[0]*((h[1]-h[0])/dz+1)
    Fl[0] = -Kp[1]*((h[1]-h[0])/dz+1)-dz/2*((theta[0]-theta0[0])/dts+WU[0]) 
    Fl[-1] = -Kp[-2]*((h[-1]-h[-2])/dz+1)-dz/2*((theta[-1]-theta0[-1])/dts+WU[-1])
    for i in range(1,(h.size-1)): 
        Fl[i]=(-Kp[i+1]*((h[i+1]-h[i])/dz+1)-Kp[i]*((h[i]-h[i-1])/dz+1))/2
    return (Fl)

#build SUM functions            #Arithmetic mean
def get_SUM(x,dz,I):
    X=np.zeros(x.size-1)
    for i in range(0,I-1):
        X[i]=(x[i+1]+x[i])/2
    return(X.sum())


# build system of equation
def build_system(h0,h,theta0, bc, q, soil_carac, ww, bcff):
    # bcff if TRUE, 
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
            if bcff=='True': # fixed pressure head 
                M[-1,-1] = 1. 
                M[-1,-2] = 0.  
                B[-1]   = 0.
            else:
                #M[-1,-2] = get_m1(I,Kp)
                #M[-1,-1] = - get_m1(I, Kp) - get_b1(I-1,C) - get_b2(I-1,theta,soil_carac)
                #B[-1] = - get_b1(I-1,C)*h[-1] - get_b2(I-1,theta,soil_carac)*h0[-1] \
                #+ (K[I-1]-K[I-2])/dz + get_b4(I-1,theta0,theta) + WU[I-1] + get_m3(I-1,Kp)*dz*(bc['top'][1][ww]/Kp[-1] +1)
                M[-1,-2] = get_m1(-2,Kp)
                M[-1,-1] = -get_b1(-1,C)/2. - get_m1(-2,Kp)
                B[-1] = -get_b1(-1,C)*h[-1]/2 + get_b4(-1,theta0,theta)/2. + Kp[-2]/dz + WU[-1]/2.+bc['top'][1][ww]/dz
                       
                
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

   # free drainage at the bottom
    #if bc['bot'][0] == 'free_drainage':
    #    M[0,1] = get_m3(0,Kp)
    #    M[0,0] = get_m2(0,Kp,C,theta) + get_m1(0,Kp)
    #    B[0] = - get_b1(0,C)*h[0] - get_b2(0,theta,soil_carac)*h0[0] \
    #    + (K[0]-K[1])/dz + get_b4(0,theta0,theta) + WU[0]
     
    # return M and B
    return(np.mat(M), np.transpose(np.mat(B)))

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
    tt=t
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
    PER = np.zeros((1))
    TA = np.zeros((1))
    VOL = np.zeros((1))
    # -- run initization
    # initialize output matrices
    S=np.mat(np.zeros(shape=(I,1))) 
    #include the time 0
    # Store initial condition
    S[:,0]= np.transpose(np.mat(h_init))
    h = h0 = h_init 
    h_list = []
    Theta = np.mat(np.zeros(shape=(I,1)))
    theta0 = get_theta(h0,soil_carac)
    Theta[:,0]=np.transpose(np.mat(theta0))
    Ka = np.mat(np.zeros(shape=(I,1))) # util para comparar con Hydrus
    K0 = get_K(h,soil_carac)
    Ka[:,0]=np.transpose(np.mat(K0))
    Flow = np.mat(np.zeros(shape=(I,1)))
    VOL[0] = (get_SUM(theta0,dz,I))*dz
    Ta=0.
    Ro=0.
    Infil=0.
    Per=0.

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
        #RO=np.hstack((RO,0.))
        #INF=np.hstack((INF,bc['top'][1][ww]*dt)) #we must choose if we want to calculate the infiltration even with fixed head, now it only calculates to fixed flow
        # Picard iteration 
        h00=h
        for m in range(PICmax):
            Iter=m
            M , B = build_system(h0, h, theta0 , bc, q, soil_carac, ww, bcff='False')
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
            
        #Infil = Kp[-1]*(((get_m1(I-1,Kp)*h1[-2]+get_b2(I-1,theta,soil_carac)*h0[-1] \
        #    -get_b3(I-2,Kp)-get_b4(I-1,theta0,theta))/(get_m3(I-1,Kp)*dz))-1)
        #Ro=(bc['top'][1][ww]-Infil[0,0])*dt
        #RO=np.hstack((RO, (bc['top'][1][ww]-Infil[0,0])*dt))
        #INF=np.hstack((INF,Infil[0,0]*dt))      
        
        theta = get_theta(h,soil_carac)
        C     = get_C(h,soil_carac)
        K     = get_K(h,soil_carac)
        Kp    = get_Kp(K)
        WU    = w_uptake(theta0, RD, q[ww], soil_carac, dz) # water update
        Ta    = get_SUM(WU,dz,I)*dts*dz+Ta # cumulative transpiration
        Flux  = get_flux(Kp,h,dz,theta,theta0,dts,WU)  # flux in each node 
        Per = Flux[0]*dts+Per # cumulative percolation
        if h[-1]>=0.:
            Infil = Flux[-1]*dts+Infil  # infiltration at the top of the soil column
            Ro = Ro+(bc['top'][1][ww]-Flux[-1])*dts # Runoff 
        else:
            Infil = bc['top'][1][ww]*dts+Infil
            Ro=Ro+0

        if t[-1]==bc['time'][ww] or t[-1]==T:
            Theta = np.c_[Theta,np.zeros(I)]
            Theta[:,-1] = np.transpose(np.mat(theta))
            Ka = np.c_[Ka,np.zeros(I)]
            Ka[:,-1] = np.transpose(np.mat(K))
            S = np.c_[S,np.zeros(I)]
            S[:,-1] = np.mat(h)
            Flow = np.c_[Flow,np.zeros(I)] # matrix, same as S
            Flow[:,-1] = np.transpose(np.mat(Flux))
            TA = np.hstack((TA, Ta)) # useful for mass balance 
            RO = np.hstack((RO, Ro)) # idem
            INF = np.hstack((INF,Infil)) # idem
            PER = np.hstack((PER, Per)) # idem
            VOL = np.hstack((VOL, (get_SUM(theta,dz,I))*dz)) # Volume of water in each cell for each 
            tt = np.hstack((tt,bc['time'][ww]))
            Ro=0.
            Infil=0.
            Per=0.
            Ta=0. 


        if t[-1]==T: #  if last time step is reached
            tt[-1] = 
            break
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
    # return simulation results

    return(t,tt,z,S,Theta,Ka,Flow,INF,RO,PER,TA,VOL)

##Mass Balance


def mass_balance(TimeLevel,t,INF,PER,TA,RO,VOL):
    if TimeLevel[0] not in t[:]:
        print('ERROR: Time levels are not correct. Choose the prescribed values in vector "t" ')
    if TimeLevel[-1] not in t[:]:
        print('ERROR: Time levels are not correct. Choose the prescribed values in vector "t" ')
    
    for i in range (0,t.size):
        if TimeLevel[0]==t[i]:
            a=i
            break
    for i in range (0,t.size): 
        if TimeLevel[-1]==t[i]:
            b=i
            break
    INFcom=sum(INF[a+1:b+1])
    PERcom=sum(PER[a+1:b+1])
    TAcom=sum(TA[a+1:b+1])
    ROcom=sum(RO[a+1:b+1])
    Erabs=VOL[b]-VOL[a]-(-INFcom+PERcom+ROcom-TAcom)
    Errel=abs(Erabs)*100/abs(VOL[b]-VOL[a])
    if Errel>0.1:
        warning='Attention. Error too high.'
    else:
        warning=''
    print('---------Mass Balance-----------')
    print('')
    print('')
    print('Interval t1:' + str(t[a]) + '  &  t2:' +str(t[b]))
    print('')
    print('Initial volume:  ' + str(VOL[a]) + '  [L]')
    print('Final volume:    ' + str(VOL[b]) + '  [L]')
    print('______')
    print('')
    print('Total infiltration:          ' + str(-INFcom) + '  [L]')
    print('Total percolation:           ' + str(-PERcom) + '  [L]')
    print('Total actual transpiration:  ' + str(TAcom) + '  [L]')
    print('Total runoff:                ' + str(-ROcom) + '  [L]')
    print('______')
    print('')
    print('Total input:     ' + str(-INFcom) + '  [L]')
    print('Total output:    ' + str(-(PERcom+ROcom-TAcom)) + '  [L]')
    print('______')
    print('')
    print('Input')
    print('Relative Rainfall: ' + str(INFcom*100/INFcom) + '  [%]')
    print('')
    print('Output')
    print('Relative percolation:            ' + str(PERcom*100/(PERcom+ROcom-TAcom)) + '  [%]')
    print('Relative actual transpiration:   ' + str(-TAcom*100/(PERcom+ROcom-TAcom)) + '  [%]')
    print('Relative runoff:                 ' + str(ROcom*100/(PERcom+ROcom-TAcom)) + '  [%]')
    print('______')
    print('')
    print('Absolute error: ' + str(Erabs) + '  [L]')
    print('Relative error: ' + str(Errel) + '  [%]' + '   ' + str(warning))









