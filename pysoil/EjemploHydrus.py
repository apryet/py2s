#K


h_vals = np.arange( -5, 1, 0.01 )
Kr_vals = [ Kr_curve(h_val,soil_carac) for h_val in h_vals ]
plt.ion()
plt.plot(-h_vals,Kr_vals,label='Kr(h)')





#  -------------- load Santa CRUZ clim data
# load clim data file
import csv
with open(askopenfilename()) as csvfile:
    csv_file = csv.reader(csvfile, delimiter=',')
    # eliminate blank rows if they exist
    rows = [row for row in csv_file if row]
    headings = rows[0] # get headings
    clim_data= {}
    for row in rows[1:]:
        # append the dataitem to the end of the dictionary entry
        # set the default value of [] if this key has not been seen
        for col_header, data_column in zip(headings, row):
            clim_data.setdefault(col_header, []).append(data_column)

rain = np.array(clim_data['Rain_mm_Tot'],dtype='float')

# --------------------------------------------------------------
# ------------------- Preamble ---------------------------------
# --------------------------------------------------------------

# Check all values below

# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35}

# time and space discretization

T = 3600*24*10. # total simulation duration [s] 
L = 0.6 # model length [m]
dt = 60*60 # time step [s]
dz = 0.01 # mesh cell size [m]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

# boundary conditions
h_bot = -10
h_top = -0.75
bc = {'top':['fixed_head',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = 0

# initial condition
h_init = np.array([h_bot]*I)

# --------------------- Simulation run --------------------------
t, z, S = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)


# --------------------------------------------------------------
# ------------------- Post-processing --------------------------
# --------------------------------------------------------------


# ----------------- Plot snapshot at t_n -----------------------
plt.ion()
plt.hold(True)
plt.plot(np.asarray(S)[:,n],z)
plt.show()


# ----------------- Plot time series at x_i -------------------
i=10
plt.ion()
plt.plot(t,np.asarray(S)[i,:])
plt.show()

n=120
theta_profile = np.array( [theta_curve(h_val,G) for h_val in np.array(S)[:,n] ] )
plt.plot(theta_profile,z)

# -------------------- Plot animated graph ---------------------

for n in range(0,N):
    plt.hold(False)     
    plt.plot(np.asarray(S)[:,n],z)
    #plt.draw()
    #sleep(0.05)   
    plt.draw()
    sleep(0.05)





#load different scenarios

# --------------------------------------------------------------
# ------------------- Preamble ---------------------------------
# --------------------------------------------------------------

# Check all values below

# Soil characteristics 
soil_carac = {'Ksat':9.44e-5, 'Ss':0, 'eta':0.287, 'theta_r' : 0.075, 'theta_s' : 0.287, 'n':4.1, 'alpha':3.3}

# time and space discretization

T = 1200. # total simulation duration [s] 
L = 0.7 # model length [m]
dt = 3 # time step [s]
dz = 0.01 # mesh cell size [m]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

# boundary conditions
h_bot = -0.315
h_top = -3.29/3600
bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = 0

# initial condition
h_init = np.array([h_bot]*I)

# --------------------- Simulation run --------------------------
t, z, S = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)


# --------------------------------------------------------------
# ------------------- Post-processing --------------------------
# --------------------------------------------------------------


# ----------------- Plot snapshot at t_n -----------------------
plt.ion()
plt.hold(True)
plt.plot(np.asarray(S)[:,n],z)
plt.show()


# ----------------- Plot time series at x_i -------------------
i=10
plt.ion()
plt.plot(t,np.asarray(S)[i,:])
plt.show()

n=120
theta_profile = np.array( [theta_curve(h_val,G) for h_val in np.array(S)[:,n] ] )
plt.plot(theta_profile,z)

# -------------------- Plot animated graph ---------------------

for n in range(0,N):
    plt.hold(False)     
    plt.plot(np.asarray(S)[:,n],z)
    #plt.draw()
    #sleep(0.05)   
    plt.draw()
    sleep(0.05)



##scenario 2

#load different scenarios

# --------------------------------------------------------------
# ------------------- Preamble ---------------------------------
# --------------------------------------------------------------

# Check all values below

# Soil characteristics 
soil_carac = {'Ksat':9.44e-5, 'Ss':0, 'eta':0.287, 'theta_r' : 0.075, 'theta_s' : 0.287, 'n':4.1, 'alpha':3.3, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 600. # total simulation duration [s] 
L = 0.7 # model length [m]
dt = 3 # time step [s]
dz = 0.01 # mesh cell size [m]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot = -0.315
h_top = -3.29/3600
bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = 0

# initial condition
h_init = np.array([h_bot]*I)

# --------------------- Simulation run --------------------------
t, z, S, INF, RO = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)



A=np.array([[1.,2.,3.,4.]])
B=np.array([[5.,6.,7.,8.]])

plt.ion()
plt.hold(True)
plt.plot(np.asarray(S)[:,N],z)
plt.show()




import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import * 
from time import *
from math import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
from scipy.sparse.linalg import cgs
from scipy.sparse import csr_matrix











####EJEMPLO PARA HYDRUS

#load different scenarios

# --------------------------------------------------------------
# ------------------- Preamble ---------------------------------
# --------------------------------------------------------------

# Check all values below

# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 3600. # total simulation duration [s] 
L = 0.7 # model length [m]
dt = 1 # time step [s]
dz = 0.01 # mesh cell size [m]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot = -.97
h_top = -0.001/3600
bc = {'top':['fixed_flow',[h_top]*N],'bot':['free_drainage']}

# source term
q = 0

# initial condition
h_init = np.array([h_bot]*I)

# --------------------- Simulation run --------------------------
t, z, S, INF, RO = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)




plt.plot(np.asarray(S)[:,N],z)
import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/compHydrus/test.csv', 'wb') as test_file:
    AA = csv.writer(test_file, delimiter=',')
    AA.writerow(np.asarray(S)[:,N]) 







# Check all values below

# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 88. # total simulation duration [s] 
L = 0.7 # model length [m]
dt = 1. # time step [s]
dz = 0.01 # mesh cell size [m]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot=-0.97
#h_bot = -0.00001/3600
h_top = -0.001/3600
bc = {'top':['fixed_flow',[h_top]*N],'bot':['free_drainage']}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_flow',[h_bot]*N]}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = 0

# initial condition
h_init = np.array([-0.97]*I)

# --------------------- Simulation run --------------------------
t, z, S, INF, RO = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)




#plt.plot(np.asarray(S)[:,N],z)
plt.plot(100.*np.asarray(S)[:,N],100.*z)


import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/compHydrus/test.csv', 'wb') as test_file:
    AA = csv.writer(test_file, delimiter=',')
    AA.writerow(100.*np.asarray(S)[:,N]) 





soil_carac = {'Ksat':1e-2, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':0.0335, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 3600. # total simulation duration [s] 
L = 70 # model length [cm]
dt = 1 # time step [s]
dz = 1 # mesh cell size [cm]

N = int(round(T/dt)) # number of time steps
I = int(round(L/dz)) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot = -0.001/3600
h_top = -0.1/3600
bc = {'top':['fixed_flow',[h_top]*N],'bot':['free_drainage']}

# source term
q = 0

# initial condition
h_init = np.array([-97.]*I)

# --------------------- Simulation run --------------------------
t, z, S, INF, RO = run_varsat(L, T, dz, dt, h_init, bc, q, soil_carac)




plt.plot(np.asarray(S)[:,N],z)


import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/compHydrus/test.csv', 'wb') as test_file:
    AA = csv.writer(test_file, delimiter=',')
    AA.writerow(np.asarray(S)[:,N]) 







# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 3600. # total simulation duration [s] 
timestep=np.array([0.,600.,1200.,1800.,2400.,3000.])
tstep=3.
L = 0.70 # model length [m]
#dt = 60 # time step [s]
dz = 0.01 # mesh cell size [m]
I = int(round(L/dz)+1) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot=-0.97
#h_bot = -0.00001/3600
h_top = -0.001/3600
bc = {'time':timestep,'top':['fixed_flow',[h_top]*len(timestep)],'bot':['free_drainage']}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_flow',[h_bot]*N]}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = np.array([0,0,0,0,0,0])

# initial condition
h_init = np.array([-0.97]*I)

# --------------------- Simulation run --------------------------
t, z, S, INF, RO, Theta = run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac)




#plt.plot(np.asarray(S)[:,N],z)
plt.plot(100.*np.asarray(S)[:,-1],100.*z)

plt.plot(100.*np.asarray(Theta)[:,-1],100.*z)



import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/compHydrus/test.csv', 'wb') as test_file:
    AA = csv.writer(test_file, delimiter=',')
    AA.writerow(100.*np.asarray(S)[:,-1]) 




##Version 3
# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.25}

# time and space discretization

T = 3600. # total simulation duration [s] 
timestep=np.array([0.,600.,1200.,1800.,2400.,3000.])
tstep=3.
L = 0.70 # model length [m]
#dt = 60 # time step [s]
dz = 0.01 # mesh cell size [m]
I = int(round(L/dz)+1) # number of cells

#Root Information
LR=0.5
type='distribution'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
h_bot=-0.97
#h_bot = -0.00001/3600
h_top = -0.001/3600
bc = {'time':timestep,'top':['fixed_flow',[h_top]*len(timestep)],'bot':['free_drainage']}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_flow',[h_bot]*N]}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = np.array([0,0,0,0,0,0])

# initial condition
h_init = np.array([-0.97]*I)

# --------------------- Simulation run --------------------------
t, tt, z, S, Theta, Ka, Flow, INF, RO, PER, TA, VOL = run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac)




#plt.plot(np.asarray(S)[:,N],z)
plt.plot(100.*np.asarray(S)[:,-1],100.*z)

TimeLevel=np.array([0.,3600.])

mass_balance(TimeLevel,tt,INF,PER,TA,RO,VOL)



##Version 3
#Ejemplo con etp
# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.36}

# time and space discretization

T = 3600. # total simulation duration [s] 
timestep=np.array([0.,600.,1200.,1800.,2400.,3000.])
tstep=60.
L = 0.70 # model length [m]
#dt = 60 # time step [s]
dz = 0.01 # mesh cell size [m]
I = int(round(L/dz)+1) # number of cells

#Root Information
LR=0.7
type='homogeneous'
#get Root Distribution
z = np.linspace(0,L,I)
RD = get_RD(L, I, LR, z, dz, type) 


# boundary conditions
#h_bot=-0.5
#h_bot = -0.00001/3600
h_top = -0.001/3600
bc = {'time':timestep,'top':['fixed_flow',[h_top]*len(timestep)],'bot':['free_drainage']}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_flow',[h_bot]*N]}
#bc = {'top':['fixed_flow',[h_top]*N],'bot':['fixed_head',[h_bot]*N]}

# source term
q = [0.0005/3600]*len(timestep)
#q = [0.]*len(timestep)

# initial condition
h_init = np.array([-0.5]*I)

# --------------------- Simulation run --------------------------
t, tt, z, S, Theta, Ka, Flow, INF, RO, PER, TA, VOL = run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac)




#plt.plot(np.asarray(S)[:,N],z)
plt.plot(100.*np.asarray(S)[:,-1],100.*z)

TimeLevel=np.array([0.,3600.])

mass_balance(TimeLevel,tt,INF,PER,TA,RO,VOL)






