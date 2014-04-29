# Numerical solution for 1D Richards equation in homogeneous medium.
# discretization: finite difference, implicit (backward)
# A. Pryet 18/04/13 - alexandre.pryet@ensegid.fr
# z-axis is vertical upward


# --------------------------------------------------------------
#------------------ Load functions ----------------------------
# --------------------------------------------------------------

execfile('varsat_1D.py')


#  -------------- load Santa CRUZ clim data
# load clim data file
import csv
with open('./data_sz/clim_vls_ago.csv', 'rb') as csvfile:
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


