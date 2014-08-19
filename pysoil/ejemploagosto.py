
#  -------------- load Santa CRUZ clim data
# load clim data file
import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/Fluxsoil/Data/int_model_ago2011copy1.csv', "rU") as csvfile:
    csv_file = csv.reader(csvfile, delimiter=',')
    # eliminate blank rows if they exist
    rows = [row for row in csv_file if row]
    taille = len(rows)
    headings = rows[0] # get headings


sup=taille
TF=np.zeros(sup)
ETP=np.zeros(sup)
rows=np.array(rows)
for i in range(1,int(sup)):
	TF[i-1]=float(rows[i][1])
	ETP[i-1]=float(rows[i][3]) #3 o 4

TF=TF/(1000*15*60)
ETP=ETP/(1000*15*60)


##Version 3
# Soil characteristics 
soil_carac = {'Ksat':1e-4, 'Ss':0, 'eta':0.368, 'theta_r' : 0.102, 'theta_s' : 0.368, 'n':2, 'alpha':3.35, 'FC1':0.2, 'FC2':0.36}

# time and space discretization
days=10
T=60.*15*4*24*days
timestep=np.linspace(0.,(T-15*60),(24*days*4))
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
h_top = -TF
bc = {'time':timestep,'top':['fixed_flow',h_top],'bot':['free_drainage']}

# source term
q = ETP

# initial condition
h_init = np.array([-0.73]*I)
#for i in range(40,71):
#	h_init[i]=-0.83



# --------------------- Simulation run --------------------------
t, tt, z, S, Theta, Ka, Flow, INF, RO, PER, TA, VOL = run_varsat(L, T, dz, tstep, h_init, bc, q, soil_carac)




#plt.plot(np.asarray(S)[:,N],z)
plt.plot(100.*np.asarray(S)[:,-1],100.*z)

TimeLevel=np.array([0.,T])

mass_balance(TimeLevel,tt,INF,PER,TA,RO,VOL)



import csv
with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/compHydrus/AGOSTO2011/test.csv', 'wb') as test_file:
    AA = csv.writer(test_file, delimiter=',')
    AA.writerow(np.asarray(S)[:,-1]) 






for n in range(0,len(timestep)):
    plt.hold(False)     
    plt.plot(np.asarray(S)[:,n],100*z)
    #plt.draw()
    #sleep(0.05)   
    plt.draw()
    sleep(0.05)



