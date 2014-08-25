####RUN INTERCEPTION MODEL FOR PEST OPTIMIZATION


# -- Load modules
import numpy as np
import scipy as sp
from time import *
from math import *

# -----------------------------------------------------------------------------
# -- Drainage functions -------------------------------------------------------
# -----------------------------------------------------------------------------

# Drainage function by Rutter 1975 y Guevara-Escobar, 2007
def drainage_func1(canopy_carac,cc):
    return(canopy_carac['Ds']*(cc)**canopy_carac['b'])

# Drainage function by Pitmann 1989, Domingo 1998 : PASTURE
def drainage_func2(canopy_carac,cc):
    return(canopy_carac['Ds']*exp(canopy_carac['b']*(cc)))


# Drainage function b  : FOREST
def drainage_func3(canopy_carac,cc):
    if cc > canopy_carac['S']:
	return(canopy_carac['Ds']*exp(canopy_carac['b']*(cc-canopy_carac['S'])))
    else : 
	return(0)

# -----------------------------------------------------------------------------
# Interception model ----------------------------------------------------------
# -----------------------------------------------------------------------------

def interception_model(clim_data, canopy_carac, compute_cwi = False, delta_t = 900, drainage_func = drainage_func2):

        # clim_data : dictionary containing climatic records (Rainfall, Wet canopy potential evaporation, Dry canopy potential evaporation, Throughfall, Stemflow, ...
	# canopy_carac : dictionary with canopy characteristics
	# compute_cwi : TRUE/FALSE, whether to compute or not cloud water interception
	# drainage_func : drainage fonction to be used
	# delta_t [s] : time interval between two measurements in clim_data

	# get number of observations
	nobs = len(clim_data['Rain'])

	# initialize output with :
	# C, simulated canopy water content
	# D, simulated canopy drainage
	# E, simulated canopy evaporation (i.e. interception)
	# TF,simulated throughfall
	# CWI, estimated cloud water interception
	clim_model={ 'C':np.zeros(nobs), 'D':np.zeros(nobs), 'E':np.zeros(nobs), 'T':np.zeros(nobs), 'TF':np.zeros(nobs), 'CWI':np.zeros(nobs) }
	
	# iterate over observations records
	for i in range(1,nobs):
	
		# initialize variables 
		cc0= clim_model['C'][i-1] # initial canopy storage 
		ee0=0  # initialize evaporation for current record
		tt0=0  # initialize transpiration for current record
		dd0=0  # initialize drainage for current record

		# define number of time steps within each delta_t
		if (cc0+clim_data['Rain'][i]) <= canopy_carac['S']: # if canopy water content is smaller than canopy storage capacity
			nstep = 1 # 1 time step of length delta_t
		else: # otherwise, use int(delta_t) time steps of 1 sec
			nstep= int(delta_t)

		# iterate over computational time steps
		for j in range(0,nstep):

			# canopy storage is updated with current rainfall
			# N.B. : input from CWI is not considered.
			cc = cc0 + clim_data['Rain'][i] / float(nstep) * (1-canopy_carac['p'])

			# evaporation and transpiration rate estimation
			if cc >= canopy_carac['S']: # if the canopy is wet, evaporation at full rate
				ee = clim_data['Evap_pot'][i]/float(nstep)*(1-canopy_carac['p'])
				tt = 0

			else : # otherwise, canopy evaporation proportional to canopy water content
				ee = clim_data['Evap_pot'][i]/float(nstep)*(1-canopy_carac['p'])*cc/canopy_carac['S']
				tt = clim_data['Tran_pot'][i]/float(nstep)*(1-canopy_carac['p'])*(1 - cc/canopy_carac['S'])
			
			# update evaporation for current record 
			ee0 = ee + ee0  
			tt0 = tt + tt0  

			# remove evaporated water from canopy
			cc = cc - ee 

			# correct storage if it goes below 0
			if cc < 0:
			    ee = ee + cc
			    cc = 0
				
			# drainage estimation
			if nstep == 1 :
			    dd = delta_t * drainage_func(canopy_carac,cc) 
			else :
			    dd = drainage_func(canopy_carac,cc)
			
			# update total drainage
			dd0 = dd0 + dd

			# remove drained water from canopy
			cc = cc - dd	
			
			# correct storage if it goes below 0
			if cc < 0:
				dd = dd + cc
				cc = 0

			# update canopy water content 
			cc0=cc

		# compute throughfall
		tf = dd0 + canopy_carac['p'] * clim_data['Rain'][i]

		# compute cloud water interception
		if compute_cwi == True : 
		    cwi = clim_data['Thfall'][i] - tf
		else : 
		    cwi = 0

		# write clim_model data #step1
		clim_model['C'][i] = cc
		clim_model['D'][i] = dd0
		clim_model['E'][i] = ee0
		clim_model['T'][i] = tt0
		clim_model['TF'][i] = tf
		clim_model['CWI'][i] = cwi

	return(clim_model)
