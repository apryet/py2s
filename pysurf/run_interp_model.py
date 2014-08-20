import csv
import sys
from os.path import expanduser
home = expanduser("~")
import numpy as np

sys.path.append( home +'/Programmes/py2s/pysurf/')

from interp_model import *

# -------------------------------------------------------------------------------------
# --- Import clim data  --------------------------------------------------
# -------------------------------------------------------------------------------------


with open('./clim_data.csv', "rU") as csvfile:
    csv_reader = csv.DictReader(csvfile, delimiter=',', dialect = 'excel')
    clim_data_raw = {}
    for row in csv_reader :
	for column, value in row.iteritems() :
	    clim_data_raw.setdefault(column, []).append(value)

#import pandas
#data = pandas.read_csv(r'..\data\data.csv')

# -------------------------------------------------------------------------------------
# --- Identify climate variables 
# -------------------------------------------------------------------------------------

##Climate Data
clim_data = {'Rain':np.array(clim_data_raw['Rain_mm_Tot'],dtype=float), \
    'Thfall': np.array(clim_data_raw['Thfall_mm'],dtype=float),
    'Evap_pot': np.array(clim_data_raw['E_mm_Tot'],dtype=float),
    'Tran_pot':np.array(clim_data_raw['T_mm_Tot'],dtype=float)
    }

# -------------------------------------------------------------------------------------
# --- Set canopy characteristics 
# -------------------------------------------------------------------------------------

##Canopy Data
p_clim=0.1
S_clim=0.5
Ds_clim=0.00000000000109265
b_clim=8.83
#Ds_clim=0.0000000223
#b_clim=7.71


canopy_carac = {'p':p_clim, 'S':S_clim, 'Ds':Ds_clim, 'b':b_clim}
# -------------------------------------------------------------------------------------
# --- Run model  --------------------------------------------------
# -------------------------------------------------------------------------------------

int_model=interception_model(clim_data,canopy_carac)

# -------------------------------------------------------------------------------------
# --- Write results  --------------------------------------------------
# -------------------------------------------------------------------------------------

with open('output.dat', 'wb') as outputfile : 
    csvwriter = csv.DictWriter(outputfile, delimiter=',', lineterminator='\n', fieldnames=int_model.keys())
    csvwriter.writeheader()
    for i in range(len(int_model.values()[0])) :
	csvwriter.writerow( { key:val for key,val in zip( int_model.keys(), [ val[i] for val in int_model.values() ] ) } )


