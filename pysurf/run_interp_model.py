
##Read clim data...

import csv
#with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/Interception/Data/clim_data2.csv', "rU") as csvfile:
with open('clim_data2.csv', "rU") as csvfile:
    csv_file = csv.reader(csvfile, delimiter=',')
    # eliminate blank rows if they exist
    rows = [row for row in csv_file if row]
    taille = len(rows)
    headings = rows[0] # get headings


#sup=20544
sup0=1	#inicio en la garua
#sup0=20545	#inicio en invierno
#sup=20544	#termina en garua
sup=taille-1	#termina en invierno
Rain_clim=np.zeros(sup)
Thfall_clim=np.zeros(sup)
Evap_clim=np.zeros(sup)
Tran_clim=np.zeros(sup)
Season_clim=np.zeros(sup)
Fogpr_clim=np.zeros(sup)

for i in range(sup0,int(sup)):
	Rain_clim[i-1]=float(rows[i][7])
	Thfall_clim[i-1]=float(rows[i][9])
	Evap_clim[i-1]=float(rows[i][10])
	Tran_clim[i-1]=float(rows[i][11])
	Season_clim[i-1]=float(rows[i][12])
	Fogpr_clim[i-1]=float(rows[i][13])

##Climate Data
clim_data = {'Rain':Rain_clim, 'Thfall': Thfall_clim, 'Evap':Evap_clim, 'Tran':Tran_clim}

##Canopy Data
p_clim=0.1
S_clim=0.5
Ds_clim=0.00000000000109265
b_clim=8.83
#Ds_clim=0.0000000223
#b_clim=7.71

canopy_carac = {'p':p_clim, 'S':S_clim, 'Ds':Ds_clim, 'b':b_clim}

int_model=interception_model(clim_data,canopy_carac, drain_func=2)


#for i in range(0, len(int_model['TF'])):
#    if int_model['C'][i]>=0.1*canopy_carac['S']:
#        clim_data['Tran'][i]=float(0)


#for i in range(0, len(int_model['TF'])):
#    if int_model['C'][i]>=canopy_carac['S']:
#        clim_data['Tran'][i]=float(0)
#    elif int_model['C'][i]>=0.05 and int_model['C'][i]<canopy_carac['S']:
#        clim_data['Tran'][i]=clim_data['Tran'][i]*int_model['C'][i]/canopy_carac['S']
        
        

#for i in range(0,len(int_model['TF'])):
#	if int_model['C'][i]>=0.03:
#		clim_data['Tran'][i]=float(0)

#with open('/Users/redhotkenny/Documents/Documentos/Doctorado/SoilModel/Interception/Data/output.dat', 'w') as outputfile :
with open('output.dat', 'w') as outputfile :
	for i in range(0,len(int_model['TF'])):
		outputfile.write( str(round(clim_data['Rain'][i],5)) + ',' + str(round(clim_data['Thfall'][i],5)) + ',' +  str(round(int_model['E'][i],5)) + ',' + str(round(clim_data['Tran'][i],5)) + ','  + str(round(int_model['C'][i],5)) + ',' + str(round(int_model['D'][i],5)) + ',' + str(round(int_model['TF'][i],5)) + ',' + str(round(int_model['CWI'][i],5))  + '\n' )



#with open('output.dat', 'w') as outputfile :
#	for i in range(0,len(int_model['TF'])):
#	    outputfile.write( str(round(int_model['TF'][i],5)) + '\n' )



#with open('obs.dat', 'w') as outputfile :
#	for i in range(0,len(clim_data['Thfall'])):
#	    outputfile.write( str(clim_data['Thfall'][i]) + '\n' )


