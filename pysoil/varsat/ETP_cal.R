# CÁLCULO DE LA EVAPOTRANSPIRACIÓN MÉTODO PENMAN-MONTEITH [FAO] MODIFICADO A NIVEL de 15min

# PASOS A SEGUIR PARA CORRER EL MODELO

#1.- Ingresar tablas con datos meteorológicos (Obtener de tabla)
	#t.media= temperatura media 15min [°C]
	#RH.media= humedad relativa media 15min [%]
	#date= fechas a nivel 15min as.POSIXct 
	#solar.rad= radiación solar media 15min [W]
	#uz= velocidad del viento a la altura [huz] desde la superficie [m/s]
	
	
	
#2.- Ingresar los siguientes datos:
	#z= elevación sobre el nivel del mar [m]
	#lat.grad= latitud en grados del lugar [grados] (hemisferio sur el valor es negativo)
	#long.grad=longitud en grados del lugar [grados](0-360)
	#huz= altura de medición del viento desde la superficie [m]
	#hRH= altura de medición de la humedad relativa [m]
	#h= altura media del cultivo o vegetación [m]
	#LAI= indice folial / leaf area index 
	#alpha= albedo
	
# # 	z=400
	# lat.grad=-0.672112	#valores positivos al norte y negativos al sur
	# long.grad=-90.323619	#valores negativos al oeste de Greenwich y positivos al este
	# huz=3.2	#no se realiza corrección de la velocidad del viento, se asume que los valores de velocidad del viento medida a esta altura es igual a la velocidad del viento a 2m.
	# hRH=2
	# h=1.2	#sacado de Fuente-Tomai, 2010
	# LAI=2.88	#Es el valor del cultivo de referencia de la FAO
	# alpha=0.23	#Es el valor del cultivo de referencia de la FAO
	# rs=100

#3.- Cargar la función etp


# FUNCIÓN PARA EL CÁLCULO DE LA EVAPOTRANSPIRACIÓN POTENCIAL - PENMAN-MONTEITH
	
etp <- function (date, t.mean, RH.mean, solar.rad, uz, option=T) 
{
# PARÁMETROS ATMOSFÉRICOS/ATMOSPHERIC PARAMETERS

#[P]:Presión Atmosférica/Atmospheric Pressure [kPa]
	#[z]: elevación sobre el nivel del mar [m]

P = 101.3 * ((293 - 0.0065*z)/293)^5.26 #[kPa]

#[gamma]:Constante Psycométrica/Psychometric constant [kPa/°C]
	#[P]:Presión Atmosférica/Atmospheric Pressure [kPa]
	#[lambda]:Calor latente de vaporización/Latent heat of vaporization [MJ kg^-1]
	lambda = 2.45 #[MJ kg^-1]
	#[cp]:Calor específico a presión constante/specific heat ar constant pressure [MJ kg^-1 °C^-1]
	cp = 1.013 * 10^-3  #[MJ kg^-1 °C^-1]
	#[epsilon]:ratio molecular weight of water vapour/dry air
	epsilon = 0.622

gamma = (cp*P)/(epsilon*lambda) #[kPa/°C]
#gamma = 0.665 * 10^-3 * P

# TEMPERATURA DEL AIRE / AIR TEMPERATURE 


#[rho]: densidad del aire media a presión constante / mean air density at constant pressure [kg/m^3]
	#[R]: constante de gas específico [kJ/kg*K]
	R= 0.287 #[kJ/kg*K]

rho= P/((1.01*(t.mean+273))* R) #[kg/m^3]


# HUMEDAD DEL AIRE/AIR HUMIDITY

#[es]:Presión de vapor de saturación media/Mean saturation vapour pressure [kPa]
	#[t.max] [t.min]: Temperatura del aire / Air Temperature [°C]
	
es = 0.6108*exp((17.27*t.mean)/(t.mean+237.3)) #[kPa]

#[m]:Pendiente curva presion de vapor de saturación/Slope of saturation vapour pressure curve [kPa °C^-1]
	#[t.mean]: Temperatura del aire / Air Temperature [°C]

m = 4098 * (0.6108 * exp ((17.27 * t.mean)/(t.mean + 237.3)))/(t.mean + 237.3)^2 #[kPa °C^-1]


#[ea]:Presión de vapor actual/Actual vapour pressure [kPa]
	#[RH.max] [RH.min]: Humedad relativa / Relative humidity [%]
	#[es]:Presión de vapor de saturación [kPa]

	ea = (es*RH.mean/100) #[kPa]


# RADIACIÓN/RADIATION 

#[Ra]:Radiacion extraterrestre/Extraterrestial radiation [MJ m^-2 hour ^-1]
	#[lat]:Latitud/Latitude [rad] 
	lat = lat.grad*pi/180 #[rad]
	#[J]:Número día del año/Number of day in the year [0-365/366]
	J = as.numeric(format(date,format="%j")) # [0-365/366]
	#[dr]:Dist. relativa inversa Tierra-Sol/Inverse relative distance Earth-Sun [rad]
	dr = 1 + 0.033*cos(2*pi*J/365) #[rad]
	#[delta]:Declinación Solar/Solar declination [rad]
	delta = 0.409*sin(2*pi*J/365 - 1.39) #[rad]
	#[Gsc]:Constante Solar/Solar Constant [MJ m^-2 min^-1]
	Gsc = 0.082 #[MJ m^-2 min^-1]
	
	#Calculo de [omega]:Ángulo horario al ocaso y amanecer/Sunset hour angle [rad]
	#correcion estacional para el tiempo solar
	cons.b<-2*pi*(J-81)/364
	Sc=0.1645*sin(2*cons.b)-0.1255*cos(cons.b)-0.025*sin(cons.b)
	omegasol = acos (-tan (lat) * tan (delta)) #[rad]
	#angulo solar en el momento que ocurre el punto medio del periodo considerado
	H = as.numeric(format(date,format="%H")) # [0-365/366]
	tstep<-0.25	#duración del perido considerado, 1 para hora, 0.5 para 30 minutos, etc
	if(long.grad>=0&long.grad<=180){
		long<-360-long.grad
	}else{
	long=-long.grad
	}
	
	omega=pi/12*((H+0.5)+0.06667*(90-long)+Sc-12)
	omega1<-omega-pi*tstep/24
	omega2<-omega+pi*tstep/24
	
Ra = (12*60/pi) * Gsc * dr * ((omega2-omega1)*sin(lat)*sin(delta) + cos (lat)*cos(delta)*(sin(omega2)-sin(omega1))) #[MJ m^-2 hour ^-1]

#[N]:Horas de luz/Daylight hours [hour]
	#[omega]:Ángulo horario al ocaso/Sunset hour angle [rad]
	
N = 24*omegasol/pi #[hour]


#[Rs]:Radiación Solar/Solar Radiation [MJ m^-2 hour^-1]
Rs = (0.0864)/(24*4)*solar.rad #[MJ m^-2 day^-1]

#[Rso]:Radiación solar a cielo despejado/Clear-sky Solar Radiation [MJ m^-2 hour^-1]
	#[Ra]:Radiacion extraterrestre/Extraterrestial radiation [MJ m^-2 hour ^-1]
	#[z]: elevación sobre el nivel del mar [m]

Rso =  (0.75 + 2*10^-5*z)*Ra #[MJ m^-2 hour^-1]

#[Rns]:Radiación onda corta neta/Net shortwave radiation [MJ m^-2 hour^-1]
	#[alpha]: albedo
	#[Rs]:Radiación Solar/Solar Radiation [MJ m^-2 hour^-1]

Rns = (1-alpha)*Rs #[MJ m^-2 hour^-1]


#[Rnl]:Radiación onda larga neta/Net longwave radiation [MJ m^-2 hour^-1]
	#[t.mean]: Temperatura del aire / Air Temperature [°C]
	#[ea]:Presión de vapor actual/Actual vapour pressure [kPa]	
	#[Rs]:Radiación Solar/Solar Radiation [MJ m^-2 hour^-1]
	#[Rso]:Radiación solar a cielo despejado/Clear-sky Solar Radiation [MJ m^-2 hour^-1]
	#Stefan-Boltzmann constant [sigma]
	sigma = (4.903*10^-9)/(24*4) #[MJ K^-4 m^-2 hour^-1]
	cons.R<-0.4	#valor que se utiliza para la noche
Rnl<-c()
for(i in 1:length(H)){
if(H[i]>=6&H[i]<=18){
	Rnl[i] = sigma*(t.mean[i] + 273.16)^4 * (0.34-0.14*ea[i]^(1/2)) * ((1.35*Rs[i]/Rso[i])-0.35) #[MJ m^-2 hour^-1]
	}else{
	Rnl[i] = sigma*(t.mean[i] + 273.16)^4 * (0.34-0.14*ea[i]^(1/2)) * ((1.35*cons.R)-0.35) #[MJ m^-2 hour^-1]
	}
}
#[Rn]:Radiación Neta/Net Radiation [MJ m^-2 hour^-1]
	#[Rns]:Radiación onda corta neta/Net shortwave radiation [MJ m^-2 hour^-1]
	#[Rnl]:Radiación onda larga neta/Net longwave radiation [MJ m^-2 hour^-1]
	
Rn = Rns - Rnl #[MJ m^-2 hour^-1]

#[G]:Flujo Calor del suelo/Soil heat flux [MJ m^-2 hour^-1]
G<-c()
for(i in 1:length(H)){
if(H[i]>=6&H[i]<=18){
	G[i]<-0.1*Rn[i]
	}else{
	G[i]<-0.5*Rn[i]
	}
}	# [MJ m^-2 hour^-1] es cero para cálculo a nivel diario
 

# VELOCIDAD DEL VIENTO [u2] : en este modelo no se necesita pero se la tiene como una referencia

#[u2]:velocidad a 2m desde la superficie [m/s]
	#[uz]:velocidad del viento a [huz] m desde la superficie [m/s]
	#[huz]:altura de medición velocidad del viento [m]

#u2 = uz / (4.87/log(67.8*huz - 5.42)) #[m s^-1]


#[ra]:RESISTENCIA AERODINÁMICA [s m^1]

	#[huz]:altura de medición velocidad del viento / height of wind measurements [m]
	#[hRH]:altura de medición de humedad / height of humidity measurements [m]
	#[uz]: velocidad del viento a la altura huz / wind speed at huz [m/s]
	
	#[d]:altura del plano de desplazamiento cero / zero plane displacement height [m]
	#d= 0.75*h #[m]
	d= 2*h/3 #[m]
	#[zom]:longitud de la rugosidad que rige la transferencia de momento / roughness length governing momentum transfer [m]
	#zom= 0.1*h #[m]
	zom= 0.123*h #[m]
	#[zoh]:longitud de la rugosidad que rige la transferencia de calor y vapor / roughness length governing transfer of heat and vapour [m]
	#zoh= 0.123*zom #[m]
	#[k]:constante de von Karman / von Karman's constant
	k= 0.4

ra= (log((huz-d)/zom))^2/(k^2*uz) #[s m^1]
#ra= (log((huz-d)/zom))^2/(k^2*u2) #[s m^1]

#[rs]:RESISTENCIA SUPERFICIAL [s m^1]

	#[r1]:resistencia stomatal de hoja bien iluminada / bulk stomatal resistance of well-illuminated leaf [s m^-1]
	#r1= 100 #[s m^-1]
	#[LAIactive]:índice folial activo / active leaf area index [LAIactive]
		#[LAI]: indice folial / leaf area index 
	#LAIactive= 0.5 * LAI

#rs= r1/LAIactive #[s m^1]
#rs=100	

#[etp]:CÁLCULO EVAPOTRANSPIRACIÓN POTENCIAL [mm hour^-1]
if (option==TRUE){
etp = (0.408 * ((m*(Rn)) + ((86400/(24*4))*rho*cp*(es-ea)/ra))) / (m + gamma*(1+(rs/ra))) #[mm hour^-1]
}else{
etp = (0.408 * ((m*(Rn)) + ((86400/(24*4))*rho*cp*(es-ea)/ra))) / (m + gamma)
}
etp[etp<0]=0

#PRESENTACIÓN DE RESULTADOS
return (etp)
}


###############################
####Read clim_data
#####Climate Data
tz <- "America/Guayaquil"
# 
clim_vls<- read.csv("clim_data300.csv")
as.POSIXct(as.character(clim_vls$date), format= "%Y-%m-%d %H:%M:%S", tz= tz)->clim_vls$date
###Load Rs parameter
#Rs<-read.csv(file.choose(),header=F)
Rs<-read.csv("param.dat",header=F)
rs=Rs[1,1]

###Load Period
#per=read.csv(file.choose(),header=F)
per=read.csv("period.dat",header=F)
sup0=per$V1[1]+1
sup=per$V1[2]


##Datos SC300
z=300
lat.grad=-0.91563	#valores positivos al norte y negativos al sur
long.grad=-89.47456	#valores negativos al oeste de Greenwich y positivos al este
h=0.13	#sacado de Fuente-Tomai, 2010
huz=h+2	#no se realiza corrección de la velocidad del viento, se asume que los valores de velocidad del viento medida a esta altura es igual a la velocidad del viento a 2m.
#LAI=4	#Es el valor del cultivo de referencia de la FAO
alpha=0.23	#Es el valor del cultivo de referencia de la FAO


###Run Model
clim_vls.p<-subset(clim_vls,date>=clim_vls$date[sup0] & date<=clim_vls$date[sup])
clim_vls.p$T_mm_Tot<-etp(clim_vls.p$date,clim_vls.p$Temp,clim_vls.p$RH,clim_vls.p$SlrW,clim_vls.p$WS_ms,option=T)

###Write 
#write.table (clim_vls.p$T_mm_Tot,file.choose(), sep=",", row.names=F)
write.table (clim_vls.p$T_mm_Tot, "transp.csv", sep=",", row.names=F)




	