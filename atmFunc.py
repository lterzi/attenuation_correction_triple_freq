#Functions to convert between different 
#atmospheric variables

import numpy as np

# This function calculates the relative
# humidity.
# 
# input:
# specific humidity [Kg/Kg]: speHum
# temperature [K]: temp
# pressure [Pa]: press
# 
# output:
# relative humidity [%]: relHum
#
def speHumi2RelHum(speHum, temp, press):
    w = speHum/(1-speHum)
    es = 611*np.exp(17.27*(temp - 273.16)/(237.3 + (temp - 273.16)))
    ws = 0.622*(es/(press-es))      
    relHum = 100 *(w/ws)
    
    return relHum

# This function calculates the vapor
# pressure.
# 
# input:
# specific humidity [Kg/Kg]: speHum
# temperature [K]: temp
# pressure [Pa]: press
# 
# output:
# vapor Pressure [Pa]: e
#
def calcVaporPress(speHum, temp, press):
    w = speHum/(1-speHum)
    es = 611*np.exp(17.27*(temp - 273.16)/(237.3 + (temp - 273.16)))
    e = w*(press - es)/0.622
    
    return e

# This function calculates the vapor
# density.
# 
# input:
# temperature [K]: temp
# vapor pressure [Pa]: press
# 
# output:
# vapor density [kg/m^3]: density
#
def calcVaporDens(vaporPress, temp):
    density = 0.0022*vaporPress/temp
    
    return density


# This function calculates the dray air
# density.
# 
# input: 
# temperature [K]: temp
# pressure [Pa]: press
# water vapor density [kg/m^3]: waterVaporDens
# 
# output:
# dry air density [kg/m^3]: dryAirDens
#
def calcDryAirDens(press, waterVaporDens, temp):
    
    molMassDryAir = 0.0289644 #Kg/mol
    idealGasCons = 8.31447 #J/mol K
    
    atmDens = (press * molMassDryAir)/(idealGasCons * temp)
    dryAirDens = atmDens - waterVaporDens
    
    return dryAirDens

# This function calculates the total
# optical depth.
# 
# input:
# two way atmospheric attenuation [dB]: attAtm2Way
# 
# output:
# total optical depth []: optDepTot 
#
def calcAtt2Opt(attAtm2Way):

    optDepTot = 1*np.log(10**(attAtm2Way/20.))
    
    return optDepTot
