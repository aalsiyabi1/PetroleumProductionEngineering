#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import math 


# In[ ]:


#Calculating solution GOR
def Rs(gSG,p,t,API):
    nume = p * (10 ** (0.0125 * API)) #numerator
    denom = 18 * (10 ** (0.00091 * t)) #denominator
    Rs = gSG * ((nume / denom) ** 1.2048)
    return Rs


#
#Calculating Density of oil / API gravity of stock tank oil
def API(oSG):
    API = (141.5 / oSG) - 131.5 #(API degree)
    return API

#Calculating Specific gravity of stock tank oil 
def oSG(API):
    oSG = 141.5 / (API + 131.5) 
    return oSG

#Calculating density of oil 
def oilDensity(Rs,oSG,gSG,t):
    '''
    Calculate density of oil using Standing correlation
    Output:
    oilDensity : density of stock tank oil (lbm/ft3)
 
    Input:
    Rs : Solution ga-oil-ratio (scf/stb)
    oSG : Specific gravity of of stock tank oil (1 for freshwater)
    gSG : Specific gravity of gas (1 for air)
    t : temperature (F)
    '''
    oilDensity = (62.4 * oSG + 0.0136 * Rs * gSG) / (0.972 + 0.000147 * (Rs * math.sqrt(gSG / oSG) + 1.25 * t) ** 1.175)
    return oilDensity

#Calcualting FVF of oil
def Bo(Rs,oSG,gSG,t):
    '''
    Calculate formation volume factor of oil
    Output:
    Bo : Formation volume factor of oil (rb/stb)
 
    Input:
    Rs : Solution ga-oil-ratio (scf/stb)
    oSG : Specific gravity of of stock tank oil (1 for freshwater)
    gSG : Specific gravity of gas (1 for air)
    t : temperature (F)
    '''
    Bo = 0.9759 + 0.00012 * (Rs * math.sqrt(gSG / oSG) + 1.25 * t) ** 1.2
    return Bo
    
#Calculating viscosity of dead oil
def viscOD(API,t):
    '''
    Calculate viscosity of dead oil using Standing's correlation
    Output:
    viscOD : Viscosity of dead oil (cp)
    
    Input:
    API : API gravity of oil (API)
    t : temperature (F)
    '''
    A = 10 ** (0.43 + (8.33 / API))
    viscOD = (0.32 + (1.8 * (10 ** 7)) / (API ** 4.53)) * ((360 / (t + 200))** A)
    return viscOD

#Calculating viscosity of saturated oil
def viscOB(viscOD,Rs):
    '''
    Calculate viscosity of saturated oil using Standing's correlation
    Output:
    viscOB : Viscosity of saturated oil (cp)
    
    Input:
    viscOD : Viscosity of dead oil (cp)
    Rs : Solution ga-oil-ratio (scf/stb)
    '''
    a = Rs * (2.2 * (10 ** -7) * Rs - 7.4 * (10 ** -4))
    c = 8.62 * Rs * (10 ** -5)
    d = 1.10 * Rs * (10 ** -3)
    e = 3.74 * Rs * (10 ** -3)
    b = (0.68 / (10 ** c)) + (0.25 / (10 ** d)) + (0.062 / (10 ** e))
    
    viscOB = (10 ** a) * (viscOD ** b)
    return viscOB

#Calculating viscosity of unsaturated oil
def viscO(viscOB,p,pb):
    '''
    Calculate viscosity of unsaturated oil using Standing's correlation
    Output:
    viscO : Viscosity of unsaturated oil (cp)
    
    Input:
    viscOB : Viscosity of saturated oil (cp)
    p : Pressure (psia)
    pb : Bubble-Point Pressure (psia)
    '''
    viscO = viscOB + 0.001 * (p - pb) * (0.024 * (viscOB ** 1.6) + 0.38 * (viscOB ** 0.56))
    return viscO


# In[ ]:


#example 2.1 (calculate oil density and viscosity)
Rs = 600
p = 4475
t = 140
pb = 2745
api = 35
gSG = 0.77

oilSG = oSG(api)
rhoO = oilDensity(Rs,oilSG,gSG,t)
muOD = viscOD(api,t)
muOB = viscOB(muOD,Rs)
muO = viscO(muOB,p,pb)

print('Density of oil = {},\n viscosity of oil = {}'.format(rhoO,muO))


# #### 
