#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Example 2.3 (calculate gas viscosity with Carr correlation)

import math

from pseudoProps import *


# In[2]:


#Add docstring
def gasVisc(gSG,N2,CO2,H2S,p,t):

  ppcCalc = ppc(gSG,N2,CO2,H2S)
  tpcCalc = tpc(gSG,N2,CO2,H2S)
  pprCalc = ppr(p,ppcCalc)
  tprCalc = tpr(t,tpcCalc)
  
  #Corrections for gas viscosity at 14.7 psia
  gasViscUncorr = 8.188 * (10 ** -3) - 6.15 * (10 ** -3) * (math.log10(gSG)) + (1.709 * (10 ** -5) - 2.062 * (10 ** -6) * gSG) * t #uncorrected gas viscosity
  gasViscN2Corr = N2 * (9.59 * (10 ** -3) + 8.48 * (10 ** -3) * (math.log10(gSG))) #N2 correction  
  gasViscCO2Corr = CO2 * (6.24 * (10 ** -3) + 9.08 * (10 ** -3) * (math.log10(gSG))) #CO2 correction 
  gasViscH2SCorr = H2S * (3.73 * (10 ** -3) + 8.49 * (10 ** -3) * (math.log10(gSG))) #H2S correction 

  #Corrected gas viscosity at 14.7 
  gasViscCorr = gasViscUncorr + gasViscN2Corr + gasViscCO2Corr + gasViscH2SCorr  

  #Calculate gas viscosity at elevated pressure using Dempsey (1965) relation
  a0,a1,a2,a3,a4 = -2.46211820, 2.97054714, -0.28626405, 0.00805420, 2.80860949
  a5,a6,a7,a8,a9 = -3.49803305, 0.36037302, -0.01044324, -0.79338568, 1.39643306
  a10,a11,a12,a13,a14,a15 = -0.14914493, 0.00441016, 0.08393872, -0.18640885, 0.02033679, -0.00060958 

  gasViscR = a0 + a1 * pprCalc + a2 * (pprCalc ** 2) + a3 * (pprCalc ** 3) + tprCalc * (a4 + a5 * pprCalc + a6 * (pprCalc ** 2) + a7 * (pprCalc ** 3)) + (tprCalc ** 2) * (a8 + a9 * pprCalc + a10 * (pprCalc ** 2) + a11 * (pprCalc ** 3)) + (tprCalc ** 3) * (a12 + a13 * pprCalc + a14 * (pprCalc ** 2) + a15 * (pprCalc ** 3))

  #gas viscosity
  gasVisc = (gasViscCorr / tprCalc) * (math.exp(gasViscR))



  return gasVisc


# In[5]:
gasVisc(0.65,0.1,0.08,0.02,10000,180)

# In[ ]:




