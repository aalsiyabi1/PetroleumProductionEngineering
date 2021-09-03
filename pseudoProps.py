#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#calculate pseudo-critical pressure
def ppc(gSG,N2,CO2,H2S):
  ppc = 678 - 50 * (gSG-0.5) - 206.7 * N2 + 440 * CO2 + 606.7 * H2S

  return ppc


#calculate pseudo-critical temperature
def tpc(gSG,N2,CO2,H2S):
  tpc = 326 + 315.7 * (gSG-0.5) -240 * N2 - 83.3 * CO2 + 133.3 * H2S

  return tpc


#Pseudo-reduced pressure
def ppr(p,ppc):
  ppr = p / ppc

  return ppr


#Pseudo-reduced temperature
def tpr(t,tpc):
  tR = t + 460 #convert Temp to Rankine
  tpr = tR / tpc 

  return tpr


# In[ ]:




