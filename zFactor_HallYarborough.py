#!/usr/bin/env python
# coding: utf-8

# In[55]:


#Add docstring

import numpy as np
import math
from scipy import optimize 
from pseudoProps import *


#More accurate z-factor value than Brill-Begs
def z_factor(p,t,gSG,N2,CO2,H2S):    
    
    pseudoP = ppc(gSG=gSG,N2=N2,H2S=H2S,CO2=CO2)
    pseudoT = tpc(gSG=gSG,N2=N2,H2S=H2S,CO2=CO2)
    
    pseudoReducedP = ppr(p=p,ppc=pseudoP)
    pseudoReducedT = tpr(t=t,tpc=pseudoT)
    
    tr = 1 / pseudoReducedT
    
    A = 0.06125 * tr * (math.exp(-1.2 * (1 - tr)**2))
    B = tr * (14.76 - 9.76 * tr + 4.58 * tr**2 )
    C = tr * (90.7 - 242.2 * tr + 42.4 * tr **2)
    D = 2.18 + 2.82 * tr
    
    fY = lambda Y: ((Y + Y**2 + Y**3 - Y**4) / (1 - Y )**3) - A * pseudoReducedP - B * Y**2 + C * Y**D
    Y_soln = optimize.fsolve(fY,0)
    fY_soln = fY(Y_soln)
    z = A * pseudoReducedP / Y_soln
    
    
    print('A: {}, B: {}, C: {}, D: {}, Y (reduced density): {}, Root: {}, z-factor: {}'.format(A,B,C,D,Y_soln,fY_soln,z))
    return z


def gasDensity(gSG,p,z,t):
    
    tR = t + 460
    rhoG = (2.7 * gSG * p) / (z * tR)
    return rhoG


# In[57]:


t = 200 
p = 2000
gSG = 0.7
N2 = 0.05
CO2 = 0.05
H2S = 0.02

z = z_factor(p,t,gSG,N2,CO2,H2S)
gDensity = gasDensity(gSG,p,z,t)

print(z,gDensity)


# In[ ]:




