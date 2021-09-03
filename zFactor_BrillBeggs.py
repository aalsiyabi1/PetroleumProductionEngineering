#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Example 2.4 (estimate z-factor)

import math
from pseudoProps import *


# In[14]:


#Add docstring
def zFactor(p,t,gSG,N2,CO2,H2S):
    
    pseudoP = ppc(gSG=gSG,N2=N2,H2S=H2S,CO2=CO2)
    pseudoT = tpc(gSG=gSG,N2=N2,H2S=H2S,CO2=CO2)

    pseudoReducedP = ppr(p=p,ppc=pseudoP)
    pseudoReducedT = tpr(t=t,tpc=pseudoT)
    
    F = 0.3106 - 0.49 * pseudoReducedT + 0.1824 * (pseudoReducedT ** 2)
    E = 9 * (pseudoReducedT - 1)
    D = 10 ** F
    C = 0.132 - 0.32 * math.log10(pseudoReducedT)
    B = (0.62 - 0.23 * pseudoReducedT) * pseudoReducedP + ((0.066/(pseudoReducedT-0.86)) - 0.037) * (pseudoReducedP ** 2) + ((0.32 * (pseudoReducedP ** 6)) / (10 ** E))
    A = 1.39 * (pseudoReducedT - 0.92) ** 0.5 - 0.36 * pseudoReducedT - 0.10
    
    z = A + ((1 - A)/(math.exp(B))) + C * (pseudoReducedP ** D)



    print(pseudoP,pseudoT,pseudoReducedP,pseudoReducedT,A,B,C,D,z)
    return z

# In[16]:


p = 5000
t = 180
gSG = 0.65
N2 = 0.1
CO2 = 0.08
H2S = 0.02


zFactor(p,t,gSG,N2,CO2,H2S)


# In[ ]:




