#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Example 3.2 solution
#get_ipython().run_line_magic('matplotlib', 'inline')
import math
import matplotlib.pyplot as plt


# In[18]:


#ADD Docstring for all functions
import numpy as np


def re(A):
    re = math.sqrt((A * 43560) / math.pi)
    return re


def plotIPR(qo,pwf,regime):
    
    plt.figure(figsize=(9,5))
    ipr = plt.plot(qo,pwf)
    
    #Determine max x and y axis for plot
    ylabel,ytick = plt.yticks()
    xlabel,xtick = plt.xticks()
    yThres = ylabel[5] - ylabel[4]
    ymax = math.ceil(pwf[-1] / yThres) * yThres
    xThres = xlabel[5] - xlabel[4]
    xmax = math.ceil(qo[0] / xThres) * xThres
    
    plt.ylim(0,ymax)
    plt.xlim(0,xmax)
    
    #Plot format 
    plt.xlabel('Production Rate ($q_{o}$), (STB/day)', size=12, style='normal')
    plt.ylabel('Flowing Bottom-Hole Pressure ($p_{wf}$), (psia)', size=12, style='normal')
    plt.title('IPR Curve of {} Flow for Undersaturated Oil'.format(regime), size=18, weight=500)
    plt.grid()
    plt.show()

    
def qoCalc(J,pb,pe):
    pwf = [pb,pe]
    qo1 = J * (pe - pb)
    qo2 = J * (pe - pe)
    qo = [qo1,qo2]
    return qo,pwf

#IPR for undersaturated oil (Single Phase Reservoir)

#Transient flow regime calculations
def transientIPR(k,h,Bo,muO,t,phi,ct,rw,pb,pe):
    a = math.log10(t * 30 * 24) + math.log10(k / (phi * muO * ct * (rw**2))) - 3.23
    J = (k * h) / (162.6 * Bo * muO * a)
    qo,pwf = qoCalc(J,pb,pe)
    plotIPR(qo,pwf,'Transient')
    return J,qo,pwf

#Steady-state flow regime calculations
def steadyStateIPR(k,h,Bo,muO,re,rw,s,pb,pe):
    a = math.log(re/rw) + s
    J = (k * h) / (141.2 * Bo * muO * a)
    qo,pwf = qoCalc(J,pb,pe)
    plotIPR(qo,pwf,'Steady-State')
    return J,qo,pwf

#Pseudo-Steady-state flow regime calculations
def pseudoSteadyStateIPR(k,h,Bo,muO,re,rw,s,pb,pe):
    a = math.log(re/rw) - (3/4) + s
    J = (k * h) / (141.2 * Bo * muO * a)
    qo,pwf = qoCalc(J,pb,pe)
    plotIPR(qo,pwf,'Pseudo-Steady-State')
    return J,qo,pwf


# In[19]:


phi = 0.19
k = 8.2
h = 53
pe = 5651
pb = 50
Bo = 1.1
muO = 1.7
ct = 0.0000129
A = 640
rw = 0.328
s = 0
t = 1 #Number of months

re = re(A)
transientIPR(k,h,Bo,muO,t,phi,ct,rw,pb,pe)

# In[20]:


steadyStateIPR(k,h,Bo,muO,re,rw,s,pb,pe)


# In[21]:


pseudoSteadyStateIPR(k,h,Bo,muO,re,rw,s,pb,pe)

# In[ ]:


