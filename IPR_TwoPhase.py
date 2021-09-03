#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Example 3.3 solution
import math
import matplotlib.pyplot as plt

# In[2]:


#ADD Docstring for all functions


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
    plt.title('IPR Curve of {} Flow for Saturated Oil'.format(regime), size=18, weight=500)
    plt.grid()
    plt.show()



#IPR for saturated oil (Two-Phase Reservoir) using Vogel's Equation
def twoPhaseIPR(k,h,Bo,muO,re,rw,s,pb,pe):
    a = math.log(re/rw) - (3/4) + s
    J = (k * h) / (141.2 * Bo * muO * a)
    qmax = (J * pe) / 1.8
    qo = []
    pwf = []
    for p in range(pe):
        q = qmax * (1 - 0.2 * (p/pe) - 0.8 * (p/pe)**2)
        qo.append(q)
        pwf.append(p)
    plotIPR(qo,pwf,'Two Phase')
    return qo,pwf


#IPR for saturated oil (Two-Phase Reservoir) using Vogel's Equation
def partialTwoPhaseIPR(k,h,Bo,muO,re,rw,s,pb,pe):
    a = math.log(re/rw) - (3/4) + s
    J = (k * h) / (141.2 * Bo * muO * a)
    qv = (J * pb) / 1.8
    qb = J * (pe - pb)
    qo = []
    pwf = []
    for p in range(pb):
        q = qb + qv * (1 - 0.2 * (p/pb) - 0.8 * (p/pb)**2)
        qo.append(q)
        pwf.append(p)
    qo.append(0)
    pwf.append(pe)
    plotIPR(qo,pwf,'Partial Two Phase')
    return qo,pwf


# In[3]:
phi = 0.19
k = 8.2
h = 53
pe = 5651
pb = 3000
Bo = 1.1
muO = 1.7
ct = 0.0000129
A = 640
rw = 0.328
s = 0
t = 1 #Number of months

Re = re(A)

qo,pwf = partialTwoPhaseIPR(k,h,Bo,muO,Re,rw,s,pb,pe)


