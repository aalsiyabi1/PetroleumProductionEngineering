#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Excercise 2.2

import pandas as pd 
import math


# In[ ]:


gas_comp = {'Compound':['C1','C2','C3','i-C4','n-C4','i-C5','n-C5','C6','C7+','N2','CO2','H2S'],
            'MWi':[16.04,30.07,44.10,58.12,58.12,72.15,72.15,86.18,114.23,28.02,44.01,34.08],
            'Pci':[673,709,618,530,551,482,485,434,361,227,1073,672],
            'Tci':[344,550,666,733,766,830,847,915,1024,492,548,1306]}


# In[ ]:


for i,feature in enumerate(gas_comp['Compound']):
    print('{}, MWi: {}, Pci: {}, Tci: {}'.format(feature,gas_comp['MWi'][i],gas_comp['Pci'][i],gas_comp['Tci'][i]))


# In[ ]:


gas_comp_calc = gas_comp
gas_comp_calc['yi'] = [0.775,0.083,0.021,0.006,0.002,0.003,0.008,0.001,0.001,0.050,0.030,0.020]

y = gas_comp_calc['yi']
MW = gas_comp_calc['MWi']
Pc = gas_comp_calc['Pci']
Tc = gas_comp_calc['Tci']
gas_comp_calc['yiMWi'] = []
gas_comp_calc['yiPci'] = []
gas_comp_calc['yiTci'] = []


for yi,MWi in zip(y,MW):
  yiMWi = yi * MWi
  gas_comp_calc['yiMWi'].append(yiMWi)

for yi,Pci in zip(y,Pc):
  yiPci = yi * Pci
  gas_comp_calc['yiPci'].append(yiPci)

for yi,Tci in zip(y,Tc):
  yiTci = yi * Tci
  gas_comp_calc['yiTci'].append(yiTci)

print('{:<15} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}'.format('Compound','MWi','Pci','Tci','yi','yiMWi','yiPci','yiTci'))

for i,feature in enumerate(gas_comp_calc['Compound']):
  print('{:<15} {:<10.2f} {:<10} {:<10} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f}'
        .format(feature,gas_comp['MWi'][i],gas_comp['Pci'][i],gas_comp['Tci'][i],gas_comp['yi'][i],gas_comp['yiMWi'][i],gas_comp['yiPci'][i],gas_comp['yiTci'][i]))


sum_yi = sum(gas_comp_calc['yi'])
MWa = sum(gas_comp_calc['yiMWi']) #Apparent Molecular Weight
Ppr = sum(gas_comp_calc['yiPci']) #Pseudo-Critical Pressure
Tpc = sum(gas_comp_calc['yiTci']) #Pseudo-Critical Temperature

gama_g = MWa / 28.97 #Specific Gravity

print('Total yi: {}'.format(sum_yi))
print('\n')
print('Gas Apparent Molecular Weight: {:<10.3f}'.format(MWa))
print('Gas Pseudo-Critical Pressure: {:<10.3f}'.format(Ppr))
print('Gas Pseudo-Critical Temperature: {:<10.3f}'.format(Tpc))
print('Gas Specific Gravity: {:<10.3f}'.format(gama_g))


# In[ ]:




