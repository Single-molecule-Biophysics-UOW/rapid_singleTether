# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 18:13:56 2021

@author: shm975
"""


from matplotlib import pyplot as plt
import pandas as pd



UvrD_data = pd.read_csv('C:/Users/shm975/Documents/UvrD_data/3\'bioblunt/integration/UvrD_synchMean.csv')
Lambda_data = pd.read_csv(r'C:\Users\shm975\Documents\tempData\211201\lambda_synchMean.csv')
phi29_data = pd.read_csv(r'C:\Users\shm975\Documents\phi29Data\averageRPAData25uM.csv')
phi29DNA_data = pd.read_csv(r'C:\Users\shm975\Documents\phi29Data\averageDNAData25uM.csv')

#%%
print(UvrD_data)

#%%
print(phi29_data.keys())
fig, ax = plt.subplots()



UvrD_zero_int = (UvrD_data['meanDNA'][(UvrD_data['seconds']>0) & (UvrD_data['seconds']<50)]).mean()
Lambda_zero_int = (Lambda_data['meanDNA'][(Lambda_data['seconds']>0) & (Lambda_data['seconds']<50)]).mean()
phi29_zero_int  = (phi29_data['mean_int'][(phi29_data['seconds']>-800) & (phi29_data['seconds']<-400)]).mean()

UvrD_data['norm'] = (UvrD_data['meanDNA']-UvrD_zero_int)/(UvrD_data['meanDNA']-UvrD_zero_int).max()
Lambda_data['norm'] = (Lambda_data['meanDNA']-Lambda_zero_int)/(Lambda_data['meanDNA']-Lambda_zero_int).max()
phi29_data['norm'] = (phi29_data['mean_int']-phi29_zero_int)/(phi29_data['mean_int']-phi29_zero_int).max()
#print(UvrD_data['norm'])
ax.plot(UvrD_data['seconds'],(1-UvrD_data['norm'])*2600,label = 'UvrD unwinding')
ax.plot(Lambda_data['seconds'],(1-Lambda_data['norm'])*2600,label = '$\lambda$ exo degradation')
ax.plot(phi29_data['seconds'],phi29_data['norm']*2600,label = 'phi29 synthesis')
ax.legend()
ax.set_xlim([-355,55])
ax.set_ylim([-0.01,2650])
fig.show()
#fig.savefig('C:/Users/shm975/Documents/OneDrive - University of Wollongong/phi29Paper/figures/allEnzymes.svg')

#%%
fig, ax = plt.subplots(2,sharex=True)

phi29DNA_zero_int  = (phi29DNA_data['mean_int'][(phi29DNA_data['seconds']>0) & (phi29_data['seconds']<400)]).mean()
print(phi29DNA_zero_int )
phi29DNA_data['norm'] = (phi29DNA_data['mean_int']-phi29DNA_zero_int)/(phi29DNA_data['mean_int']-phi29DNA_zero_int).max()


import numpy as np
diff_kernel = np.array([1,1,1,1,1,-1,-1,-1,-1,-1])
der = np.convolve(phi29_data['norm'],diff_kernel ,'same')

ax[0].plot(phi29_data['seconds'],phi29_data['norm'],label = 'RPA')
ax[0].plot(phi29DNA_data['seconds'],phi29DNA_data['norm'],label = 'DNA')
ax[1].plot(phi29_data['seconds'],der,label = 'RPA derivative')
fig.show()