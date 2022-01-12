# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 11:35:11 2021

@author: stefa
"""

#lambda exo synchronization in sequence space

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
#%% first load already normalized dataset
dataFolder = r"C:\Users\stefa\OneDrive - University of Wollongong\phi29Paper\data\lambda\211221\norm_integration/"

files = os.listdir(dataFolder)
dataframe = pd.read_csv(dataFolder+files[0])
dataframe.name = files[0]
print(dataframe.name)
print(dataframe.keys())
#%%
dataframe['rounded'] = (dataframe['nts'].round())
grouped_df = dataframe.groupby('trajectory')

traj = grouped_df.get_group(1)

plt.figure()
plt.plot(traj['seconds'],traj['nts'])
plt.plot(traj['seconds'],traj['piecewise_triple'])
plt.plot(traj['seconds'],traj['rounded'])
plt.show()

#%%




print(dataframe['rounded'])

#%%
nt_synch = pd.DataFrame()
nt_synch['seconds'] = dataframe.groupby('rounded')['seconds'].mean()
nt_synch['smooth_seconds'] = nt_synch['seconds'].rolling(500).mean()





fig,ax = plt.subplots(2,sharex=True)
ax[0].plot(nt_synch.index,nt_synch['seconds'])
ax[0].plot(nt_synch.index,nt_synch['smooth_seconds'])
ax[0].set_xlim([-20,2700])
#ax[1].plot(nt_synch.index,nt_synch['seconds'].diff()/ nt_synch.index.to_series().diff())
ax[1].plot(nt_synch.index,nt_synch['smooth_seconds'].diff(periods=10)/ nt_synch.index.to_series().diff(periods=10))
#ax[1].plot(nt_synch.index,nt_synch['smooth_seconds'].diff(periods = 10))
plt.show()




