#-*- coding: utf-8 -*-
"""
Created on Tue Nov  9 16:03:08 2021

@author: Stefan H. Mueller
"""

import pandas as pd
from synchronization import rst_tools_v2 as rst
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gamma



#%% read the data

#change this to the filepath and filename containing integrated rois
dataFolder = r"C:\Users\smueller\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unbinding/integration_used/"
fileName = "new_corr_UvrD100nM_5'biodsDNA_ATP12.5uM_OD3.3_001.nd2.csv"

dataframe = pd.read_csv(dataFolder+fileName)
dataframe.name = fileName


#%% now find all molecules that dissociate in one step:
events_all_data_ssDNA = []

events,sortOut,sortOutReason = rst.find_unbinding_v2(dataframe, 'Intensity', 0.7, 60000)
events.name = dataframe.name

#%% and now get the time off dissociation for every molecule:

timing = rst.unbinding_timing(events,groupName = 'trajectory',column = 'seconds')
timing.name = events.name


#%% plot a histogram and fit with gamma distribution

#the time to the first disapearance event is governed by the flow kinetics.
#It's resonable  to set this time to zero for the gamma fit with floc=
gammaFit = gamma.fit(timing['seconds'],floc=timing['seconds'].min()-0.01)
x=np.arange(timing['seconds'].min(),timing['seconds'].max(),0.1)
plt.figure()
plt.hist(timing, edgecolor = 'black',density=True)
plt.plot(x, gamma.pdf(x,*gammaFit))
plt.show()


#%%plot number of molecules visible:

timepoints = events['seconds'].unique()
timepoints =timepoints[timepoints>=0]
N_molecules = []
#loop through all timepoints and count how many disTimes are bigger than this time point:
for t in timepoints:
    molecules = timing[timing['seconds']>t]
    N_molecules.append(len(molecules))
    print(len(timepoints))
    print(len(N_molecules))
N_molecules = pd.DataFrame({'N': N_molecules, 'time': timepoints})
N_molecules.name = events.name

plt.figure()
plt.plot(N_molecules)
plt.show()    


