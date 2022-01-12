# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 15:09:06 2022

@author: stefa
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from synchronization import rst_tools_v2 as rst
import os


folder = r"C:\Users\stefa\OneDrive\Documents\tempData\UvrdData\UvrD_data\template_comparisson/"

files = os.listdir(folder)

print(files)
data = []
for i in files:
    df = pd.read_csv(folder+i)
    df.name = i
    data.append(df)
    
#%% now find all molecules that dissociate in one step:
events_all_data = []
N = len(data)
n = 0
for df in data:
    print(len(df['trajectory'].unique()))
    events,sortOut,sortOutReason = rst.find_unbinding_v2(df, 'Intensity', 0.7, 60000)
    events.name = df.name
    events_all_data.append(events)
    print('finished {} out of {}. Found {} events'.format(n+1,N, len(events['trajectory'].unique())))
    n+=1


#%%first lets see if the DNA intensity is reproducable between experiments. To do so look at the intensities at slice 0!
import itertools
slice0_Intensity = pd.DataFrame()
fig,ax = plt.subplots(3,3,sharex=True)
indeces = list(itertools.product([0,1,2],repeat =2))
print(len(indeces))
n=0
for df in data:
    
    doc_string = df.name.split('_')
    print(doc_string)
    conc = "".join(filter(str.isdigit, doc_string[4]))
    #RPA =  "".join(filter(str.isdigit, doc_string[4]))
    #OD =   "".join(filter(str.isdigit, doc_string[6]))
    preinc = "".join(filter(str.isdigit, doc_string[3]))
    #DATE = "".join(filter(str.isdigit, doc_string[7]))
    mean = (df['Intensity'][df['slice']==0]).mean()
    
    ax[indeces[n]].hist(df['Intensity'][df['slice'] ==3],edgecolor='black',label = 'conc:'+conc+'mean: {:.2E}'.format(mean))
    ax[indeces[n]].legend()
    n+=1
    if n ==9:
        break
plt.show()

#%% and now get the time off dissociation for every molecule:
timing_all_data = []
n=0
for events in events_all_data:
    timing = rst.unbinding_timing(events,groupName = 'trajectory',column = 'seconds')
    timing.name = events.name
    print(timing.name)
    timing_all_data.append(timing)    
    print('finished {} out of {}'.format(n+1,N))    
    n+=1

#%%plot number of molecules visible:

#loop through all timepoints and count how many disTimes are bigger than this time point:
N_molecules_all_data = []

count = 0
for events,timing in zip(events_all_data,timing_all_data):
    
    
    timepoints = events['seconds'].unique()
    
    timepoints =timepoints[timepoints>=0]
    N_molecules = []
    count+=1
    for t in timepoints:
        molecules = timing[timing['seconds']>t]
        N_molecules.append(len(molecules))
    
    N_molecules = pd.DataFrame({'N': N_molecules, 'time': timepoints})
    N_molecules.name = events.name
    print(N_molecules.name)
    N_molecules_all_data.append(N_molecules)
    
#%%
for i in N_molecules_all_data:
    print(i.name)

#%%

n=0
fig,ax = plt.subplots(figsize = (4,8))
done=False
N0 = [1144,755,1273]
conc = ["'5'blunt","'5'ssDNA","3'bioblunt"]
styles = ['-', '--', '-.']

print(len(N_molecules_all_data))

for data,conc in zip(N_molecules_all_data,conc):
    print(data.name)
    
    #if data.name not in ['new_corr_UvrD100nM_nopreinc_5\'BiossDNA_ATP50uM_OD3.3_001.csv','new_corr_UvrD100nM_nopreinc_5\'BiossDNA_ATP0uM_OD3.3_009.csv']:
    #    continue

    #N0 = len(all_data_dsDNA[n]['trajectory'].unique())
    #print("Nmax:",N0)
    #ax.plot(N_molecules_all_ssDNA[n]['time'],N_molecules_all_ssDNA[n]['N']/N_molecules_all_ssDNA[n].max(),label=conc)
    #if conc in ['50','10']:
     #   if done:
      #      pass
    N = (data['N'].max())
    N0_real = N0[n]
    #dataN = N0 - data['N'].diff().fillna(0).abs().cumsum()
    ax.plot(data['time'],(data['N']+(N0_real-N))/(N0_real),label = 'N={}, template:{}'.format(N0_real,conc),color = 'black',ls = styles[n])
    ax.set_xlim([0,120])
    #ax.plot(data['time'],(data['N']+N0[n])/N0[n],label = 'N=691,file:{} '.format(data.name))
    ax.legend()
    n+=1

fig.savefig('unbinding_templateComp.svg')

#%%
fig,ax = plt.subplots()
for df in N_molecules_all_data:
    ax.plot(df['time'],df['N'])
plt.show()

