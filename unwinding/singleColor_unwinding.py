# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:39:59 2021

@author: shm975
"""



import os
import pandas as pd
import rst_tools as rst
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import expon
from scipy.stats import gamma
from scipy.stats import laplace
from scipy.stats import gumbel_r
from scipy.stats import lomax
from scipy.stats import rv_histogram
#%%
files = os.listdir()
print(files)
#%%
all_data=[]
for f in files:
    dataframe = pd.read_csv(f)
    dataframe.name = f
    all_data.append( dataframe)

#%%first lets see if the DNA intensity is reproducable between experiments. To do so look at the intensities at slice 0!
import itertools
slice0_Intensity = pd.DataFrame()

fig,ax = plt.subplots(3,3)
indeces = list(itertools.product([0,1,2],repeat =2))
n=0
data = all_data[0]
for frame in [0,40,80,120,160,300,400,600,]:
    doc_string = data.name.split('_')
    print(doc_string)
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    mean = (data['Intensity'][data['slice']==frame]).mean()
    ax[indeces[n]].hist(data['Intensity'][data['slice'] ==frame],edgecolor='black',label = 'mean: {:.2E}, frame ={}'.format(mean,frame))
    ax[indeces[n]].legend()
    n+=1

plt.show()


#%%
print(data.keys())
#%% now plot the trajectories to have a look:
plt.ioff()
eventFolder = 'C:/Users/shm975/Documents/UvrD_data/3\'bioblunt/events/'
#data.reset_index(inplace=True)

#normalize first:
for data in all_data:
    rst.normalize_all_trajs(data,'Intensity')
    rst.smooth_all_trajectories(data,'norm_Intensity',15,result_prefix='savgol',polyorder = 0,method = 'savgol')

    #rst.smooth_all_trajectories(data,'norm_Intensity',15,result_prefix='window',method = 'window')
#%%save data containing the smoothed trajectories:
print(data.keys())
for data in all_data:
    data.to_csv('smooth_'+data.name)
#%%
segmentFiles = os.listdir('../segments/')
all_segments =[]
for f in segmentFiles:
    dataframe = pd.read_csv('../segments/'+f)
    dataframe.name = f
    all_segments.append( dataframe)


#%%


#rst.plot_all_trajs(data, eventFolder,df_segment=segment_data,segments=True,xcolumn = 'seconds',ycolumn = ['norm_Intensity','savgol_norm_Intensity'])
#%%
n=0
fig,ax = plt.subplots(2,sharex=True)
for segment_data in all_segments:
    print(segment_data.head(10))
    #include last data of each trajectory:
    for name, group in segment_data.groupby('trajectory'):
        print(group)
        group.drop(group.tail(1).index,inplace=True)
        print(group)
    #print(segment_data.head(10))
    #print (n)
    #segment_data['rates'] =( (segment_data['y2']-segment_data['y1'])/segment_data['Duration'])
    #segment_data = segment_data[segment_data['rates']<0]
    #low_threshold = segment_data['rates'].mean()-np.abs(segment_data['rates'].std())
    #segment_data = segment_data[segment_data['rates']>low_threshold]
    #print(segment_data['rates'].head(10))
    #ax[n].hist(segment_data['rates']*2600,edgecolor='black',label='n={}'.format(len(segment_data)),bins = 100)
    #ax[n].legend()
    #n+=1
    
fig.show()
#%%plot as boxplot:
fig,ax = plt.subplots(1,sharex=True)
for segment,i in zip(all_segments,range(len(all_segments))):
    bp = segment.boxplot(column='rates',positions = [i+1],labels=str(i+1) ,grid=False)

    y = segment['rates']
    # Add some random "jitter" to the x-axis
    x = np.random.normal(i, 0.04, size=len(y))
    ax.plot(x, y, 'r.', alpha=0.2)


fig.show()

#std = segment_data['rates'].std()
#threshold = segment_data['rates'].mean()+(std)
#low_threshold = segment_data['rates'].mean()-(std)
#%%
#segment_data = segment_data [segment_data['rates'] < 0] 
#segment_data = segment_data [segment_data['rates'] > low_threshold] 
#segments_trunc = segments_trunc [segments_trunc['rates'] > -50]

#%%

ax.hist(segment_data['rates']*2600,edgecolor='black',label='n={}'.format(len(segment_data)),bins = 100)
ax.legend()
fig.show()