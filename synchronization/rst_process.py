# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 17:43:29 2021

@author: stefa
"""

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

#load files

folder = 'C:/Users/shm975/Documents/tempData/211005/phi29/integrations/exp1_10uM/'
files=os.listdir(folder)

preReactions = [f for f in files if 'prereaction' in f]
reactions = [f for f in files if '_reaction' in f]

print(preReactions)
print(reactions)

#read files and merge color channels.

C1_reaction = pd.read_csv(folder+reactions[0]) 
C1_preReaction = pd.read_csv(folder+preReactions[0])
C2_reaction = pd.read_csv(folder+reactions[1]) 
C2_preReaction = pd.read_csv(folder+preReactions[1])

df_reaction = pd.merge(C1_reaction,C2_reaction, on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])
df_preReaction = pd.merge(C1_preReaction,C2_preReaction, on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])

def calibrate(merged,conv):
    """
    add column with seconds
    conv is seconds/frame
    """
    merged['seconds'] = merged['slice'].mul(conv)   
    return merged
#%%
def subtract_bg(preReaction,reaction):    
    #now subtract the bg from the RPA intensity:
    #group background by trajectory. Here we make the assumption that the trajectory numbering is consistent!
    groupedPreReaction = preReaction.groupby('trajectory')
    groupedReaction = reaction.groupby('trajectory')
    
    corr_Intensity_RPA = []
    for name,group in groupedPreReaction:
        #get baseline:
        baseline = group['Intensity_RPA'].mean()
        corr_Intensity_RPA.append(groupedReaction.get_group(name)['Intensity_RPA'] - baseline)
    corrIntColumn = pd.concat(corr_Intensity_RPA)
    reaction['corr_Intensity_RPA'] = corrIntColumn
    
calibrate(df_reaction,2)
calibrate(df_preReaction,2)
import rst_tools as rst
#st.plot_allBig(df_reaction,'C:/Users/shm975/Documents/tempData/211005/phi29/all_events/')
#subtract_bg(df_preReaction,df_reaction)
#%%
def rolling_avg(df):
    df["smoothRPA"] = df['Intensity_RPA'].rolling(3,min_periods=1).mean()
    df["smoothDNA"] = df['Intensity_DNA'].rolling(3,min_periods=1).mean()
df_reaction.reset_index(inplace=True)
dfFiltered,sortOutDNA = rst.findEvents_v2(df_reaction,column = 'Intensity_DNA',Rscore=0.65, return_events=True)
#eventList = [3,9,15,17,48,61,62,68,71,77,89,91,106,111,115,117,121,122,123,129,132,137,148,157,158,161,164,167,177,
#             195,214,215,216,218,220,240,243,244,246,248,252,257,286,296,299,303,322,325,326,330,336,342,344,346,349,350,355,357,366,367,375,
#             382,386,387,390,399,404,408,418,422,426,427,437,442,447,451,455,458,464,474,475,479,481,485,492,501,503]
#dfFiltered = rst.selectEvents(df_reaction, eventList)

rolling_avg(dfFiltered)




#%%plot some trajectories
import itertools

def plot_all_trajs(df,out):
    fig,ax = plt.subplots(3,3)
    indeces = list(itertools.product([0,1,2],repeat =2))
    n=0
    for name,group in df.groupby("trajectory"):                
        rpa = ((group["Intensity_RPA"]/group["Intensity_RPA"].max())*group["Intensity_DNA"].max()).rolling(4, min_periods=1).mean()
        ax[indeces[n]].plot(group["seconds"],rpa,color='magenta')
        ax[indeces[n]].plot(group["seconds"],group["Intensity_DNA"],color='black')
        ax[indeces[n]].plot(group["seconds"],group["regression"],color='red')
        #ax[indeces[n]].plot(group["seconds"],-group["derivative"],color='blue')
        ax[indeces[n]].set_title('trajectory {}'.format(name))
        #find the corresponding segments:
        #segment = df_segment.loc[df_segment['trajectory'] == name]
        #segmentPlotter(segment,ax[indeces[n]])
        n+=1
        if n == 9:
            #save figure and close it:
            fig.savefig(out+'trajectory_'+str(name))
            plt.close(fig)
            #make new one:
            fig,ax = plt.subplots(3,3)
            n= 0

rst.plot_allBig(dfFiltered, 'C:/Users/shm975/Documents/tempData/211005/phi29/events_10uM/')

#%%
def alignTime(df):
    regr = np.array(df)
    oldTime = df.index.get_level_values(1)
    changePoint = np.where(np.gradient(regr)!=0)[0][0]     
    newTime = oldTime -changePoint
    return newTime


def calcShift(df):
    """calculate the shift for aligment as in alignTime, but only return the shift"""
    regr = np.array(df)
    oldTime = df.index.get_level_values(1)
    changePoint = np.where(np.gradient(regr)!=0)[0][0]     
    return changePoint
def alignData(df,timeConversion,traj_column = 'trajectory',slice_column = 'slice',returnShift = False):
    df = df.set_index([traj_column,slice_column]).sort_index()
    df['alignedTime'] = df['regression'].groupby(traj_column).transform(alignTime)
    if returnShift:
        cpList = []
        for name,groups in df['regression'].groupby(traj_column):
            changePoint = calcShift(groups)
            cpList.append(changePoint)
#convert time from frames to seconds:
    df['seconds'] = df['alignedTime'].mul(timeConversion)
    df['raw time'] = df.index.get_level_values(1)*(timeConversion)
    if returnShift:
        return df,cpList
    return df

dfFiltered.reset_index(inplace=True)
rst.regression_analysis(dfFiltered)
#%%
dfFiltered.reset_index(inplace=True)
print(dfFiltered.keys())
#%%



#%%
df_synchronyzed = alignData(dfFiltered,1.)

#df_synchronyzed.reset_index(inplace=True)
def rolling_avg(df):
    df["smoothRPA"] = df['Intensity_RPA'].rolling(3,min_periods=1).mean()
    df["smoothDNA"] = df['Intensity_DNA'].rolling(3,min_periods=1).mean()
def trunc_trajs(df,xcolumn,lower,upper):
    truncatedLower = df.groupby('trajectory').apply(
        lambda x: x.loc[x[xcolumn] >lower ]
        )
    truncatedUpper = truncatedLower.groupby('trajectory').apply(
        lambda x: x.loc[x[xcolumn] <upper ]
        )
    return truncatedUpper

    


    

#rolling_avg(df_synchronyzed)    
df_synchronyzed = trunc_trajs(df_synchronyzed,'alignedTime',-100,100)

df_synchronyzed = rst.normalize_all_trajs(df_synchronyzed,'smoothRPA')

df_synchronyzed = trunc_trajs(df_synchronyzed,'alignedTime',-80,100)

df_synchronyzed.to_csv('C:/Users/shm975/Documents/tempData/211005/phi29/selectedEvents10uM.csv')


#%%

plot_all_trajs(df_synchronyzed, 'C:/Users/shm975/Documents/tempData/211005/phi29/man_events/')

print(df_reaction['seconds'])

#%%make mean
mean_df = pd.DataFrame()
mean_df['medianRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].median()
mean_df['raw_medianRPA'] = df_synchronyzed.groupby('seconds')['Intensity_RPA'].median()
mean_df['medianDNA'] = df_synchronyzed.groupby('seconds')['Intensity_DNA'].median()
mean_df['meanRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].mean()
mean_df['meanDNA'] = df_synchronyzed.groupby('seconds')['Intensity_DNA'].mean()
mean_df['semRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].std()

#%%plot it
fig, ax = plt.subplots()
ax.plot(mean_df['medianRPA'],label='median',color='magenta',ls=':')
ax.plot(mean_df['raw_medianRPA'],label='median',color='magenta')

ax.plot(mean_df['medianDNA'],label='median',color='black')

ax.plot(mean_df['meanRPA'],label='mean')

ax.fill_between(mean_df.index,mean_df['meanRPA']-mean_df['semRPA'],mean_df['meanRPA']+mean_df['semRPA'],label='mean',alpha=0.2)

ax.plot(mean_df['meanDNA'],label='mean')
ax.legend()
fig.show()

#%% now the whole thing again with start synch:
import start_synch as s

df_synchronyzed.reset_index(inplace=True)
synchRPA = s.start_synch(df_synchronyzed,column = 'Intensity_RPA')

#%%make mean of start synchronized data
mean_df = pd.DataFrame()
mean_df['medianRPA'] = synchRPA.groupby('seconds_new')['smoothRPA'].median()
mean_df['raw_medianRPA'] = synchRPA.groupby('seconds_new')['Intensity_RPA'].median()
mean_df['medianDNA'] = synchRPA.groupby('seconds_new')['Intensity_DNA'].median()
mean_df['meanRPA'] = synchRPA.groupby('seconds_new')['smoothRPA'].mean()
mean_df['meanDNA'] = synchRPA.groupby('seconds_new')['Intensity_DNA'].mean()
mean_df['semRPA'] = synchRPA.groupby('seconds_new')['smoothRPA'].std()

#%%
print(mean_df)
mean_df.reset_index(inplace=True)
print(mean_df['seconds_new'])
#%%
mean_trunc = mean_df#[mean_df['seconds']>-200. ]
#mean_trunc = mean_trunc[mean_trunc['seconds']<200. ]


#%%plot it
fig, ax = plt.subplots()

ax.plot(mean_trunc['seconds_new'], mean_trunc ['medianRPA']/mean_trunc ['medianRPA'].max(),label='median',color='magenta',ls=':')
ax.plot(mean_trunc['seconds_new'],mean_trunc ['raw_medianRPA']/mean_trunc ['raw_medianRPA'].max(),label='median',color='magenta')

ax.plot(mean_trunc['seconds_new'],mean_trunc ['medianDNA']/mean_trunc ['medianDNA'].max(),label='median',color='black')

ax.plot(mean_trunc['seconds_new'],mean_trunc ['meanRPA']/mean_trunc ['meanRPA'].max(),label='mean')

normMean = mean_trunc ['meanRPA']/mean_trunc ['meanRPA'].max()
ax.fill_between(mean_trunc['seconds_new'],normMean-mean_trunc ['semRPA']/mean_trunc ['meanRPA'].max(),
                normMean+mean_trunc ['semRPA']/mean_trunc ['meanRPA'].max(),label='mean',alpha=0.2)

ax.plot(mean_trunc['seconds_new'],mean_trunc ['meanDNA']/mean_trunc ['meanDNA'].max(),label='mean')
ax.legend()
fig.show()