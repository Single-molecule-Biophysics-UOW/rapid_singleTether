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

folder = 'C:/Users/stefa/OneDrive/Documents/tempData/210922/integrations/'
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
subtract_bg(df_preReaction,df_reaction)
#%%
import STA as sta
dfFiltered,sortOutDNA = sta.findEvents_v2(df_reaction,column = 'corr_Intensity_RPA',Rscore=0.7, return_events=True)
#%%plot some trajectories
import itertools

def plot_all_trajs(df,out):
    fig,ax = plt.subplots(3,3)
    indeces = list(itertools.product([0,1,2],repeat =2))
    n=0
    for name,group in df.groupby("trajectory"):                 
        ax[indeces[n]].plot(group["seconds"],group["corr_Intensity_RPA"]/group["corr_Intensity_RPA"].max(),color='magenta')
        ax[indeces[n]].plot(group["seconds"],group["Intensity_DNA"]/group["Intensity_DNA"].max(),color='black')
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
    
plot_all_trajs(df_reaction, 'C:/Users/stefa/OneDrive/Documents/tempData/210922/events/')

#%%
print(df_reaction['seconds'])
