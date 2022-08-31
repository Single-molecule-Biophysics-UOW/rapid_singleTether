# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 17:43:29 2021

@author: stefa
"""

import os
import pandas as pd

from matplotlib import pyplot as plt
import sys

sys.path.append("C:/Users/StefanMueller/Documents/rapid_singleThether/")
from synchronization import rst_tools_v2 as rst

#%%
#load files

folder = r"D:\NAR_upload\Mueller_et_al_NAR\phi29\trajectory/"
files=os.listdir(folder)

#we need 4 files. reaction and prereaction each in two colors
C1_reaction = pd.read_csv(folder+"C1_corrected_reaction_phi29_1uMdNTPs_OD3.6_001.nd2_series_1.csv") 
C1_preReaction = pd.read_csv(folder+"C1_corrected_prereaction_phi29_1uMdNTPs_OD3.6_001.nd2_series_1.csv")
C2_reaction = pd.read_csv(folder+"C2_corrected_reaction_phi29_1uMdNTPs_OD3.6_001.nd2_series_1.csv") 
C2_preReaction = pd.read_csv(folder+"C2_corrected_prereaction_phi29_1uMdNTPs_OD3.6_001.nd2_series_1.csv")

#frame rate in s/frame
timeConversion = 1.0 
#read files and merge color channels.
df_reaction = C1_reaction.merge(C2_reaction, how='outer', on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])


#%%
df_preReaction = pd.merge(C1_preReaction,C2_preReaction,how='outer', on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])

rst.subtract_bg(df_preReaction,df_reaction)


df_reaction= df_reaction.dropna()
    
dfFiltered,sortOut,sortOutReason = rst.find_unbinding_v2(df_reaction,'Intensity_DNA', 0.65, 50000)
dfFiltered["smoothRPA"] = rst.smooth(dfFiltered['Intensity_RPA'])
dfFiltered["smoothDNA"] = rst.smooth(dfFiltered['Intensity_DNA'])

df_synchronyzed = rst.alignData(dfFiltered, timeConversion)

#%%
#the further away from the synchronized point, the less dataPoints potentially.
# it can make sense to truncate around this point:
df_synchronyzed = rst.trunc_trajs(df_synchronyzed,'alignedTime',-100,100)




#%%make synchronized mean and plot
mean_df = pd.DataFrame()
# mean_df['medianRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].median()
mean_df['raw_medianRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].mean()
mean_df['medianDNA'] = df_synchronyzed.groupby('seconds')['Intensity_DNA'].mean()
mean_df['meanRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].mean()
mean_df['meanDNA'] = df_synchronyzed.groupby('seconds')['Intensity_DNA'].mean()
mean_df['semRPA'] = df_synchronyzed.groupby('seconds')['smoothRPA'].sem()

fig, ax = plt.subplots()
# ax.plot(mean_df['medianRPA'],label='median',color='magenta',ls=':')
# ax.plot(mean_df['raw_medianRPA'],label='median',color='magenta')

ax.plot(mean_df['meanDNA'],label='mean',color='black')

ax.plot(mean_df['meanRPA'],label='mean')

ax.fill_between(mean_df.index,
                mean_df['meanRPA']-mean_df['semRPA'],
                mean_df['meanRPA']+mean_df['semRPA'],
                label='mean',alpha=0.2)

ax.plot(mean_df['meanDNA'],label='mean')
ax.legend()
fig.show()
