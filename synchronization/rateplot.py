# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 12:27:14 2021

@author: shm975
"""

import pandas as pd
from matplotlib import pyplot as plt
import itertools

folder = 'C:/Users/shm975/Documents/tempData/211005/phi29/'

events = pd.read_csv(folder + 'selectedEvents10uM.csv')
segments = pd.read_csv(folder + 'selectedSegments10uM_.01.csv')

plt.ioff()

def segmentPlotter(seg_data,ax):
    #steps = True:
    for index, row in seg_data.iterrows():
        ax.plot([row["x1"],row["x2"]],[row["y1"],row["y2"]],color='red')
    
def plot_all_trajs_big(df_integration,df_segment,xcolumn,ycolumn,out,segments=True): 
    for name,group in df_integration.groupby("trajectory"):                 
        fig,ax = plt.subplots()
        ax.plot(group[xcolumn],group[ycolumn])
        ax.set_title('trajectory {}'.format(name))
        if segments:
            #find the corresponding segments:
            segment = df_segment.loc[df_segment['trajectory'] == name]
            segmentPlotter(segment,ax)
        fig.savefig(out+'trajectory_'+str(name))
        plt.close(fig)
    
def plot_all_trajs(df_integration,df_segment,xcolumn,ycolumn, trajs_per_figure,out,segments=True):
    fig,ax = plt.subplots(3,3)
    ax_positions = itertools.permutations((0,1,2),2)
    indeces = list(itertools.product([0,1,2],repeat =2))
    n=0
    for name,group in df_integration.groupby("trajectory"):                 
        ax[indeces[n]].plot(group[xcolumn],group[ycolumn])
        ax[indeces[n]].set_title('trajectory {}'.format(name))
        if segments:
            #find the corresponding segments:
            segment = df_segment.loc[df_segment['trajectory'] == name]
            segmentPlotter(segment,ax[indeces[n]])
        n+=1
        if n == 9:
            #save figure and close it:
            fig.savefig(out+'trajectory_'+str(name))
            plt.close(fig)
            #make new one:
            fig,ax = plt.subplots(3,3)
            n= 0
#%%
plot_all_trajs_big(events,segments,'alignedTime','norm_smoothRPA',folder+'changePoint10uM/')  
#%%plot distributions



segments['rates'] =( (segments['y2']-segments['y1'])/segments['Duration'])*2600

#%%
segments_trunc = segments [segments['rates'] < 200] 
segments_trunc = segments_trunc [segments_trunc['rates'] > -50]

#%%
fig,ax = plt.subplots()
ax.hist(segments_trunc['rates'],bins=28,weights = 1/(segments_trunc['sigma_A']+segments_trunc['sigma_B']),edgecolor='black',label='n={}'.format(len(segments_trunc)))
ax.legend()
fig.show()
          