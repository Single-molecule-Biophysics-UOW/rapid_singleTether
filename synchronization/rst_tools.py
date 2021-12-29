# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 14:04:54 2021

@author: shm975
"""


import itertools
from matplotlib import pyplot as plt
import numpy as np
import STA as sta
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
import pandas as pd

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
def alignData(df,timeConversion,traj_column = 'trajectory',slice_column = 'slice',changePoint_column = 'regression',returnShift = False):
    df = df.set_index([traj_column,slice_column]).sort_index()
    df['alignedTime'] = df[changePoint_column].groupby(traj_column).transform(alignTime)
    if returnShift:
        cpList = []
        for name,groups in df[changePoint_column].groupby(traj_column):
            changePoint = calcShift(groups)
            cpList.append(changePoint)
#convert time from frames to seconds:
    df['seconds'] = df['alignedTime'].mul(timeConversion)
    df['raw time'] = df.index.get_level_values(1)*(timeConversion)
    if returnShift:
        return df,cpList
    return df
def regression(df,maxDepth=1,output = 'predict'):

    x,y = np.array(df.index.get_level_values(1)).reshape(-1,1),np.array(df).reshape(-1,1)
    
    regr = DecisionTreeRegressor(max_depth=maxDepth)
    fit = regr.fit(x,y)
    score = fit.score(x, y)  
    fity = fit.predict(x)
    #d = {'fity': fity, 'score': np.ones(fity.shape)*score}
    #regrdf = pd.DataFrame(data=d)
        
    if output == 'predict':
        return fity
    elif output == 'score':
        return score
    return score
def linearRegression(df, output ='Score'):
    x,y = np.array(df.index.get_level_values(1)).reshape(-1,1),np.array(df).reshape(-1,1)
    
    regr = LinearRegression()
    fit = regr.fit(x,y)
    score = fit.score(x, y)  
    fity = fit.predict(x)
    #d = {'fity': fity, 'score': np.ones(fity.shape)*score}
    #regrdf = pd.DataFrame(data=d)   
    if output == 'predict':
        return fity
    elif output == 'score':
        return score
    return score

def selectEvents(df,eventList):
    """select events manually
    """       
    df.set_index(['trajectory','slice'],inplace=True)
    
    sortOut = []
    for names,groups in df.groupby('trajectory'):        
        if names not in eventList:                  
            sortOut.append(names)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    final.reset_index(inplace=True)
    return final

def regression_analysis(df, column = 'Intensity_DNA'):
    df.set_index(['trajectory','slice'],inplace=True)
    df.sort_index(inplace=True)
    df['regression'] = df[column].groupby('trajectory').transform(regression,output='predict')
    df['score'] = df[column].groupby('trajectory').transform(regression,output='score')
    df['score_null'] = df[column].groupby('trajectory').transform(linearRegression,output='score')
    
    
def findEvents_v2(df,Rscore,return_events=True, column = 'Intensity_DNA'):
    """find trajectories where DNA is clearly unbinding
    by using a regression analysis and subsequent threshold.
    V2: compare this to another regression with max_depth=0, meaning a single line fit
    # Two lines have to fit better than the null-hypothesis. Idealy this should be a 
    # maximum likelyhood test, see changePoint! MAybe even just use changePoint!
    Rscore is the combined R**2 of linear regressions
    if return_events=False, the trajectories not showing this will be returned
    """       
    df.set_index(['trajectory','slice'],inplace=True)
    df.sort_index(inplace=True)

    df['regression'] = df[column].groupby('trajectory').transform(regression,output='predict')
    df['score'] = df[column].groupby('trajectory').transform(regression,output='score')
    df['score_null'] = df[column].groupby('trajectory').transform(linearRegression,output='score')
    sortOut = []
    for names,groups in df.groupby('trajectory'):        
        if groups['score'].mean()>Rscore:
            #if groups['regression'].tail(5).mean() < groups['regression'].head(5).mean():
            values = groups['regression'].unique() #contains two values
            print(values)
            if 2*list(values)[1] < list(values)[0]:
                #print(groups['score'].mean(),groups['score_null'].mean())
                if groups['score'].mean()>groups['score_null'].mean():
                #print(names,'looks promising')
                    continue       
        #print('filter pout:',names)
        sortOut.append(names)
            #print(sortOut)
        #print('drop',i)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    final.reset_index(inplace=True)
    return final,sortOut   



# def plot_all_trajs(df,out):
#     fig,ax = plt.subplots(3,3)
#     indeces = list(itertools.product([0,1,2],repeat =2))
#     n=0
#     for name,group in df.groupby("trajectory"):                
#         rpa = ((group["Intensity_RPA"]/group["Intensity_RPA"].max())*group["Intensity_DNA"].max()).rolling(4, min_periods=1).mean()
#         ax[indeces[n]].plot(group["seconds"],rpa,color='magenta')
#         ax[indeces[n]].plot(group["seconds"],group["Intensity_DNA"],color='black')
#         #ax[indeces[n]].plot(group["seconds"],group["regression"],color='red')
#         #ax[indeces[n]].plot(group["seconds"],-group["derivative"],color='blue')
#         ax[indeces[n]].set_title('trajectory {}'.format(name))
#         #find the corresponding segments:
#         #segment = df_segment.loc[df_segment['trajectory'] == name]
#         #segmentPlotter(segment,ax[indeces[n]])
#         n+=1
#         if n == 9:
#             #save figure and close it:
#             fig.savefig(out+'trajectory_'+str(name))
#             plt.close(fig)
#             #make new one:
#             fig,ax = plt.subplots(3,3)
#             n= 0

def segmentPlotter(seg_data,ax):
    #steps = True:
    for index, row in seg_data.iterrows():
        ax.plot([row["x1"],row["x2"]],[row["y1"],row["y2"]],color='red')
        
def plot_all_trajs(df_integration,out,df_segment=None,segments=False, xcolumn = 'slice', ycolumn = 'Intensity', group = 'trajectory'):
    fig,ax = plt.subplots(3,3)
    ax_positions = itertools.permutations((0,1,2),2)
    indeces = list(itertools.product([0,1,2],repeat =2))
    n=0
    
    if isinstance(ycolumn,str):        
        for name,group in df_integration.groupby(group):                 
            ax[indeces[n]].plot(group[xcolumn],group[ycolumn])
            ax[indeces[n]].set_title('trajectory {}'.format(name))
            if segments:
                #find the corresponding segments:
                segment = df_segment.loc[df_segment[group] == name]
                segmentPlotter(segment,ax[indeces[n]])
            n+=1
    if isinstance(ycolumn,list):
        for name,group in df_integration.groupby(group):                 
            for column in ycolumn:
                ax[indeces[n]].plot(group[xcolumn],group[column])
                ax[indeces[n]].set_title('trajectory {}'.format(name))
                if segments:
                    #find the corresponding segments:
                    segment = df_segment.loc[df_segment[group] == name]
                    segmentPlotter(segment,ax[indeces[n]])
            n+=1
    
    if n == 9:
        #save figure and close it:
        fig.savefig(out+'trajectory_'+str(name))
        plt.close(fig)
        #make new one:
        fig,ax = plt.subplots(3,3)
        n= 0
            
def normalize_all_trajs(df,ycolumn,head=0,tail=0):
    df['norm_'+ycolumn] = df.groupby('trajectory')[ycolumn].transform(norm,head=head,tail=tail)
    return df

def norm(traj,head=100,tail=180):
    #traj_temp = traj.drop(traj.index[0:head],inplace=False)
    #traj_temp.drop(traj.index[traj.index[-1]-tail:],inplace=True)    
    normalizationValue = traj.max()  
    traj = (traj/normalizationValue)*1.
    return traj             

def smooth(df,wd):
    if isinstance(df,pd.DataFrame):   
        smoothed = df.rolling(wd).mean().fillna(method='ffill')
        return smoothed
    else:
        df = pd.Series(df)
        smoothed = df.rolling(wd).mean()
        return smoothed

            
def plot_allBig(df,out):
    
    #indeces = list(itertools.product([0,1,2],repeat =2))
    #n=0
    for name,group in df.groupby("trajectory"):                        
        fig,ax = plt.subplots(2,sharex=True)
        ax[0].plot(group["seconds"],group["Intensity_DNA"],color='black')
        rpa = (group["Intensity_RPA"]).rolling(4, min_periods=1).mean()
        ax[1].plot(group["seconds"],rpa,color='magenta')
        #ax[1].plot(group["seconds"],group["regression"],color='red')
        #ax[indeces[n]].plot(group["seconds"],-group["derivative"],color='blue')
        ax[0].set_title('trajectory {}'.format(name))
        #find the corresponding segments:
        #segment = df_segment.loc[df_segment['trajectory'] == name]
        #segmentPlotter(segment,ax[indeces[n]])
        fig.savefig(out+'trajectory_'+str(name))
        plt.close(fig)
        #n+=1
        # if n == 9:
        #     #save figure and close it:
        #     fig.savefig(out+'trajectory_'+str(name))
        #     plt.close(fig)
        #     #make new one:
        #     fig,ax = plt.subplots(3,3)
        #     n= 0