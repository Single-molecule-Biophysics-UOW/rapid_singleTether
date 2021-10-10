# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 21:04:57 2021

@author: shm975
"""


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
def start_synch(df, column = 'corr_Intensity_RPA'):
    """
    atempt to find the begin of strand displacement by thresholding.
    Intensity of 1 RPA is 6500 => set threshold to 2 RPAs binding
    """
    n=0
    
    
    
    
    thresholded = []
    
    df['smoothRPA'] = df[column].rolling(4,min_periods=1).mean()
    
    
    df['thresholded'] = df['smoothRPA'].clip(upper=13000)
    #df.reset_index(inplace=True)
    synchRPA = alignData(df)
        #first do a rolling avg
    return synchRPA


def alignTime(df):
    thresholded = np.array(df)
    oldTime = df.index.get_level_values(1)
    changePoint = thresholded.argmax()#np.where(np.gradient(regr)!=0)[0][0]     
    newTime = oldTime -changePoint
    return newTime
def calcShift(df):
    """calculate the shift for aligment as in alignTime, but only return the shift"""
    thresholded = np.array(df)
    oldTime = df.index.get_level_values(1)
    changePoint = thresholded.argmax()   
    return changePoint


def alignData(df,timeConversion=1.,traj_column = 'trajectory',slice_column = 'slice',returnShift = False):
    df = df.set_index([traj_column,slice_column]).sort_index()
    df['alignedTimeRPA'] = df['thresholded'].groupby(traj_column).transform(alignTime)
    if returnShift:
        cpList = []
        for name,groups in df['thresholded'].groupby(traj_column):
            changePoint = calcShift(groups)
            cpList.append(changePoint)
#convert time from frames to seconds:
    df['seconds_new'] = df['alignedTimeRPA'].mul(timeConversion)
    
    if returnShift:
        return df,cpList
    return df
    # n=0    
    # for name,group in df.groupby('trajectory'):
    #     fig = plt.figure()
    #     plt.plot(df['seconds'],df['thresholded'])
    #     plt.plot(df['seconds'],df['smoothRPA'])
    #     fig.show()
    #     if n>1:
    #         break
    #     n+=1
        #now threshold:
        
    
    