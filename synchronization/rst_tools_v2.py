# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 14:04:54 2021

@author: shm975
"""

import traceback
import itertools
from matplotlib import pyplot as plt
import numpy as np
from synchronization import STA as sta
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

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
        
def plot_all_trajs(df_integration,out,df_segment=None,segments=False, xcolumn = 'slice', ycolumn = 'Intensity', groupName = 'trajectory'):
    """
    

    Parameters
    ----------
    df_integration : pd.DataFrame
        datframe containing all trajectories to plot.
    out : String
        path to output folder
    df_segment : pd.DataFrame, optional
        dataframe containing segments from changePoint
    segments : bool, optional
        Segments will be plotted if True. The default is False.
    xcolumn : String, optional
        xcolumn name. The default is 'slice'.
    ycolumn :  String, optional
        ycolumn name. The default is 'Intensity'.
    groupName : String, optional
        Name to group data by. The default is 'trajectory'.
        Necessary if data contains more than one trajectory

    Returns
    -------
    None.

    """
    fig,ax = plt.subplots(3,3)
    ax_positions = itertools.permutations((0,1,2),2)
    indeces = list(itertools.product([0,1,2],repeat =2))
    n=0
    
    if isinstance(ycolumn,str):        
        for name,group in df_integration.groupby(groupName):                             
            ax[indeces[n]].plot(group[xcolumn],group[ycolumn])
            ax[indeces[n]].set_title('trajectory {}'.format(name))
            if segments:
                #find the corresponding segments:
                segment = df_segment.loc[df_segment[groupName] == name]
                segmentPlotter(segment,ax[indeces[n]])
            n+=1
            if n == 9:
                #save figure and close it:
                fig.savefig(out+'trajectory_'+str(name)+'.jpg')
                plt.close(fig)
                #make new one:
                fig,ax = plt.subplots(3,3)
                n= 0
            
    if isinstance(ycolumn,list):
        for name,group in df_integration.groupby(groupName):
            if 'rate' in ycolumn:
                label = '{:.2f} bp/s'.format((group['rate'].unique())[0])
            else:
                label = None                 
            for column in ycolumn:
                #print ('n:{}, column:{}'.format(n,column))
                try:
                    t = (indeces[n])
                except:
                    print("indexing error")
                try:
                    t=(group[column])
                except:
                    print('ycolumn problem')
                line, =ax[indeces[n]].plot(group[xcolumn],group[column])
                
                
                #ax[indeces[n]].plot([0,group[xcolumn].max()],[group[column].mean()+10*group[column].std(),group[column].mean()+10*group[column].std()])
                ax[indeces[n]].set_title('trajectory {}'.format(name))
                if segments:
                    #find the corresponding segments:                    
                    segment = df_segment.loc[df_segment[groupName] == name]
                    segmentPlotter(segment,ax[indeces[n]])  
            line.set_label(label)
            ax[indeces[n]].legend()
            
            n+=1
            if n == 9:
                #save figure and close it:
                fig.savefig(out+'trajectory_'+str(name)+'.png')
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

def smooth(df,wd=5):
    if isinstance(df,pd.DataFrame):   
        smoothed = df.rolling(wd).mean().fillna(method='ffill')
        return smoothed
    else:
        df = pd.Series(df)
        smoothed = df.rolling(wd).mean()
        return smoothed
def savgol_filter_smooth(df,wd=2,polyorder=2):
    ydata = np.array(df)
    smoothed = savgol_filter(ydata,wd,polyorder)
    return smoothed
    
def smooth_all_trajectories(df,column,wd,groupName = 'trajectory',result_prefix = 'smooth',polyorder = 3,method = 'savgol'):
    print(df.keys())
    groupedDF = df.groupby(groupName)
    print(groupedDF)
    
    if method == 'savgol':
        df[result_prefix + '_'+column] = groupedDF[column].transform(savgol_filter_smooth,wd=wd,polyorder = polyorder)
    if method == 'window':
        print('window mean')
        df[result_prefix + '_'+column] = groupedDF[column].transform(smooth,wd=wd)
    #return df


def regression(df,maxDepth=1,output = 'predict'):
    x,y = np.array(df.index.get_level_values(0)).reshape(-1,1),np.array(df).reshape(-1,1)
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
          
def find_unbinding(df, column, Rscore, threshold):
    #first do regression:
    df['regression'] = df.groupby('trajectory')[column].transform(regression,output='predict')
    df['score'] = df.groupby('trajectory')[column].transform(regression,output='score')
    sortOut = []
    sortOutReason = []
    for names,groups in df.groupby('trajectory'):        
        if groups['score'].mean()>Rscore:
            #if groups['regression'].tail(5).mean() < groups['regression'].head(5).mean():
            values = groups['regression'].unique() #contains two values
            
            stepSize = values[0]-values[1]
            
            if stepSize >= threshold:
                continue
            else:
                sortOut.append(names)
                sortOutReason.append('smallStep')
        else:
            sortOut.append(names)
            sortOutReason.append('lowScore')
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    final.reset_index(inplace=True)
    return final, sortOut, sortOutReason
def derivative(df):
    return df.diff().fillna(0).abs()

def find_unbinding_v2(df, column, Rscore, threshold):
    """
    This will find dissociation events based on the first derivative.
    The regression technique is good at finding the timepoint of dissociation, but not 
    so much for deciding if it is one step or not!

    Parameters
    ----------
    df : pd.DataFrame
        input data
    column : String
        column name of y data to analyse.
    Rscore : TYPE
        DESCRIPTION.
    threshold : TYPE
        DESCRIPTION.

    Returns
    -------
    final : TYPE
        DESCRIPTION.
    sortOut : TYPE
        DESCRIPTION.
    sortOutReason : TYPE
        DESCRIPTION.

    """
    #find first derivative:
    df['derivative'] = df.groupby('trajectory')[column].transform(derivative)
    #the derivative should have a very clear peak (a minimum?) at the time of dissociation!
    #find that time by the regression method:
    df['regression'] = df.groupby('trajectory')[column].transform(regression,output='predict')
    df['score'] = df.groupby('trajectory')[column].transform(regression,output='score')
    sortOut = []
    sortOutReason = []
    for name,group in df.groupby('trajectory'):        
        timing = group[group['regression'].diff().fillna(0) !=0]    #this is where the dissociation happends
        #print(len(timing))
        #now see if the derivative at that time is above/below a threshold:
        #to avoid +-1 indexing errors take derivative
        if timing['derivative'].mean() < threshold:
            sortOut.append(name)
            sortOutReason.append('smallDerivative')
            continue
        #now make sure the stepsize is large too:
        values = group['regression'].unique() #contains two values
        stepSize = values[0]-values[1]
        if stepSize < threshold:
            sortOut.append(name)
            sortOutReason.append('smallStep')
            continue
        
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    final.reset_index(inplace=True)
    return final, sortOut, sortOutReason        
        
        
        
        
def sortOut_trajectories(df, sortOut, groupName = 'trajectory'):
    """
    removes trajectories (groups) from dataframe based on list

    Parameters
    ----------
    df : pd.DataFrame
        dataFrame that can be grouped into trajectories
    sortOut : array-like
        array/list containing trajectory names to sort out

    Returns
    -------
    filtered pd.DataFrame

    """
    print(df)   
    filtered = df[df.groupby(groupName).filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    final.reset_index(inplace=True)
    return final
        
    #     if groups['score'].mean()>Rscore:
    #         #if groups['regression'].tail(5).mean() < groups['regression'].head(5).mean():
    #         values = groups['regression'].unique() #contains two values
            
    #         stepSize = values[0]-values[1]
            
    #         if stepSize >= threshold:
    #             continue
    #         else:
    #             sortOut.append(names)
    #             sortOutReason.append('smallStep')
    #     else:
    #         sortOut.append(names)
    #         sortOutReason.append('lowScore')
    # filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    # final = filtered.dropna()
    # final.reset_index(inplace=True)
    # return final, sortOut, sortOutReason



def unbinding_timing(df,column, groupName,time_corr = False):
    #loop through all trajectories
    # and append the time of unbinding (i.e regression.diff()!=0 to a list)
    all_timing = []
    for name,group in df.groupby(groupName):
        timing = group[group['regression'].diff().fillna(0) !=0]
        all_timing.append(timing[column].iloc[0])
        #if time_corr is True apply a corrected time axis that starts when the first molecule unbinds:
        
    timingDF = pd.DataFrame({column:all_timing})
    return timingDF
    
    
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
    
def piecewise_linear(x,m,a,t):
    """
    \
     \
      \
       \_________

    """
    f1 = lambda x: m*x+t
    f2 = lambda x: m*a+t
    piecewise = np.piecewise(x,[x<=a],[f1,f2])
    return piecewise

def piecewise_linear_triple(x,m,a,b,t):
    """
    ________
            \
             \
              \
               \_________

    """
    f0 = lambda x: m*a+t
    f1 = lambda x: m*x+t
    f2 = lambda x: m*b+t
    condlist = (x<a,((x>=a)&(x<=b)),x>b)
    
    
    piecewise = np.piecewise(x,condlist,[f0,f1,f2])
    return piecewise

def lin(x,m,t):
    return m*x+t

def complex_norm(df,xy = None,sem=None,initial_guess =[-200,30,300,100000], decrease =  True):
    if xy:
        t = np.array(df[xy[0]])
        Intensity = np.array(df[xy[1]])
        print(Intensity)
    else:
        t = df.index
        Intensity = np.array(df)
        print(Intensity)
    if isinstance(sem,str):
        semData = np.array(df[sem])
        # print('norm sem')
    

    # print(t)
    #first do the fit and return popt:
    try:
        # print('trying to fit')
        popt,pcov= curve_fit(piecewise_linear_triple, t,Intensity,p0=initial_guess,maxfev=10000)
    except ValueError:
        traceback.print_exc()
        raise ValueError      
    except Exception as inst:
        print(type(inst))
        # print(t)
        # print(Intensity)
        # plt.plot(t,Intensity)
        # plt.plot(t,piecewise_linear_triple(t,-20,30,300,100000))
        # plt.show()
        norm_Intensity = np.ones(len(t))*np.nan
        return norm_Intensity
    # intensity when time <popt[1]: 2600 bp
    # intensity when time >popt[2]: 0 bp
    if decrease:
        maxValue = popt[0]*popt[1]+popt[3]
        offsetValue = popt[0]*popt[2]+popt[3]
        #now normalize:
        norm_Intensity = (Intensity - offsetValue)/(maxValue- offsetValue)
        if isinstance(sem,str):
            norm_sem = (semData)/(maxValue)
    else:
        maxValue = popt[0]*popt[2]+popt[3]
        offsetValue = popt[0]*popt[1]+popt[3]
        #now normalize:
        norm_Intensity = (Intensity - offsetValue)/(maxValue-offsetValue)
        if isinstance(sem,str):
            norm_sem = (semData)/(maxValue)

    
    if isinstance(sem,str):
        # print('normed sem')
        return [norm_Intensity, norm_sem]
    else:
        return norm_Intensity