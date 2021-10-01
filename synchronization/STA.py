# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:13:48 2020

@author: shm975
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
import toolbox as tb
import os
import seaborn as sns


from tkinter import filedialog
from tkinter import Tk



#from skimage import filters





def importerGui(allData=True,n=0):
    """
    This methods provides a user-interface to load the data. It asks for folder, which needs to contain
    the following subfolders:
        1.: /integrations/ the trajectory integrations are saved in here
        2.: /meta/ the meta data for every traj-integration file is saved here

    Returns
    -------
    1.: loaded data as pandas dataframe. if allData = True a list of dataFrames is returned
    2.: loaded preReactionData as pandas dataframe. if allData = True a list of dataFrames is returned
    3.: timeConversion factor calculated from metadata. if allData = True a list of of floats is returned

    """
    print('this function was changed, it returns the folder and the chosen filename now!')
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    #ask for directory containing the data
    folder = filedialog.askdirectory(title='choose directory')
    subfolders = list(os.walk(folder))[0][1]
    
    
    if 'meta' not in subfolders:
        print('folder \'meta\' not found!')
        return None
    if 'integrations' not in subfolders:
        print('folder \'integrations\' not found. all trajectories should be in a folder of this name!')
        return None

    C1_preReactions = []
    C2_preReactions = []
    C1_reactions = []
    C2_reactions = []

    for name in os.listdir(folder+'/integrations'):
        if 'C1' in name:
            if 'prereaction' in name:
                C1_preReactions.append(name)
            else:
                C1_reactions.append(name)
        elif 'C2' in name:
            if 'prereaction' in name:
                C2_preReactions.append(name)
            else:
                C2_reactions.append(name)
    metaFileNames = []            
    for name in os.listdir(folder+'/meta'):
        metaFileNames.append(name)
    #do a sanity check:
    lengths = np.array([len(C1_preReactions),len(C1_reactions),len(C2_reactions)])
    if not np.all(lengths == lengths[0]):
        print('not all same length',lengths)
        return None
    if len(metaFileNames) != len(C1_preReactions)+len(C1_reactions):
        print('some meta data seems to be missing!')
        if len(metaFileNames) > len(C1_preReactions)+len(C1_reactions):
            print ('there is more meta data than trajectory integrations')
        elif len(metaFileNames) < len(C1_preReactions)+len(C1_reactions):
            print ('there is more trajectory integrations data than meta data')
        elif len(C1_preReactions) > len(C1_reactions):
            print('more preReaction files than reaction files')
        elif len(C1_preReactions) < len(C1_reactions):
            print('more reaction files than prereaction files')
        return None
    print('sucess!')

    dataList = []
    preReactionList = []
    timeConversionList = []
    baseFileNameList = []
    #loop over all data
    #if allData = True n=0 means looping over whole list, if
    #allData is false and n is provided loop over datafile of interest to the end and only return first file
    
    print((C1_reactions))
    print(list(range(n,len(C1_reactions))))
    for data in range(n,len(C1_reactions)):
        print(data)
        print('start loading files:')
        baseFileName = C1_reactions[data].replace('reaction_','').replace('C1-corrected_','').replace('.nd2_series_1.csv','').replace(
            '_series_1','').replace('.csv','')
        print(baseFileName)
        print(metaFileNames)
        print(C1_reactions)
        for file in C1_reactions:
            if baseFileName in file:
                C1_reactionName = file
            else:
                print('C1_reactionFile not found!',C1_reactions)
                print('baseFileName:',baseFileName)
        for file in C1_preReactions:
            if baseFileName in file:
                C1_preReactionName = file
        for file in C2_reactions:
            if baseFileName in file:
                C2_reactionName = file
        for file in C2_preReactions:
            if baseFileName in file:
                C2_preReactionName = file
        for file in metaFileNames:
            print(baseFileName)
            print(file)
            print('is the base File name in the list?')
            if baseFileName in file:
                metaFileName = file    
        try:
            print('C1 reation name:'+C1_reactionName)
            C1_reaction = pd.read_csv(folder+'/integrations/'+C1_reactionName)
        except Exception:
            print('C1 reaction file not found!', C1_reactionName)
        try:        
            C2_reaction = pd.read_csv(folder+'/integrations/'+C2_reactionName)
        except Exception:
            print('C2 reaction file not found!', C1_reactionName)
        try:        
            C1_preReaction = pd.read_csv(folder+'/integrations/'+C1_reactionName)
        except Exception:
            print('C1 preReaction file not found!', C1_reactionName)
        try:        
            C2_preReaction = pd.read_csv(folder+'/integrations/'+C2_reactionName)
        except Exception:
            print('C2 preReaction file not found!', C1_reactionName)
        try:
            timeConversion = findTimeCov(folder+'/meta/'+metaFileName)
        except Exception:
            print('meta file not found!',metaFileName)
        print('files loaded!')
        df = pd.merge(C1_reaction,C2_reaction, on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])
        df.name = baseFileName
        dataList.append(df)
        preReactionList.append(C2_preReaction)
        timeConversionList.append(timeConversion)
        baseFileNameList.append(baseFileName)
        #if len(dataList) == 1:
        if allData == False:
            return folder,dataList[0],preReactionList[0],timeConversionList[0],baseFileName
            print('all data is false')
            #print('return single file')
        
    return folder, dataList,preReactionList,timeConversionList,baseFileNameList
        #else: 

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
def findEvents(df,Rscore,return_events=True,column = 'Intensity_DNA'):
    """find trajectories where DNA is clearly unbinding
    by using a regression analysis and subsequent threshold. Rscore is the combined R**2 of linear regressions
    if return_events=False, the trajectories not showing this will be returned
    """
    
    
    df.set_index(['trajectory','slice'],inplace=True)
    df.sort_index(inplace=True)
    #print(df.keys())
    df['regression'] = df[column].groupby('trajectory').transform(regression,output='predict')
    df['score'] = df[column].groupby('trajectory').transform(regression,output='score')
    
    sortOut = []
    for names,groups in df.groupby('trajectory'):
        #print(groups)
        #print(groups['score'].head(1))
        #print(groups['score'].mean())        
        if groups['score'].mean()>Rscore:
            if groups['regression'].tail(5).mean() < groups['regression'].head(5).mean():
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
    #print(df.keys())
    df['regression'] = df[column].groupby('trajectory').transform(regression,output='predict')
    # df['regression_null'] = df['Intensity_DNA'].groupby('trajectory').transform(regression,output='predict',max_depth=0)
    df['score'] = df[column].groupby('trajectory').transform(regression,output='score')
    df['score_null'] = df[column].groupby('trajectory').transform(linearRegression,output='score')
    sortOut = []
    for names,groups in df.groupby('trajectory'):        
        if groups['score'].mean()>Rscore:
            if groups['regression'].tail(5).mean() < groups['regression'].head(5).mean():
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
    
    
    
    # if return_events == True:
    #     events = df[df['score']>Rscore]
    #     events.reset_index(inplace=True)
    #     grouped = events.groupby(['trajectory'])    
    #     #make sure intensity of DNA decreases
    #     sort_out = lambda x:x['Intensity_DNA'].tail(int(len(x)*0.2)).mean() < x['Intensity_DNA'].head(int(len(x)*0.2)).mean()
    #     sortOut = []
    #     for name,groups in grouped:
    #         if sort_out(groups) ==True:                
    #             sortOut.append(name)
    #     filteredEvents = (grouped.filter(lambda x:x['Intensity_DNA'].tail(int(len(x)*0.2)).mean() < x['Intensity_DNA'].head(int(len(x)*0.2)).mean()))      
    #     #print('sort out because of DNA: ',sortOut)
    #     return filteredEvents,sortOut
    
def process(merged,Rscore,complete=True):
    grouped = merged.groupby('trajectory')
    newM = merged.set_index(['trajectory','slice']).sort_index()
    newM['regression'] = newM['Intensity_DNA'].groupby('trajectory').transform(regression,output='predict')
    newM['score'] = newM['Intensity_DNA'].groupby('trajectory').transform(regression,output='score')
    #add number of RPAs column  by dividing by intensity
    newM['n_RPA'] = newM['Intensity_RPA'].div(4000.)
    #add number of nucleotides synthesized:
    newM['nts'] = newM['n_RPA'].mul(30.)
    
    if complete == True:
        newF = newM[newM['score']>Rscore]
        data = newF.reset_index()
        grouped = data.groupby(['trajectory'])    
        #make sure intensity of DNA decreases
    
        filteredData = (grouped.filter(lambda x:x['Intensity_DNA'].tail(int(len(x)*0.2)).mean() < x['Intensity_DNA'].head(int(len(x)*0.2)).mean()))      
        return filteredData
    else:
        newF = newM[newM['score']<Rscore]
        data = newF.reset_index()       
        return data

def RPAfilter(df):
    #make sure intensity of RPA increases
    grouped = df.groupby('trajectory')
    sortOut = []
    
    for names,groups in grouped:
        #print(groups['regression'].gradient())
        changePoint = np.where(np.gradient(groups['regression'])!=0)[0][0] 
        RPA = np.array(groups['Intensity_RPA'])
        #print(len(groups))
        #print(np.array(groups['Intensity_RPA'])[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))])
        #if np.mean(smooth(RPA[0:int(0.1*traj_len)],15)) <= 0.2 * np.mean(smooth(RPA[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))],15)):
        if np.mean(smooth(RPA[changePoint-10:changePoint+10],15)) <= np.mean(RPA[0:10]):
            sortOut.append(names)              
            #print(sortOut)
        #print('drop',i)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    return(final,sortOut)

def RPAfilter_v2(df,minimumCP=15):
    """
    do something smarter than just a threshold!

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    grouped = df.groupby('trajectory')
    sortOut = []
    
    for names,groups in grouped:
        #print(groups['regression'].gradient())
        changePoint = np.where(np.gradient(groups['regression'])!=0)[0][0] 
        RPA = np.array(groups['Intensity_RPA'])
        DNA = np.array(groups['Intensity_DNA'])
        
        if changePoint<minimumCP:
            #that means the dissociation happened inthe first 15 datapoints, therefore most of the event might be missed
            print('WARNING: event is sorted out because dissociation happened too early. You might want to adjust that if you are looking at something fast!')
            sortOut.append(names)
        
        #print(len(groups))
        #print(np.array(groups['Intensity_RPA'])[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))])
        #if np.mean(smooth(RPA[0:int(0.1*traj_len)],15)) <= 0.2 * np.mean(smooth(RPA[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))],15)):
        #the window lengths are onlhy different to accomodate for the short amount of data looked at! the window hsa to be shorter than the input vector
        if np.nanmean(smooth(RPA[changePoint-minimumCP:changePoint+10],15)) < np.nanmean(smooth(RPA[0:15],7)+20000):
            
            sortOut.append(names)
            continue
        if np.mean(smooth(RPA[changePoint-minimumCP:changePoint+10],15)) < 20000:
            sortOut.append(names)
            continue
        if np.mean(smooth(RPA[0:20],15)) > 30000:
            #if the starting Intensity is unreasonably high, there might be an aggregate etc. sort out.
            sortOut.append(names)
        # else:
        #     plt.figure()
        #     print('weird event')
        #     print((RPA[0:10]))
        #     plt.plot(smooth(RPA,15),label = 'RPA')
        #     plt.plot(DNA)
        #     plt.plot([changePoint,changePoint],[0,20000])
        #     plt.title('mean around CP:{}, mean at start:{}'.format(np.mean(smooth(RPA[changePoint-10:changePoint+10],15)),np.nanmean(smooth(RPA[0:15],7))))
        #     plt.legend()
            #raise Exception
            #print(sortOut)
        #print('drop',i)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    return(final,sortOut)


def RPAfilter_v3(df,minimumCP=15):
    """
    run regression analysis on rpa curve to categorise all data points in two 
    categories. Then take the delta between the two lines. It has to be above a certain 
    threshold!

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    df['RPA_regression'] = df['Intensity_RPA'].groupby('trajectory').transform(regression,output='predict')
    
    
    
    grouped = df.groupby('trajectory')
    sortOut = []
    
    for names,groups in grouped:
        
        high,low = groups['RPA_regression'].unique()
        # print(groups.keys())
        # plt.plot(groups['seconds'],groups['RPA_regression'],label = 'high={},low={}'.format(high,low))
        # plt.plot(groups['seconds'],groups['Intensity_RPA'],label = 'high={},low={}'.format(high,low))
        # plt.legend()
        
        if np.abs(high-low) < 100000.:
            sortOut.append(names)
        #print(groups['regression'].gradient())
        # changePoint = np.where(np.gradient(groups['regression'])!=0)[0][0] 
        # RPA = np.array(groups['Intensity_RPA'])
        # DNA = np.array(groups['Intensity_DNA'])
        
        # if changePoint<minimumCP:
        #     #that means the dissociation happened inthe first 15 datapoints, therefore most of the event might be missed
        #     print('WARNING: event is sorted out because dissociation happened too early. You might want to adjust that if you are looking at something fast!')
        #     sortOut.append(names)
        
        # #print(len(groups))
        # #print(np.array(groups['Intensity_RPA'])[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))])
        # #if np.mean(smooth(RPA[0:int(0.1*traj_len)],15)) <= 0.2 * np.mean(smooth(RPA[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))],15)):
        # #the window lengths are onlhy different to accomodate for the short amount of data looked at! the window hsa to be shorter than the input vector
        # if np.nanmean(smooth(RPA[changePoint-minimumCP:changePoint+10],15)) < np.nanmean(smooth(RPA[0:15],7)+20000):
            
        #     sortOut.append(names)
        #     continue
        # if np.mean(smooth(RPA[changePoint-minimumCP:changePoint+10],15)) < 20000:
        #     sortOut.append(names)
        #     continue
        # if np.mean(smooth(RPA[0:20],15)) > 30000:
        #     #if the starting Intensity is unreasonably high, there might be an aggregate etc. sort out.
        #     sortOut.append(names)
        # else:
        #     plt.figure()
        #     print('weird event')
        #     print((RPA[0:10]))
        #     plt.plot(smooth(RPA,15),label = 'RPA')
        #     plt.plot(DNA)
        #     plt.plot([changePoint,changePoint],[0,20000])
        #     plt.title('mean around CP:{}, mean at start:{}'.format(np.mean(smooth(RPA[changePoint-10:changePoint+10],15)),np.nanmean(smooth(RPA[0:15],7))))
        #     plt.legend()
            #raise Exception
            #print(sortOut)
        #print('drop',i)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()
    return(final,sortOut)



    
def RPAfilter2(df):
    #make sure intensity of RPA increases
    grouped = df.groupby('trajectory')
    sortOut = []
    
    for names,groups in grouped:
        #print(groups['regression'].gradient())
        #changePoint = np.where(np.gradient(groups['regression'])!=0)[0][0] 
        RPA = groups['Intensity_RPA']
        #print(len(groups))
        #print(np.array(groups['Intensity_RPA'])[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))])
        #if np.mean(smooth(RPA[0:int(0.1*traj_len)],15)) <= 0.2 * np.mean(smooth(RPA[changePoint-int(0.1*(len(groups))):changePoint+int(0.1*(len(groups)))],15)):
        if np.mean(RPA.tail(25)) <= 20000.+np.mean(RPA.head(25)):
            sortOut.append(names)              
            #print(sortOut)
        #print('drop',i)
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOut).astype('bool')]
    final = filtered.dropna()

    return(final,sortOut)

def Manualfilter(df,sortOutList):    
    filtered = df[df.groupby('trajectory').filter(lambda x: x.name not in sortOutList).astype('bool')]
    final = filtered.dropna()
    return(final)

    #filteredData = (grouped.filter(lambda x:x['Intensity_RPA'].tail(int(len(x)*0.2)).mean() > x['Intensity_DNA'].head(int(len(x)*0.2)).mean()))  
    #return df
 
def int_cal(df, RPA_int,RPA_footprint=30.,assume_complete=False, template_length= 2600):
    df['n_RPA'] = df['Intensity_RPA'].div(RPA_int)
    if not assume_complete:
        df['nucleotides'] = df['n_RPA'].mul(RPA_footprint)
    else:
        #assume the reaction goes to completion, determine nucleotides
        #by setting max intensity to template length
        df['nucleotides'] = df['Intensity_RPA']
        fun = lambda x: x*template_length/x.max()
        df['nucleotides'] = df.groupby('trajectory')['nucleotides'].transform(fun)
        #print(df.groupby('trajectory')['nucleotides'].transform(lambda x: x*template_length/x.max()))
        #df.trajectory.map(df.groupby('trajectory')['Intensity_RPA'].sum())
        # for names,groups in df.groupby('trajectory'):
        #     groups['nucleotides'] = groups['Intensity_RPA'].div(groups['Intensity_RPA'].max()).mul(template_length)
        #df['nucleotides'] = df['Intensity_RPA'].div(df['Intensity_RPA'].max()).mul(template_length)
    return(df)
def truncate(df,cutoff):
    df = df[df.seconds <= cutoff]
    #fun = lambda g: np.where(g.seconds <= cutoff)
    #print(fun(df.groupby('trajectory')))
    #df.groupby('trajectory').filter(lambda g: g.seconds <= cutoff)#[lambda x: x<=cutoff]
    #for names,groups 
    
    #
    #df.groupby(['Date','Advertiser']).ID.count()[lambda x: x >= 500]
    #print(df.groupby('trajectory').filter(lambda x: x[x['seconds']>=cutoff]))
    return df
def calibrate(merged,conv):
    #grouped = merged.groupby('trajectory')
    #newM = merged.set_index(['trajectory','slice']).sort_index()
    #add time column
    merged['time'] = merged['slice'].mul(conv)   
    return merged
# def filterRPAIntensity(df):
#     regr = np.array(df)
#     changePoint = np.where(np.gradient(regr)!=0)[0][0]
    #print(changePoint)
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

        # for segment in temp:
        # plt.plot(segment[2],segment[4],color='magenta')
        # plt.plot([segment[2],segment[3]],[segment[4],segment[5]],color='magenta')
        #plt.plot(segment[3],segment[5])
    
    # plt.savefig(trajFolder+'traj_{}'.format(names))
    # plt.close()
    
def traj_plotter(trajData,trajNumber,ycolumn='nucleotides',xcolumn = 'seconds',segments = False,segmentData=None,normalize = True):
#plt.switch_backend('Agg')
    if normalize:
        trajgroups = trajData.groupby('trajectory')
        
        #get data for group:
        traj = trajgroups.get_group(trajNumber)
        
        #print(traj['nucleotides'].max())
        fig,ax = plt.subplots()
        ax.plot(traj[xcolumn],traj[ycolumn]/traj[ycolumn].max(),color = 'grey',alpha=0.5)
        #ax.plot(traj[xcolumn],tb.window((traj['Intensity_DNA']/traj['Intensity_DNA'].max()),window_len=5),color = 'grey',alpha=0.5)
        ax.set_ylabel('intensity [arb Units]')
        #ax2=plt.twinx(ax)
        normedRPA = traj[ycolumn]/traj[ycolumn].max()
        ax.plot(traj[xcolumn],normedRPA,color='magenta',alpha = 0.5)
        #ax2.plot(traj[xcolumn],tb.window(traj[ycolumn],window_len=5),color='magenta')
        if segments:
            seggroups = segmentData.groupby('trajectory')
            seg = seggroups.get_group(trajNumber)
            ax.plot([seg['x1'],seg['x2']],[seg['y1'],seg['y2']],color='black',linestyle='--')
        return fig,ax
    else:
        trajgroups = trajData.groupby('trajectory')
        
        #get data for group:
        traj = trajgroups.get_group(trajNumber)
        
        #print(traj['nucleotides'].max())
        fig,ax = plt.subplots()
        ax.plot(traj[xcolumn],traj[ycolumn],color = 'grey',alpha=0.5)
        #ax.plot(traj[xcolumn],tb.window((traj['Intensity_DNA']/traj['Intensity_DNA'].max()),window_len=5),color = 'grey',alpha=0.5)
        ax.set_ylabel('intensity [arb Units]')
        #ax2=plt.twinx(ax)
        normedRPA = traj[ycolumn]
        ax.plot(traj[xcolumn],normedRPA,color='magenta',alpha = 0.5)
        #ax2.plot(traj[xcolumn],tb.window(traj[ycolumn],window_len=5),color='magenta')
        if segments:
            seggroups = segmentData.groupby('trajectory')
            seg = seggroups.get_group(trajNumber)
            ax.plot([seg['x1'],seg['x2']],[seg['y1'],seg['y2']],color='black',linestyle='--')
        return fig,ax
        
def twocolor_traj_plotter(trajData,trajNumber,ycolumn_1='nucleotides',ycolumn_2 = 'DNA',xcolumn = 'seconds',segments = False,segmentData=None,normalize = True):
#plt.switch_backend('Agg')
    if normalize:
        trajgroups = trajData.groupby('trajectory')
        
        #get data for group:
        traj = trajgroups.get_group(trajNumber)
        
        #print(traj['nucleotides'].max())
        fig,ax = plt.subplots()
        ax.plot(traj[xcolumn],traj['Intensity_DNA']/traj['Intensity_DNA'].max(),color = 'grey',alpha=0.5)
        #ax.plot(traj[xcolumn],tb.window((traj['Intensity_DNA']/traj['Intensity_DNA'].max()),window_len=5),color = 'grey',alpha=0.5)
        ax.set_ylabel('intensity [arb Units]')
        #ax2=plt.twinx(ax)
        normedRPA = traj[ycolumn_1]/traj[ycolumn_1].max()
        ax.plot(traj[xcolumn],normedRPA,color='magenta',alpha = 0.5)
        normedDNA = traj[ycolumn_2]/traj[ycolumn_2].max()
        ax.plot(traj[xcolumn],normedDNA,color='black',alpha = 0.5)
        #ax2.plot(traj[xcolumn],tb.window(traj[ycolumn],window_len=5),color='magenta')
        if segments:
            seggroups = segmentData.groupby('trajectory')
            seg = seggroups.get_group(trajNumber)
            ax.plot([seg['x1'],seg['x2']],[seg['y1'],seg['y2']],color='black',linestyle='--')
        return fig,ax
    else:
        trajgroups = trajData.groupby('trajectory')
        
        #get data for group:
        traj = trajgroups.get_group(trajNumber)
        
        #print(traj['nucleotides'].max())
        fig,ax = plt.subplots()
        ax.plot(traj[xcolumn],traj['Intensity_DNA'],color = 'grey',alpha=0.5)
        #ax.plot(traj[xcolumn],tb.window((traj['Intensity_DNA']/traj['Intensity_DNA'].max()),window_len=5),color = 'grey',alpha=0.5)
        ax.set_ylabel('intensity [arb Units]')
        #ax2=plt.twinx(ax)
        normedRPA = traj[ycolumn_1]
        ax.plot(traj[xcolumn],normedRPA,color='magenta',alpha = 0.5)
        normedDNA = traj[ycolumn_2]
        ax.plot(traj[xcolumn],normedDNA,color='black',alpha = 0.5)
        #ax2.plot(traj[xcolumn],tb.window(traj[ycolumn],window_len=5),color='magenta')
        if segments:
            seggroups = segmentData.groupby('trajectory')
            seg = seggroups.get_group(trajNumber)
            ax.plot([seg['x1'],seg['x2']],[seg['y1'],seg['y2']],color='black',linestyle='--')
        return fig,ax
    # for names,groups in new.groupby('trajectory'):
    # #print(groups.keys())
    # plt.figure()
    # end = np.max(np.where(groups['seconds']==1))
    # rpaInt = np.array(tb.window(groups['Intensity_RPA'],window_len=5))[:end+20]    
    # dnaInt = np.array((tb.window(groups['Intensity_DNA'],window_len=5)))[:end+20]
    # time = np.array(groups['seconds'])[:end+20]
    # traj = np.array(groups['trajectory'])[:end+20]
    # plt.plot(time,groups['Intensity_RPA'][:end+20],alpha=0.5,color='magenta')
    # plt.plot(time,groups['Intensity_DNA'][:end+20],alpha=0.5,color='green')
    # plt.plot(time,dnaInt,color='magenta')
    # plt.plot(time,rpaInt,color='green')
    # array = np.column_stack((traj,time,rpaInt))
    # start = t.time()
    # temp = cp.binary_search(array,0, len(rpaInt),5000, 0.99)
    # cptime.append(t.time()-start)
    
    # #print(temp[0])
    # #print('x1:',temp[0][2],'x2:',temp[0][3],'y1:',temp[0][4],'y2:',temp[0][5])
    # if temp==0:
    #     continue
    # for segment in temp:
    #     plt.plot(segment[2],segment[4],color='magenta')
    #     plt.plot([segment[2],segment[3]],[segment[4],segment[5]],color='magenta')
    #     #plt.plot(segment[3],segment[5])
    
    # plt.savefig(trajFolder+'traj_{}'.format(names))
    # plt.close()
    # if temp != 0:
    #     cptable.append(temp)
    # del temp
    

# def makeMean(newfiltered,templateLength=2600):
#     meanDNA = newfiltered.groupby('seconds')['Intensity_DNA'].mean()
#     meanDNA = meanDNA.div(meanDNA.max())
#     stdDNA = newfiltered.groupby('seconds')['Intensity_DNA'].std(ddof=0)
#     stdDNA =stdDNA.div(meanDNA.max())
#     meanRPA = newfiltered.groupby('seconds')['Intensity_RPA'].mean()
#     meanRPA = meanRPA.div(meanRPA.max()).mul(templateLength)
#     stdRPA = newfiltered.groupby('seconds')['Intensity_RPA'].std(ddof=0)
#     stdRPA = stdRPA.div(meanRPA.max())
#     return meanDNA,stdDNA, meanRPA, stdRPA

def makeMean(newfiltered,templateLength=2600,xname='seconds',yDNA = 'Intensity_DNA',yRPA='Intensity_RPA'):
    meanDNA = newfiltered.groupby(xname)[yDNA].mean()
    meanDNA = meanDNA.div(meanDNA.max())
    stdDNA = newfiltered.groupby(xname)[yDNA].std(ddof=0)
    stdDNA =stdDNA.div(meanDNA.max())
    meanRPA = newfiltered.groupby(xname)[yRPA].mean()    
    stdRPA = newfiltered.groupby(xname)[yRPA].std(ddof=0)    
    maxRPA = meanRPA.max()
    meanRPA = meanRPA.div(maxRPA)#.mul(templateLength)
    stdRPA = stdRPA.div(maxRPA)#.mul(templateLength)
    print(np.max(meanRPA))
    print(np.max(stdRPA))
    return meanDNA,stdDNA, meanRPA, stdRPA

def makeMeanSingleColor(newfiltered,templateLength=2600,xname='seconds',yname = 'Intensity'):
    meanDNA = newfiltered.groupby(xname)[yname].mean()
    meanDNA = meanDNA.div(meanDNA.max())
    stdDNA = newfiltered.groupby(xname)[yname].std(ddof=0)
    stdDNA =stdDNA.div(meanDNA.max())
    return meanDNA,stdDNA



def findTimeCov(path):
    timestamps = []
    with open(path) as meta:
        for line in meta:
            if 'timestamp' in line:
                timestamps.append(line.split(',')[1])
    times = np.array(timestamps).astype(float)
    conv = np.round(np.mean(np.diff(times)),2)
    return conv

def smooth(df,wd):
    if isinstance(df,pd.DataFrame):   
        smoothed = df.rolling(wd).mean().fillna(method='ffill')
        return smoothed
    else:
        df = pd.Series(df)
        smoothed = df.rolling(wd).mean()
        return smoothed

def makeMeanNonorm(newfiltered):
    meanDNA = newfiltered.groupby('seconds')['Intensity_DNA'].mean()
    stdDNA = newfiltered.groupby('seconds')['Intensity_DNA'].std(ddof=0)
    meanRPA = newfiltered.groupby('seconds')['Intensity_RPA'].mean()
    stdRPA = newfiltered.groupby('seconds')['Intensity_RPA'].std(ddof=0)
    return meanDNA,stdDNA, meanRPA, stdRPA

def getMeanTraj(df,x='seconds',y1='Intensity_DNA',y2='nucleotides'):
    meanDNA = df.groupby(x)[y1].mean()
    stdDNA = df.groupby(x)[y1].std(ddof=0)
    meanRPA = df.groupby(x)[y2].mean()
    stdRPA = df.groupby(x)[y2].std(ddof=0)
    return meanDNA,stdDNA, meanRPA, stdRPA

def analyse(name, file_path,changePointAnalysis = False, sigma = 1.):
    plt.switch_backend('Agg')
    pathparse = file_path.split('/')
    print(pathparse)
    Intfolder = pathparse[-2]
    filename = pathparse[-1][3:-4]
    print(filename)
    baseFolder = pathparse[-3]
    basePath = os.path.dirname(os.path.dirname(file_path))
    segFolder = basePath + '/segments'
    #assume metafile is in meta folder, same root-dir and called file_meta.csv:
    #remove the _series and the corrected bit!
    metaPath = basePath+'/meta/' + filename[10:-6] + '_meta.csv'
    print(metaPath)
    timeConversion = findTimeCov(metaPath)


    dataDNA = basePath + '/' + Intfolder + '/C1_' + filename + '.csv'
    dataRPA = basePath + '/' + Intfolder + '/C2_' + filename + '.csv'
    
    print('DNA file:',dataDNA)
    print('RPA file:',dataRPA)
    
    #read the data
    dfDNA = pd.read_csv(dataDNA)
    dfRPA = pd.read_csv(dataRPA)    
    #merge the two color channels
    df = pd.merge(dfDNA,dfRPA, on=['trajectory','slice','x','y'],suffixes = ["_DNA","_RPA"])
    print("merged data")
    #find events by DNA intensity
    dfFiltered,sortOutDNA = findEvents(df,Rscore=0.7, return_events=True)
    print('filtered by DNA int')
    #align trajectories by DNA unbinding:
    dfFilteredAligned = alignData(dfFiltered, timeConversion)
    print('align trajectories')
    #truncate: cut off data 50s after DNA unbinding:
    dfFilteredAligned = truncate(dfFilteredAligned,60)
    print('truncated')
    #filter out trajectories where RPA intensity does not increase:
    dfFilteredRPA,sortOutRPA = RPAfilter2(dfFilteredAligned)
    print('filtered by RPA')
    #calibrate trajectories: intensity_RPA -> number of RPAs -> nucleotides
    int_cal(dfFilteredRPA, 5000, assume_complete=True)
    print('calibrated')
    dfFilteredAligned = dfFilteredRPA

    

    trajFolder = basePath +'/trajectories/'+filename+'/'
    if not os.path.exists(trajFolder):
        os.makedirs(trajFolder)
    
  
    
    
    avgFolder = basePath +'/avgTraj/'
    if not os.path.exists(avgFolder):
        os.makedirs(avgFolder)
    mean_DNA, std_DNA, mean_RPA, std_RPA = getMeanTraj(dfFilteredAligned)

    
    plt.figure()
    color = ['green','blue','magenta','orange','red','black','grey']
    ax = sns.lineplot(data=mean_RPA,alpha = 0.2,color = color[0])
    ax2=ax.twinx()
    sns.lineplot(data=mean_RPA,color = color[0],ax=ax2)

    sns.lineplot(data=mean_DNA,ax=ax,color = 'black')

    ax.set_xlabel('aligned time [s]')
    ax.set_ylabel('nts synthesized [nt/s]')
    N = len(dfFilteredAligned['trajectory'].unique())
    ax.set_title('average trajectory, N={}'.format(N))
    lineRPA = ax.lines[0]
    x = lineRPA.get_xydata()[:,0]
    std1 = mean_RPA-std_RPA#ax.lines[1]
    std2 = mean_RPA+std_RPA
    ax2.fill_between(x,std1,std2, color='green', alpha=0.3)
    plt.savefig(avgFolder+filename+'.png')

    meanDNAData = pd.concat([mean_DNA,std_DNA],axis =1)#,columns = ['mean','std'])
    meanDNAData.columns = ['mean_int','std']
    meanRPAData = pd.concat([mean_RPA,std_RPA],axis =1)#,columns = ['mean','std'])
    meanRPAData.columns = ['mean_int','std']

#print(mean_DNA)
#print(mean_DNA.keys())
# meanDNAData = pd.merge(mean_DNA,std_DNA,on=['seconds'])
# meanRPAData = pd.merge(mean_RPA,std_RPA,on=['seconds'])
    meanDNAData.to_csv(avgFolder+filename+'_DNA.csv')
    meanRPAData.to_csv(avgFolder+filename+'_RPA.csv')
    return None
