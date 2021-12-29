# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 16:04:40 2021

@author: shm975
"""




import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from synchronization import rst_tools_v2 as rst
#%%

dataFolder = r"C:\Users\shm975\Documents\tempData\211221\lambda/integration/"

files = os.listdir(dataFolder)
print(files)
#%%
dataframe = pd.read_csv(dataFolder+files[0])
dataframe.name = files[0]

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
    piecewise = np.piecewise(x,(x<a,(x>=a)&(x<=b),x>b),[f0,f1,f2])
    return piecewise

def lin(x,m,t):
    return m*x+t
x = np.arange(0,500,0.1)

testData = dataframe[dataframe['trajectory']==1]
popt,pcov= curve_fit(piecewise_linear_triple, np.array(testData['seconds']),np.array(testData['Intensity']),p0=[-1000,30,300,400000])

plt.figure()
plt.plot(testData['seconds'],testData['Intensity'])
plt.plot(x,piecewise_linear_triple(x,*popt))
plt.show()
                                      

#%% fit every trajectory:

for name,group in dataframe.groupby('trajectory'):
    print(group['Intensity'])
    print(group['seconds'])
    break
#%%

def piecewise_fit(df,output = 'predict',function = piecewise_linear):
    #print(df.head(20))
    if function == piecewise_linear:
        p0 = [-300,300,100000]
    else:
        p0 = p0=[-1000,30,300,400000]
    x,y = np.array(df.index.get_level_values(0)),np.array(df)
    try:
        popt,pcov = curve_fit(function,x,y,p0=p0)
    except RuntimeError:
        if function == piecewise_linear_triple:
            popt = [1,0,300,100000]
    fity = function(x,*popt)
    #d = {'fity': fity, 'score': np.ones(fity.shape)*score}
    #regrdf = pd.DataFrame(data=d)

    if output == 'predict':
        return fity
    if output == 'score':
        return pcov
    if output == 'slope':        
        return popt[0]
#dataframe.reset_index(inplace=True)    
dataframe.set_index('seconds',inplace=True)

dataframe['piecewise'] = dataframe.groupby('trajectory')['Intensity'].transform(piecewise_fit,output='predict')
dataframe['piecewise_triple'] = dataframe.groupby('trajectory')['Intensity'].transform(piecewise_fit,function=piecewise_linear_triple,output='predict')
#%%
#dataframe.reset_index(inplace=True)    
rates = []
rates_old = []
for name,group in dataframe.groupby('trajectory'):
    group.set_index('seconds',inplace=True)
    slope = (piecewise_fit(group['Intensity'],function=piecewise_linear_triple,output='slope'))
    slope2 = (piecewise_fit(group['Intensity'],function=piecewise_linear,output='slope'))
    rates.append(-slope)
    rates_old.append(-slope2)
print(rates)
#%%
rates = np.array(rates)
rates = rates[rates>0]
rates = rates[rates<5000]

rates2 = np.array(rates_old)
rates2 = rates2[rates2>0]
rates2 = rates2[rates2<5000]

plt.figure()
plt.hist(rates,edgecolor = 'black',label='n={}'.format(len(rates)),bins = int(np.sqrt(len(rates))))
plt.hist(rates2,facecolor='gold',alpha=0.5,edgecolor = 'black',label='n={}'.format(len(rates2)),bins = int(np.sqrt(len(rates2))))
plt.xlabel('time (s)')
plt.ylabel('count')
plt.legend()
plt.show()

#%%
path=r"C:\Users\shm975\Documents\tempData\211221\lambda/events/"
dataframe.reset_index(inplace=True)
eventFolder = os.path.isdir(path)
if not eventFolder:
    os.makedirs(path)
rst.plot_all_trajs(dataframe, path,xcolumn='seconds',ycolumn = ['Intensity','piecewise','piecewise_triple'])

#%%
n=0
rates = []
confidence = []
sortOutList = []
for name,group in dataframe.groupby('trajectory'):
    
    #ax[n].plot(group['slice'],group['Intensity'])    
    popt_uncal,pcov_uncal = curve_fit(piecewise_linear,np.array(group['slice']),np.array(group['Intensity']),p0=[300,300,100000])
    print(popt_uncal)
    #assume it goes the whole way, then 2600bp equals f(0) - f(a)
    f0 =piecewise_linear(0,*popt_uncal)
    fa = piecewise_linear(popt_uncal[1],*popt_uncal)
    calib_intensity = np.array((group['Intensity']*2600 / (f0-fa)))
    calib_time = np.array(group['slice']*1)
    
    
    
    popt,pcov = curve_fit(piecewise_linear,calib_time,np.array(calib_intensity),p0=[300,30,10000])
    print(np.sqrt(pcov[0,0]))
    if np.sqrt(pcov[0,0]) > 100:
        sortOutList.append(name)
        continue
    
    #x = np.linspace(0,group['slice'].max(),10000)
    #print(x)
    #ax[n].plot(x,piecewise_linear(x,*popt), label = 'slope: {:.2f}, offset: {:.2f}, a={:.2f}'.format(popt[0],popt[1],popt[2]))
    #ax[n].legend()
    confidence.append(np.sqrt(pcov[0,0]))
    rates.append(popt[0])
    n+=1