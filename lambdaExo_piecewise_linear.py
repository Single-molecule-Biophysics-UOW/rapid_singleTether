# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 14:27:49 2021

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


dataframe = pd.read_csv(dataFolder+files[2])
dataframe.name = files[2]





def piecewise_linear(x,m,a,t):
    f1 = lambda x: m*x+t
    f2 = lambda x: m*a+t
    piecewise = np.piecewise(x,[x<=a],[f1,f2])
    return piecewise

def lin(x,m,t):
    return m*x+t
x = np.arange(0,100,0.1)

data = np.random.randn(1000) + lin(x,0.5,5)


fig, ax = plt.subplots(5)

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
    
    
#%%
#print(sortOutList)
dataframe_new = rst.sortOut_trajectories(dataframe, sortOutList)
for name,group in dataframe_new.groupby('trajectory'):
    print(group)
#%%
from synchronization import rst_tools_v2 as rst

plt.ioff()
def fit_traj(df,output='fity'):    
    y = np.array(df)
    x = np.arange(0,len(y),1)
    print(df)
    popt,pcov = curve_fit(piecewise_linear,x,y,p0=[300,300,100000])
    #print(x)
    new_y = piecewise_linear(x,*popt)
    changePoint = np.zeros(new_y.shape)
    
    #print(len(changePoint))
    #print(len(new_y))
    if int(np.round(popt[1])) >= len(changePoint):
        changePoint[-1] = 1
    else:
        changePoint[int(np.round(popt[1]))] = 1
    if output == 'fity':
        return np.array(new_y)
    elif output == 'changePoint':
        return changePoint
    

dataframe_new['regression'] = dataframe_new.groupby('trajectory').Intensity.transform(fit_traj,output = 'fity')

dataframe_new['changePoint'] = dataframe_new.groupby('trajectory').Intensity.transform(fit_traj,output = 'changePoint')


#%%
df_synchronized = rst.alignData(dataframe_new,1.,traj_column = 'trajectory',slice_column = 'slice',changePoint_column = 'changePoint',returnShift = False)

#%%
cp=[]
for name,group in df_synchronized.groupby('trajectory'):
    index = group['changePoint'][group['changePoint']>0].index[0][1]
    cp.append(index)
    
#%%
rates = np.array(rates)
rates=rates[rates>-25]
rates=rates[rates<0]
fig,ax = plt.subplots()
#ax.hist(cp,edgecolor='black')
ax.hist(rates,edgecolor='black',bins=30)
fig.show()
print(rates)
#%%

#for name,group in dataframe.groupby('trajectory'):
#    print(group[['Intensity','regression']])

#%%
df_synchronized.reset_index(inplace=True)
#%%make mean
mean_df = pd.DataFrame()


mean_df['medianDNA'] = df_synchronized.groupby('seconds')['Intensity'].median()
mean_df['meanDNA'] = df_synchronized.groupby('seconds')['Intensity'].mean()
mean_df['semDNA'] = df_synchronized.groupby('seconds')['Intensity'].sem()


mean_df = mean_df[(mean_df.index>-800) & (mean_df.index<0)]

offset = (mean_df['meanDNA'][(mean_df.index>-100) & (mean_df.index<0)]).mean()

mean_df['normMean'] = ((mean_df['meanDNA']- offset)/(mean_df['meanDNA']- offset).max())*2600


linRegion = mean_df[(mean_df.index > -800) & (mean_df.index <= 0)]
print(linRegion)
popt_mean,pcov_mean = curve_fit(lin,mean_df['normMean'].index*0.2,mean_df['normMean'])




#%%plot it
fig, ax = plt.subplots()


#mean_df.reset_index(inplace=True)
x_fit = np.arange(-160,-20,0.1)

# ax[0].fill_between(mean_df['seconds']*0.2,mean_df['meanDNA']-mean_df['semDNA'],mean_df['meanDNA']+mean_df['semDNA'],label='sem',alpha=0.2)
#ax.plot(mean_df['normMean'],label='median',color='black')
ax.plot(mean_df.index*0.2,mean_df['normMean'],label='median',color='black')
#ax[0].plot(linRegion,label='lin')
#ax.plot(x_fit,lin(x_fit,*popt_mean),label = 'linear fit, m={}'.format(popt[0]))
ax.legend()
#print(mean_df['seconds']*0.2)
#ax[0].legend()
#print(mean_df['meanDNA'].diff().fillna(0))
#ax[1].plot(mean_df['seconds']*0.2,mean_df['meanDNA'].diff().fillna(0))

fig.show()


#%%
#%%
print(mean_df)
#%% save data:

mean_df.to_csv(dataFolder+'UvrD_synchMean.csv')


#%%

directory = r"K:\Tirf\Lisannes_microscope\211215\Lambda\events\\" + dataframe.name +"\\"
if not os.path.exists(directory):
    os.makedirs(directory)
    
        #print(events['score'])
rst.plot_all_trajs(dataframe_new,directory,xcolumn = 'seconds',ycolumn = ['Intensity','regression'])    
#%%
print(dataframe['regression'])

    
#%%
fig,ax = plt.subplots()
print(rates)
rates = np.abs(np.array(rates)/5.)
#rates=rates[rates<100]
confidence  = np.array(confidence)
confidence = confidence[confidence< 10]
print(len(confidence))
#ax.hist(np.abs(confidence),edgecolor='black',bins =25)

ax.hist(np.abs(rates),edgecolor='black',label = 'n={}'.format(len(rates)),bins=25)
ax.set_xlabel('Rate (nt\s)')
ax.set_ylabel('# of occurance')
#ax.set_xlim([0,35])
ax.legend()
fig.show()
fig.savefig('UvrD_rate.svg')

    
