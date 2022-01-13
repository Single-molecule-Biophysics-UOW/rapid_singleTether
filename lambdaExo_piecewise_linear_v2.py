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
from scipy.signal import savgol_filter
#%%

dataFolder = r"C:\Users\shm975\Documents\tempData\211216\lambda\integration/"


files = os.listdir(dataFolder)
print(files)


#%% load data
dataframe = pd.read_csv(dataFolder+files[1])
dataframe.name = files[1]

print(dataframe.name)

#%%
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

dataframe.reset_index(inplace=True)


 #%% normalization to basepairs
#Start intensity is 2600 bp dsDNA, end is 2600 ssDNA
#use fitted values to do normalization


def complex_norm(df,smooth = False):
    t = df.index
    Intensity = df
    if smooth:
        Intensity = savgol_filter(Intensity,7,3)

    # print(t)
    #first do the fit and return popt:
    try:
        popt,pcov= curve_fit(piecewise_linear_triple, t,Intensity,p0=[-200,30,300,np.max(Intensity)-1],bounds=((-np.inf,np.min(t),np.min(t),0),(np.inf,np.max(t),np.max(t),np.max(Intensity))),maxfev=10000)
    
    except Exception:
        print(Exception)
        # print(t)
        # print(Intensity)
        # plt.plot(t,Intensity)
        # plt.plot(t,piecewise_linear_triple(t,-20,30,300,100000))
        # plt.show()
        norm_Intensity = np.ones(len(t))*np.nan
        return norm_Intensity
    # intensity when time <popt[1]: 2600 bp
    # intensity when time >popt[2]: 0 bp
    maxValue = popt[0]*popt[1]+popt[3]
    offsetValue = popt[0]*popt[2]+popt[3]


    #now normalize:
    norm_Intensity = (Intensity - offsetValue)/(maxValue-offsetValue)
    
    return norm_Intensity

#first make sure seconds are passed to norm function by setting them as index:
try:
    dataframe.set_index('seconds',inplace=True)
except KeyError:
    pass
dataframe['norm_Intensity'] = dataframe.groupby('trajectory')['Intensity'].transform(complex_norm,smooth = True)
dataframe['nts'] = dataframe['norm_Intensity']*2620
try:
    dataframe.reset_index(inplace=True)
except KeyError:
    #already reset.
    pass




#%% fit every trajectory again to get rates:


def piecewise_fit(df,output = 'predict',function = piecewise_linear):
    #print(df.head(20))
    #check for nans first:
    x,y = np.array(df.index.get_level_values(0)),np.array(df)
    if np.isnan(np.sum(y)):
        if output == 'predict':
            return np.ones(len(x))*np.nan
        if output == 'score':
            return [np.nan,np.nan,np.nan]
        if output == 'slope':        
            return np.ones(len(x))*np.nan
    if function == piecewise_linear:
        p0 = [-10,10,0]
    else:
        p0 = p0=[-10,10,200,0]
    
    try:
        popt,pcov = curve_fit(function,x,y,p0=p0)
    except RuntimeError:
        if function == piecewise_linear_triple:
            popt = [np.nan,np.nan,np.nan,np.nan]
        else:
            popt = [np.nan,np.nan,np.nan]
        fity = np.ones(len(x))*np.nan
    fity = function(x,*popt)
    #d = {'fity': fity, 'score': np.ones(fity.shape)*score}
    #regrdf = pd.DataFrame(data=d)

    if output == 'predict':
        return fity
    if output == 'score':
        return pcov
    if output == 'slope':        
        return np.ones(len(x))*popt[0]
    if output == 'a':        
        return np.ones(len(x))*popt[1]
    if output == 'b':        
        return np.ones(len(x))*popt[2]
#dataframe.reset_index(inplace=True)    
try:
    dataframe.set_index('seconds',inplace=True)
except KeyError:
    pass

dataframe['piecewise'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,output='predict')
dataframe['piecewise_triple'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=piecewise_linear_triple,output='predict')
dataframe['rate'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=piecewise_linear_triple,output='slope')
dataframe['a'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=piecewise_linear_triple,output='a')
dataframe['b'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=piecewise_linear_triple,output='b')


#%% plot all events:
plt.ioff()
try:
    dataframe.reset_index(inplace=True)
except ValueError:
    #already reset index
    pass
eventFolder = r"C:\Users\shm975\Documents\tempData\211216\lambda\events_006/"

rst.plot_all_trajs(dataframe, eventFolder, xcolumn = 'seconds',ycolumn = ['nts','piecewise_triple','rate'])

#%% save the normalized data:
norm_dataFolder = r"C:\Users\shm975\Documents\tempData\211216\lambda\norm_integration/"

if os.path.isdir(norm_dataFolder) == False:
    os.makedirs(norm_dataFolder)
try:
    dataframe.reset_index(inplace=True)
    print(dataframe.keys())
except Exception:
    pass
dataframe.to_csv(norm_dataFolder+'norm_'+dataframe.name, index=False)


#%%
try:
    dataframe.reset_index(inplace = True)
except ValueError:
    pass
print(dataframe.keys())
#%%
rates = []
for name,group in dataframe.groupby('trajectory'):
    maxT = group['seconds'].max()
    if (group['a'].unique()[0] < 0) or  (group['b'].unique()[0] < 0):
        print('sort out traj {}, negative a'.format(name))
        continue
    if (group['a'].unique()[0] > maxT) or  (group['b'].unique()[0] > maxT):
        print('sort out traj {}, too big params'.format(name))
    rates.append(-group['rate'].unique()[0])
#rates = (-dataframe['rate'].unique())
rates = np.array(rates)
#%%
rates = rates[rates > 0]
rates = rates[rates<100]
plt.figure()
plt.hist(rates,edgecolor = 'black',label='n={}'.format(len(rates)),bins = int(np.sqrt(len(rates))))
#plt.hist(rates2,facecolor='gold',alpha=0.5,edgecolor = 'black',label='n={}'.format(len(rates2)),bins = int(np.sqrt(len(rates2))))
plt.xlabel('time (s)')
plt.ylabel('count')
plt.legend()
plt.show()
#%% save the rates:
print(dataframe.name)
#%%
    
df_rates = pd.DataFrame({'rate':rates})

folder = r"C:\Users\shm975\Documents\OneDrive - University of Wollongong\phi29Paper\data\lambda\all_rates/"

df_rates.to_csv(folder+"new_rates_006.csv",index=False)

#%%
# #dataframe.reset_index(inplace=True)    
# rates = []
# rates_old = []
# for name,group in dataframe.groupby('trajectory'):
#    # group.set_index('seconds',inplace=True)
#     slope = (piecewise_fit(group['Intensity'],function=piecewise_linear_triple,output='slope'))
#     slope2 = (piecewise_fit(group['Intensity'],function=piecewise_linear,output='slope'))
#     rates.append(-slope)
#     rates_old.append(-slope2)
# print(rates)
# #%%
# rates = np.array(rates)
# rates = rates[rates>0]
# #rates = rates[rates<5000]

# rates2 = np.array(rates_old)
# rates2 = rates2[rates2>0]
# #rates2 = rates2[rates2<5000]

# plt.figure()
# plt.hist(rates,edgecolor = 'black',label='n={}'.format(len(rates)),bins = int(np.sqrt(len(rates))))
# plt.hist(rates2,facecolor='gold',alpha=0.5,edgecolor = 'black',label='n={}'.format(len(rates2)),bins = int(np.sqrt(len(rates2))))
# plt.xlabel('time (s)')
# plt.ylabel('count')
# plt.legend()
# plt.show()

# #%%
# path=r"C:\Users\shm975\Documents\tempData\211221\lambda/events/"
# dataframe.reset_index(inplace=True)
# eventFolder = os.path.isdir(path)
# if not eventFolder:
#     os.makedirs(path)
# rst.plot_all_trajs(dataframe, path,xcolumn='seconds',ycolumn = ['Intensity','piecewise','piecewise_triple'])

# #%%
# n=0
# rates = []
# confidence = []
# sortOutList = []
# for name,group in dataframe.groupby('trajectory'):
    
#     #ax[n].plot(group['slice'],group['Intensity'])    
#     popt_uncal,pcov_uncal = curve_fit(piecewise_linear,np.array(group['slice']),np.array(group['Intensity']),p0=[300,300,100000])
#     print(popt_uncal)
#     #assume it goes the whole way, then 2600bp equals f(0) - f(a)
#     f0 =piecewise_linear(0,*popt_uncal)
#     fa = piecewise_linear(popt_uncal[1],*popt_uncal)
#     calib_intensity = np.array((group['Intensity']*2600 / (f0-fa)))
#     calib_time = np.array(group['slice']*1)
    
    
    
#     popt,pcov = curve_fit(piecewise_linear,calib_time,np.array(calib_intensity),p0=[300,30,10000])
#     print(np.sqrt(pcov[0,0]))
#     if np.sqrt(pcov[0,0]) > 100:
#         sortOutList.append(name)
#         continue
    
#     #x = np.linspace(0,group['slice'].max(),10000)
#     #print(x)
#     #ax[n].plot(x,piecewise_linear(x,*popt), label = 'slope: {:.2f}, offset: {:.2f}, a={:.2f}'.format(popt[0],popt[1],popt[2]))
#     #ax[n].legend()
#     confidence.append(np.sqrt(pcov[0,0]))
#     rates.append(popt[0])
#     n+=1