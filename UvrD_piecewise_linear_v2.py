# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 13:05:32 2022

@author: stefa
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from synchronization import rst_tools_v2 as rst
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter




#%% C:\Users\stefa\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unwinding\211212\integration
folder = r'C:\Users\shm975\Documents\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unwinding\211214/integration/'
#folder = r'C:/Users\stefa\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unwinding\211214\integration/'

files = os.listdir(folder)
files = [x for x in files if 'csv' in x]
print(files)
data = []
for i in files:
    df = pd.read_csv(folder+i)
    df.name = i
    data.append(df)
dataframe = data[0]    
#%%#%% normalization to basepairs
#Start intensity is 2620 bp dsDNA, end is 2620 ssDNA
#use fitted values to do normalization


def complex_norm(df, smooth = False):
    t = df.index
    Intensity = df
    if smooth:
        Intensity = savgol_filter(Intensity,7,3)
    # print(t)
    #first do the fit and return popt:
    try:
        popt,pcov= curve_fit(rst.piecewise_linear_triple, t,Intensity,p0=[-200,30,300,np.max(Intensity-1)],
                             bounds=((-np.inf,np.min(t),np.min(t),0),(np.inf,np.max(t),np.max(t),np.max(Intensity))),maxfev=10000)
    
    except RuntimeError:
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
    norm_Intensity = (Intensity - offsetValue)/(maxValue- offsetValue)

    #now normalize:
    norm_Intensity = (Intensity - offsetValue)/(maxValue-offsetValue)
    return norm_Intensity

#first make sure seconds are passed to norm function by setting them as index:
try:
    dataframe.set_index('seconds',inplace=True)
except KeyError:
    pass
dataframe['norm_Intensity'] = dataframe.groupby('trajectory')['Intensity'].transform(complex_norm,smooth=True)
dataframe['nts'] = dataframe['norm_Intensity']*2620
try:
    dataframe.reset_index(inplace=True)
except KeyError:
    #already reset.
    pass
#%% plot some examples nicley

options = (dataframe['trajectory'].unique())    




grouped = dataframe.groupby('trajectory')
good = [1009]#[158,1009,406]
choices = good#list(np.random.choice(options,size=10))


n=0

for i in choices:
    print('choice:',i)
    traj = grouped.get_group(i)
    maxI = traj['Intensity'].max()
    maxT, minT = traj['seconds'].max(), traj['seconds'].min()
    popt,pcov= curve_fit(rst.piecewise_linear_triple, np.array(traj['seconds']),np.array(traj['Intensity']),p0=[-200,30,300,maxI],
                             bounds=((-np.inf,minT,minT,0),(np.inf,maxT,maxT,maxI)),maxfev=10000)
    fig,ax = plt.subplots()
    x = np.arange(-10,200,0.1)
    
    ax.plot(traj['seconds'], traj['Intensity'], label = 'traj {}'.format(i),color = 'black')
    ax.plot(traj['seconds'],rst.piecewise_linear_triple(np.array(traj['seconds']),*popt),label = 'm={:.2f},a={:.2f},b={:.2f},o={:.2f}'.format(*popt),color ='red')
    ax.plot(x,(popt[0]*x)+popt[3],label = 'f$_1$(x) = mx+o'.format(popt[0],popt[3]),color ='red',ls = '--')
    ax.plot(x,(popt[0]*popt[1]*np.ones(x.shape))+popt[3],label = 'f$_2$(x) = ma+o'.format(popt[0],popt[3]),color ='red',ls = '-.')
    ax.plot(x,(popt[0]*popt[2]*np.ones(x.shape))+popt[3],label = 'f$_3$(x) = mb+o'.format(popt[0],popt[3]),color ='red',ls = ':')
    #ax.plot(traj['raw time'], traj['Intensity_DNA'], label = 'traj {}'.format(i))
    ax.legend()
    ax.set_xlim([0,200])
    n+=1
    
    #fig.savefig('UwrD_unwind_normalisation_figure.svg')
plt.show()


#%% fit every trajectory again to get rates:


def piecewise_fit(df,output = 'predict',function = rst.piecewise_linear):
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
    if function == rst.piecewise_linear:
        p0 = [-10,10,0]
    else:
        p0 = p0=[-10,10,200,0]
    
    try:
        popt,pcov = curve_fit(function,x,y,p0=p0)
    except RuntimeError:
        if function == rst.piecewise_linear_triple:
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


dataframe['piecewise_triple'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=rst.piecewise_linear_triple,output='predict')
print('finished piecewisetriple')
dataframe['rate'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=rst.piecewise_linear_triple,output='slope')
print('finished rate')
dataframe['a'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=rst.piecewise_linear_triple,output='a')
dataframe['b'] = dataframe.groupby('trajectory')['nts'].transform(piecewise_fit,function=rst.piecewise_linear_triple,output='b')



#%% plot all events:
plt.ioff()
try:
    dataframe.reset_index(inplace=True)
except ValueError:
    #already reset index
    pass

eventFolder = r'C:\Users\shm975\Documents\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unwinding\211212/events_complete/'
rst.plot_complete_trajs(dataframe, eventFolder, xcolumn = 'seconds',ycolumn = ['nts','piecewise_triple','rate','a','b'])

#%% save the normalized data:
norm_dataFolder = r'C:/Users\stefa\OneDrive - University of Wollongong\phi29Paper\data\UvrD_unwinding\211214\norm_integration/'

if os.path.isdir(norm_dataFolder) == False:
    print('make folder')
    os.makedirs(norm_dataFolder)
try:
    dataframe.reset_index(inplace=True)
    print(dataframe.keys())
except Exception:
    pass
dataframe.to_csv(norm_dataFolder+'norm_'+dataframe.name, index=False)


#%%
try:
    dataframe.reset_index(inplace=True)
except ValueError:
    pass
rates = []
for name,group in dataframe.groupby('trajectory'):
    
    maxT = group['seconds'].max()
    if (group['a'].unique()[0] < 0) or  (group['b'].unique()[0] < 0):
        print('sort out traj {}, negative a'.format(name))
        continue
    if (group['a'].unique()[0] > maxT) or  (group['b'].unique()[0] > maxT):
        print('sort out traj {}, too big params'.format(name))
        continue
    rates.append(-group['rate'].unique()[0])
#rates = (-dataframe['rate'].unique())
rates = np.array(rates)

#%%
rates = rates[rates > 0]
rates = rates[rates<80]
ratesToFit = rates[rates<80]
fig,ax = plt.subplots()
ax.hist(rates,edgecolor = 'black',facecolor='grey',label='n={}'.format(len(rates)),density=True,bins=25)
#plt.hist(rates2,facecolor='gold',alpha=0.5,edgecolor = 'black',label='n={}'.format(len(rates2)),bins = int(np.sqrt(len(rates2))))
plt.xlabel('time (s)')
plt.ylabel('count')
plt.legend()
plt.show()
#fig.savefig('UvrD_unwinding_hist_50uM.svg')


from sklearn.mixture import GaussianMixture

from sklearn.mixture import BayesianGaussianMixture

#model = GaussianMixture(2,tol=10e-5,covariance_type='diag',max_iter=1000, n_init = 10,init_params='random', weights_init=[0.1,0.9]).fit(np.array(mixed_random).reshape(-1,1))
model = GaussianMixture(2,tol=10e-5,covariance_type='full',weights_init=[0.3,0.7],means_init = np.array([10,40]).reshape(-1,1)).fit(np.array(ratesToFit).reshape(-1,1))



# plot resulting fit
x_range = np.arange(0,70,0.1)
pdf = np.exp(model.score_samples(x_range.reshape(-1, 1)))
responsibilities = model.predict_proba(x_range.reshape(-1, 1))
pdf_individual = responsibilities * pdf[:, np.newaxis]
print(model.means_)
print(model.covariances_.shape)

cov = model.covariances_

print(np.diag(model.covariances_[1]))
print(model.predict(np.array(ratesToFit).reshape(-1,1)))
print([ np.sqrt(  np.trace(cov[i])/2) for i in range(0,2) ])
print('sem1 = {}'.format(np.sqrt(model.covariances_[0,0])))
print('sem2 = {}'.format(np.sqrt(model.covariances_[1,0])))

print(model.weights_)
#plt.hist(logged_counts, bins='auto', density=True, histtype='stepfilled', alpha=0.5)
ax.plot(x_range, pdf, label='Mixture diag',color='black')
ax.plot(x_range, pdf_individual,ls='--', label='Components',color='black')

# ax.plot(x_range, pdf_b, label='Mixture diag',color='orange')
# ax.plot(x_range, pdf_individual_b,ls='--', label='Components',color='orange')


# ax.plot(x_range, pdf2, label='Mixture diag',color='black')
# ax.plot(x_range, pdf_individual2,ls='--', label='Components',color='black')


# fig.savefig('dg_uvrD_unwinding.svg')
plt.show()
