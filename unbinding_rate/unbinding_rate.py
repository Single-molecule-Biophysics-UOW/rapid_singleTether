# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 16:03:08 2021

@author: shm975
"""


import os
import pandas as pd
from synchronization import rst_tools_v2 as rst
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import expon
from scipy.stats import gamma
from scipy.stats import laplace
from scipy.stats import gumbel_r
from scipy.stats import lomax
from scipy.stats import rv_histogram
#%%
dataFolder = r"C:\Users\stefa\OneDrive\Documents\tempData\UvrdData\UvrD_data\5'BioBlunt\integration\\"
files = os.listdir(dataFolder)
print(files)
#%%
all_data=[]
for f in files:
    if f[-3:]!='csv':
        continue
    dataframe = pd.read_csv(dataFolder+f)
    dataframe.name = f
    all_data.append( dataframe)

#%%first lets see if the DNA intensity is reproducable between experiments. To do so look at the intensities at slice 0!
import itertools
slice0_Intensity = pd.DataFrame()
fig,ax = plt.subplots(3,3,sharex=True)
indeces = list(itertools.product([0,1,2],repeat =2))
n=0
for data in all_data:
    print(data)
    doc_string = data.name.split('_')
    print(doc_string)
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    mean = (data['Intensity'][data['slice']==0]).mean()
    ax[indeces[n]].hist(data['Intensity'][data['slice'] ==3],edgecolor='black',label = 'OD:'+OD+' Date:'+ DATE +'mean: {:.2E}'.format(mean))
    ax[indeces[n]].legend()
    n+=1

plt.show()
#%% now find all molecules that dissociate in one step:
events_all_data = []
N = len(all_data)
n = 0
for data in all_data:
    events,sortOut,sortOutReason = rst.find_unbinding_v2(data, 'Intensity', 0.7, 80000)
    events.name = data.name
    events_all_data.append(events)
    print('finished {} out of {}'.format(n+1,N))
    n+=1
    
#%%plot events:
plt.ioff()
for events in events_all_data:
    directory = r"C:\Users\stefa\OneDrive\Documents\tempData\UvrdData\UvrD_data\5'biossDNA\events\\" + events.name +"\\"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
        #print(events['score'])
    rst.plot_all_trajs(events,directory,xcolumn = 'seconds',ycolumn = ['Intensity','regression','derivative'])

#%%synchronize to the first unbinding event


def corr_time(df):
    #find first molecule that dissociates and set that tio time-point zero:
    earliest_unbinding = df['seconds'].max()
    
    for name, group in df.groupby('trajectory'):
        timing = group[group['regression'].diff().fillna(0) !=0]
        print('timing:',timing['seconds'].iloc[0])
        print('earliest unbinidng:',earliest_unbinding)
        if timing['seconds'].iloc[0] < earliest_unbinding:
            earliest_unbinding = timing['seconds'].iloc[0]
            print('new earliest:',timing['seconds'].iloc[0])
    #global correction: seconds - earliest unbiding:
    df['synch_time'] = df['seconds'] - earliest_unbinding

for events in events_all_data:    
    corr_time(events)
#%%
print(events_all_data[0][['seconds','synch_time']])


#%% and now get the time off dissociation for every molecule:
timing_all_data = []
for events in events_all_data:
    timing = rst.unbinding_timing(events,groupName = 'trajectory',column = 'synch_time')
    timing.name = events.name
    timing_all_data.append(timing)


print(timing_all_data[1].min())


    
timing_all_data_old = []
for events in events_all_data:
    timing = rst.unbinding_timing(events,groupName = 'trajectory',column = 'seconds')
    timing.name = events.name
    timing_all_data_old.append(timing)
print(timing_all_data_old[1].min())    
#%%
print(timing_all_data[0].keys())

#%% plot efficiencies as bar plot:
efficiency = []
xlabels = []
n = 0



for events,data in zip(events_all_data,all_data): 
   # data.reset_index(inplace=True)
    #events.reset_index(inplace=True)
    efficiency.append(len(events['trajectory'].unique())/len(data['trajectory'].unique()))
    doc_string = data.name.split('_')
    
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    print(conc)
    print(RPA)
    xlabels.append(str(conc)+'_'+str(n))
    n+=1

fig,ax = plt.subplots()
ax.bar(xlabels,efficiency,color = ['grey'])

print(len(data['trajectory'].unique()))
print(len(events['trajectory'].unique()))
#%% now concatenate data with same conditions:
all_timing_corrected = {'5uM':[],'8uM': [], '10uM': [], '20uM': [], '50uM': [], '5000uM': []}

for time in timing_all_data:
    #time['corr_sec'] = time['seconds']-time['seconds'].min()
    doc_string = time.name.split('_')    
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    
    key = conc + 'uM'
    
    data = all_timing_corrected[key]
    print(time['synch_time'].min())
    
    new_data = data + list(time['synch_time'])
    all_timing_corrected.update({key:new_data})

    
#%%plot number of molecules visible:


   

print (events_all_data[0].keys())

#%%
#loop through all timepoints and count how many disTimes are bigger than this time point:
N_molecules_all = []
c=[]
count = 0
for events,timing in zip(events_all_data,timing_all_data):
    #we need the corrected time for all events:
    doc_string = events.name.split('_')
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    c.append(conc)
    print(events.keys())
    timepoints = events['synch_time'].unique()
    timepoints =timepoints[timepoints>=0]
    print(timing.keys())
    N_molecules = []
    print('N0:',len(timing))
    N_molecules.append(len(timing))
    print('N_molecules:',N_molecules)
    print('loop:',count)
    count+=1
    for t in timepoints:
        molecules = timing[timing['synch_time']>t]
        N_molecules.append(len(molecules))
    N_molecules_all.append(np.array(N_molecules))
#%%
print(N_molecules_all[0])
    
fig,ax = plt.subplots()
for n,conc in zip(range(len(N_molecules_all)),c):
    ax.plot(N_molecules_all[n]/N_molecules_all[n].max(),label=conc)
ax.legend()
fig.show()
    
    
    
#%%
from unbinding_rate.hypoexp import hypoexp
from unbinding_rate.statsmodel_hypoExp import hypoexpon
from scipy.optimize import curve_fit 

x = np.arange(1,250,0.1)
print(all_timing_corrected.keys())
fig2,ax2 = plt.subplots(1)
ax2 = [ax2]
n=0
for key in ['20uM']:
    try:
        m,bins,patches = ax2[n].hist(all_timing_corrected[key],weights=all_timing_corrected[key],edgecolor='black',label = '{}, n={}'.format(key, len(all_timing_corrected[key])),density=True,cumulative = False,bins = int(np.sqrt(len(all_timing_corrected[key]))))
    except ValueError:
        continue
    
    
    mean = np.mean(all_timing_corrected[key])
    std = np.std(all_timing_corrected[key])
    cV = std/mean
    
    #rate1 = (2/mean)* (1+np.sqrt(1+2*(cV**2 -1)))**(-1)
    #rate2 = (2/mean)* (1-np.sqrt(1+2*(cV**2 -1)))**(-1)
    #try:
    #X2 = np.sort(all_timing_corrected[key])
    #F2 = np.array(range(len(all_timing_corrected[key])))/float(len(all_timing_corrected[key]))
    print(np.min(all_timing_corrected[key]))
    params = hypoexp.fit(np.array(all_timing_corrected[key])+1,0.1,0.02,floc=0, fscale = 1)
    model = hypoexpon(all_timing_corrected[key])
    cons = ({'type': 'ineq',
         'fun' : lambda x: np.all(x)},
            {'type': 'ineq',
         'fun' : lambda x: x[0]-x[1]})
    

    p1 = 0.108
    p2 = 0.016
    #results = model.fit(start_params= [p1,p2],maxiter=100000, maxfunc= 10000)
    params = hypoexp.fit(np.array(all_timing_corrected[key])+0.1,p1,p2,floc=0, fscale = 1)
    print('succesfull at: {} and {}'.format(p1,p2))
        
    results = model.fit(start_params= [p1,p2],maxiter=100000, maxfunc= 100000, constraints=cons)
    print('params:',results.params)

        
    #print(params)
    
    #popt,pcov = curve_fit(hypoExp,X2,F2,bounds=[0.001,1.],method='dogbox')
    #print(popt)
    #print(pcov)
    ax2[n].plot(x,hypoexp.pdf(x,*params),label = 'scipy, p1 ={:.2}, p2={:.2}'.format(params[0],params[1]))
    ax2[n].plot(x,hypoexp.pdf(x,0.019,0.108),label = 'scipy, p1 ={:.2}, p2={:.2}'.format(0.019,0.108))
    ax2[n].plot(x,hypoexp.pdf(x,0.016,0.108),label = 'scipy, p1 ={:.2}, p2={:.2}'.format(0.016,0.108))
    ax2[n].plot(x,model._pdf(x,results.params[0],results.params[1]),label = 'statsmodells, p1 ={:.2}, p2={:.2}'.format(results.params[0],results.params[1]))
    ax2[n].legend()
    #except:
    #    pass
    #print('{} : rate1:{:.2}, rate2: {:.2}, c$_V$={}'.format(key,rate1,rate2,cV))
    n+=1
fig2.show()
    
#%%



#%%    
for key in['5uM','8uM','10uM','20uM','50uM']:
    plt.figure()
    plt.hist(all_timing_corrected[key],edgecolor='black',label = 'c = {}, n={}'.format(key,len(all_timing_corrected[key])))
    plt.legend()
    plt.show()


#%% and now plot histograms of the unbinding timing and apply hmm
from scipy.optimize import curve_fit
import toolbox as tb
from scipy.special import gamma as gf
def hypoExp(x,l1,l2):
    #if l1 ==l2:
    #    l2 = l2 +1E-8
    if l1 == l2:
        #in this case we look at a gamma-distribution with k=2!
        a=2
        b = l1
        nom = b**a * x *np.exp(-b*x)
        denom = gf(a)       #gamma function
        f = nom/denom
    else:
        C = l1/(l2-l1)
        f = C *( np.exp(-l1*x) - np.exp(-l2*x) )
    return f


fig,ax = plt.subplots(3,3,sharex=True)
indeces = list(itertools.product([0,1,2],repeat =2))
fig2,ax2 = plt.subplots(5,1,sharex = True)

n=0
for time in timing_all_data:
    time['corr_sec'] = time['seconds']-time['seconds'].min()
    doc_string = time.name.split('_')    
    conc = "".join(filter(str.isdigit, doc_string[5]))
    RPA =  "".join(filter(str.isdigit, doc_string[4]))
    OD =   "".join(filter(str.isdigit, doc_string[6]))
    DATE = "".join(filter(str.isdigit, doc_string[7]))
    
    # if conc != '5':
        # continue
    
    
    # histDataX,histDataY = tb.hist(time['seconds']-time['seconds'].min(),bins = 20,density = True)
    # try:
    #     popt, pcov = curve_fit(hypoExp,histDataX,histDataY,p0=[0.001,0.05],max_nfev=10000000)
    # except:
    #     pass
    #     #popt = [0.02,0.5]
    x = np.arange(1,250,0.2)
    params = gamma.fit(time['corr_sec'],f0=2)
    print(params)
    ax[indeces[n]].plot(x,gamma.pdf(x,*params))
    ax[indeces[n]].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
    equation = r" f(x) = $e^{-(e^{\frac{(-x-\mu)}{\beta}})}$"
    #ax2.plot(x,hypoExp(x,*popt))
    #ax2.plot(x,hypoExp(x, 0.001,0.01),color = 'red',label = '0.001, 0.01')
    #ax2.plot(x,hypoExp(x, 0.001,0.05),color = 'blue',label = '0.0007, 0.05')
    #ax2.plot(x,hypoExp(x, 0.001,0.09),color = 'yellow',label = '0.0005, 0.09')
    #ax2.legend()
    #ax2.hist(time['seconds']-time['seconds'].min(),bins = 30,edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' $\mu$={:.1}, $\beta$={:.1}'.format(popt[0],popt[1]),density=True)
    #ax[indeces[n]].plot(x,gumbel_r.pdf(x,*params),color='red',label = equation)
    ax[indeces[n]].legend()
    n+=1
    if conc == '5':
        ax2[0].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
    if conc == '8':
        ax2[1].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
    if conc == '10':
        ax2[2].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
    if conc == '20':
        ax2[3].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
    if conc == '50':
        ax2[4].hist(time['corr_sec'],edgecolor='black',label = '{}$\mu$M ATP, n={},\n'.format(conc,len(time))+r' a={}, loc={:.1}m, scale={:.1}'.format(params[0],params[1],params[2]),density=True)
fig2.show()
#%%
print(xlabels)
    
#%%

data.reset_index(inplace=True)
print(len(data['trajectory'].unique()))
print(len(events['trajectory'].unique()))


#%%plot sortet out events:
out = 'C:/Users/shm975/Documents/UvrD_data/events/'
rst.plot_selection(data,sortOut,out,segments=False, regression =True)

intsortOut = sortOutReason.count('Intensity')
RscoreSortOut = sortOutReason.count('Rscore')
plt.bar(['R$^2$','Intensity'],[RscoreSortOut,intsortOut])
#%%
plt.show()

#%%
timing = rst.unbinding_timing(events)
#%%

#print(timing)
import toolbox as tb
from toolbox import double_expon
from toolbox import double_gaussian

hist = np.histogram(timing['slice'],bins=len(timing))
hist_dist = rv_histogram(hist)
plt.close(fig)
fig,ax = plt.subplots()

x = np.arange(1,800,0.2)
params = expon.fit(timing['slice'])
paramsG = gamma.fit(timing['slice'])
paramsLap = laplace.fit(timing['slice'])
paramsLomax = lomax.fit(timing['slice'])

#print(params)
ax.hist(timing['slice'],bins = 24,edgecolor = 'black',density=True,label='n={}'.format(len(timing)),cumulative=False)

paramsSE = tb.fit_hist(timing['slice'], tb.exponential,p0=[0.008])[0]
paramsDE = tb.fit_hist(timing['slice'], tb.double_gaussian,p0 = [0.1,0.008,15])[0]
plt.plot(x,tb.exponential(x,*paramsSE))
plt.plot(x,tb.double_expon(x,*paramsDE))
plt.show()

#%%
#ax.plot(x,hist_dist.pdf(x))
ax.plot(x,expon.pdf(x,*params))

ax.plot(x,laplace.pdf(x, *paramsLap))
#ax.plot(x,gamma.pdf(x,*paramsG))
ax.legend()
fig.show()