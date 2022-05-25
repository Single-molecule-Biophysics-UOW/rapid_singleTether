# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 16:04:40 2021

__author__ = Stefan H. Mueller
__email__ = smueller@uow.edu.au

"""



#import useful packages
import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import sys
#provide path to directory containing rst_tools_v2.py

sys.path.append("C:/Users/smueller/Documents/GitHub/rapid_singleThether/")
from synchronization import rst_tools_v2 as rst

#%% read the data

#change this to the filepath and filename containing integrated rois
dataFolder = r"D:\Tirf\Lisannes_microscope\211215\Lambda\integration/"
fileName = "new_lambda_5'biossDNA_ATP0uM_OD3.3_003"

dataframe = pd.read_csv(dataFolder+fileName)
dataframe.name = fileName



#%% normalization to basepairs
#Start intensity is 2600 bp dsDNA, end is 2600 nts ssDNA and 0bp dsDNA
#use fitted values to do normalization



#first make sure seconds are passed to norm function by setting them as index:
try:
    dataframe.set_index('seconds',inplace=True)
except KeyError:
    pass

#apply normalization to the whole dataset trajectory by trajectory
#by using pandas.transform. the normalization is implemented in rst_tools_v2.py
#and described in the methods section
dataframe['norm_Intensity'] = dataframe.groupby('trajectory')['Intensity'].transform(rst.complex_norm)
dataframe['nts'] = dataframe['norm_Intensity']*2620
try:
    dataframe.reset_index(inplace=True)
except KeyError:
    #already reset.
    pass

#%%plot random trajectory:
    
dataframe.sort_values('slice',inplace=True)
    
grouped = dataframe.groupby('trajectory')
# traj = grouped.get_group(107)
groups = list(grouped.groups.keys())
choice = np.random.choice(groups,size=1)[0] #52
traj = grouped.get_group(choice)


fig, ax = plt.subplots()
ax.plot(traj['slice'], traj['nts'])
# ax.plot(traj['slice'], traj['mapped_fit'])


plt.show()



#%% calculate rates:
# now the data is normalized to basepairs. the slope of the fit theoretically
# contains the rate of the enzyme but the unit would be arb. Unit/frame. To get
# the right value it is easiest to just fit again. The function piecewise_fit
# implemented in rst_tools_v2.py does almost the same as the complex_norm function
# but it returns the individual fit parameter rather a normalized Intensity
# In future versions of this script the two functions might be merged.

   
try:
   # just make sure the calibrated 'seconds' column is the index
   dataframe.set_index('seconds',inplace=True)
except KeyError:
    pass

# apply the function 4 times returning the actual fit y-values or the parameters
# of the fit respectively. In future versions of the script this will be done
# in one go. Somehow.

dataframe['piecewise_triple'] = dataframe.groupby('trajectory')['nts'].transform(
    rst.piecewise_fit,function=rst.piecewise_linear_triple,output='predict')
dataframe['rate'] = dataframe.groupby('trajectory')['nts'].transform(
    rst.piecewise_fit,function=rst.piecewise_linear_triple,output='slope')
dataframe['a'] = dataframe.groupby('trajectory')['nts'].transform(
    rst.piecewise_fit,function=rst.piecewise_linear_triple,output='a')
dataframe['b'] = dataframe.groupby('trajectory')['nts'].transform(
    rst.piecewise_fit,function=rst.piecewise_linear_triple,output='b')


#%% plot all events:
plt.ioff()
try:
    dataframe.reset_index(inplace=True)
except ValueError:
    #already reset index
    pass
# enter the path to a folder where you want the plots (1 png file per 9 trajectories) to be saved
# Be warned if files with the same names do exist in the same folder they will
# be overwritten.
eventFolder = r"D:\Tirf\Lisannes_microscope\211215\Lambda\plots/"
#if the folder doesnt exist we create it
if os.path.isdir(eventFolder) == False:
    os.makedirs(eventFolder)

# the actual plot function is implemented in rst_tools_v2.py.
# xcolumn specifies the name of the x column to be plotted. ycolumn is the name
# of the ycolumn. Multiple y-columns are allowed in the form of a list.
# if rate is added to ycolumn it will be added as label.
rst.plot_all_trajs(dataframe, eventFolder, xcolumn = 'seconds',ycolumn = ['nts','piecewise_triple', 'rate'])

#%% save the normalized data:
norm_dataFolder = r"D:\Tirf\Lisannes_microscope\211215\Lambda\norm_data"

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
#%% make rate histogram:

# to get a rate histogram we loop through all trajectories. We check if
# the fit-parameters are sensible (see methods section of paper). If they are
# realistic we add the value to the list

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
#%% plot histogram

# a last step sorting out events that escaped the previous selection:
# sort out negative rates.
rates = rates[rates > 0]
# and sort out outliers. Either set to sensible value manually or exclude
# values based on the standard deviation of the distribution:
# warning: running this cell multiple times will filter out more and more values
# since you always sort out whats below the new standard deviation!
threshold = np.mean(rates) + 5* np.std(rates)
rates = rates[rates<threshold]
plt.figure()
plt.hist(rates,edgecolor = 'black',label='n={}'.format(len(rates)),bins = int(np.sqrt(len(rates))))

plt.xlabel('time (s)')
plt.ylabel('count')
plt.legend()
plt.show()

# if you want to save the figure:
# the format is automatically inferred from the fileending, supported are standard
# image filetypes like png, jpg,...
# available vector graphic formats are svg pdf eps...
filetype = 'svg'
fig.savefig('rate_hist.{}'.format(filetype))

#%% save the rates to csv file:

    
df_rates = pd.DataFrame({'rate':rates})

folder = r"C:\Users\shm975\Documents\OneDrive - University of Wollongong\phi29Paper\data\lambda\all_rates/"

df_rates.to_csv(folder+"new_rates_006.csv",index=False)

