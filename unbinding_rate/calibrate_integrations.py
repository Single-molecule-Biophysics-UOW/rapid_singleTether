# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:43:50 2021

@author: shm975
"""


import os
import pandas as pd
from synchronization import rst_tools_v2 as rst
import numpy as np
from matplotlib import pyplot as plt

meta_folder = "C:/Users/stefa/OneDrive/Documents/tempData/UvrdData/UvrD_data/5'bioBlunt/meta/"
uncalib_integration_folder = "C:/Users/stefa/OneDrive/Documents/tempData/UvrdData/UvrD_data/5'bioBlunt/uncalibrated_integrations/"
calib_integration_folder = "C:/Users/stefa/OneDrive/Documents/tempData/UvrdData/UvrD_data/5'bioBlunt/integration/"

meta_files = os.listdir(meta_folder)
integration_files = os.listdir(uncalib_integration_folder)
print(integration_files)
#%%
def findTimeCov(path):
    timestamps = []
    with open(path) as meta:
        for line in meta:
            if 'timestamp' in line:
                timestamps.append(line.split(',')[1])
    times = np.array(timestamps).astype(float)
    conv = np.round(np.mean(np.diff(times)),2)
    return conv

for f in integration_files:
    
    #we assume that the meta_file has the same name, apart from the _meta.csv!
    #find the right meta file:
    for m in meta_files:
        if '8uM' in f:
            print(f[5:-4])
            if f[5:-4] in m:
                print('f:',f)
            #print('meta:',m)
            #print('found metadata!')
                meta = m
                break
        # else:
        #     print('meta data nmot found for:')
        #     print(f)
    try:
        #print('meta:',meta)
        #print('integration:',f)
        pass
    except NameError:
        continue
    
    conv = findTimeCov(meta_folder+meta)
    print(conv)
    #now load the actual file as pandas, add the calibrated time column save and close:
    data = pd.read_csv(uncalib_integration_folder+f)
    
    data['seconds'] = data['slice']*conv
    
    #carefull! this overwrites the data. don't make a mistake!
    data.to_csv(calib_integration_folder+'new_'+f)
