# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:43:50 2021

@author: shm975
"""


import os
import pandas as pd
import rst_tools as rst
import numpy as np
from matplotlib import pyplot as plt

meta_folder = "C:/Users/shm975/Documents/UvrD_data/3'Bioblunt/meta/"
uncalib_integration_folder = "C:/Users/shm975/Documents/UvrD_data/3'Bioblunt/uncalibrated/"
calib_integration_folder = "C:/Users/shm975/Documents/UvrD_data/3'Bioblunt/integration/"

meta_files = os.listdir(meta_folder)
integration_files = os.listdir(uncalib_integration_folder)

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
        if f[5:-4] in m:
            meta = m
            break    
    conv = findTimeCov(meta_folder+meta)
    #now load the actual file as pandas, add the calibrated time column save and close:
    data = pd.read_csv(uncalib_integration_folder+f)
    
    data['seconds'] = data['slice']*conv
    
    #carefull! this overwrites the data. don't make a mistake!
    data.to_csv(calib_integration_folder+f)
