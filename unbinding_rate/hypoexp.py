# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:27:47 2021

@author: stefa
"""

from scipy.stats import rv_continuous
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import expon
from scipy.stats import gamma

class hypoexp(rv_continuous): 
    "hypo-exponential distribution"           
    def _pdf(self,x,l1,l2):
        C = (l1*l2/(l1-l2)) 
        return (C* (np.exp(-l2*x)-np.exp(-l1*x)))
    def _argcheck(self, l1,l2):
        return l1!=l2 and l1>0 and l2>0
    
hypoexp = hypoexp(name = 'hypoexp')