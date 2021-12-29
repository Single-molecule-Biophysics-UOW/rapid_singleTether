# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 11:19:39 2021

@author: stefa
"""



from __future__ import division
import statsmodels
from matplotlib import  pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
from statsmodels.base.model import GenericLikelihoodModel



class hypoexpon(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        if exog is None:
            exog = np.zeros_like(endog)
            endog_names = ['l1,l2']
            
        super(hypoexpon, self).__init__(endog, exog, **kwds)
    def _argcheck(self, l1,l2):
        return l1!=l2 and l1>0 and l2>0
    def _pdf(self,x,l1,l2):
        C = (l1*l2/(l1-l2))
        return (C* (np.exp(-l2*x)-np.exp(-l1*x)))
    def nloglikeobs(self, params):
        l1 = params[0]
        l2 = params[1]

        return -np.log(self._pdf(self.endog, l1, l2))
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params is None:
            l1_start = 0.2#self.endog.mean()
            l2_start = 0.02
            #excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)
            
            start_params = np.array([l1_start, l2_start])
            
        return super(hypoexpon, self).fit(start_params=start_params,
                                                    maxiter=maxiter, maxfun=maxfun, **kwds)
    
class hypoexpon_triple(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        if exog is None:
            exog = np.zeros_like(endog)
            endog_names = ['l1,l2,l3']
            
        super(hypoexpon_triple, self).__init__(endog, exog, **kwds)
    def _argcheck(self, l1,l2,l3):
        return l1!=l2 and l1>0 and l2>0 and l3>0
    def _pdf(self,x,l1,l2,l3):
        C = (l1*l2*l3/(l1-l2-l3))
        return (C* (np.exp(-l3*x)-np.exp(-l2*x)-np.exp(-l1*x)))
    def nloglikeobs(self, params):
        l1 = params[0]
        l2 = params[1]
        l3 = params[2]

        return -np.log(self._pdf(self.endog, l1, l2,l3))
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params is None:
            l1_start = 0.2#self.endog.mean()
            l2_start = 0.02
            l3_start = 0.04
            #excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)
            
            start_params = np.array([l1_start, l2_start, l3_start])
            
        return super(hypoexpon_triple, self).fit(start_params=start_params,
                                                    maxiter=maxiter, maxfun=maxfun, **kwds)    

    
class exphypoexpon(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        if exog is None:
            exog = np.zeros_like(endog)
            endog_names = ['l1,l2']
            
        super(exphypoexpon, self).__init__(endog, exog, **kwds)
    
    def _pdf(self,x,l1,l2):
        C = (l1*l2/(l1-l2))
        return (C* (np.exp(-l2*x)-np.exp(-l1*x)))
    def nloglikeobs(self, params):
        l1 = params[0]
        l2 = params[1]

        return -np.log(self._pdf(self.endog, l1, l2))
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params is None:
            lambda_start = self.endog.mean()
            excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)
            
            start_params = np.array([excess_zeros, lambda_start])
            
        return super(exphypoexpon, self).fit(start_params=start_params,
                                                    maxiter=maxiter, maxfun=maxfun, **kwds)