# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 17:31:30 2021

@author: shm975
"""


import numpy as np
from matplotlib import pyplot as plt

def simulate(p1,p2,N=10000, timesteps = 200):
    
    
    states = np.ones(N)+1 #all particles have state 2 in hte beginning
    
    time = np.arange(1,timesteps+1,1)    
    nativeList = []
    decayList = []
    rate1 = 1-p1
    rate2 = 1-p2
    intermediateList = []
    timeOfDis = []
        
    for t in time:
        #Count occurances of states at this particular time:
        native = list(states).count(2)
        intermediate = list(states).count(1)
        decay = list(states).count(0)
        nativeList.append(native)
        decayList.append(decay)
        intermediateList.append(intermediate)
        #go through particles and decide if they decay:
        for particle in range(N):
            if states[particle] == 2:
                r = np.random.random()    
                if r > rate1:
                    states[particle] -= 1
            if states[particle] == 1:
                r = np.random.random()
                if r > rate2:
                    states[particle] -= 1
                    timeOfDis.append(t)
                #change state of this particle:
            if states[particle] == 0:
                    pass
    return np.array(nativeList), np.array(decayList), np.array(intermediateList), np.array(timeOfDis)

nativeList = []
decayList = []
intermediateList = []
timeOfDis = []
scenarios = [(0.08,0.02),(0.08,0.1),(0.08,0.9)]
for s in scenarios:
    native, decay, intermediate,time = simulate(*s)
    nativeList.append(native)
    decayList.append(decay)
    intermediateList.append(intermediate)
    timeOfDis.append(time)

                
    
#%%    
from scipy.stats import rv_continuous
from scipy.special import gamma as gf   
from scipy.optimize import curve_fit 
from unbinding_rate.statsmodel_hypoExp import hypoexpon 
    #print(np.random.exponential(0))
def hypoExp(x,l1,l2, cumulative=True):
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
    #normalize:
    f = f/ np.sum(f)
    # for i in f:
    #     print(i)
    #     print(np.sum(f))
    #     print(i/np.sum(f))
    
    
    
    if cumulative == True:
        f = np.cumsum(f)
    else:
        f=(f*np.sum(f))/np.max(f)
    return f

class hypoexp(rv_continuous): 
    "hypo-exponential distribution"           
    def _pdf(self,x,l1,l2):
        
        C = (l1*l2/(l1-l2)) 
        return (C* (np.exp(-l2*x)-np.exp(-l1*x)))
    #def _logpdf(self,x,l1,l2)
    # def _logpdf(self, x, l1, l2):
    #     logC = np.log(1/(l1-l2))
    #     logpdf = logC - l2*x + np.log(1+np.exp(x*(l2-l2)))
    #     return logpdf
    def _argcheck(self, l1,l2):
        return l1!=l2 and 0<l1 and 0<l2

hypoexp = hypoexp(name = 'hypoexp')
#fig,[ax1,ax2] = plt.subplots(1,2)
plt.figure()
ax1 = plt.subplot(1,2,1)
#ax1.set_yscale('log')
#ax1.margins(0.05) 
ax2 = plt.subplot(3,2,2)
#ax2.set_xscale('log')
#ax2.margins(2, 2)
ax3 = plt.subplot(3,2,4)
#ax3.set_xscale('log')
#ax3.margins(x=0, y=-0.25) 
ax4 = plt.subplot(3,2,6)
#ax4.set_xscale('log')
#ax4.margins(x=0, y=-0.25) 
#ax5 = plt.subplot(3,3,3)
#ax6 = plt.subplot(3,3,6)
#ax7 = plt.subplot(3,3,9)
axes = [ax1,ax2,ax3,ax4]
x = np.arange(0.1,200,0.1)
for i in range(3):
    
    ax1.plot(nativeList[i], label = 'native state, p1={:.3},p2={:.3}'.format(scenarios[i][0],scenarios[i][1]),ls=':')
    ax1.plot(decayList[i],label='disappeared, p1={:.3},p2={:.3}'.format(scenarios[i][0],scenarios[i][1]),ls = '-.')
    ax1.plot(intermediateList[i],label='intermediate, p1={:.3},p2={:.3}'.format(scenarios[i][0],scenarios[i][1]),ls = '--')
    ax1.plot(nativeList[i]+intermediateList[i],label = 'observed, p1={:.3},p2={:.3}'.format(scenarios[i][0],scenarios[i][1]))
    #ax1.legend()
    
    cons = ({'type': 'ineq',
         'fun' : lambda x: np.all(x)})

    
    model = hypoexpon(timeOfDis[i])
    results = model.fit(start_params= [0.01,10],maxiter=100000, maxfunc= 10000, constraints=cons)
    print(results.summary())
    #params = hypoexp.fit(timeOfDis[i],0.2,'MLE',floc = 0,fscale = 1)
    params = results.params
    #axes[i+1].hist(timeOfDis[i],edgecolor='black',bins=50,density = True,label = 'p1={:.3}, p2={:.3}'.format(scenarios[i][0],scenarios[i][1]),alpha = 0.5,cumulative = True)
    axes[i+1].plot(x,model._pdf(x, *params),label = 'fit: p1 = {:.2}, p2 = {:.2}'.format(params[0],params[1]),alpha = 0.7,color = 'red')
    
    axes[i+1].hist(timeOfDis[i],edgecolor='black',bins=50,density = True,label = 'p1={:.3}, p2={:.3}'.format(scenarios[i][0],scenarios[i][1]),alpha = 0.5,cumulative = False)
    mean = np.mean(timeOfDis[i])
    std = np.std(timeOfDis[i])
    cV = std/mean
    print('mean:{}'.format(mean))
    
    
        #rate1 = (2/mean)* (1+np.sqrt(1+2*(cV**2 -1)))**(-1)
        #rate2 = (2/mean)* (1-np.sqrt(1+2*(cV**2 -1)))**(-1)
        #print('rate1:{:.2}, rate2: {:.2}, c$_V$={}'.format(rate1,rate2,cV))
    
    print(len(timeOfDis[i]))


    #binStep = bins[1]-bins[0]
    
    #cumData = np.cumsum(hist)*binStep    
    
    #popt,pcov = curve_fit(hypoExp,X2,F2,p0=[0.08,0.2])
    #axes[i+1].plot(x,hypoExp(x, scenarios[i][0], scenarios[i][1]))
    #axes[i+4].plot(x,hypoExp(x, scenarios[i][0], scenarios[i][1],cumulative=False))
    #axes[i+4].plot(x,hypoexp.pdf(x, *params),label = 'fit: p1 = {:.2}, p2 = {:.2}'.format(params[0],params[1]),alpha = 0.7,color = 'red')
    axes[i+1].plot(x,hypoexp.pdf(x, scenarios[i][0],scenarios[i][1]),label = 'theoretical',color = 'green',lw=3)
    axes[i+1].legend()
    axes[i+1].legend()
    print(scenarios[i][0],scenarios[i][1])
    #print(params)
#ax2.plot(time,hypoExp(time,1-rate1,1-rate2))




plt.show()




#the main time iteration:
#for t in time:
    #we start with N particles in state 1.