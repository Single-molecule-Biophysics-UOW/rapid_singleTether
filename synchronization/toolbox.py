# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:07:53 2018

@author: Stefan Mueller
This program implements often used functions, for processing data, fitting and plotting
"""

import numpy as np
from math import factorial
#from matplotlib import pyplot as plt
from scipy import optimize as opt
import os
from scipy import stats
import pandas as pd



import zipfile as z
#import poisson_MLE as pmle
from statsmodels.base.model import GenericLikelihoodModel
from matplotlib.collections import LineCollection
from matplotlib import pyplot as plt
import shutil as sh
import traceback
from collections import Iterable

import re
import pandas as pd
import sys
#from matplotlib import collections  as mc
#from matplotlib.patches import ConnectionPatch
from scipy.signal import savgol_filter as sgf
import re
# from read_roi import read_roi_zip

#%%


"""
define some basic functions for fitting curves, by passing number_of_parameters to the function, it returns the numbver of parameters it needs instead of the function
"""
def getSeriesEvents(filenames):
    series = []
    event = []
    for i in filenames:
        filenames_parse = re.split('_|\.',i)
        series.append(filenames_parse[1])
        event.append(filenames_parse[3])
    return series,event

def cureTIP(TIP,filenames,sortOutList):
    proc = []
    for j in range(len(filenames)):
        if filenames[j] in sortOutList:
            print(filenames[j])
            continue
        proc.append(np.nanmax(TIP[2][j])-np.nanmin(TIP[2][j]))
    return proc

def parsecureData(data,filenames,sortOutList):
    Time = []
    Intensity = []
    Position =[]
    j=0
    for i in range(data.shape[0]):
        #print (j)
        
        if i%3==0:
            if filenames[j] in sortOutList:
#                print("Time for filename"+str(filenames[j])+"was sorted out")
                continue
            Time.append(data[i])
        if i%3==1:
            if filenames[j] in sortOutList:
#                print("Intensity for filename"+str(filenames[j])+"was sorted out")
                continue
            Intensity.append(data[i])
        if i%3==2:
            if filenames[j] in sortOutList:
#                print("Position for filename"+str(filenames[j])+"was sorted out")
                j+=1
                continue
            else:
                j+=1
            Position.append(data[i])
    return Time,Intensity,Position
def parseData(data):
    Time = []
    Intensity = []
    Position =[]
    for i in range(data.shape[0]):
        if i%3==0:           
            Time.append(data[i])
        if i%3==1:
            Intensity.append(data[i])
        if i%3==2:
            Position.append(data[i])
    return Time,Intensity,Position
def cureFilenames(filenames,sortOutList):
    curFilenames = []
    for i in filenames:
        if i in sortOutList:
            continue
        else:
            curFilenames.append(i)
    return curFilenames
def fancy_plot(line_ends):
    """
    make a fancy plot for length of DNA pieces or replication events.
    Input: Endpoints of the lines as list.
    """
    sort_ends = np.sort(line_ends,0)
    lines = []
    d = 0
    for i in sort_ends:
        lines.append([(0,d), (0+i,d)])
        d+=2
    #gs1 = gridspec.GridSpec(2, 1)
    #gs1.update(hspace=0.0) # set the spacing between axes. 
    #print lines
    lc = LineCollection(lines,color = 'gray', linewidths=3,zorder=2,alpha = 0.8)
    fig, ax0 = plt.subplots()
    ax0.add_collection(lc)
    ax0.autoscale()
    ax0.margins(0.1)
    return fig,ax0
def fancy_color_plot(line_ends,colors):
    """
    make a fancy plot for length of DNA pieces or replication events.
    Input: Endpoints of the lines as list.
    """
    sort_ends = np.sort(line_ends,0)
    lines = []
    d = 0
    for i in sort_ends:
        lines.append([(0,d), (0+i,d)])
        d+=2
    #gs1 = gridspec.GridSpec(2, 1)
    #gs1.update(hspace=0.0) # set the spacing between axes. 
    #print lines
    lc = LineCollection(lines,color = colors,cmap=plt.cm.get_cmap('RdBu') ,linewidths=3,zorder=2,alpha = 0.8)
    fig, ax0 = plt.subplots()
    ax0.add_collection(lc)
    ax0.autoscale()
    ax0.margins(0.1)
    return fig,ax0

# def poisson_MLE(data):
#     result = pmle.minimize_loglikelyhood(data)
#     return result

def load_file_names(folder):
    """
    load all Filenames of Files in Folder  and return as list of strings
    input: absolute path to folder
    output: List
    """
    file_names = []
    for file in os.listdir(folder):
        file_names.append(file)
    return file_names

def flatten(lis):
    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:        
            yield item 
def read_zip(folder,extract_folder):
    file_names = []
    for file in os.listdir(folder):
        file_names.append(file)                   
    #files = []
    data = []
    folders = []
    for i in file_names:
        with z.ZipFile(folder+i,'r') as myzip:
            data.append(myzip.extractall(path = extract_folder+"/"+i))
            folders.append(extract_folder+i+'/')
    return data, folders
    
def load_csv(folder):
    names = load_file_names(folder)
    data = []
    for i in names:
        data.append(np.transpose(np.loadtxt(folder+i,delimiter=',',skiprows=1)))
    return data
def load_single_csv(filepath):
    data = np.transpose(np.loadtxt(filepath,delimiter=',',skiprows=1))
    return data  



def load_files(folder,delimiter=None,skiprows=0,usecols = None,dataend1=None,dataend2=None,namesonly = False, dtype = None):
    """
    reads all data in a folder and returns filenames and the data\n
    the path to the folder must be passed as argument \n
    Arguments:
        folder: path to folder were data is saved in
    Keyword arguments: 
        delimiter: is passed on to numpy.loadtxt, default = None \n
        skiprows: is passed on to numpy.loadtxt, default = 0 \n
        dataend: if specified dataname must contain "dataend", default = None \n
        good to filter for .csv data only or certain names\n
        namesonly: if True only a list of filenames is returned
    Output:
        filenames: list of read files\n
        data: list of arrays contains the actual data\n
    Code:\n

    file_names = []\n
    data = []\n
    if dataend!=None:\n       
        for file in os.listdir(folder):\n
            if dataend in file:\n
                file_names.append(file)\n
    else:\n
        for file in os.listdir(folder):\n
            file_names.append(file)\n
    for file in file_names:\n
        data.append(np.transpose(np.loadtxt(folder+file,delimiter=delimiter,skiprows=skiprows))) \n
    return(file_names,data)\n
    
    """
    file_names = []
    data = []
    if dataend1 == None and dataend2 == None:
        for file in os.listdir(folder):
            file_names.append(file)            
            signal = np.loadtxt(folder+file,delimiter=delimiter,skiprows=skiprows, usecols = usecols,dtype=dtype)
            data.append(np.transpose(signal))  
                
        if namesonly == False:
            return(file_names,data)
        elif namesonly == True:
            return file_names
    if dataend1!=None:
        if dataend2 == None:
            for file in os.listdir(folder):
                
                if dataend1 in file:
                    file_names.append(file)
            if namesonly == True:
                return file_names
            else:
                for file in file_names:
                    data.append(np.transpose(np.loadtxt(folder+file,delimiter=delimiter,skiprows=skiprows, usecols = usecols,dtype=dtype))) 
                return(file_names,data)
        
        elif dataend2 != None:
            for file in os.listdir(folder):
                if dataend1 in file and dataend2 in file:
                    file_names.append(file)
            if namesonly == True:
                #print ('returning names')
                return file_names
    elif namesonly == True and dataend1 == None and dataend2 ==None:
        for file in os.listdir(folder):
            file_names.append(file)
        return(file_names)
    else:
        #print('doing else')
        for file in os.listdir(folder):
            file_names.append(file)
            for file in file_names:
                data.append(np.transpose(np.loadtxt(folder+file,delimiter=delimiter,skiprows=skiprows, usecols = usecols,dtype=dtype))) 
                return(file_names,data)

def stage_pos_checker(folder,img1,img2):
    #print(img1[-7:-4])
    columns = ['key','value1','value2','value3']

    filenames = load_files(folder,dataend1=img1[-26:-23], dataend2= 'Metadata',namesonly = True)

    data = []
    for i in filenames:
        data.append(pd.read_csv(folder+i,names=columns, engine='python'))#,delimiter=',')
        for i in range(len(data[0]['key'])):
            if data[0]['key'][i] == 'dXPos':
                xpos = i
            if data[0]['key'][i] == 'dYPos':
                ypos = i
        #print data[0]['key'][100]
        for i in data:   
            X1 = i['value1'][xpos]
    
            Y1 = i['value1'][ypos]
    
    filenames2 = load_files(folder,dataend1=img2[-26:-23], dataend2= 'Metadata',namesonly = True) 
    data2 = []
    for i in filenames2:
        data2.append(pd.read_csv(folder+i,names=columns, engine='python'))#,delimiter=',')
        for i in range(len(data2[0]['key'])):
            if data2[0]['key'][i] == 'dXPos':
                xpos2 = i
            if data2[0]['key'][i] == 'dYPos':
                ypos2 = i
            
        for i in data:
            X2 = i['value1'][xpos2]
            Y2 = i['value1'][ypos2]
    if (X1,Y1) == (X2,Y2):
        return True
    else: return False
    

def fit(function, data,x=None,start_values=None):
    loc,scale = function.fit(data)#,weights=np.ones(len(data)) / len(data)
    try:
        if x == None:
            x=np.linspace(0,np.max(data)+np.max(data)/2,5000)
    except:
        print ("error")
    y = function.pdf(x,loc=loc, scale=scale)
    #y = stats.norm.pdf(loc) / scale
    return loc, scale, x,y 
def fit_gamma(data,start_values=None):
    a, loc, scale = stats.gamma.fit(data,fa=2)
    x=np.linspace(0,np.max(data)+np.max(data)/2,5000)
    y = stats.gamma.pdf(x,a,loc=loc, scale=scale)
    return a,loc,scale,x,y
def fit_expon(data):
    loc,scale = stats.expon.fit(data)
    x=np.linspace(0,np.max(data)+np.max(data)/2,5000)
    y = stats.expon.pdf(x,loc=loc, scale=scale)
    return loc, scale,x,y

def fit_expon_fixedLoc(data,fixLoc):
    loc,scale = stats.expon.fit(data,floc=fixLoc)
    x=np.linspace(0,np.max(data)+np.max(data)/2,5000)
    y = stats.expon.pdf(x,loc=loc, scale=scale)
    return loc, scale,x,y

def fit_expon_loc_scale(data,fixLoc):
    loc,scale = stats.expon.fit_loc_scale(data)
    x=np.linspace(0,np.max(data)+np.max(data)/2,5000)
    y = stats.expon.pdf(x,loc=loc, scale=scale)
    return loc, scale,x,y

def gaussian(x=0,a1=0,b1=0,c1=0, number_of_params=0): 
    """
    gaussian function: g=a1*np.exp(-(x-b1)**2/2*c1**2) \n
    by passing number_of_parameters to the function, it returns the numbver of parameters it needs instead of the function \n
    Code:\n
    if number_of_params!=0:\n
        return 3\n
    else:\n
        g=a1*np.exp(-(x-b1)**2/2*c1**2)\n
        return g\n
    
    """       
    if number_of_params!=0:
        return 3
    else:
        g=a1*np.exp(-(x-b1)**2/2*c1**2)
        return g

def double_gaussian(x=0,a1=0,b1=0,c1=0, a2=0,b2=0,c2=0, number_of_params=0):
    """
    double gaussian function: g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2) \n
    by passing number_of_parameters to the function, it returns the numbver of parameters it needs instead of the function \n
    Code:\n
    if number_of_params!=0:\n
        return 6\n
    else:\n
        g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2)\n
        return g\n
    
    """           
    if number_of_params!=0:
        return 6
    else:
        g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2)
        return g

def triple_gaussian(x=0,a1=0,b1=0,c1=0,a2=0,b2=0,c2=0,a3=0,b3=0,c3=0,number_of_params=0):
    """
    triple gaussian function: g=a1*np.exp(-(x-b1)**2/2*c1**2) \n
    by passing number_of_parameters to the function, it returns the numbver of parameters it needs instead of the function \n
    Code:\n
    if number_of_params!=0:\n
        return 9\n
    else:\n
        g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2)+a3*np.exp(-(x-b3)**2/2*c3**2)\n
        return g\n
    
    """  
    if number_of_params!=0:
        return 9
    else:
        g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2)+a3*np.exp(-(x-b3)**2/2*c3**2)
        return g
    
    
def poisson(l=0,k=0, number_of_params=0):
    """
    triple gaussian function: g=a1*np.exp(-(x-b1)**2/2*c1**2) \n
    by passing number_of_parameters to the function, it returns the numbver of parameters it needs instead of the function \n
    Code:\n
    if number_of_params!=0:\n
        return 9\n
    else:\n
        g=a1*np.exp(-(x-b1)**2/2*c1**2)+a2*np.exp(-(x-b2)**2/2*c2**2)+a3*np.exp(-(x-b3)**2/2*c3**2)\n
        return g\n
    """
    if number_of_params!=0:
        return 2
    else:
        g=np.exp(-l)*((l**k)/factorial(int(k)))
        g=g.astype(float)
        return g

def exponential(x,l):
    f=l*np.exp(-l*x)
    return f

import numpy
def window(y,window_len=7,pol_order = 1):
    """
    

    Parameters
    ----------
    y : array-like
        1-d column vector containing signal to smooth
    window_len : int, optional
        length of the applied window
    pol_order : int, optional
        polynomial order used to fit sample, default is 1 which equals mean

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return sgf(y,window_len,pol_order)
def smooth(x,window_len=11,window='hanning'):
    print('smooth is now depricated. Use scipys savitzki golay filtering directly!')
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
#    print y
    
    y_new = y[int(np.round((window_len/2-1))):-int(np.round((window_len/2)))]
    if len(y_new) != len(x):
        y_trunc = y_new[0:len(x)]
        return y_trunc
    else:
        return y[int(np.round((window_len/2-1))):-int(np.round((window_len/2)))]


def excl(data,lower_limit, upper_limit):
    """ 
    /n
    exclude points from data set below lower and above upper_limit\n
    code:
       data=data[(data>lower_limit) & (data< upper_limit)]\n
       return data \n
    """
    #print data[(data>=lower_limit) & (data<= upper_limit)]
    matrix = np.where((data>lower_limit) & (data< upper_limit))
    #print matrix
    #print data[75]
    #thresholded=data[(data>lower_limit) & (data< upper_limit)]
    return matrix



def exclude(data,lower_limit, upper_limit):
    """ 
    deprecated! use excl instead, it also returns the indeces that were not excluded!
    /n
    exclude points from data set below lower and above upper_limit\n
    code:
       data=data[(data>lower_limit) & (data< upper_limit)]\n
       return data \n
    """
    #print data[(data>=lower_limit) & (data<= upper_limit)]
    thresholded=data[(data>lower_limit) & (data< upper_limit)]
    return thresholded


def bin_mid(bin_edges):
    """
    calculates middles of bins created with numpy.histogram
    """
    #print(len(bin_edges))
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
    #print (len(bin_middles))
    return bin_middles

def hist(data, bins = None, density = False):
    """
    returns bins, counts of data \n
    bins is a list of the middle points of all bins, counts is a list of the corresponding counts\n
    where the number of bins is the sqrt of the amount of points in the data set \n
    the numpy function histogram is used to do that
    """
    if bins:
        counts, bin_edges = np.histogram(data, bins=bins, density = density)
    else:
        counts, bin_edges = np.histogram(data, density = density)
    bins = bin_mid(bin_edges)
    print (bins, counts)
    return bins, counts



def fit_hist(data, function, p0 = [], bounds =((),())):
    """
    Definition: fit_hist(data, function, p0 = [], bounds =((),())) \n
    fits the **data** with the given **function** \n
    with **p0** initial guesses can be passed to scipy.optimize \n
    must be a list of parameters\n
    with **bounds**, boundary for the fits can be passed to scipy.optimize\n
    must be a tuple of tuples: ((lower_boundaries,...),(upper_boundaries,...)) \n
    returns fit data as fit_data=popt, std_error, pcov, function \n
    **popt**: estimated optimal parameters \n
    **std_error**: standart error of parmeters, calculated as sqrt of diagonal elements of pcov \n
    **pcov**: covariance matrix calculated by scipy optimize
    """
    #p0 = list(np.zeros((1,function(number_of_params=1)))[0])
    #bounds = tuple(np.ndarray.flatten((-np.inf*np.ones((function(number_of_params=1),1))))),tuple(np.ndarray.flatten((+np.inf*np.ones((function(number_of_params=1),1)))))
    bins, counts = hist(data)
    #print(bins,counts)
    try:
        popt, pcov = opt.curve_fit(function, bins, counts,p0=p0, bounds=bounds, maxfev=100000)
        std_error = np.sqrt(np.diag(pcov))
        fit_data = popt, std_error, pcov, function
        return fit_data
    except Exception:
        exc_info = sys.exc_info()
        traceback.print_exception(*exc_info)
        print("Error - curve_fit failed")

def fit_linear(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return slope, intercept, r_value, p_value, std_err
# def plot_fit_linear(fit_data,x_lim=0):
#     """Plot a line from slope and intercept"""
#     axes = plt.gca()
#     if x_lim == 0:
#         x_vals = np.array(axes.get_xlim())
#         y_vals = fit_data[1] + fit_data[0] * x_vals
#     else:
#         x_vals = np.linspace(0,x_lim,1000)
#         y_vals = fit_data[1] + fit_data[0] * x_vals
#     plt.plot(x_vals, y_vals, '--',label='linear fit: y={0:.2f}*x+{1:.2f}'.format(fit_data[0],fit_data[1]))
#     plt.legend()

# def fit_plot(data,fit_data):
#     if fit_data == None:
#         print('Nothing to plot')
#     else:
#         xfit = np.linspace(0,len(data),len(data))#np.linspace(0,300,300)
#         yfit = fit_data[3](xfit,*fit_data[0])
#         #print ('yfit='+str(yfit))
#         plt.plot(xfit,yfit)
        
#         #plt.plot(x,triple_gaussian(x,*popt),color='b')
        
# def hist_plot(data):
#     plt.hist(data,bins=10)#int(np.sqrt(len(data))))
    
    
def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm    


def line_prepender(filename):
    print(filename)
    with open(filename, 'r+') as f:
        line = f.readline()
        content = f.read()
        f.seek(0, 0)            
        f.write( 'n'+ line.rstrip('\r\n') + '\n' + content)
        f.close()
        return (line.rstrip('\r\n') + '\n' + content)




class Poisson_MLE(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        if exog is None:
            exog = np.zeros_like(endog)
            
        super(Poisson_MLE, self).__init__(endog, exog, **kwds)
    
    def nloglikeobs(self, params):
        #pi = params[0]
        lambda_ = params[0]

        return -np.log(stats.poisson.pmf(self.endog, lambda_))
    
    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        if start_params is None:
            lambda_start = self.endog.mean()
            excess_zeros = (self.endog == 0).mean() - stats.poisson.pmf(0, lambda_start)            
            start_params = np.array([excess_zeros, lambda_start])
            
        return super(Poisson_MLE, self).fit(start_params=start_params,
                                                    maxiter=maxiter, maxfun=maxfun, **kwds)
# def linepatch(line):
#     """function for draw linerois on images \n
#        input: lineROI ibject
#        output: ConnectionPatch object, can be added to plot by ax.add_artist(patch)
#        """
#     coordsA = "data"
#     coordsB = "data"
#     con = ConnectionPatch((line.x1,line.y1), (line.x2,line.y2), coordsA,coordsB, linewidth=2, color= 'red')
#     return con  


def polyLinePlot(image,polyLine):
    """
    returns pixel values of image along a polyLine \n
    input: line: polyLineROI object, image: numpy nd array
    output: values along line as list
    """ 
    print(polyLine)
    ycoords = polyLine.y_coords
    xcoords = polyLine.x_coords
    ncoords = polyLine.nPoints
    #iterate all line segments:
    all_xcoords = []
    all_ycoords = []
    for i in range(ncoords-1):
        yseg = list(range(np.round(ycoords[i]),np.round(ycoords[i+1])))

        if np.round(xcoords[i]) == np.round(xcoords[i+1]):
            xseg = (np.round(xcoords[i])*np.ones(len(yseg))).astype(np.int)
        elif xcoords[i] < xcoords[i+1] :
            xseg = list(range(int(np.round(xcoords[i])),int(np.round(xcoords[i+1]))))
#        print(int(np.round(line.x2)),int(np.round(line.x1)))
            xseg = np.round(np.linspace(xseg[0],xseg[-1],len(yseg))).astype(np.int)
        elif xcoords[i] > xcoords[i+1]:
            xseg = list(range(int(np.round(xcoords[i+1])),int(np.round(xcoords[i]))))
            xseg = np.flip(np.round(np.linspace(xseg[0],xseg[-1],len(yseg))).astype(np.int),axis=0)
        all_xcoords.append(xseg)            
        all_ycoords.append(yseg)
    coordlist = list(map(list, zip(all_xcoords,all_ycoords)))
    values = []
    for i in coordlist:
        print (i)
        values.append(brightness(image,i))
    return values
    
def multi_polyLinePlot(image,rois):
    """
    multiplexed version of polyLinePlot
    depends on polyLinePlot.
    Output: np Array containing lists
    """    
    lineplots = []    
    for i in rois:
        lineplots.append(polyLinePlot(i,image))
    return np.array(lineplots)        



def lineplot(line, image):
    """
    returns pixel values of image along a line \n
    input: line: lineROI object, image: numpy nd array
    output: values along line as list
    """
    ycoords = list(range(int(np.round(line.y1)),int(np.round(line.y2))))
    if np.round(line.x1) == np.round(line.x2):
        xcoords = (np.round(line.x1)*np.ones(len(ycoords))).astype(np.int)
    elif line.x1 < line.x2 :
        
        xcoords = list(range(int(np.round(line.x1)),int(np.round(line.x2))))
#        print(int(np.round(line.x2)),int(np.round(line.x1)))
        
        xcoords = np.round(np.linspace(xcoords[0],xcoords[-1],len(ycoords))).astype(np.int)
    elif line.x1 > line.x2:
        xcoords = list(range(int(np.round(line.x2)),int(np.round(line.x1))))
        xcoords = np.flip(np.round(np.linspace(xcoords[0],xcoords[-1],len(ycoords))).astype(np.int),axis=0)
    
    coordlist = list(map(list, zip(xcoords,ycoords)))
    #lineplot = img(coordlist)
    #print(coordlist)
    values = []
    
    for i in coordlist:
#        print (i)
        values.append(brightness(image,i))
    return values
    #print(image[coordlist[0][0 1]])
    #print(image[[410,94]])
    #print (lineplot)
    
    #return coordlist
def multi_lineplot(image,rois):
    """
    multiplexed version of lineplot
    depends on lineplot.
    Output: np Array containing lists
    """    
    lineplots = []    
    for i in rois:
        lineplots.append(lineplot(i,image))
    return np.array(lineplots)

def brightness(img, choord):
    """
    returns pixel value of image at choordinates
    Input choord as [x,y]
    Returns brightness of image as input datatype
    """
    # print('find brightness:')
    # print(img.shape)
    # print(choord)
    #print(choord)
    x = choord[0]
    y = choord[1]
    value = img[y,x]
    return value    

def DNA_length(lineplot,return_bounds = False):
    """
    determines length of DNA molecules by calculating the first derivative
    the length is typically argmax(dx/dy)-argmin(dx/dy)
    """
    x = np.arange(0,len(lineplot),1)
    try:
        y = sgf(np.array(lineplot),7,1)
    except ValueError:
        print("encountered ValueError while smoothing, try decreasing window size")
        try:
            y = sgf(np.array(lineplot),3,1)
        except ValueError:
            print("encountered ValueError multiple times...")
            return None
    dy = np.zeros(y.shape,np.float)
    dy[0:-1] = np.diff(y)/np.diff(x)
    dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
    DNA_length = np.argmin(dy)-np.argmax(dy)
    if return_bounds==True:
        return DNA_length, [np.argmax(dy),np.argmin(dy)]
    else:
        return DNA_length
    

def makeID(fileName):
    date = (re.findall(r'(\d{6})', fileName))[0]
    event = (re.findall(r'event_(\d{3})', fileName))[0]
    series = ((re.findall(r'series_(\d{1,2})', fileName))[0]).zfill(3)
    ID = (date+series+event)
    return ID