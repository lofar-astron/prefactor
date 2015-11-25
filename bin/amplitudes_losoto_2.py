#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
#matplotlib.use('GTK')
import lofar.parmdb
import sys
import os
import scipy
import time
import numpy
import math
import pyrap.tables
import scipy.signal
import scipy.ndimage
import scipy.interpolate
from pylab import *
#import lofar.expion.fitting as fitting
pi = numpy.pi
from scipy.interpolate import splprep, splev
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
import sys, os, glob, re
import numpy as np
import shutil
import progressbar
import logging
import pyrap.tables as pt
import lofar.parmdb
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter, solFetcher
import multiprocessing as mp

args = sys.argv

########################
###### USER INPUT ######

globaldbname = args[1] # input h5 parm file
calsource    = args[2] # name for writing outputfiles
n_chan       = args[3] # number of channels solved for per subband (i.e., the number of solutions along the frequencies axis of the MS) 
if len(args) > 4:
    bad_sblist_str = args[4]
    bad_sblist = [int(SB) for SB in bad_sblist_str.strip('\"\'').split(';')]
else:
    bad_sblist   = [] # run first with all subbands, 
                      # then determine which subbands are bad based on the plots and set "bad_sblist" 
                      # accordingly and re-run (some trial and error will be needed)

print "bad SBs:",bad_sblist

#### END USER INPUT ####
########################


n_chan = numpy.int(n_chan) # convert to integer, just in case
def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)
    
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]
	 
def rebin2d(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions

    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array

    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated

    See Also
    --------
    resize : Return a new array with the specified shape.

    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])

    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])

    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return numpy.repeat(numpy.repeat(a, m/M, axis=0), n/N, axis=1)


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
def median_window_filter(ampl, half_window, threshold):
    ampl_tot_copy = numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]

    #fix oct 2012
    median_array  = scipy.signal.medfilt(sol,half_window*2.-1)

    sol_flag = numpy.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.
        if abs(sol[i] - median) > (threshold * q):
            sol_flag[i] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
           ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy
    
    
def running_median(ampl,half_window) :

    ampl_tot_copy = numpy.copy(ampl)

    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl
    std = numpy.zeros(len(ampl))

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]
    
    for i in range(len(ampl)):
      #print i, i+half_window
      std[i] =  numpy.median(sol[i:i+(2*half_window)])  

    return std


ionmodel = h5parm(globaldbname ,readonly=True)
amptab = ionmodel.getSoltab('sol000','amplitude000')
# save for quicker access
#numpy.save('amplitude_array.npy',ionmodel.amplitudes)

#amplitude_arraytmp = numpy.load('amplitude_array.npy')
amplitude_arraytmp = amptab.val[:]

print numpy.shape(amplitude_arraytmp)

print len(amptab.ant[:]),  len(amptab.freq[:]) ,  len(amptab.time[:])


# in reinout version 62 is the number of freqs, 75 is time, 2912/8 is averaging in freq
#amplitude_array = rebin(amplitude_arraytmp, (75, 2912/8, 62, 1, 2))
amplitude_array = rebin(amplitude_arraytmp, (2, 1, len(amptab.ant[:]),  len(amptab.freq[:])/n_chan ,  len(amptab.time[:])))
#amplitude_array = rebin(amplitude_arraytmp, (2, 1, 61, 244, 4754))

# in reinout version averaging of 2912/8 bins to one value persubband.
#freqsarray = rebin(numpy.copy(ionmodel.freqs),(2912/8,))
freqsarray = rebin(numpy.copy(amptab.freq[:]),( len(amptab.freq[:])/n_chan ,))
#freqsarray = rebin(numpy.copy(amptab.freq[:]),(244,))

#print numpy.shape(amplitude_array)
#print numpy.shape(freqsarray)


source_id  = 0
show_plot  = True
freqs_tmp = freqsarray



goodfreq_el = range(len(amptab.freq[:])/n_chan)




# remove the badd subbands given by the user
for bad in bad_sblist:
  print 'removing subband: ',  bad
  goodfreq_el.remove(bad)


print 'Using subbands ', goodfreq_el

freqs = numpy.copy(freqs_tmp[goodfreq_el])
times = numpy.arange(0,len(amptab.time[:]),1)
freqs_new  = numpy.arange(numpy.min(freqs)+1e3,numpy.max(freqs)-1e3, 195.3125e3)
amps_array = numpy.zeros( (len(amptab.ant[:]),len(amptab.time[:]),len(freqs_new),2), dtype='float')



os.system('rm -f ' + 'freqs_for_amplitude_array.npy')
numpy.save('freqs_for_amplitude_array.npy',freqs_new)


print numpy.shape(freqs)
print numpy.shape(times)
make_matrixplot = True
if make_matrixplot:
        Nr = int(np.sqrt(len(amptab.ant[:])))
        Nc = int((len(amptab.ant[:]))/Nr)+1
        fp,axp = matplotlib.pyplot.subplots(Nr,Nc,sharex=True,sharey=True,figsize=(16,12))
	fp.subplots_adjust(hspace = 0.4)
        axsp = axp.reshape((Nr*Nc,1))
        for antenna_id in range(0,len(amptab.ant[:])):
                amp_xx_tmp = numpy.copy(amplitude_array[0,source_id,antenna_id,:,:])# array(time, freq)
                amp_yy_tmp = numpy.copy(amplitude_array[1,source_id,antenna_id,:,:]) # array(time, freq)
        
                amp_xx = numpy.copy(amp_xx_tmp[goodfreq_el,:])
                amp_yy = numpy.copy(amp_xx_tmp[goodfreq_el,:])
                
                try:
                        axsp[antenna_id][0].imshow(numpy.transpose(amp_xx), vmax=maxscale, vmin=minscale, aspect='auto')
                except NameError:
                        minscale = numpy.median(amp_xx)*0.3
                        maxscale = numpy.median(amp_xx)*2.0
                        axsp[antenna_id][0].imshow(numpy.transpose(amp_xx), vmax=maxscale, vmin=minscale, aspect='auto')
                axsp[antenna_id][0].set_title(amptab.ant[antenna_id], fontsize=8)
                axsp[antenna_id][0].set_ylabel('time', fontsize=8)
                axsp[antenna_id][0].set_xlabel('subband', fontsize=8)
		
        fp.savefig('matrix_xx.png')

        fp,axp = matplotlib.pyplot.subplots(Nr,Nc,sharex=True,sharey=True,figsize=(16,12))
	fp.subplots_adjust(hspace = 0.4)
        axsp = axp.reshape((Nr*Nc,1))
        for antenna_id in range(0,len(amptab.ant[:])):
                amp_xx_tmp = numpy.copy(amplitude_array[0,source_id,antenna_id,:,:])# array(time, freq)
                amp_yy_tmp = numpy.copy(amplitude_array[1,source_id,antenna_id,:,:]) # array(time, freq)
        
                amp_xx = numpy.copy(amp_xx_tmp[goodfreq_el,:])
                amp_yy = numpy.copy(amp_yy_tmp[goodfreq_el,:])
                
                try:
                        axsp[antenna_id][0].imshow(numpy.transpose(amp_yy), vmax=maxscale, vmin=minscale, aspect='auto')
                except NameError:
                        minscale = numpy.median(amp_yy)*0.3
                        maxscale = numpy.median(amp_yy)*2.0
                        axsp[antenna_id][0].imshow(numpy.transpose(amp_yy), vmax=maxscale, vmin=minscale, aspect='auto')
                axsp[antenna_id][0].set_title(amptab.ant[antenna_id],  fontsize=8)
                axsp[antenna_id][0].set_ylabel('time', fontsize=8)
                axsp[antenna_id][0].set_xlabel('subband', fontsize=8)
      
        fp.savefig('matrix_yy.png')
        
       
for antenna_id in range(0,len(amptab.ant[:])):
 amp_xx_tmp = numpy.copy(amplitude_array[0,source_id,antenna_id,:,:])# array(time, freq)
 amp_yy_tmp = numpy.copy(amplitude_array[1,source_id,antenna_id,:,:]) # array(time, freq)
        
 amp_xx = numpy.copy(amp_xx_tmp[goodfreq_el,:])
 amp_yy = numpy.copy(amp_yy_tmp[goodfreq_el,:])
 
 print 'Doing', amptab.ant[antenna_id]
         
 if show_plot:
   minscale = numpy.median(amp_xx)*0.3
   maxscale = numpy.median(amp_xx)*2.0
   subplots_adjust(wspace = 0.6)

   matplotlib.pyplot.subplot(121)
   matplotlib.pyplot.imshow(numpy.transpose(amp_xx), vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(amptab.ant[antenna_id]+ ' XX ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   #matplotlib.pyplot.tight_layout()
   matplotlib.pyplot.subplot(122)
   matplotlib.pyplot.imshow(numpy.transpose(amp_yy), vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(amptab.ant[antenna_id]+ ' YY ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   #matplotlib.pyplot.tight_layout()
   matplotlib.pyplot.savefig('%s_ampmat.pdf'%(amptab.ant[antenna_id]))
   matplotlib.pyplot.close()
   matplotlib.pyplot.cla()

 #print np.shape(amp_xx),'asdfdf'
 amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (3,3))
 amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (1,7))
 amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (3,3))
 amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (1,7))
 #print numpy.shape(amp_xx)
 
 if show_plot:
   subplots_adjust(wspace = 0.6)

   matplotlib.pyplot.subplot(121)
   matplotlib.pyplot.imshow(numpy.transpose(amp_xx), vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(amptab.ant[antenna_id]+' Smoothed XX ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   matplotlib.pyplot.subplot(122)
   matplotlib.pyplot.imshow(numpy.transpose(amp_yy), vmax=maxscale, vmin=minscale, aspect='auto')
   matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
   matplotlib.pyplot.ylabel('time')
   matplotlib.pyplot.title(amptab.ant[antenna_id]+' Smoothed YY ampl')
   matplotlib.pyplot.colorbar()#(orientation='horizontal')
   matplotlib.pyplot.savefig('%s_ampmat_smooth.pdf'%(amptab.ant[antenna_id]))
   matplotlib.pyplot.close()
   matplotlib.pyplot.cla()

  #2D interpol XX
 for time in range(0,len(amptab.time[:])):
  
   amps      = numpy.copy(amp_xx[:,time])
   fl = scipy.interpolate.interp1d(freqs,amps,kind='slinear',bounds_error=False)
   ynew= fl(freqs_new)
   #print 'AMPSARRAY: ',amps_array
   amps_array[antenna_id,time,:,0] = numpy.copy(savitzky_golay(ynew, 17, 2))
   #if show_plot:
      #print numpy.shape(freqs_new), numpy.shape(amps_array)
      #matplotlib.pyplot.plot(freqs_new/1e6,amps_array[antenna_id,time,:,0])
 if show_plot:
   matplotlib.pyplot.plot(freqs_new/1e6,numpy.median(amps_array[antenna_id,:,:,0], axis=0))
   matplotlib.pyplot.xlabel('freq [MHz]')
   matplotlib.pyplot.ylabel('ampl')
   matplotlib.pyplot.savefig('%s_profileXX.pdf'%(amptab.ant[antenna_id]))
   matplotlib.pyplot.close()
   matplotlib.pyplot.cla()

  #2D interpol YY
 for time in range(0,len(amptab.time[:])):
   amps      = numpy.copy(amp_yy[:,time])
   fl = scipy.interpolate.interp1d(freqs,amps,kind='slinear',bounds_error=False)
   ynew= fl(freqs_new)
   amps_array[antenna_id,time,:,1] = numpy.copy(savitzky_golay(ynew, 17, 2))
   #if show_plot:
     
 if show_plot:
   matplotlib.pyplot.plot(freqs_new/1e6,numpy.median(amps_array[antenna_id,:,:,1], axis=0))
   matplotlib.pyplot.xlabel('freq [MHz]')
   matplotlib.pyplot.ylabel('ampl')
   matplotlib.pyplot.savefig('%s_profileYY.pdf'%(amptab.ant[antenna_id]))
   matplotlib.pyplot.close()
   matplotlib.pyplot.cla()





   #sys.exit()
 #for SB in range(len(ionmodel.freqs[:])):
 # print 'Doing SB', SB
 # amp_xx[:,SB] = median_window_filter(amp_xx[:,SB], 200,5.0) 
 # amp_xx[:,SB] = median_window_filter (amp_xx[:,SB], 50, 3.0)
 # amp_xx[:,SB] = median_window_filter(amp_xx[:,SB], 50, 2.5)

 #matplotlib.pyplot.imshow(amp_xx, vmax=numpy.median(amp_xx)*2.0, vmin=numpy.median(amp_xx)*0.3, aspect='auto')
 #matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
 #matplotlib.pyplot.ylabel('time')
 #matplotlib.pyplot.colorbar()
 #matplotlib.pyplot.show()
os.system('rm -f ' + calsource + '_amplitude_array.npy') 
numpy.save(calsource + '_amplitude_array.npy',amps_array)
print numpy.shape(amps_array)
