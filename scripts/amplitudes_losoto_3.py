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

make_matrixplot = True
source_id  = 0
show_plot  = False

print "bad SBs:",bad_sblist

#### END USER INPUT ####
########################


n_chan = np.int(n_chan) # convert to integer, just in case
def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)
    
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]


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
    ampl_tot_copy = np.copy(ampl)
    ndata = len(ampl)
    flags = np.zeros(ndata, dtype=bool)
    sol = np.zeros(ndata+2*half_window)
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

    sol_flag = np.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = np.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = np.median(window_masked)
        q = 1.4826 * np.median(np.abs(window_masked - median))

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

    ampl_tot_copy = np.copy(ampl)
    
    ndata = len(ampl)
    flags = np.zeros(ndata, dtype=bool)
    sol = np.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl
    std = np.zeros(len(ampl))

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]
    
    for i in range(len(ampl)):
        #print i, i+half_window
        std[i] =  np.median(sol[i:i+(2*half_window)])  

    return std


ionmodel = h5parm(globaldbname ,readonly=True)
amptab = ionmodel.getSoltab('sol000','amplitude000')

amplitude_arraytmp = amptab.val[:]

print "Shape of amplitudes array:",np.shape(amplitude_arraytmp)

nfreqs = len(amptab.freq[:])
ntimes = len(amptab.time[:])
nants = len(amptab.ant[:])

print "Number of Antennas:",len(amptab.ant[:]),  " of Frequnecies:",nfreqs ,  " of Times:", ntimes

# Convert all frequencies to subband-frequencies 
SBfreqs = rebin(np.copy(amptab.freq[:]),( len(amptab.freq[:])/n_chan ,))






  
freqidx = np.arange(n_chan/2,nfreqs,n_chan)
freqs = amptab.freq[freqidx]
timeidx = np.arange(ntimes)
freqgrid = np.arange(nfreqs)
SBgrid = np.floor(freqgrid/n_chan)
SBvals = freqgrid/n_chan

freqs_new  = np.arange(np.min(freqs),np.max(freqs)+100e3, 195.3125e3)
amps_array_flagged = np.zeros( (nants,ntimes,len(freqs_new),2), dtype='float')
amps_array = np.zeros( (nants,ntimes,len(freqs_new),2), dtype='float')
minscale = np.zeros( nants )
maxscale = np.zeros( nants )

if len(freqs_new) < 20:
    print "Frequency span is less than 20 subbands! The filtering will not work!"
    print "Please run the calibrator pipeline on the full calibrator bandwidth." 
    raise ValueError("Frequency span is less than 20 subbands! Amplitude filtering will not work!")

# remove the badd subbands given by the user
print "Have",max(SBgrid),"subbands."
for bad in bad_sblist:
    print 'removing subband: ',  bad

try:
    os.remove('freqs_for_amplitude_array.npy')
except OSError:
    pass
np.save('freqs_for_amplitude_array.npy',freqs_new)

if show_plot:
    for antenna_id in range(0,len(amptab.ant[:])):
        amp_xx_raw = np.copy(amplitude_arraytmp[0,source_id,antenna_id,freqidx,:])
        amp_yy_raw = np.copy(amplitude_arraytmp[1,source_id,antenna_id,freqidx,:])
        minscale[antenna_id] = np.median(amp_xx_raw)*0.3
        maxscale[antenna_id] = np.median(amp_xx_raw)*2.0
        subplots_adjust(wspace = 0.6)
        matplotlib.pyplot.subplot(121)
        matplotlib.pyplot.imshow(np.transpose(amp_xx_raw), vmax=maxscale[antenna_id], vmin=minscale[antenna_id], aspect='auto')
        matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
        matplotlib.pyplot.ylabel('time')
        matplotlib.pyplot.title(amptab.ant[antenna_id]+ ' XX ampl')
        matplotlib.pyplot.colorbar()#(orientation='horizontal')
        matplotlib.pyplot.subplot(122)
        matplotlib.pyplot.imshow(np.transpose(amp_yy_raw), vmax=maxscale[antenna_id], vmin=minscale[antenna_id], aspect='auto')
        matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
        matplotlib.pyplot.ylabel('time')
        matplotlib.pyplot.title(amptab.ant[antenna_id]+ ' YY ampl')
        matplotlib.pyplot.colorbar()#(orientation='horizontal')
        matplotlib.pyplot.savefig('%s_ampmat.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()

# Flag bad data
# I'm too lazy right now to think up a more time-efficient way 
# The way I read and write the amplitudes transposes the time-frequency matrix.
for antenna_id in range(0,len(amptab.ant[:])):
    for time in range(0,len(amptab.time[:])):
        amp_xx_tmp = np.copy(amplitude_arraytmp[0,source_id,antenna_id,:,time])
        amp_yy_tmp = np.copy(amplitude_arraytmp[1,source_id,antenna_id,:,time])
        freq_tmp = amptab.freq[:]
        assert len(amp_xx_tmp[:]) == len(freq_tmp[:])
        mask_xx = np.not_equal(amp_xx_tmp,1.)
        for bad in bad_sblist:
            mask_xx = np.logical_and(SBgrid!=bad,mask_xx)
        if np.sum(mask_xx)>2:
            amps_xx_tointer = amp_xx_tmp[mask_xx]
            freq_xx_tointer = freq_tmp[mask_xx]
            #fl_xx = scipy.interpolate.interp1d(freq_xx_tointer,amps_xx_tointer,kind='linear',bounds_error=False)
            #amps_array_flagged[antenna_id,time,:,0] = fl_xx(freqs_new)
            amps_array_flagged[antenna_id,time,:,0] = np.interp(freqs_new,freq_xx_tointer,amps_xx_tointer)
        elif time>0:
            amps_array_flagged[antenna_id,time,:,0] = amps_array_flagged[antenna_id,(time-1),:,0]
        mask_yy = np.not_equal(amp_yy_tmp,1.)
        for bad in bad_sblist:
            mask_yy = np.logical_and(SBgrid!=bad,mask_yy)
        if np.sum(mask_yy)>2:
            amps_yy_tointer = amp_yy_tmp[mask_yy]
            freq_yy_tointer = freq_tmp[mask_yy]
            #fl_yy = scipy.interpolate.interp1d(freq_yy_tointer,amps_yy_tointer,kind='linear',bounds_error=False)
            #amps_array_flagged[antenna_id,time,:,1] = fl_yy(freqs_new)
            amps_array_flagged[antenna_id,time,:,1] = np.interp(freqs_new,freq_yy_tointer,amps_yy_tointer)
        elif time>0:
            amps_array_flagged[antenna_id,time,:,1] = amps_array_flagged[antenna_id,(time-1),:,1]

if make_matrixplot:
    print "Going to make Matrix-Plots."
    Nr = int(np.sqrt(len(amptab.ant[:])))
    Nc = int((len(amptab.ant[:]))/Nr)+1
    fp_xx,axp_xx = matplotlib.pyplot.subplots(Nr,Nc,sharex=True,sharey=True,figsize=(16,12))
    fp_xx.subplots_adjust(hspace = 0.4)
    axsp_xx = axp_xx.reshape((Nr*Nc,1))
    fp_yy,axp_yy = matplotlib.pyplot.subplots(Nr,Nc,sharex=True,sharey=True,figsize=(16,12))
    fp_yy.subplots_adjust(hspace = 0.4)
    axsp_yy = axp_yy.reshape((Nr*Nc,1))   
    for antenna_id in range(0,len(amptab.ant[:])):        
        amp_xx = amps_array_flagged[antenna_id,:,:,0]
        amp_yy = amps_array_flagged[antenna_id,:,:,1]
        try:
            axsp_xx[antenna_id][0].imshow(amp_xx, vmax=maxscale_mat, vmin=minscale_mat, aspect='auto')
        except NameError:
            minscale_mat = np.median(amp_xx)*0.3
            maxscale_mat = np.median(amp_xx)*2.0
            axsp_xx[antenna_id][0].imshow(amp_xx, vmax=maxscale_mat, vmin=minscale_mat, aspect='auto')
        axsp_xx[antenna_id][0].set_title(amptab.ant[antenna_id], fontsize=8)
        axsp_xx[antenna_id][0].set_ylabel('time', fontsize=8)
        axsp_xx[antenna_id][0].set_xlabel('subband', fontsize=8)
        try:
            axsp_yy[antenna_id][0].imshow(amp_yy, vmax=maxscale_mat, vmin=minscale_mat, aspect='auto')
        except NameError:
            minscale_mat = np.median(amp_yy)*0.3
            maxscale_mat = np.median(amp_yy)*2.0
            axsp_yy[antenna_id][0].imshow(amp_yy, vmax=maxscale_mat, vmin=minscale_mat, aspect='auto')
        axsp_yy[antenna_id][0].set_title(amptab.ant[antenna_id], fontsize=8)
        axsp_yy[antenna_id][0].set_ylabel('time', fontsize=8)
        axsp_yy[antenna_id][0].set_xlabel('subband', fontsize=8)
    fp_xx.savefig('matrix_xx.png')
    fp_yy.savefig('matrix_yy.png')
    fp_xx.clf()
    fp_yy.clf()
  

# Smooth the data further
ampsoutfile = open(calsource + '_amplitude_array.txt','w')
ampsoutfile.write('# Antenna name, Antenna ID, subband, XXamp, YYamp, frequency\n')
for antenna_id in range(0,len(amptab.ant[:])):
    amp_xx = np.copy(amps_array_flagged[antenna_id,:,:,0])
    amp_yy = np.copy(amps_array_flagged[antenna_id,:,:,1])

    if show_plot:
        matplotlib.pyplot.plot(np.median(amp_xx, axis=0),'b+')
        matplotlib.pyplot.grid(b=True,which='major',axis='x')
        matplotlib.pyplot.xlabel('Subband')
        matplotlib.pyplot.ylabel('ampl')
        matplotlib.pyplot.savefig('%s_flaggedXX.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()
        
        matplotlib.pyplot.plot(np.median(amp_yy, axis=0),'b+')
        matplotlib.pyplot.grid(b=True,which='major',axis='x')
        matplotlib.pyplot.xlabel('Subband')
        matplotlib.pyplot.ylabel('ampl')
        matplotlib.pyplot.savefig('%s_flaggedYY.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()

    amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (3,3))
    amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (7,1))
    amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (3,3))
    amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (7,1))



    for i in range(0,len(freqs_new)):
        ampsoutfile.write('%s %s %s %s %s %s\n'%(amptab.ant[antenna_id], antenna_id, i, np.median(amp_xx[:,i], axis=0), np.median(amp_yy[:,i], axis=0), freqs_new[i]) )

    for time in range(0,len(amptab.time[:])):
        amps_array[antenna_id,time,:,0] = np.copy(savitzky_golay(amp_xx[time,:], 17, 2))
        amps_array[antenna_id,time,:,1] = np.copy(savitzky_golay(amp_yy[time,:], 17, 2))
    
    if show_plot:
        subplots_adjust(wspace = 0.6)
        matplotlib.pyplot.subplot(121)
        matplotlib.pyplot.imshow(amps_array[antenna_id,:,:,0], vmax=maxscale[antenna_id], vmin=minscale[antenna_id], aspect='auto')
        matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
        matplotlib.pyplot.ylabel('time')
        matplotlib.pyplot.title(amptab.ant[antenna_id]+' Smooth XX ampl')
        matplotlib.pyplot.colorbar()#(orientation='horizontal')
        matplotlib.pyplot.subplot(122)
        matplotlib.pyplot.imshow(amps_array[antenna_id,:,:,1], vmax=maxscale[antenna_id], vmin=minscale[antenna_id], aspect='auto')
        matplotlib.pyplot.xlabel('calibrator SB (incrasing freq)')
        matplotlib.pyplot.ylabel('time')
        matplotlib.pyplot.title(amptab.ant[antenna_id]+' Smooth YY ampl')
        matplotlib.pyplot.colorbar()#(orientation='horizontal')
        matplotlib.pyplot.savefig('%s_ampmat_smooth.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()

        matplotlib.pyplot.plot(freqs_new/1e6,np.median(amps_array[antenna_id,:,:,0], axis=0))
        matplotlib.pyplot.xlabel('freq [MHz]')
        matplotlib.pyplot.ylabel('ampl')
        matplotlib.pyplot.savefig('%s_profileXX.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()

        matplotlib.pyplot.plot(freqs_new/1e6,np.median(amps_array[antenna_id,:,:,1], axis=0))
        matplotlib.pyplot.xlabel('freq [MHz]')
        matplotlib.pyplot.ylabel('ampl')
        matplotlib.pyplot.savefig('%s_profileYY.pdf'%(amptab.ant[antenna_id]))
        matplotlib.pyplot.close()
        matplotlib.pyplot.cla()

try:
    os.remove(calsource + '_amplitude_array.npy')
except OSError:
    pass
np.save(calsource + '_amplitude_array.npy',amps_array)
