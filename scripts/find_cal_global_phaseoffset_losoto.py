#!/usr/bin/env python

# NOTE: THIS finds the median phase offset between XX and YY per station
#
# Written by Wendy Williams, 9 Jan 2015     (original version)
#         by Andreas Horneffer, 12 Oct 2015 (losoto, pythonplugin)

import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import scipy.signal as s
import sys
import os
from losoto.h5parm import h5parm, solWriter, solFetcher

import pylab as pl

def normalize(phase):
    """
    Normalize phase to the range [-pi, pi].
    """

    # Convert to range [-2*pi, 2*pi].
    out = np.fmod(phase, 2.0 * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2.0 * np.pi
    out[out > np.pi] -= 2.0 * np.pi

    return out

def smooth(x,window_len=11,window='hanning'):
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
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def smooth_array(x, axis=0,window_len=11,window='hanning'):
    nx, ny = x.shape
    y = np.zeros((nx,ny))
    d = window_len/2
    for i in range(ny):
        xi = x[:,i]
        yi = smooth(xi,window_len=window_len,window='hanning')
        yi = yi[d:-d]
        y[:,i] = yi
    return y

def main(losotoname, store_basename, refstationID=2, sourceID=0):
    # station-ID 2 is more likely to be on the superterp
    # I don't imagine someone would use a sourceID != 0
    inh5parm = h5parm(losotoname ,readonly=True)
    phasetab = inh5parm.getSoltab('sol000','phase000')
    phases_tmp = np.copy(phasetab.val)
    freqs = np.copy(phasetab.freq)
    stationsnames = [ stat for stat in phasetab.ant]
    # this gets the subband number to any given frequency in HBA-low
    subbands = np.unique(np.round(freqs/195.3125e3-512.))
    nsubbands = len(subbands)
    nchan = len(freqs)/nsubbands
    if nsubbands*nchan != len(freqs):
        print "find_cal_global_phaseoffset_losoto.py: irregular number of ch/SB detected! Bailing out!"
        print "  nchan %d, nSB: %d, nfreq: %d" % (nchan, nsubbands, len(freqs))
        sys.exit(1)
    tmpfreqs = freqs.reshape([nsubbands,nchan])
    freq_per_sb = np.mean(tmpfreqs,axis=1)
    nstations = len(stationsnames)
    refphases = phases_tmp[:,sourceID,refstationID,:,:]
    # loop over all stations:
    for istat in xrange(nstations):
        phases_00 = phases_tmp[0,sourceID,istat,:,:]-refphases[0,:,:]
        phases_11 = phases_tmp[1,sourceID,istat,:,:]-refphases[1,:,:]
        phases_diff = normalize(phases_00-phases_11)
        tmp_phases_diff = np.median(phases_diff,axis=1)
        med_phases_diff = np.median(tmp_phases_diff.reshape([nsubbands,nchan]),axis=1)
        if istat == 0:
            global_stat_offsets = med_phases_diff
        else:
            global_stat_offsets = np.vstack( (global_stat_offsets, med_phases_diff) )
    global_stat_offsets_smoothed = np.zeros([nsubbands,nstations])
    for istat in xrange(nstations):
        global_stat_offsets_smoothed[:,istat] = s.medfilt(global_stat_offsets[istat,:], kernel_size=15)

    np.save('freqs_for_phase_array.npy', freq_per_sb)
    np.save(store_basename.strip() + '_phase_array.npy', global_stat_offsets_smoothed)
    np.save(store_basename.strip() + '_station_names.npy', stationsnames)

    # do the plotting!
    Nr = int(np.sqrt(nstations))
    if Nr*(Nr+1) > nstations:
        Nc = Nr+1
    else:
        Nc = Nr+1
        Nr = Nc
    f, axs = pl.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(16,12))
    ax = axs.reshape((Nr*Nc,1))
    for istat, stat in enumerate(stationsnames):
        ax[istat][0].set_title(stat)
        ax[istat][0].plot(subbands, global_stat_offsets[istat,:],'b')
        ax[istat][0].plot(subbands, global_stat_offsets_smoothed[:,istat],'g')
        ax[istat][0].set_ylim(-3.2,3.2)
        ax[istat][0].set_xlim(np.min(subbands), np.max(subbands))

    f.savefig('phase_xx_yy_offset.png')

