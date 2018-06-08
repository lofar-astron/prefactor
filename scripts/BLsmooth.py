#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2015 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Usage: BLavg.py vis.MS
# Load a MS, average visibilities according to the baseline lenght,
# i.e. shorter BLs are averaged more, and write a new MS

import os, sys
import optparse, itertools
import logging
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import pyrap.tables as pt
logging.basicConfig(level=logging.DEBUG)

def addcol(ms, incol, outcol):
    if outcol not in ms.colnames():
        logging.info('Adding column: '+outcol)
        coldmi = ms.getdminfo(incol)
        coldmi['NAME'] = outcol
        ms.addcols(pt.makecoldesc(outcol, ms.getcoldesc(incol)), coldmi)
    if outcol != incol:
        # copy columns val
        logging.info('Set '+outcol+'='+incol)
        pt.taql("update $ms set "+outcol+"="+incol)

opt = optparse.OptionParser(usage="%prog [options] MS", version="%prog 0.1")
opt.add_option('-f', '--ionfactor', help='Gives an indication on how strong is the ionosphere [default: 0.2]', type='float', default=0.2)
opt.add_option('-s', '--bscalefactor', help='Gives an indication on how the smoothing varies with BL-lenght [default: 0.5]', type='float', default=0.5)
opt.add_option('-i', '--incol', help='Column name to smooth [default: DATA]', type='string', default='DATA')
opt.add_option('-o', '--outcol', help='Output column [default: SMOOTHED_DATA]', type="string", default='SMOOTHED_DATA')
opt.add_option('-w', '--weight', help='Save the newly computed WEIGHT_SPECTRUM, this action permanently modify the MS! [default: False]', action="store_true", default=False)
opt.add_option('-r', '--restore', help='If WEIGHT_SPECTRUM_ORIG exists then restore it before smoothing [default: False]', action="store_true", default=False)
opt.add_option('-b', '--nobackup', help='Do not backup the old WEIGHT_SPECTRUM in WEIGHT_SPECTRUM_ORIG [default: do backup if -w]', action="store_true", default=False)
opt.add_option('-a', '--onlyamp', help='Smooth only amplitudes [default: smooth real/imag]', action="store_true", default=False)
opt.add_option('-c', '--onlycopy', help='Copies a column instead of smoothing', type="string", default=True)
(options, msfile) = opt.parse_args()

if msfile == []:
    opt.print_help()
    sys.exit(0)
msfile = msfile[0]
print msfile

if not os.path.exists(msfile):
    logging.error("Cannot find MS file.")
    sys.exit(1)

# open input/output MS
ms = pt.table(msfile, readonly=False, ack=False)

freqtab = pt.table(msfile + '/SPECTRAL_WINDOW', ack=False)
freq = freqtab.getcol('REF_FREQUENCY')
freqtab.close()
wav = 299792458. / freq
timepersample = ms.getcell('INTERVAL',0)
all_time = ms.getcol('TIME_CENTROID')

# check if ms is time-ordered
if not all(all_time[i] <= all_time[i+1] for i in xrange(len(all_time)-1)):
    logging.critical('This code cannot handle MS that are not time-sorted.')
    sys.exit(1)

# create column to smooth
addcol(ms, options.incol, options.outcol)
# if smoothing should not be performed
if options.onlycopy == 'False':
    sys.exit(0)

# retore WEIGHT_SPECTRUM
if 'WEIGHT_SPECTRUM_ORIG' in ms.colnames() and options.restore:
    addcol(ms, 'WEIGHT_SPECTRUM_ORIG', 'WEIGHT_SPECTRUM')
# backup WEIGHT_SPECTRUM
elif options.weight and not options.nobackup:
    addcol(ms, 'WEIGHT_SPECTRUM', 'WEIGHT_SPECTRUM_ORIG')

ant1 = ms.getcol('ANTENNA1')
ant2 = ms.getcol('ANTENNA2')

all_uvw = ms.getcol('UVW')
# iteration on baseline combination
# do for loop twice to remove all_uvw and save memory
stddevs = np.zeros( len(set(ant1)) * len(set(ant2)) )
for i, ant in enumerate(itertools.product(set(ant1), set(ant2))):

    if ant[0] >= ant[1]: continue
    sel = (ant1 == ant[0]) & (ant2 == ant[1])

    # compute the FWHM
    uvw = all_uvw[sel,:]
    uvw_dist = np.sqrt(uvw[:, 0]**2 + uvw[:, 1]**2 + uvw[:, 2]**2)
    dist = np.mean(uvw_dist) / 1.e3
    if np.isnan(dist): continue # fix for missing anstennas

    stddev = options.ionfactor * (25.e3 / dist)**options.bscalefactor * (freq / 60.e6) # in sec
    stddev = stddev/timepersample # in samples
    logging.debug("For BL %i - %i (dist = %.1f km): sigma=%.2f samples." % (ant[0], ant[1], dist, stddev))
    stddevs[i] = stddev

del all_uvw

nchans = len(pt.table(msfile+"/SPECTRAL_WINDOW",ack=False)[0]["CHAN_FREQ"])
pols = [0,1,2,3] # full-pol smoothing
for pol in pols:
    logging.debug("Workign on pol %i" % pol)
    all_data = ms.getcolslice(options.outcol, [0,pol], [nchans-1,pol])
    all_weights = ms.getcolslice('WEIGHT_SPECTRUM', [0,pol], [nchans-1,pol])
    all_flags = ms.getcolslice('FLAG', [0,pol], [nchans-1,pol])

    all_flags[ np.isnan(all_data) ] = True # flag NaNs
    all_weights[all_flags] = 0 # set weight of flagged data to 0
    del all_flags

    # iteration on baseline combination
    logging.info('Smoothing baselines')
    for i, ant in enumerate(itertools.product(set(ant1), set(ant2))):

        if ant[0] >= ant[1]: continue
        if stddevs[i] == 0: continue # fix for missing anstennas

        sel = (ant1 == ant[0]) & (ant2 == ant[1])

        #Multiply every element of the data by the weights, convolve both the scaled data and the weights, and then
        #divide the convolved data by the convolved weights (translating flagged data into weight=0). That's basically the equivalent of a
        #running weighted average with a Gaussian window function.

        # get cycle values
        weights = all_weights[sel]
        data = all_data[sel]

        # set bad data to 0 so nans do not propagate
        data = np.nan_to_num(data*weights)

        # smear weighted data and weights
        if options.onlyamp:
            dataAMP = gfilter(np.abs(data), stddevs[i], axis=0)
            dataPH = np.angle(data)
        else:
            dataR = gfilter(np.real(data), stddevs[i], axis=0)#, truncate=4.)
            dataI = gfilter(np.imag(data), stddevs[i], axis=0)#, truncate=4.)

        weights = gfilter(weights, stddevs[i], axis=0)#, truncate=4.)

        # re-create data
        if options.onlyamp:
            data = dataAMP * ( np.cos(dataPH) + 1j*np.sin(dataPH) )
        else:
            data = (dataR + 1j * dataI)
        data[(weights != 0)] /= weights[(weights != 0)] # avoid divbyzero
        all_data[sel] = data
        all_weights[sel] = weights

        #print np.count_nonzero(data[~flags]), np.count_nonzero(data[flags]), 100*np.count_nonzero(data[flags])/np.count_nonzero(data)
        #print "NANs in flagged data: ", np.count_nonzero(np.isnan(data[flags]))
        #print "NANs in unflagged data: ", np.count_nonzero(np.isnan(data[~flags]))
        #print "NANs in weights: ", np.count_nonzero(np.isnan(weights))

    logging.warning('Writing %s column.' % options.outcol)
    ms.putcolslice(options.outcol, all_data, [0,pol], [nchans-1,pol])
    del all_data

    if options.weight:
        logging.warning('Writing WEIGHT_SPECTRUM column.')
        ms.putcol('WEIGHT_SPECTRUM', all_weights, [0,pol], [nchans-1,pol])
        del all_weights

ms.close()
logging.info("Done.")
