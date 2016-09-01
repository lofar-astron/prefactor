#!/usr/bin/env python

#################################################################################
#                                                                               #
# Written by Wendy Williams, 26 May 2014                                        #
#                                                                               #
#                                                                               #
#################################################################################

import lofar.parmdb as lp
import numpy as np
import sys, os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse

mpl.rc('font',size =8 )
mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95 )




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

def solplot_phaseonly(parmdb, imageroot, refstationi, plot_international=False):
    parmdbmtable = lp.parmdb(parmdb)
    soldict = parmdbmtable.getValuesGrid('*')
    names = parmdbmtable.getNames()

    'Gain:1:1:Phase:RS508HBA'
    stationsnames = np.array([name.split(':')[-1] for name in names])
    stationsnames = np.unique(stationsnames)
    if not plot_international:
        stationsnames = np.array([name for name in stationsnames if name[0] in ['C','R'] ])
    Nstat = len(stationsnames)

    refstation = stationsnames[refstationi]
    phase11_ref = soldict['Gain:1:1:Phase:{s}'.format(s=refstation)]['values']
    phase00_ref = soldict['Gain:0:0:Phase:{s}'.format(s=refstation)]['values']
    times= soldict['Gain:1:1:Phase:{s}'.format(s=refstation)]['times']

    #Nr = int(np.ceil(np.sqrt(Nstat)))
    #Nc = int(np.ceil(np.float(Nstat)/Nr))
    Nr = int(Nstat)
    Nc = 1
    f, ax = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(12,72))
    axs = ax.reshape((Nr*Nc,1))
    for istat, station in enumerate(stationsnames):
        phase11 = soldict['Gain:1:1:Phase:{s}'.format(s=station)]['values']
        phase00 = soldict['Gain:0:0:Phase:{s}'.format(s=station)]['values']
        
        # don't plot flagged phases
        phase00 = np.ma.masked_where(phase00==0, phase00)
        phase11 = np.ma.masked_where(phase11==0, phase11)
        
        #try:
            #if len(phase11) > 1000:
                #fmt = ','
            #else:
                #fmt = '.'
        #except TypeError:
            #print "no phases"
        fmt = ','
        
        axs[istat][0].plot(times, normalize(phase00-phase00_ref), color='b',  marker=fmt, label='Gain:0:0:Phase')
        axs[istat][0].plot(times, normalize(phase11-phase11_ref), color='g',  marker=fmt, label='Gain:1:1:Phase')
        axs[istat][0].set_ylim(-3.2, 3.2)
        axs[istat][0].set_xlim(times.min(), times.max())
        axs[istat][0].set_title(station)


    f.savefig(imageroot+"_phase.png",dpi=100)
    return


def solplot_ampphase(parmdb, imageroot, refstationi, norm_amp_lim=False, median_amp=False, plot_international=False):

    parmdbmtable = lp.parmdb(parmdb)
    soldict = parmdbmtable.getValuesGrid('*')
    names = parmdbmtable.getNames()

    'Gain:1:1:Phase:RS508HBA'
    stationsnames = np.array([name.split(':')[-1] for name in names])
    stationsnames = np.unique(stationsnames)
    if not plot_international:
        stationsnames = np.array([name for name in stationsnames if name[0] in ['C','R'] ])
    Nstat = len(stationsnames)

    refstation = stationsnames[refstationi]
    times = soldict['Gain:1:1:Real:{s}'.format(s=refstation)]['times']
    
    real11_ref = soldict['Gain:1:1:Real:{s}'.format(s=refstation)]['values']
    real00_ref = soldict['Gain:0:0:Real:{s}'.format(s=refstation)]['values']
    imag11_ref = soldict['Gain:1:1:Imag:{s}'.format(s=refstation)]['values']
    imag00_ref = soldict['Gain:0:0:Imag:{s}'.format(s=refstation)]['values']

    valscorr00 = real00_ref +1.j*imag00_ref
    valscorr11 = real11_ref +1.j*imag11_ref


    amp00_ref = np.abs(valscorr00)
    phase00_ref = np.angle(valscorr00)
    amp11_ref = np.abs(valscorr11)
    phase11_ref = np.angle(valscorr11)

    Nr = int(np.ceil(np.sqrt(Nstat)))
    Nc = int(np.ceil(np.float(Nstat)/Nr))
    fp, axp = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(16,12))
    axsp = axp.reshape((Nr*Nc,1))
    fa, axa = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=(16,12))
    axsa = axa.reshape((Nr*Nc,1))
    ymin = 2
    ymax = 0
    for istat, station in enumerate(stationsnames):
        
        real11 = soldict['Gain:1:1:Real:{s}'.format(s=station)]['values']
        real00 = soldict['Gain:0:0:Real:{s}'.format(s=station)]['values']
        imag11 = soldict['Gain:1:1:Imag:{s}'.format(s=station)]['values']
        imag00 = soldict['Gain:0:0:Imag:{s}'.format(s=station)]['values']
        
        valscorr00 = real00 +1.j*imag00
        valscorr11 = real11 +1.j*imag11
        
        if len(valscorr11) > 1000:
            fmt = ','
        else:
            fmt = '.'
        amp00 = np.abs(valscorr00)
        phase00 = np.angle(valscorr00)
        amp11 = np.abs(valscorr11)
        phase11 = np.angle(valscorr11)
        #phase11 = soldict['Gain:1:1:Phase:{s}'.format(s=station)]
        #phase00 = soldict['Gain:0:0:Phase:{s}'.format(s=station)]
        
        ## for y scale: check max and min values
        amp00m = np.ma.masked_where(amp00==1, amp00).compressed()
        amp11m = np.ma.masked_where(amp11==1, amp11).compressed()
        
        if len(amp00m) > 0:
            ymax = max(np.max(amp00m),ymax)
        if len(amp11m) > 0:
            ymax = max(np.max(amp11m),ymax)
        if len(amp00m) > 0:
            ymin = min(np.min(amp00m),ymin)
        if len(amp11m) > 0:
            ymin = min(np.min(amp11m),ymin)
            
        # don't plot flagged amplitudes
        amp00 = np.ma.masked_where(amp00==1, amp00)
        amp11 = np.ma.masked_where(amp11==1, amp11)
        
        # don't plot flagged phases
        phase00 = np.ma.masked_where(phase00==0, phase00)
        phase11 = np.ma.masked_where(phase11==0, phase11)
        
        axsp[istat][0].plot(times, normalize(phase00-phase00_ref), color='b',  marker=fmt, label='Gain:0:0:Phase')
        axsp[istat][0].plot(times, normalize(phase11-phase11_ref), color='g',  marker=fmt, label='Gain:1:1:Phase')
        axsp[istat][0].set_ylim(-3.2, 3.2)
        axsp[istat][0].set_xlim(times.min(), times.max())
        axsp[istat][0].set_title(station)
        
        axsa[istat][0].plot(times, amp00, color='b', marker=fmt, label='Gain:0:0:Amp')
        axsa[istat][0].plot(times, amp11, color='g', marker=fmt, label='Gain:1:1:Amp')
        if median_amp:
            median_amp00 = np.median(amp00)
            median_amp11 = np.median(amp11)
            #if abs(median_amp00) < 1.0e-10:
                #if np.sum(amp00==median_amp00) != len(amp00):
                    #amp00 = np.ma.masked_where(amp00==1, amp00).compressed()
                    #median_amp00 = np.median(amp00)
            #if abs(median_amp11) < 1.0e-10:
                #if np.sum(amp11==median_amp11) != len(amp11):
                    #amp11 = np.ma.masked_where(amp11==1, amp11).compressed()
                    #median_amp11 = np.median(amp11)
            #print median_amp00,median_amp11
            
            axsa[istat][0].plot([times[0], times[-1]], [median_amp00,median_amp00], color='b', label='<Gain:0:0:Amp>')
            axsa[istat][0].plot([times[0], times[-1]], [median_amp11,median_amp11], color='g', label='<Gain:1:1:Amp>')
            
        if norm_amp_lim:
            axsa[istat][0].set_ylim(0,2 )
        else:
            axsa[istat][0].set_ylim(ymin,ymax)
        #print istat, station,ymin, ymax
        axsa[istat][0].set_xlim(times.min(), times.max())
        axsa[istat][0].set_title(station)


    fp.savefig(imageroot+"_phase.png",dpi=100)
    
    
    fa.savefig(imageroot+"_amp.png",dpi=100)
    return


if __name__ == "__main__":   
    parser = argparse.ArgumentParser() #prog='plot_solutions_all_stations.py',usage='[options] <parmdb> <imageroot> ')
    parser.add_argument('-p', '--phase-only', dest='phaseonly', action="store_true", default=False, help="plot phase-only solutions")
    parser.add_argument('-n', '--norm-amplitude-limits', dest='norm_amp_lim', action="store_true", default=False, help="plot amps between 0 and 2")
    parser.add_argument('-m', '--plot-median-amplitude', dest='median_amp', action="store_true", default=False, help="plot median amplitudes")
    parser.add_argument('-i', '--plot-international-stations', dest='plot_international', action="store_true", default=False, help="plot international stations")
    parser.add_argument('-r', '--refstation', dest='refstation', default=0, help="given reference station (integer)")
    parser.add_argument('parmdb', help="Name of solution parmdb")
    parser.add_argument('imageroot', help="Root name for output images")
    
    args = parser.parse_args()
    
    parmdb = args.parmdb
    imageroot = args.imageroot
    norm_amp_lim = args.norm_amp_lim
    median_amp = args.median_amp
    plot_international = args.plot_international
    refstation = int(args.refstation)

    #if len(args) <= 2:
        #print 'Insufficient arguments'
        #print 'Usage: python %prog [options] <parmdbname> <imagename>'

    #filename = args[0]
    #imagename = args[1]
    
    if args.phaseonly:
        solplot_phaseonly(parmdb, imageroot, refstation, plot_international=plot_international)
    else:
        solplot_ampphase(parmdb, imageroot, refstation, norm_amp_lim=norm_amp_lim, median_amp=median_amp, plot_international=plot_international)
