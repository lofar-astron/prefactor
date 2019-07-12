#!/usr/bin/env python
import argparse
import glob
import signal
import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
import numpy
import casacore.tables as pt
from matplotlib import cm


def input2strlist_nomapfile(invar):
   """
   from bin/download_IONEX.py
   give the list of MSs from the list provided as a string
   """
   str_list = None
   if type(invar) is str:
       if invar.startswith('[') and invar.endswith(']'):
           str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
       else:
           str_list = [invar.strip(' \'\"')]
   elif type(invar) is list:
       str_list = [str(f).strip(' \'\"') for f in invar]
   else:
       raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
   return str_list


def string2bool(instring):
    if isinstance(instring, bool):
        return instring
    if instring.upper() == 'TRUE' or instring == '1':
        return True
    elif instring.upper() == 'FALSE' or instring == '0':
        return False
    else:
        raise ValueError('string2bool: Cannot convert string "'+instring+'" to boolean!')


def main(input, output, title='uv coverage', limits=',,,', timeslots='0,10,0', antennas='-1',
         wideband=True):

        MSlist = input2strlist_nomapfile(input)
        if len(MSlist) == 0:
            print('Error: You must specify at least one MS name.')
            sys.exit(1)
        plottitle = title
        fileformat = output.split('.')[-1]
        if fileformat not in ['png','pdf','eps','ps']:
            print('Error: Unknown file extension')
            sys.exit(1)
        axlimits = limits.strip().split(',')
        if len(axlimits) == 4:
            xmin,xmax,ymin,ymax = axlimits
        else:
            print('Error: You must specify four axis limits')
            sys.exit(1)
        timeslots = timeslots.split(',')
        if len(timeslots) != 3:
            print('Error: Timeslots format is start,skip,end')
            sys.exit(1)
        for i in range(len(timeslots)):
            timeslots[i] = int(timeslots[i])
            if timeslots[i] < 0:
                print('Error: timeslots values must not be negative')
                sys.exit(1)
        antToPlotSpl = antennas.split(',')
        antToPlot = []
        for i in range(len(antToPlotSpl)):
            tmpspl = antToPlotSpl[i].split('..')
            if len(tmpspl) == 1:
                antToPlot.append(int(antToPlotSpl[i]))
            elif len(tmpspl) == 2:
                for j in range(int(tmpspl[0]),int(tmpspl[1])+1):
                    antToPlot.append(j)
            else:
                print('Error: Could not understand antenna list.')
                sys.exit(1)
        wideband = string2bool(wideband)

        outdir = os.path.dirname(output)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        badval = 0.0
        xaxisvals = []
        yaxisvals = []
        flagvals = []
        savex = numpy.array([])
        savey = numpy.array([])
        prev_firstTime = None
        prev_lastTime = None
        prev_intTime = None
        for inputMS in MSlist:
            # open the main table and print some info about the MS
            t = pt.table(inputMS, readonly=True, ack=False)
            tfreq = pt.table(t.getkeyword('SPECTRAL_WINDOW'),readonly=True,ack=False)
            ref_freq = tfreq.getcol('REF_FREQUENCY',nrow=1)[0]
            ch_freq = tfreq.getcol('CHAN_FREQ',nrow=1)[0]
            if wideband:
                ref_wavelength = 2.99792458e8/ch_freq
            else:
                ref_wavelength = [2.99792458e8/ref_freq]

            firstTime = t.getcell("TIME", 0)
            lastTime = t.getcell("TIME", t.nrows()-1)
            intTime = t.getcell("INTERVAL", 0)
            nTimeslots = (lastTime - firstTime) / intTime
            if timeslots[1] == 0:
                    timeskip = 1
            else:
                    timeskip = int(timeslots[1])
            if timeslots[2] == 0:
                    timeslots[2] = nTimeslots
            # open the antenna subtable
            tant = pt.table(t.getkeyword('ANTENNA'), readonly=True, ack=False)

            # Station names
            antList = tant.getcol('NAME')
            if len(antToPlot)==1 and antToPlot[0]==-1:
                    antToPlot = list(range(len(antList)))

            # select by time from the beginning, and only use specified antennas
            tsel = t.query('TIME >= %f AND TIME <= %f AND ANTENNA1 IN %s AND ANTENNA2 IN %s' % (firstTime+timeslots[0]*intTime,firstTime+timeslots[2]*intTime,str(antToPlot),str(antToPlot)), columns='ANTENNA1,ANTENNA2,UVW,FLAG_ROW')

            # Now we loop through the baselines
            i = 0
            nb = (len(antToPlot)*(len(antToPlot)-1))/2
            for tpart in tsel.iter(["ANTENNA1","ANTENNA2"]):
                ant1 = tpart.getcell("ANTENNA1", 0)
                ant2 = tpart.getcell("ANTENNA2", 0)
                if ant1 not in antToPlot or ant2 not in antToPlot: continue
                if ant1 == ant2: continue
                i += 1
                # Get the values to plot (only need to read again if this ms has new times)
                if prev_firstTime != firstTime or prev_lastTime != lastTime or prev_intTime != intTime:
                    uvw = tpart.getcol('UVW', rowincr=timeskip)
                    savex = numpy.append(savex,[uvw[:,0],-uvw[:,0]])
                    savey = numpy.append(savey,[uvw[:,1],-uvw[:,1]])

                # Get the flags
                flags = tpart.getcol('FLAG_ROW', rowincr=timeskip)
                for w in ref_wavelength:
                    flagvals.extend(flags.tolist())
                    flagvals.extend(flags.tolist())
            for w in ref_wavelength:
                xaxisvals.extend((savex/w/1000.).tolist())
                yaxisvals.extend((savey/w/1000.).tolist())
            prev_firstTime = firstTime
            prev_lastTime = lastTime
            prev_intTime = intTime

        # Plot the uv coverage
        xaxisvals = numpy.array(xaxisvals)
        yaxisvals = numpy.array(yaxisvals)
        zvals = numpy.zeros(xaxisvals.shape)
        flagvals = numpy.array(flagvals)
        flagged_ind = numpy.where(flagvals)
        zvals[flagged_ind] = 1.0
        uvmax = max(xaxisvals.max(),yaxisvals.max())
        uvmin = min(xaxisvals.min(),yaxisvals.min())
        uvuplim = 0.02*(uvmax-uvmin)+uvmax
        uvlolim = uvmin-0.02*(uvmax-uvmin)
        if xmin == '':
            minx = uvlolim
        else:
            minx = float(xmin)
        if xmax == '':
            maxx = uvuplim
        else:
            maxx = float(xmax)
        if ymin == '':
            miny = uvlolim
        else:
            miny = float(ymin)
        if ymax == '':
            maxy = uvuplim
        else:
            maxy = float(ymax)
        if minx == maxx:
            minx = -1.0
            maxx = 1.0
        if miny == maxy:
            miny = -1.0
            maxy = 1.0
        plt.xlabel(r'u [k$\lambda$]')
        plt.ylabel(r'v [k$\lambda$]')
        plt.xlim([minx,maxx])
        plt.ylim([miny,maxy])
        gridsize = 300
        plt.hexbin(xaxisvals, yaxisvals, C=zvals, gridsize=gridsize, cmap=cm.jet, bins=None)
        plt.title(plottitle)
        plt.axes().set_aspect('equal')
        plt.grid(True)
        cb = plt.colorbar()
        cb.set_label('flagged fraction')
        plt.savefig(output)

        # Plot histogram of the flagged fraction vs. uv distance
        uvdist = numpy.sqrt(xaxisvals**2 + yaxisvals**2)
        sorted_ind = numpy.argsort(uvdist)
        uvdist = uvdist[sorted_ind]
        zvals = 1.0-zvals[sorted_ind]
        hist, bins = numpy.histogram(uvdist, 50)
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)]
        inds = []
        for b in bins[:-1]:
            inds.append(numpy.where(uvdist > b)[0][0])
        inds.append(len(bins))
        flagged_mean = [numpy.mean(zvals[inds[i]:inds[i+1]]) for i in range(len(inds)-1)]
        plt.clf()
        plt.bar(bins_mean, flagged_mean, (bins[1]-bins[0])/2.0)
        plt.xlabel(r'uv distance [k$\lambda$]')
        plt.ylabel('unflagged fraction')
        plt.savefig('{0}_uvdist{1}'.format(*os.path.splitext(output)))




