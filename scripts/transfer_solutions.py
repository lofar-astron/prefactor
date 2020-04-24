#!/usr/bin/env python
# Written for prefactor by Alexander Drabent (alex@tls-tautenburg.de), 19th July 2019.

from losoto.h5parm import h5parm
from losoto.lib_operations import *

import logging
import os
import re
import numpy
import subprocess
########################################################################
def main(h5parmdb, refh5parm = '/data/solutions/3C48.h5', insolset='sol000', outsolset='sol000', insoltab='amplitude000', outsoltab='amplitude000', antenna = '[FUSPID].*', trusted_sources = ['3C48', '3C147'], parset = None):

    
    logging.info('Transferring solutions from ' +  str(refh5parm) + ' to ' + str(h5parmdb) + '.')
    logging.info('Solutions will be transferred from soltab ' + str(insoltab) + ' to ' + str(outsoltab) + '.')
    
    ### Open up the h5parm, get an example value
    data           = h5parm(h5parmdb, readonly=False)
    refdata        = h5parm(refh5parm, readonly=True)
    insolset       = refdata.getSolset(insolset)
    outsolset      = data.getSolset(outsolset)
    insoltab       = insolset.getSoltab(insoltab)
    outsoltab      = outsolset.getSoltab(outsoltab, useCache = True)
    
    ### Get the calibrator information from the target solset and check whether it is trusted or not.
    calibrators = []
    sources     = outsolset.obj._f_get_child('source')
    for i in numpy.arange(len(sources)):
        calibrators.append(sources[i][0])
    calibrator = numpy.unique(calibrators)
    calibrator = [c.decode() for c in calibrator]
    if len(calibrator) > 1:
        logging.error('There is more than one calibrator used in the target solution set: ' + str(calibrator) + '. No solutions will be transferred.')
        data.close()
        refdata.close()
        return(1)
    if calibrator[0].upper() in trusted_sources:
        logging.info(calibrator[0] + ' is a trusted calibrator source. No solutions from a reference solution set will be transferred!')
        data.close()
        refdata.close()
        return(0)
    else:
        logging.info(calibrator[0] + ' is not a trusted calibrator source. Solutions from a reference solution set will be transferred!')
    
  
    ### check for matching antennas
    station_names      = outsoltab.ant
    stations_to_transfer = [ station_name for station_name in station_names if re.match(antenna, station_name) ]           
    if len(stations_to_transfer) == 0:
        logging.warning('No stations found matching the regular expression: ' + antenna)
        data.close()
        refdata.close()
        return(0)
    
    ### Properly define axes order
    out_axes          = outsoltab.getAxesNames()
    out_axes.insert(0, out_axes.pop(out_axes.index('time'))) ## put time table to the first place
    out_axes.insert(1, out_axes.pop(out_axes.index('freq'))) ## put frequency table to the second place
    out_axes.insert(2, out_axes.pop(out_axes.index('ant')))  ## put antenna table to the third place
    in_axes           = insoltab.getAxesNames()
    in_axes.insert(0, in_axes.pop(in_axes.index('freq'))) ## put frequency table to the second place
    in_axes.insert(1, in_axes.pop(in_axes.index('ant')))  ## put antenna table to the third place
    
    if 'time' in in_axes:
        logging.warning('Reference solution table provides a time axis. The bandpass is assumed to be constant in time. Thus, only the first time step will be transferred.')
        in_axes.insert(2, in_axes.pop(in_axes.index('time')))  ## put time table to the third place

    ### resort the val and weights array according to the new axes order
    for outvals, outweights, coord, selection in outsoltab.getValuesIter(returnAxes=out_axes, weight=True):
        outvals = reorderAxes( outvals, outsoltab.getAxesNames(), out_axes)
        outweights = reorderAxes( outweights,  outsoltab.getAxesNames(), out_axes )
    
    ### resort the val and weights array according to the new axes order
    for invals, inweights, coord, selection in insoltab.getValuesIter(returnAxes=in_axes, weight=True):
        invals = reorderAxes( invals, insoltab.getAxesNames(), in_axes)
        inweights = reorderAxes( inweights,  insoltab.getAxesNames(), in_axes )
    
    ### check for equistance in frequency
    freq_resolution = numpy.unique(numpy.diff(outsoltab.freq))
    if len(freq_resolution) > 1:
        logging.error('Frequency axis is not equidistant!')
        data.close()
        refdata.close()
        return(1)
    
    ### look for nearest neighbours in frequency and write them into a dictionary
    nearest_freq = {}
    freq_idx = []
    for i, frequency in enumerate(numpy.array(outsoltab.freq)):
        nearest_freq[i] = numpy.abs(frequency - numpy.array(insoltab.freq)).argmin()
        freq_idx.append(nearest_freq[i])

    ### do the job!
    for station_to_transfer in stations_to_transfer:
        logging.info('Transferring solutions for antenna: ' + str(station_to_transfer))
        unique_elements, counts_elements = numpy.unique(freq_idx, return_counts = True)
        ant_index     = list(outsoltab.ant).index(station_to_transfer)
        try:
            ref_ant_index = list(insoltab.ant).index(station_to_transfer)
        except ValueError:
            logging.warning('Station ' + station_to_transfer + ' does not exist in reference soltab. No solutions will be transferred!')
            continue
        for j, frequency in enumerate(numpy.array(outsoltab.freq)):
            if 'time' in in_axes:
                weight = inweights[nearest_freq[j], ref_ant_index, 0]
            else:
                weight = inweights[nearest_freq[j], ref_ant_index]
            if counts_elements[list(unique_elements).index(nearest_freq[j])] > 1:
                if frequency < insoltab.freq[0] or frequency > insoltab.freq[-1]:
                    counts_elements[list(unique_elements).index(nearest_freq[j])] -= 1
                    logging.warning('No entry for frequency ' + str(frequency) + '. Solution will be flagged.')
                    weight = 0
            for i, time in enumerate(outsoltab.time):
                outweights[i,j,ant_index] = weight
                if 'time' in in_axes:
                    outvals[i,j,ant_index] = invals[nearest_freq[j], ref_ant_index, 0]
                else:
                    outvals[i,j,ant_index] = invals[nearest_freq[j], ref_ant_index]
                

    ### writing to the soltab
    outsoltab.setValues(outvals, selection)
    outsoltab.setValues(outweights, selection, weight = True)
    outsoltab.addHistory('Transferred solutions from ' + os.path.basename(refh5parm) + ' for the stations ' + str(stations_to_transfer).lstrip('[').rstrip(']') )
    outsoltab.flush()
    
    ### plotting new tables if provided
    if parset:
        if os.path.exists(parset):
            download  = subprocess.Popen(['losoto',  '-v', h5parmdb, parset], stdout=subprocess.PIPE)
            errorcode = download.wait()
            if errorcode != 0:
                logging.error('An error has occured while plotting.')
                data.close()
                refdata.close()
                return(1)
        else:
            logging.error('Parset file ' + parset + ' has not been found.')
            data.close()
            refdata.close()
            return(1)

    data.close()
    refdata.close()
    return(0)
    
########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Transfer solutions from a solution table of a reference h5parm to a solution table of a target h5parm.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which the solutions should be transferred.')
    parser.add_argument('--refh5parm', '--outh5parm', type=str,
                        help='Name of the h5parm from which the solutions should be transferred.')
    parser.add_argument('--insolset', '--insolset', type=str, default='sol000',
                        help='Name of the input h5parm solution set (default: sol000)')
    parser.add_argument('--outsolset', '--outsolset', type=str, default='sol000',
                        help='Name of the output h5parm solution set (default: sol000)')
    parser.add_argument('--insoltab', '--insoltab', type=str, default='amplitude000',
                        help='Name of the input h5parm solution set (default: amplitude000)')
    parser.add_argument('--outsoltab', '--outsoltab', type=str, default='amplitude000',
                        help='Name of the output h5parm solution set (default: amplitude000)')
    parser.add_argument('--antenna', '--antenna', type=str, default='[FUSPID].*',
                        help='Regular expression of antenna solutions to be transferred (default: [FUSPID].*)')
    parser.add_argument('--parset', '--parset', type=str, default=None,
                        help='Parset for plotting diagnostic plots after transfer with LoSoTo.')

    args = parser.parse_args()
    
    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    main(args.h5parm, refh5parm = args.refh5parm, insolset=args.insolset, outsolset=args.outsolset, insoltab=args.insoltab, outsoltab=args.outsoltab, antenna = args.antenna, parset = args.parset)
