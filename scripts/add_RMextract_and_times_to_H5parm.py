#!/usr/bin/env python
import os
import numpy as np
import sys
import glob
import pyrap.tables as pt
import logging
from losoto.h5parm import h5parm
from losoto.lib_operations import *

########################################################################
def get_RM_vals(MSinfo, server, prefix, ionexPath, timestep):
    """
    Call getRM() from RMextract to get the RM values for the opbservation

    Parameters
    ----------
    MSinfo : ReadMs object
       the ReadMs object for the MS for which we compute the rotation values
    server : str or None
       URL of the server where we can find the IONEX files, or None if no download is wanted
    prefix : str
       prefix of the IONEX files
    path : str
       path where we can find or store the IONEX files
    timestep : float
        timestep to use
    """
    if ionexPath[-1] != '/':
        print "get_RM_vals: ionexPath doesn't end in \"/\", adding that character."
        ionexPath += '/'

    from RMextract import getRM
    rmdict = getRM.getRM(MSinfo.msname, server=server, prefix=prefix, ionexPath=ionexPath,
                         timestep=timestep)
    if not rmdict:
        if not server:
            raise ValueError("One or more IONEX files is not found on disk and download is disabled!\n"
                             "(You can run \"bin/download_IONEX.py\" outside the pipeline if needed.)")
        else:
            raise ValueError("Couldn't get RM information from RMextract! (But I don't know why.)")

    return rmdict


###Reading in the the parameters of target data with PYRAP and putting them into directories for further use###############
class ReadMs:
    def __init__(self, ms, skiptime=False):
        self.timepara={'start':0, 'end':0, 'step':0, 'cent':0}
        self.freqpara={'start':0, 'end':0, 'step':0, 'cent':0}
        self.msname = ms
        print ms
        if not os.path.isdir(ms): sys.exit('INPUT MS DOES NOT EXIST!')
        t = pt.table(ms, readonly=True, ack=False)
        if not skiptime:
                ##########Getting Time parameters first#############
            t1 = t.sort ('unique desc TIME')
            self.timepara['step'] = t.getcell('EXPOSURE',0)
            self.timepara['start'] =  np.min(t.getcol('TIME'))-self.timepara['step']/2.
            self.timepara['end'] =  np.max(t.getcol('TIME'))+self.timepara['step']/2.
            self.timepara['cent'] = self.timepara['start']+(self.timepara['end']-self.timepara['start'])/2.
            self.mstimevalues = t1.getcol('TIME')[::-1]
            t1.close()
            ##########Getting Frequency Parameters###################
        freq=pt.table(t.getkeyword("SPECTRAL_WINDOW"), readonly=True, ack=False)
        self.fullband = freq.getcell('TOTAL_BANDWIDTH', 0)
        self.freqpara['cent'] = freq.getcell('REF_FREQUENCY', 0)
        self.freqpara['step'] = freq.getcell('CHAN_WIDTH', 0)[0]
        self.msfreqvalues = freq.getcell('CHAN_FREQ', 0)
        self.freqpara['start'] = self.msfreqvalues[0]-self.freqpara['step']/2.
        self.freqpara['end'] = self.msfreqvalues[-1]+self.freqpara['step']/2.
        freq.close()
        ##########Getting Station Names###################
        antennas = pt.table(t.getkeyword("ANTENNA"), readonly=True, ack=False)
        self.stations = antennas.getcol('NAME')
        antennas.close()
        t.close()

    def GetTimepara(self, p=''):
        if p != '': return self.timepara[p]
        else: return self.timepara
    def GetFreqpara(self, p=''):
        if p != '': return self.freqpara[p]
        else: return self.freqpara
    def GetMSNamepara(self): return self.msname

########################################################################
def main(MSfiles, h5parmdb, ionex_server="ftp://ftp.unibe.ch/aiub/CODE/", ionex_prefix='CODG', ionexPath="IONEXdata/", solsetName='sol000'):


    mslist = MSfiles.lstrip('[').rstrip(']').replace(' ','').replace("'","").split(',')

    if len(mslist) == 0:
        raise ValueError("Did not find any existing directory in input MS list!")


    data          = h5parm(h5parmdb, readonly=False)
    solset        = data.getSolset(solsetName)

    station_names = solset.getAnt().keys()


    msinfo = ReadMs(mslist[0]) # process first MS
    ## this is the same for all antennas
    starttime = msinfo.timepara['start']
    endtime   = msinfo.timepara['end']
    timestep  = msinfo.timepara['step'] * int(300.0 / msinfo.timepara['step']) # aim for ~ 300 seconds
    freqstep  = msinfo.GetFreqpara('step')

    for ms in mslist[1:]:
        msinfo = ReadMs(ms)
        # this is the same for all antennas
        assert starttime == msinfo.timepara['start']
        assert endtime   == msinfo.timepara['end']
        assert freqstep  == msinfo.GetFreqpara('step')

    if ionex_server.strip(' []\'\"').lower() == 'none':
        ionex_server = None
    rmdict = get_RM_vals(msinfo, ionex_server, ionex_prefix, ionexPath, timestep)
    rm_list = []
    for antenna in msinfo.stations:
        rm_list.append(rmdict['RM'][antenna])
    rm_vals = np.array(rm_list)

    # check if we have data for all target stations
    for antenna in msinfo.stations:
        if antenna not in station_names:
            if antenna[:2] == 'CS' or antenna[:2] == 'RS':
                # fail if it is a Dutch station
                logging.error("Station %s not found in list of calibrator data!"%(antenna))
                raise ValueError("Station "+antenna+" missing in calibrator data!")
            else:
                # just print a warning for international stations
                logging.warning("No calibratior data for station %s, but international stations will be flagged anyhow."%(antenna))

    for antenna_id, antenna in enumerate(station_names):
        if antenna not in msinfo.stations:
            logging.warning("Station %s not found in target observation, skipping generation of calibration values."%(antenna))
            continue

    logging.info('Adding rotation measure values to: ' + solsetName + ' of ' + h5parmdb)
    try:
        new_soltab = solset.getSoltab('RMextract')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='rotationmeasure', soltabName='RMextract',
                             axesNames=['ant', 'time'], axesVals=[msinfo.stations, rmdict['times']],
                             vals=rm_vals,
                             weights=np.ones_like(rm_vals, dtype=np.float16))
    new_soltab.addHistory('CREATE (by add_RMextract_and_times_to_H5parm.py script)')

    logging.info('Adjusting times of clock offsets in: ' + solsetName + ' of ' + h5parmdb)
    # For simplicity, we use the same time axis as for the RMextract values
    clock_soltab = solset.getSoltab('clock000')
    clock_single = clock_soltab.val[0, :]
    clock = np.resize(clock_single, (len(rmdict['times']), len(clock_soltab.ant[:])))
    try:
        new_soltab = solset.getSoltab('clock')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='clock', soltabName='clock',
                             axesNames=['time', 'ant'], axesVals=[rmdict['times'], clock_soltab.ant[:]],
                             vals=clock,
                             weights=np.ones_like(clock, dtype=np.float16))
    new_soltab.addHistory('CREATE (by add_RMextract_and_times_to_H5parm.py script)')

    logging.info('Adjusting times of bandpass in: ' + solsetName + ' of ' + h5parmdb)
    # For simplicity, we use the same time axis as for the RMextract values
    bp_soltab = solset.getSoltab('bandpass')
    try:
        new_soltab = solset.getSoltab('bandpass_notimes')
        new_soltab.delete()
    except:
        pass
    bp_soltab.rename('bandpass_notimes')
    bp = bp_soltab.val[:]
    bp = np.resize(bp, (len(rmdict['times']), len(bp_soltab.ant[:]), len(bp_soltab.freq[:]), len(bp_soltab.pol[:])))
    try:
        new_soltab = solset.getSoltab('bandpass')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='amplitude', soltabName='bandpass',
                             axesNames=['time', 'ant', 'freq', 'pol'], axesVals=[rmdict['times'],
                             bp_soltab.ant[:], bp_soltab.freq[:],bp_soltab.pol[:]],
                             vals=bp,
                             weights=np.ones_like(bp, dtype=np.float16))
    new_soltab.addHistory('CREATE (by add_RMextract_and_times_to_H5parm.py script)')

    logging.info('Adjusting times of XYoffsets in: ' + solsetName + ' of ' + h5parmdb)
    # For simplicity, we use the same time axis as for the RMextract values
    xy_soltab = solset.getSoltab('XYoffset')
    try:
        new_soltab = solset.getSoltab('XYoffset_notimes')
        new_soltab.delete()
    except:
        pass
    xy_soltab.rename('XYoffset_notimes')
    xy = xy_soltab.val[:]
    xy = np.resize(xy, (len(rmdict['times']), len(xy_soltab.freq[:]), len(xy_soltab.ant[:]), len(xy_soltab.pol[:])))
    try:
        new_soltab = solset.getSoltab('XYoffset')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='phase', soltabName='XYoffset',
                             axesNames=['time', 'freq', 'ant', 'pol'], axesVals=[rmdict['times'],
                             xy_soltab.freq[:], xy_soltab.ant[:], xy_soltab.pol[:]],
                             vals=xy,
                             weights=np.ones_like(xy, dtype=np.float16))
    new_soltab.addHistory('CREATE (by add_RMextract_and_times_to_H5parm.py script)')

    return 0


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Adds CommonRotationAngle to an H5parm from RMextract.')

    parser.add_argument('MSfiles', type=str,
                        help='One or more MSs for which the parmdb should be created.')
    parser.add_argument('h5parm', type=str,
                        help='H5parm to which the results of the CommonRotationAngle is added.')
    parser.add_argument('--server', type=str, default='None',
                        help='URL of the server to use. (default: None)')
    parser.add_argument('--prefix', type=str, default='CODG',
                        help='Prefix of the IONEX files. (default: \"CODG\")')
    parser.add_argument('--ionexpath', '--path', type=str, default='IONEXdata/',
                        help='Directory where to store the IONEX files. (default: \"IONEXdata/\")')
    parser.add_argument('--solsetName', '--solset', type=str, default='sol000',
                        help='Name of the h5parm solution set (default: sol000)')


    args = parser.parse_args()

    MS = args.MSfiles
    h5parmdb = args.h5parm
    logging.info("Working on:", MS, h5parmdb)
    main(MS, h5parmdb, ionex_server=args.server, ionex_prefix=args.prefix, ionexPath=args.ionexpath, solsetName=args.solsetName)
