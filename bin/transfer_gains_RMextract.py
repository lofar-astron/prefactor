#!/usr/bin/env python
import os
import numpy as np
import sys
import glob
import pyrap.tables as pt
import lofar.parmdb as pdb

###Reading in the the parameters of target data with PYRAP and putting them into directories for further use###############
class ReadMs:
    def __init__(self, ms):
        self.timepara={'start':0, 'end':0, 'step':0, 'cent':0}
	self.freqpara={'start':0, 'end':0, 'step':0, 'cent':0}
	self.msname = ms
	if not os.path.isdir(ms): sys.exit('INPUT MS DOES NOT EXIST!')
        ##########Getting Time parameters first#############
	t = pt.table(ms, readonly=True, ack=False)
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

# Make an empty parmDB with only the defaults and return the parmdb object
def make_empty_parmdb(outname):
    myParmdb=pdb.parmdb(outname,create=True)
    myParmdb.addDefValues("Gain:0:0:Ampl",1.)
    myParmdb.addDefValues("Gain:1:1:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Ampl",1.)
    myParmdb.addDefValues("Gain:0:0:Real",1.)
    myParmdb.addDefValues("Gain:1:1:Real",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Real",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Real",1.)
    myParmdb.addDefValues("AntennaOrientation",5.497787144)
    myParmdb.addDefValues("RotationMeasure",1e-6)
    return myParmdb

def get_COMMONROTATION_vals(MSinfo, server, prefix, ionexPath):
    """
    Call getRM() from RMextract to get the RM values for the opbservation,
    convert this to rotation values and write COMMONROTATION to the parmdb

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
    """
    from RMextract import getRM
    rmdict = getRM.getRM(MSinfo.msname,server=server,prefix=prefix,ionexPath=ionexPath)

    return rmdict

def main(msname, store_basename='caldata_transfer', store_directory='.', newparmdbext='-instrument_amp_clock_offset', 
         ionex_server="ftp://ftp.unibe.ch/aiub/CODE/", ionex_prefix='CODG', ionexPath="IONEXdata/"):

    # name (path) for parmdb to be written
    newparmDB = msname+newparmdbext

    # load the numpy arrays written by the previous scripts
    # (filenames constructed in the same way as in these scripts)
    from os.path import join
    freqs_ampl    = np.load(join(store_directory,'freqs_for_amplitude_array.npy'))
    amps_array    = np.load(join(store_directory,store_basename + '_amplitude_array.npy'))
    clock_array   = np.load(join(store_directory,'fitted_data_dclock_' + store_basename + '_1st.npy'))
    freqs_phase   = np.load(join(store_directory,'freqs_for_phase_array.npy'))
    phases_array  = np.load(join(store_directory,store_basename + '_phase_array.npy'))
    station_names = np.load(join(store_directory,store_basename + '_station_names.npy'))

    #print "phases shape:",np.shape(phases_array)
    #print "amps shape:",np.shape(amps_array)
    #print "clock shape:",np.shape(clock_array)

    #for ms in mslist: #this script works only on one MS!
    msinfo = ReadMs(msname)
    # this is the same for all antennas
    starttime = msinfo.timepara['start']
    endtime   = msinfo.timepara['end']
    startfreqs = msinfo.msfreqvalues-msinfo.GetFreqpara('step')/2.
    endfreqs   = msinfo.msfreqvalues+msinfo.GetFreqpara('step')/2.
    ntimes  = 1
    nfreqs  = len(startfreqs)


    if ionex_server.strip(' []\'\"').lower() == 'none':
        ionex_server = None
    rmdict = get_COMMONROTATION_vals(msinfo, ionex_server, ionex_prefix, ionexPath)

    c = 299792458.0
    lambdaSquared = (c/msinfo.msfreqvalues)**2
    # get an array wwith the same size as rmdict['times'] but filled with rmdict['timestep']
    timesteps = np.full_like(rmdict['times'],rmdict['timestep'])
    # same for frequencies
    freqsteps = np.full_like(msinfo.msfreqvalues,msinfo.freqpara['step'])
        

    outDB = make_empty_parmdb(newparmDB)

    # Now do the interpolating
    for antenna_id, antenna in enumerate(station_names):
        if antenna not in msinfo.stations:
            pass

        # form median of amplitudes along the time axis, for both polarizations
        amp_cal_00_all = np.median(amps_array[antenna_id,:,:,0],axis=0)
        amp_cal_11_all = np.median(amps_array[antenna_id,:,:,1],axis=0)
        # interpolate to target frequencies
        amp_cal_00 = np.interp(msinfo.msfreqvalues, freqs_ampl, amp_cal_00_all)
        amp_cal_11 = np.interp(msinfo.msfreqvalues, freqs_ampl, amp_cal_11_all)
        # interpolate phases
        phase_cal_00   = 0.
        phase_cal_11   = np.interp(msinfo.msfreqvalues, freqs_phase, phases_array[:,antenna_id])

        # convert to real and imaginary
        real_00 = amp_cal_00*np.cos(phase_cal_00)
        imag_00 = amp_cal_00*np.sin(phase_cal_00)
        real_11 = amp_cal_11*np.cos(-1.*phase_cal_11)
        imag_11 = amp_cal_11*np.sin(-1.*phase_cal_11)

        real_00_pdb = real_00.reshape( (ntimes,nfreqs) )
        imag_00_pdb = imag_00.reshape( (ntimes,nfreqs) )
        real_11_pdb = real_11.reshape( (ntimes,nfreqs) )
        imag_11_pdb = imag_11.reshape( (ntimes,nfreqs) )

        # generate parmDB entries
        ValueHolder = outDB.makeValue(values=real_00_pdb,
                                      sfreq=startfreqs, efreq=endfreqs,
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Gain:0:0:Real:'+antenna,ValueHolder)
        ValueHolder = outDB.makeValue(values=imag_00_pdb,
                                      sfreq=startfreqs, efreq=endfreqs,
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Gain:0:0:Imag:'+antenna,ValueHolder)
        ValueHolder = outDB.makeValue(values=real_11_pdb,
                                      sfreq=startfreqs, efreq=endfreqs,
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Gain:1:1:Real:'+antenna,ValueHolder)
        ValueHolder = outDB.makeValue(values=imag_11_pdb,
                                      sfreq=startfreqs, efreq=endfreqs,
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Gain:1:1:Imag:'+antenna,ValueHolder)

        #now handle the clock-value (no fancy interpolating needed)
        clock_pdb = np.array( np.median(clock_array[:,antenna_id]) ,ndmin=2)
        ValueHolder = outDB.makeValue(values=clock_pdb,
                                      sfreq=startfreqs[0], efreq=endfreqs[-1],
                                      stime=starttime, etime=endtime, asStartEnd=True)
        outDB.addValues('Clock:'+antenna,ValueHolder)
        
        rotation_angles = np.outer(rmdict['RM'][antenna],lambdaSquared)
        newValue = outDB.makeValue(values=rotation_angles, sfreq=msinfo.msfreqvalues, efreq=freqsteps, stime=rmdict['times'], etime=timesteps, asStartEnd=False)
        outDB.addValues('CommonRotationAngle:'+antenna,newValue)


    outDB.flush()
    outDB = False
    return {'transfer_parmDB': newparmDB }


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create a parmDB with values from the calibrator and rotaton values from RMextract.')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One or more MSs for which the IONEX data should be downloaded.')
    parser.add_argument('--server', type=str, default='None',               
                        help='URL of the server to use. (default: None)')
    parser.add_argument('--prefix', type=str, default='CODG', 
                        help='Prefix of the IONEX files. (default: \"CODG\")')
    parser.add_argument('--ionexpath', '--path', type=str, default='IONEXdata/',
                        help='Directory where to store the IONEX files. (default: \"IONEXdata/\")')
    parser.add_argument('--basename', type=str, default='caldata_transfer',
                        help='Base-name of the numpy files with the calibrator values. (default: \"caldata_transfer\")')
    parser.add_argument('--storedir', type=str, default='.',
                        help='Directory of the numpy files with the calibrator values. (default: \".\")')
    parser.add_argument('--extension', type=str, default='-instrument_amp_clock_offset',
                        help='Extension to the MS-name to get the name of the parmDB. (default: \"-instrument_amp_clock_offset\")')


    args = parser.parse_args()

    for MS in args.MSfile:
        print "Working on:", MS
        #main(MS, args.basename, args.extension, server=args.server, prefix=args.prefix, ionexPath=args.ionexpath)
        #msinfo = ReadMs(MS)
        #newparmDB = MS+args.extension
        #outDB = make_empty_parmdb(newparmDB)
        #add_COMMONROTATION_vals(outDB, msinfo, args.server, args.prefix, args.ionexpath)
        main(MS, store_basename=args.basename, store_directory=args.storedir, newparmdbext=args.extension, 
         ionex_server=args.server, ionex_prefix=args.prefix, ionexPath=args.ionexpath)
