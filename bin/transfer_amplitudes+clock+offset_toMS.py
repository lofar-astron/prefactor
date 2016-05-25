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
	freq=pt.table(t.getkeyword("SPECTRAL_WINDOW"))
	self.fullband = freq.getcell('TOTAL_BANDWIDTH', 0)
	self.freqpara['cent'] = freq.getcell('REF_FREQUENCY', 0)
	self.freqpara['step'] = freq.getcell('CHAN_WIDTH', 0)[0]
        self.msfreqvalues = freq.getcell('CHAN_FREQ', 0)
	self.freqpara['start'] = self.msfreqvalues[0]-self.freqpara['step']/2.
	self.freqpara['end'] = self.msfreqvalues[-1]+self.freqpara['step']/2.
	freq.close()
        ##########Getting Station Names###################
        antennas = pt.table(t.getkeyword("ANTENNA"))
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

def main(msname, store_basename, newparmdbext='-instrument_amp_clock_offset'):

    # name (path) for parmdb to be written
    newparmDB = msname+newparmdbext

    # load the numpy arrays written by the previous scripts
    # (filenames constructed in the same way as in these scripts)
    freqs_ampl = np.load('freqs_for_amplitude_array.npy')
    amps_array = np.load(store_basename + '_amplitude_array.npy')
    clock_array = np.load('fitted_data_dclock_' + store_basename + '_1st.npy')
    freqs_phase = np.load('freqs_for_phase_array.npy')
    phases_array  = np.load(store_basename + '_phase_array.npy')
    station_names = np.load(store_basename + '_station_names.npy')

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

    outDB = False
    return {'transfer_parmDB': newparmDB }
