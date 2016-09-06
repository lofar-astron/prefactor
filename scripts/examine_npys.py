#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import os,sys
import pylab

args=sys.argv
########################
###### USER INPUT ######

calsource    = args[1] # name for writing outputfiles

#### END USER INPUT ####
########################


amparray   = np.load(calsource + '_amplitude_array.npy')
clockarray = np.load('fitted_data_dclock_' + calsource + '_1st.sm.npy')
dtecarray  = np.load('fitted_data_dTEC_'   + calsource + '_1st.sm.npy')
numants = len(dtecarray[0,:])
#numfreqs = 242
#numtimes = 75


#print np.shape(dtecarray)

#sys.exit()


for i in range(0,numants):
    pylab.plot(dtecarray[:,i])
pylab.xlabel('Time')
pylab.ylabel('dTEC [$10^{16}$ m$^{-2}$]')
pylab.savefig('dtec_allsols.png')
pylab.close()
pylab.cla()

for i in range(0,numants):
    pylab.plot(1e9*clockarray[:,i])
pylab.xlabel('Time')
pylab.ylabel('dClock [ns]')
pylab.savefig('dclock_allsols.png')
pylab.close()
pylab.cla()


for i in range(0,numants):
  pylab.plot(np.median(amparray[i,:,:,0], axis=0))
  pylab.plot(np.median(amparray[i,:,:,1], axis=0))
pylab.xlabel('Subband number')
pylab.ylabel('Amplitude')
pylab.ylim(0,2.*np.median(amparray))
pylab.savefig('amp_allsols.png')
pylab.close()
pylab.cla()


#sys.exit(0)
