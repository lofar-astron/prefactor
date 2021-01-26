#!/usr/bin/env python

## changelog
# W.Williams 2014/11/03  add - to give input/output statistics per channel
# W.Williams 2014/11/03  fix - statistics per correlation
# A.Drabent 2019/07/24   write fraction of flagged data into output file (for prefactor3)

import numpy
import pyrap.tables as pt
import os
import sys

msname = str(sys.argv[1])

cliplevelhba = 5.0
cliplevellba = 50.0

t = pt.table(msname, readonly=False)
data = t.getcol('MODEL_DATA')
flag = t.getcol('FLAG')
freq_tab= pt.table(msname + '::SPECTRAL_WINDOW')
freq    = freq_tab.getcol('REF_FREQUENCY')

if freq[0] > 100e6:
 cliplevel = cliplevelhba
if freq[0] < 100e6:
 cliplevel = cliplevellba


print('------------------------------')
print('SB Frequency [MHz]', freq[0]/1e6)
for chan in range(0,numpy.size(data[0,:,0])):
  print('chan %i : %.5f%% input XX flagged' %( chan, 100.*numpy.sum(flag[:,chan,0] == True)/numpy.size(flag[:,chan,0]) ))
  print('chan %i : %.5f%% input YY flagged' %( chan, 100.*numpy.sum(flag[:,chan,3] == True)/numpy.size(flag[:,chan,3]) ))
input_flags_xx = 100. * numpy.sum(flag[:,:,0] == True)/numpy.size(flag[:,:,0])
input_flags_yy = 100. * numpy.sum(flag[:,:,3] == True)/numpy.size(flag[:,:,3])
print('Total : %.5f%% input XX flagged' %(   input_flags_xx ))
print('Total : %.5f%% input YY flagged' %(   input_flags_yy ))
print('')
print('Cliplevel used [Jy]', cliplevel)
print('\n\n')

for pol in range(0,numpy.size(data[0,0,:])):
 for chan in range(0,numpy.size(data[0,:,0])):
  print('Doing polarization,chan', pol, chan)
  idx = numpy.where(abs(data[:,chan,pol]) > cliplevel)
  flag[idx,chan,0] = True
  flag[idx,chan,1] = True
  flag[idx,chan,2] = True
  flag[idx,chan,3] = True 

print('')
for chan in range(0,numpy.size(data[0,:,0])):
  print('chan %i : %.5f%% output XX flagged' %( chan, 100.*numpy.sum(flag[:,chan,0] == True)/numpy.size(flag[:,chan,0]) ))
  print('chan %i : %.5f%% output YY flagged' %( chan, 100.*numpy.sum(flag[:,chan,3] == True)/numpy.size(flag[:,chan,3]) ))
output_flags_xx = 100. * numpy.sum(flag[:,:,0] == True)/numpy.size(flag[:,:,0])
output_flags_yy = 100. * numpy.sum(flag[:,:,3] == True)/numpy.size(flag[:,:,3])
print('Total : %.5f%% output XX flagged' %(   output_flags_xx ))
print('Total : %.5f%% output YY flagged' %(   output_flags_yy ))
print('')
t.putcol('FLAG', flag)
t.close()
freq_tab.close()

os.system('echo ' + str(freq[0]) + ' ' + str(output_flags_xx - input_flags_xx) + ' ' + str(output_flags_yy - input_flags_yy) + ' >> Ateamclipper.txt')
