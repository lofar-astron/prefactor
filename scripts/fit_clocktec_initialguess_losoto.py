#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import lofar.parmdb
import sys
import os
import scipy
import time
import numpy
import math
import pyrap.tables
import pp
import scipy.signal
from pylab import *
import sys, os, glob, re
import numpy as np
import shutil
import progressbar
import logging
import pyrap.tables as pt
import lofar.parmdb
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter, solFetcher
import multiprocessing as mp
#import lofar.expion.fitting as fitting

args = sys.argv
globaldbname=args[1]
calsource=args[2]
ncpus=int(args[3])
########################
###### USER INPUT ######

#globaldbname = 'L128487.h' # input h5 parm file
#calsource    = '3C295' # name for writing outputfiles
#ncpus        = 24 # number of CPUs avaulable for parallel fitting

#### END USER INPUT ####
########################


pi = numpy.pi
c  = 2.99792458e8
ionmodel = h5parm(globaldbname ,readonly=True)
solsetNames = ionmodel.getSolsets()
for solsetName in solsetNames:
    print solsetName
solset = ionmodel.getSolset(solsetName)
# solset name seems to always be sol000
soltabs = ionmodel.getSoltabs('sol000')

amptab = ionmodel.getSoltab('sol000','amplitude000')
phasetab = ionmodel.getSoltab('sol000','phase000')
rottab =  ionmodel.getSoltab('sol000','rotation000')
anttab = ionmodel.getAnt('sol000')
source_id     = 0  # source ID in global_db (usually 0)


def median_window_filter(ampl, half_window, threshold):
    ampl_tot_copy = numpy.copy(ampl)
    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]

    #fix oct 2012
    median_array  = scipy.signal.medfilt(sol,int(half_window*2.-1))

    sol_flag = numpy.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = numpy.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = numpy.median(window_masked)
        q = 1.4826 * numpy.median(numpy.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.
        if abs(sol[i] - median) > (threshold * q):
            sol_flag[i] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
           ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy
    
    
def running_median(ampl,half_window) :

    ampl_tot_copy = numpy.copy(ampl)

    ndata = len(ampl)
    flags = numpy.zeros(ndata, dtype=bool)
    sol = numpy.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl
    std = numpy.zeros(len(ampl))

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]
    
    for i in range(len(ampl)):
      #print i, i+half_window
      std[i] =  numpy.median(sol[i:i+(2*half_window)])  

    return std

def fit_dTEC_dclock_dFR(phases_rr, phases_ll, amp_rr, amp_ll, freq, distance_station):
      c = 2.99792458e8
      freq_old = numpy.copy(freq)
      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
      par3complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))
      par2complex     = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y)) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))
      rmwavcomplex    = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y)) + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))

      # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll) 
      par3complex_w   = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.cos(y))*(freq/1e5) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - numpy.sin(y))*(freq/1e5)
      par2complex_w   = lambda p, freq, y: abs(numpy.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.cos(y))*(freq/1e5) + abs(numpy.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - numpy.sin(y))*(freq/1e5)
      rmwavcomplex_w  = lambda RM, wav, y: abs(numpy.cos(2.*RM[0]*wav*wav)  - numpy.cos(y))/wav + abs(numpy.sin(2.*RM[0]*wav*wav)  - numpy.sin(y))/wav

      plotrm          = lambda RM, wav: numpy.mod( (2.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
      fitfuncfastplot    = lambda p, freq: numpy.mod((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 


      idx_rr    = numpy.where(amp_rr != 1.0)
      if len(idx_rr) != 0:
       freq       = freq[:][idx_rr]
       phases_rr  = phases_rr[:][idx_rr]
       phases_ll  = phases_ll[:][idx_rr]
       amp_ll     = amp_ll[:][idx_rr] 

      idx_ll     = numpy.where(amp_ll != 1.0)
      if len(idx_ll) != 0:
       freq       = freq[:][idx_ll]
       phases_rr  = phases_rr[:][idx_ll]
       phases_ll  = phases_ll[:][idx_ll]
   
      if len(freq) < len(freq_old):
        print 'Number of filtered out data points:',  len(freq_old)-len(freq)
      # ---------- end filter bad data -----------

      if len(freq) != 0: # prepare and make arrays if there is valid data 
        freq      = (freq[0:len(freq_old)])
        phases_ll = (phases_ll[0: len(freq_old)])
        phases_rr = (phases_rr[0: len(freq_old)])
      
        phase       = (phases_rr + phases_ll)      # not divide by 2, then later fix this
        #phase       = (phases_ll + phases_ll)  # temp, just fit RR for now 
        phase_diff  = (phases_rr - phases_ll)      # not divide by 2, then later fix this
      
        wav      = c/freq
        pi = numpy.pi
        chi_old=1e9 
      else:
        phase = (phases_rr[:] + phases_ll[:])
	phase_diff = (phases_rr[:] - phases_ll[:])
      
      if len(freq) > 10: 
      
	# FIND INTIAL GUESS
	for dTEC in numpy.arange(-1.0,1.0, 0.01):
	 for dclock in numpy.arange(-200e-9,200e-9,5e-9):
            phase_model = numpy.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC]
	      #print 'Better fit', dclock, dTEC

	fitguess_1 = numpy.copy(fitguess)
	#print 'iter 1', fitguess

	for dTEC in numpy.arange(fitguess_1[1]-0.02,fitguess_1[1]+0.02, 0.002):
	 for dclock in numpy.arange(fitguess_1[0]-8e-9,fitguess_1[0]+ 8e-9,1e-9):
            phase_model = numpy.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = numpy.mod (phase, 2*pi)
	    angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	    chi = numpy.sum(angle)

	    if chi < chi_old:
	      chi_old = chi
	      fitguess     = [dclock,dTEC]
	      #print 'Better fit', dclock, dTEC


	#print 'iter 2', fitguess

	chi_old = 1e9
	for dFR in numpy.arange(-0.1,0.1,2e-4):
          phase_model = numpy.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
	  phase_data  = numpy.mod (phase_diff, 2*pi)
          angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	  chi = numpy.sum(angle)

	  if chi < chi_old:
	    chi_old = chi
	    fitrmguess     = dFR
	    #print 'Better fit', fitrmguess

	fitrmguess_1 = numpy.copy(fitrmguess)
	for dFR in numpy.arange(fitrmguess_1-5e-4,fitrmguess_1+5e-4,0.5e-5):
          phase_model = numpy.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
	  phase_data  = numpy.mod (phase_diff, 2*pi)
          angle       = pi - numpy.abs(numpy.abs(phase_model - phase_data) - pi)
	  chi = numpy.sum(angle)

	  if chi < chi_old:
	    chi_old = chi
	    fitrmguess     = dFR
	    #print 'Better fit', fitrmguess

        # DO THE FITTING 
        # SOLVE Clock-TEC anticorrelation problem on short baselines             
        freq = freq.astype(numpy.float64)
	phase= phase.astype(numpy.float64) #epsfcn=1e-7
      
        if distance_station < 0. :   #15.0*1e3: DOES NOT WORK, NEED 3 par FIT
          fitresult, success  = scipy.optimize.leastsq(par2complex, fitguess,    args=(freq, phase))
	  #fitresult = fitguess
        else:
          fitresult, success   = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1], 0.0], args=(freq, phase),maxfev=10000)
	  #fitresult = [fitguess[0], fitguess[1], 0.0]
          #print fitresult, success
        fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))

      else:
         print 'No valid data found'
         fitresult = [0.0,0.0, 0.0]
	 fitresultrm_wav= 0.0
            
      
      
   
      show_plot = False
      if show_plot:  
        
	#if len(fitresult ==2): 
	#  fitresult =[fitresult[0], fitresult[1],0]
	#  print 'Here'
	#  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]
	
	#fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]
	
        matplotlib.pyplot.plot(freq, numpy.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        matplotlib.pyplot.plot(freq, numpy.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )			   

        TEC   = numpy.mod((-8.44797245e9*(2.*fitresult[1])/freq)+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
        Clock = numpy.mod((2.*numpy.pi*   2.*fitresult[0]*freq )+numpy.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
    
        phase_total = (2.*numpy.pi*2.*fitresult[0]*freq)+(-8.44797245e9*(2.*fitresult[1])/freq)+fitresult[2]
        residual    = numpy.mod(phase-phase_total+pi,2.*pi)-pi
        matplotlib.pyplot.plot(freq, residual, '.', color='yellow')
    
        idxl = int(min(freq_old)/1e4) 
	idxh = int(max(freq_old)/1e4) 
        bigfreqaxistmp = range(idxl, idxh)
	bigfreqaxis    =  numpy.array([float(i) for i in bigfreqaxistmp])
	bigfreqaxis    = bigfreqaxis*1e4
	
	matplotlib.pyplot.plot (bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")	
	matplotlib.pyplot.plot (bigfreqaxis, plotrm(fitresultrm_wav, c/bigfreqaxis[:]), "-", color='purple')
	#matplotlib.pyplot.plot (freq, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")	
	
	matplotlib.pyplot.plot(freq, Clock,  ',g') 
	matplotlib.pyplot.plot(freq, TEC,  ',b') 
	matplotlib.pyplot.xlabel('freq')
	matplotlib.pyplot.ylabel('phase')

	matplotlib.pyplot.show()
      
      return [fitresult[0], fitresult[1],fitresult[2], fitresultrm_wav]
   





time_id       = 100
pol_id        = 0
#antenna_id    = 29
refantenna_id = 0
source_id     = 0
goodstartidx  = 0

CScorrect     = False

if CScorrect:
  csclockvals = numpy.load('../CS_clocks.npy')

A = numpy.zeros((len(amptab.freq[:]), 2), dtype = float)
A[:,0] = amptab.freq[:]*2*pi
A[:,1] = -8.44797245e9/amptab.freq[:]
sol = numpy.zeros((len(amptab.ant), 2))

clockfit = 0.*numpy.copy(amptab.time)
TECfit   = 0.*numpy.copy(amptab.time)
RMfit   = 0.*numpy.copy(amptab.time)
phaseoffset = 0.*numpy.copy(amptab.time)

clockarray = numpy.zeros([len(amptab.time),len(amptab.ant)])
tecarray   = numpy.zeros([len(amptab.time),len(amptab.ant)])
rmarray    = numpy.zeros([len(amptab.time),len(amptab.ant)])
phaseoffsetarray = numpy.zeros([len(amptab.time),len(amptab.ant)])

print 'REF STATION:', amptab.ant[refantenna_id]
print '# TIMESLOTS  ', len(amptab.time)
print '# FREQUENCIES', len(amptab.freq)

ppservers = ()
# Creates jobserver with ncpus workers
mysecret = 'fit_clocktec'+str(os.getpid())
job_server = pp.Server(ncpus, ppservers=ppservers, secret=mysecret)
print "Starting pp with", job_server.get_ncpus(), "workers"


phases_all = numpy.copy(phasetab.val)
#phases_all = numpy.load('../phases_3C196.npy')


start_time_id = 0
stop_time_id  = len(phasetab.val[0,source_id,0,0,:])

if CScorrect:
 # CORRECT CS CLOCKS
 for time_id in range(start_time_id,stop_time_id):
   print time_id
   phases_rr = numpy.copy(phases_all[0,source_id,:,:,time_id])
   phases_ll = numpy.copy(phases_all[1,source_id,:,:,time_id])

   #RR  correct
   for ss in range(0,len(phasetab.ant[:])):
         sol[ss,0] = (csclockvals[ss,0]) # clock RR
         sol[ss,1] = 0.0  # TEC
	 #print 'RR', sol[ss,0]
   phases_rr = phases_rr - dot(A, sol.T)

   #LL  correct
   for ss in range(0,len(phasetab.ant[:])):
       sol[ss,0] = (csclockvals[ss,1]) # clock LL
       sol[ss,1] = 0.0  # TEC 
       #print 'LL', sol[ss,0]  
   phases_ll = phases_ll - dot(A, sol.T) 

   phases_all[time_id,:,:,source_id,0] = numpy.copy(phases_rr)
   phases_all[time_id,:,:,source_id,1] = numpy.copy(phases_ll)


# count number of RS stations
N_RS = 0
for ss in range(0,len(phasetab.ant[:])):
  sname = phasetab.ant[ss]
  if sname[0:2] == 'RS':
    N_RS=N_RS+1

#for antenna_id in range(3,  4):
#for antenna_id in range(len(ionmodel.stations[:])-N_RS,  len(ionmodel.stations[:])) :

freq       = numpy.copy(phasetab.freq)

for antenna_id in range(1,  len(phasetab.ant[:])) :



 if antenna_id != refantenna_id:
   stationspos      =  anttab[phasetab.ant[refantenna_id]] - anttab[phasetab.ant[antenna_id]]
   distance_station = numpy.sqrt(stationspos[0]**2 + stationspos[1]**2 + stationspos[2]**2)
   print 'Distance stations to reference station', distance_station/1e3, ' km'

   jobs = []

   for time_id in range(start_time_id,stop_time_id):
   #for time_id in range(3050,3051):
     
          
     
     phases_rr = numpy.copy((phases_all[0,source_id,antenna_id,:,time_id] \
	  -phases_all[0,source_id,refantenna_id,:,time_id])) 
     phases_ll = numpy.copy((phases_all[1,source_id,antenna_id,:,time_id] \
	  -phases_all[1,source_id,refantenna_id,:,time_id]))

     # -------- filter bad data: where amplitudes equal 1.0 -----------
     amp_rr = numpy.copy(amptab.val[0,source_id,antenna_id,:,time_id]) 
     amp_ll = numpy.copy(amptab.val[0,source_id,antenna_id,:,time_id])
                 
     #print phases_rr,phases_ll
     #print fit_dTEC_dclock_dFR(phases_rr, phases_ll, amp_rr, amp_ll, freq, distance_station),'testing non parallel'
     jobs.append(job_server.submit(fit_dTEC_dclock_dFR,(phases_rr, phases_ll, amp_rr, amp_ll, freq, distance_station),(), ("numpy","scipy.optimize","matplotlib.pyplot",)))
     print 'Submitting job #', time_id

   i= 0
   for job in jobs:
     fitresult = job()
     print fitresult,'sdfasfasfasdf'
     clockfit[start_time_id+i]  = fitresult[0]
     TECfit[start_time_id+i]    = fitresult[1]
     RMfit[start_time_id+i]     = fitresult[3]
     phaseoffset[start_time_id+i]  =  fitresult[2]

     clockarray[start_time_id+i,antenna_id] = fitresult[0]
     tecarray[start_time_id+i,antenna_id]   = fitresult[1]
     rmarray[start_time_id+i,antenna_id]    = fitresult[3]
     phaseoffsetarray[start_time_id+i,antenna_id] = fitresult[2] 
     print 'TIME_ID', start_time_id+i, 'FIT (dclock, dTEC, offset, dRM)', fitresult
     i = i + 1



os.system('rm -f ' + 'fitted_data_dclock_' + calsource + '_1st.npy')
os.system('rm -f ' + 'fitted_data_dTEC_'   + calsource + '_1st.npy')
numpy.save('fitted_data_dclock_' + calsource + '_1st.npy', clockarray)
numpy.save('fitted_data_dTEC_'   + calsource + '_1st.npy', tecarray)



for antenna_id in range(0, len(amptab.ant[:])):
   print 'Cleaning up Clock and TEC values for: ', amptab.ant[antenna_id]
   clockfit = clockarray[:,antenna_id]
   TECfit   = tecarray[:,antenna_id]
 
   clockfit = median_window_filter(clockfit, 7, 5)
   TECfit   = median_window_filter(TECfit, 5, 5)

   clockfit = median_window_filter(clockfit, 5, 3)
   TECfit   = median_window_filter(TECfit, 5, 3)
   
   clockfit = median_window_filter(clockfit, 3, 3)
   TECfit   = median_window_filter(TECfit, 3, 3)

   clockfit = running_median(clockfit, 3)
   TECfit = running_median(TECfit, 3)

   
   clockarray[:,antenna_id] = clockfit
   tecarray[:,antenna_id]   = TECfit
   #matplotlib.pyplot.plot(clockarray[:,antenna_id])
   

os.system('rm -f ' + 'fitted_data_dclock_' + calsource + '_1st.sm.npy')
os.system('rm -f ' + 'fitted_data_dTEC_'   + calsource + '_1st.sm.npy')
numpy.save('fitted_data_dclock_' + calsource + '_1st.sm.npy', clockarray)
numpy.save('fitted_data_dTEC_'   + calsource + '_1st.sm.npy', tecarray)
      
 
#matplotlib.pyplot.savefig('fit_clocktec.png')
                          

