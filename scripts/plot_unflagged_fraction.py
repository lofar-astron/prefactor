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
import numpy as np
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
   return(str_list)


def main(ms_list, frac_list, outfile='unflagged_fraction.png'):

    ms_list = input2strlist_nomapfile(ms_list)
    frac_list = input2strlist_nomapfile(frac_list)
    frac_list = np.array([float(f) for f in frac_list])
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Get frequencies
    freq_list = []
    for ms in ms_list:
        # open the main table and print some info about the MS
        t = pt.table(ms, readonly=True, ack=False)
        tfreq = pt.table(t.getkeyword('SPECTRAL_WINDOW'),readonly=True,ack=False)
        ref_freq = tfreq.getcol('REF_FREQUENCY',nrow=1)[0]
        freq_list.append(ref_freq)
    freq_list = np.array(freq_list) / 1e6  # MHz

    # Plot the unflagged fraction vs. frequency
    plt.scatter(freq_list, frac_list)
    plt.xlabel('frequency [MHz]')
    plt.ylabel('unflagged fraction')
    plt.savefig(outfile)
