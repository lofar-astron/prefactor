#!/usr/bin/env python
# -* coding: utf-8 -*-

"""
Adds phases (=0) and amplitudes (=1) to any missing station if they appear in an h5parm, but not in a particular soltab. 

Created on Tue Jul 24 2019

@author: Alexander Drabent
"""

import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
import numpy

def main(txtfile = 'Ateamclipper.txt', outfile = 'Ateamclipper.png'):

    frac_list_xx = []
    frac_list_yy = []
    freq_list    = []
    with open(txtfile, 'r') as infile:
        for line in infile:
            freq_list.append(float(line.split()[0]))
            frac_list_xx.append(float(line.split()[1]))
            frac_list_yy.append(float(line.split()[2]))
    
    # Plot the amount of clipped data vs. frequency potentially contaminated by the A-team
    plt.scatter(numpy.array(freq_list) / 1e6, numpy.array(frac_list_xx), marker = '.', s = 10)
    plt.xlabel('frequency [MHz]')
    plt.ylabel('A-team clipping fraction [%]')
    plt.savefig(outfile)
    return(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Adds phases and amplitudes to any missing station if they appear in an h5parm, but not in a particular soltab.')

    parser.add_argument('txtfile', type=str,
                        help='Input text file containing frequency and flag fraction of the XX and YY polarization.')
    parser.add_argument('outfile', type=str,
                        help='Input text file containing frequency and flag fraction of the XX and YY polarization.')


    args = parser.parse_args()

    main(txtfile = args.txtfile, outfile = args.outfile)

