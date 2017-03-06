#! /usr/bin/env python
"""
Script to sort a list of MSs by into frequency groups by time-stamp
"""
import pyrap.tables as pt
import sys, os
import re

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    descriptiontext = ("Add station-subband number to the file names.\n"
                       "It will be added as \"SSBNNN\" before the subband-number \"SBMMM\".\n"
                       "(The station-subband is a measure for the observed frequency)\n")
    
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_files', nargs='+', help='list of ms files to rename')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Go Vebose! (default=False)')

    args = parser.parse_args()

    # inverse of: "_SB\d{3}"
    subReg2 = re.compile(r'\d{3}BS_')
    for ms in args.ms_files:
        sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
        freq = sw.col('REF_FREQUENCY')[0]
        sw.close()

        station_subband = int(freq/100e6*512.)%512
        ssb_string = "_SSB%03d"%(station_subband)

        matchobj = subReg2.search(ms[::-1])
        if matchobj:
            splitpos = -matchobj.end()
        else:
            if args.verbose:
                print "Couldn\'t find subband number in filename \""+ms+"\", attaching station-subband number to end."
            splitpos = len(ms)
        newname = ms[:splitpos]+ssb_string+ms[splitpos:]
        if args.verbose:
            print "Renaming: \""+ms+"\" to: \""+newname+"\""
        os.rename(ms,newname)
