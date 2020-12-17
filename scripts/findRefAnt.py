#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import argparse
import logging

###############################################################################
def find_flagged_fraction(ms_file):
    
   outputs  = os.popen('DPPP msin=' + ms_file + ' msout=. steps=[count] count.type=counter count.warnperc=0.0000001 | grep NOTE').readlines()
   fraction_flagged = { output.split('(')[-1].rstrip(')\n'):output.split('%')[0][-5:].strip() for output in outputs if 'station' in output }
   return(fraction_flagged)

###############################################################################
def main(ms_file):

	## derive the fraction of flagged data of the entire observation
	logging.info('Reading data.')
	flagged_fraction_dict = find_flagged_fraction(ms_file)
	
	return(flagged_fraction_dict)

if __name__ == '__main__':
    descriptiontext = "Check the flagged fraction of stations in a MS.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inputms', help='name of the input MS')
    args = parser.parse_args()

    erg = main(ms_file = args.inputms)
