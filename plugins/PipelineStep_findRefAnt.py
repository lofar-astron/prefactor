#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import argparse
import re
import multiprocessing
import logging
from lofarpipe.support.data_map import DataMap, DataProduct

###############################################################################
def find_flagged_fraction(ms_file):
    
   outputs  = os.popen('DPPP msin=' + ms_file + ' msout=. steps=[count] count.type=counter count.warnperc=0.0000001 | grep NOTE').readlines()
   fraction_flagged = { output.split('(')[-1].rstrip(')\n'):output.split('%')[0][-5:].strip() for output in outputs if 'station' in output }
   return(fraction_flagged)

###############################################################################
def plugin_main(args, **kwargs):
    
	mapfile_in     = kwargs['mapfile_in']
	station_filter = kwargs['station_filter']
	data           = DataMap.load(mapfile_in)
	mslist         = [data[i].file for i in range(len(data))]
	
	## derive the fraction of flagged data of the entire observation
	print('Reading data.')
	logging.info('Reading data.')
	pool = multiprocessing.Pool(processes = multiprocessing.cpu_count())
	flagged_fraction_dict = pool.map(find_flagged_fraction, mslist)
	
	print('Apply station filter ' + str(station_filter))
	logging.info('Apply station filter ' + str(station_filter))
	flagged_fraction_data = {}
	for entry in flagged_fraction_dict:
		antennas = entry.keys()
		selected_stations = [ station_name for station_name in antennas if re.match(station_filter, station_name) ]
		if len(selected_stations) == 0:
			logging.error('No stations left after filtering.')
			return(1)
		for antenna in selected_stations:
			try:
				flagged_fraction_data[antenna].append(float(entry[antenna]))
			except KeyError:
				flagged_fraction_data[antenna] = [float(entry[antenna])]

	flagged_fraction_list = []
	sorted_stations = sorted(flagged_fraction_data.keys())
	for antenna in sorted_stations:
		flagged_fraction = sum(flagged_fraction_data[antenna]) / len(flagged_fraction_data[antenna])
		flagged_fraction_list.append(flagged_fraction)
		try:
			flagged_fraction_data[flagged_fraction].append(antenna)
		except KeyError:
			flagged_fraction_data[flagged_fraction] = [antenna]
		
		
	min_flagged_fraction = min(flagged_fraction_list)
	refant = flagged_fraction_data[min_flagged_fraction][0]
	logging.info('Selected station ' + str(refant) + ' as reference antenna. Fraction of flagged data is ' + '{:>3}'.format('{:.1f}'.format(min_flagged_fraction) + '%'))
	print('Selected station ' + str(refant) + ' as reference antenna. Fraction of flagged data is ' + '{:>3}'.format('{:.1f}'.format(min_flagged_fraction) + '%'))
    
	## return results
	result = {'refant':str(refant)}
	return(result)
