#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import json

###############################################################################
def main(pipeline = 'prefactor', run_type = 'calibrator', filter = '[CR]S*&', bad_antennas = '[CR]S*&', output_fname = 'summary.json', structure_file = None, Ateam_separation_file = None):
	"""
	Creates summary of a given prefactor3-CWL run
	
	Parameters
	----------
	filter: antenna string which is pre-selected (by user or pipeline)
	bad_anntenas: pre-selected antennas and removed stations (separated by ";")
	"""
	# location of logfile
	print('Summary logfile is written to ' + output_fname)
	
	bad_antennas_list = bad_antennas.lstrip(filter).replace('!','').replace('*','').replace('&','').split(';')

	## define contents of JSON file
	json_output = { 'metrics': { pipeline : { } } }
	json_output['metrics'][pipeline]['run_type'] = run_type

	## get diffractive_scale info
	if structure_file:
		diffractive_scale = { 'unit' : 'km'}
		with open(structure_file, 'r') as infile:
			for line in infile:
				diffractive_scale['XX'] = float(line.split()[1].replace('*',''))
				diffractive_scale['YY'] = float(line.split()[2])
				break
		json_output['metrics'][pipeline]['diffractive_scale'] = diffractive_scale

	## get Ateam_separation information
	if Ateam_separation_file:
		f = open(Ateam_separation_file, 'r')
		json_output['metrics'][pipeline]['close_sources'] = json.load(f)
		if len(json_output['metrics'][pipeline]['close_sources']) > 0:
			Ateam_list = ''
			for i in range(len(json_output['metrics'][pipeline]['close_sources'])):
				Ateam_list += i['source'] + ','
		else:
			Ateam_list = 'NONE'
		print('A-Team sources close to the phase reference center: ' + Ateam_list.rstrip(','))
        
	## write JSON file
	with open(output_fname, 'w') as fp:
		json.dump(json_output, fp)

	print('Summary has been created.')
	return(0)


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description='Creates summary of a given prefactor3-CWL run.') 
	parser.add_argument('--pipeline', type=str, default='prefactor', help='Name of the pipeline.')
	parser.add_argument('--run_type', type=str, default='calibrator', help='Type of the pipeline')
	parser.add_argument('--filter', type=str, default='[CR]S*&', help='Filter these antenna string from the processing.')
	parser.add_argument('--bad_antennas', type=str, default='[CR]S*&', help='Antenna string to be processed')
	parser.add_argument('--output_fname', '--output_fname', type=str, default='summary.json', help='Name of the output filename (default: summary.json)')
	parser.add_argument('--structure_file', type=str, default=None, help='Location of the structure function logfile')
	parser.add_argument('--Ateam_separation_file', type=str, default=None, help='Location of the Ateam_separation JSON file')
	
	args = parser.parse_args()
	
	# start running script
	main(args.pipeline, args.run_type, args.filter, args.bad_antennas, args.output_fname, args.structure_file)
	
	sys.exit(0)
