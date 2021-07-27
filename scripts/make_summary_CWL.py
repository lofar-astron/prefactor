#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import json

###############################################################################
def main(flagFiles = None, pipeline = 'prefactor', run_type = 'calibrator', filtered_antennas = '[CR]S*&', bad_antennas = '[CR]S*&', output_fname = 'summary.json', structure_file = None, Ateam_separation_file = None):
	"""
	Creates summary of a given prefactor3-CWL run
	
	Parameters
	----------
	pipeline: name of the pipeline
	filter: antenna string which is pre-selected (by user or pipeline)
	bad_anntenas: pre-selected antennas and removed stations (separated by ";")
	
	"""
	# location of logfile
	header_string = '*** ' + pipeline + ' ' + run_type + ' pipeline summary ***'
	print('*' * len(header_string) + '\n' + header_string + '\n' + '*' * len(header_string) + '\n')
	print('Summary logfile is written to ' + output_fname + '\n')
	
	## define contents of JSON file
	json_output = { 'metrics': { pipeline : { } } }
	json_output['metrics'][pipeline]['run_type'] = run_type
	
	## print antennas removed from the data
	bad_antennas_list = list(filter(None, set(bad_antennas.replace(filtered_antennas,'').replace('!','').replace('*','').replace('&','').split(';'))))
	if bad_antennas_list == []:
		print('Antennas removed from the data: NONE')
	else:
		print('Antennas removed from the data: ' + ', '.join(bad_antennas_list))
	json_output['metrics'][pipeline]['stations'] = []
	for bad_antenna in bad_antennas_list:
		json_output['metrics'][pipeline]['stations'].append({'station' : bad_antenna, 'removed' : 'yes'})
	
	## get Ateam_separation info
	if Ateam_separation_file:
		f = open(Ateam_separation_file, 'r')
		json_output['metrics'][pipeline]['close_sources'] = json.load(f)
		f.close()
		if len(json_output['metrics'][pipeline]['close_sources']) > 0:
			Ateam_list = ''
			for i in range(len(json_output['metrics'][pipeline]['close_sources'])):
				Ateam_list += i['source'] + ','
		else:
			Ateam_list = 'NONE'
		print('A-Team sources close to the phase reference center: ' + Ateam_list.rstrip(',') + '\n')

	## get diffractive_scale info
	if structure_file:
		diffractive_scale = { 'unit' : 'km'}
		with open(structure_file, 'r') as infile:
			for line in infile:
				diffractive_scale['XX'] = float(line.split()[1].replace('*','')) / 1000.
				diffractive_scale['YY'] = float(line.split()[2])                 / 1000.
				print('XX diffractive scale: %3.1f km'%(diffractive_scale['XX']))
				print('YY diffractive scale: %3.1f km'%(diffractive_scale['YY']) + '\n')
				break
		json_output['metrics'][pipeline]['diffractive_scale'] = diffractive_scale

	## get flag statistics from flag files
	states = []
	for flagFile in flagFiles:
		f = open(flagFile, 'r')
		flagged_fraction_antenna = json.load(f)
		state = flagged_fraction_antenna['state']
		states.append(state)
		station_statistics = json_output['metrics'][pipeline]['stations']
		for antenna in flagged_fraction_antenna.keys():
			if antenna in bad_antennas_list or antenna == 'state':
				continue
			try:
				index = [i for (i, item) in enumerate(station_statistics) if item['station'] == antenna][0]
				json_output['metrics'][pipeline]['stations'][index]['percentage_flagged'][state] = flagged_fraction_antenna[antenna]
			except IndexError:
				json_output['metrics'][pipeline]['stations'].append({'station' : antenna, 'removed' : 'no', 'percentage_flagged' : {state : flagged_fraction_antenna[antenna]}})
		f.close()

	## printing results human readable
	antennas = sorted([antenna for antenna in flagged_fraction_antenna.keys() if antenna != 'state'])
	antenna_len  = '{:<' + str(max([ len(antenna)     for antenna in flagged_fraction_antenna.keys()])) + '}'
	state_len    = '{:^' + str(max([ len(state)       for state   in states                         ])) + '}'
	state_names  = ' '.join([ state_len.format(state) for state   in states                         ])
	print('Amount of flagged data per station and steps:')
	print(antenna_len.format('Station') + ' ' + state_names)
	for antenna in antennas:
		if antenna in bad_antennas_list:
			continue
		values_to_print = []
		index = next(i for (i, item) in enumerate(station_statistics) if item['station'] == antenna)
		for state in states:
			try:
				values_to_print.append(state_len.format('{:6.2f}'.format(100 * station_statistics[index]['percentage_flagged'][state]) + '%'))
			except KeyError:
				values_to_print.append(state_len.format(' '))
		values_to_print = ' '.join(values_to_print)
		print(antenna_len.format(antenna) + ' ' + values_to_print)

	## write JSON file
	with open(output_fname, 'w') as fp:
		json.dump(json_output, fp)

	print('\nSummary has been created.')
	return(0)


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description='Creates summary of a given prefactor3-CWL run.')
	parser.add_argument('flagFiles', nargs='+', help='List of flag files in JSON format')
	parser.add_argument('--pipeline', type=str, default='prefactor', help='Name of the pipeline.')
	parser.add_argument('--run_type', type=str, default='calibrator', help='Type of the pipeline')
	parser.add_argument('--filtered_antennas', type=str, default='[CR]S*&', help='Filter these antenna string from the processing.')
	parser.add_argument('--bad_antennas', type=str, default='[CR]S*&', help='Antenna string to be processed')
	parser.add_argument('--output_fname', '--output_fname', type=str, default='summary.json', help='Name of the output filename (default: summary.json)')
	parser.add_argument('--structure_file', type=str, default=None, help='Location of the structure function logfile')
	parser.add_argument('--Ateam_separation_file', type=str, default=None, help='Location of the Ateam_separation JSON file')
	
	args = parser.parse_args()
	
	# start running script
	main(args.flagFiles, args.pipeline, args.run_type, args.filtered_antennas, args.bad_antennas, args.output_fname, args.structure_file, args.Ateam_separation_file)
	
	sys.exit(0)
