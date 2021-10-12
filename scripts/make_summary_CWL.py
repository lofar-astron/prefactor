#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import argparse
import json
import numpy

###############################################################################
def main(flagFiles = None, pipeline = 'prefactor', run_type = 'calibrator', filtered_antennas = '[CR]S*&', bad_antennas = '[CR]S*&', output_fname = None, structure_file = None, Ateam_separation_file = None, solutions = None, clip_sources = '', demix_sources = '', demix = False, removed_bands = '', min_unflagged = 0.5, refant = ''):
	"""
	Creates summary of a given prefactor3-CWL run
	
	Parameters
	----------
	pipeline: name of the pipeline
	filter: antenna string which is pre-selected (by user or pipeline)
	bad_anntenas: pre-selected antennas and removed stations (separated by ";")
	
	"""
	## header
	header_string = '*** ' + pipeline + ' ' + run_type + ' pipeline summary ***'
	print('*' * len(header_string) + '\n' + header_string + '\n' + '*' * len(header_string) + '\n')
	
	## define contents of JSON file
	json_output = { 'metrics': { pipeline : { } } }
	json_output['metrics'][pipeline]['run_type'] = run_type
	json_output['metrics'][pipeline]['stations'] = []
	
	## get bad antennas
	bad_antennas_list = list(filter(None, set(bad_antennas.replace(filtered_antennas,'').replace('!','').replace('*','').replace('&','').split(';'))))

	## read solset
	source = ''
	values_to_print = ''
	soltab_names = []
	if solutions:
		flagged_solutions = {}
		from losoto.h5parm import h5parm
		import losoto.lib_operations as losoto
		data     = h5parm(solutions, readonly = True)
		solset   = data.getSolset(run_type)
		antennas = sorted(solset.getAnt().keys())
		for antenna in antennas:
			json_output['metrics'][pipeline]['stations'].append({'station' : antenna, 'removed' : 'no', 'percentage_flagged' : {}})
		source   = solset.obj._f_get_child('source')[0][0].decode('utf-8')
		print('Field name: ' + source + '\n')
		json_output['metrics'][pipeline]['field_name'] = source
		for solset_name in data.getSolsetNames():
			solset   = data.getSolset(solset_name)
			soltabs  = list([soltab.name for soltab in solset.getSoltabs()])
			for soltab_name in soltabs:
				soltab_names.append(soltab_name)
				soltab  = solset.getSoltab(soltab_name)
				history = soltab.getHistory()
				ants    = soltab.ant
				axes    = soltab.getAxesNames()
				axes.insert(0, axes.pop(axes.index('ant')))
				flagged_solutions[soltab_name] = {}
				if not history.strip('\n') == '':
					values_to_print += history + '\n'
				for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axes, weight=True):
					weights = losoto.reorderAxes(weights, soltab.getAxesNames(), axes)
				for i, ant in enumerate(ants):
					if ant not in antennas:
						continue
					flagged_solutions[soltab_name][ant] = 1. - float(numpy.mean(weights[i]))
					if 'percentage_flagged_solutions' not in json_output['metrics'][pipeline]['stations'][antennas.index(ant)].keys():
						json_output['metrics'][pipeline]['stations'][antennas.index(ant)]['percentage_flagged_solutions'] = {}
					json_output['metrics'][pipeline]['stations'][antennas.index(ant)]['percentage_flagged_solutions'][soltab_name] = flagged_solutions[soltab_name][ant] * 100

	## print antennas removed from the data
	print('User-specified baseline filter: ' + filtered_antennas)
	json_output['metrics'][pipeline]['filtered_baselines_by_user'] = filtered_antennas
	if bad_antennas_list == []:
		print('Additional antennas removed from the data: NONE')
	else:
		print('Additional antennas removed from the data: ' + ', '.join(bad_antennas_list))
	
	if refant != '':
		print('Selected reference antenna: ' + refant)
	
	## get Ateam_separation info
	if Ateam_separation_file:
		f = open(Ateam_separation_file, 'r')
		json_output['metrics'][pipeline]['close_sources'] = json.load(f)
		f.close()
		if len(json_output['metrics'][pipeline]['close_sources']) > 0:
			Ateam_list = ''
			clip_list  = []
			demix_list = []
			for i in range(len(json_output['metrics'][pipeline]['close_sources'])):
				Ateam_name = json_output['metrics'][pipeline]['close_sources'][i]['source']
				if demix:
					if Ateam_name in demix_sources:
						json_output['metrics'][pipeline]['close_sources'][i]['mitigation'] = 'demix'
						demix_list.append(Ateam_name)
				elif Ateam_name in clip_sources:
					json_output['metrics'][pipeline]['close_sources'][i]['mitigation'] = 'clip'
					clip_list.append(Ateam_name)
				else:
					json_output['metrics'][pipeline]['close_sources'][i]['mitigation'] = 'none'
				Ateam_list += json_output['metrics'][pipeline]['close_sources'][i]['source'] + ','
			if clip_list != []:
				clip_list = ', '.join(clip_list)
				print(clip_list)
			else:
				clip_list = 'NONE'
			if demix_list != []:
				demix_list = ', '.join(demix_list)
			else:
				demix_list = 'NONE'
			Ateam_list = Ateam_list.rstrip(',') + '\n\tOf which were demixed: ' + demix_list + '\n\tOf which were clipped: ' + clip_list
		else:
			Ateam_list = 'NONE'
		print('A-Team sources close to the phase reference center: ' + Ateam_list + '\n')
	
	if removed_bands != '':
		print('Removed bands due to an unsufficient fraction of unflagged data (< ' + str(min_unflagged * 100) + '%): ' + ', '.join(list(filter(lambda a: a != 'None', removed_bands.replace('out_','').split(',')))) + '\n')
	
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
		antennas = sorted(flagged_fraction_antenna.keys())
		antennas.remove('state')
		station_statistics = json_output['metrics'][pipeline]['stations']
		if station_statistics == []:
			for antenna in antennas:
				json_output['metrics'][pipeline]['stations'].append({'station' : antenna, 'removed' : 'no'})
		for antenna in antennas:
			index = [i for (i, item) in enumerate(station_statistics) if item['station'] == antenna][0]
			if 'percentage_flagged' not in json_output['metrics'][pipeline]['stations'][index].keys():
				json_output['metrics'][pipeline]['stations'][index]['percentage_flagged'] = {}
			json_output['metrics'][pipeline]['stations'][index]['percentage_flagged'][state] = flagged_fraction_antenna[antenna] * 100
		f.close()

	## printing results human readable
	antennas = sorted([item['station'] for (i, item) in enumerate(station_statistics)])
	antenna_len      = '{:<' + str(max([ len(antenna) for antenna in antennas])) + '}'
	if solutions:
		## print modifications from solutions file
		if values_to_print != '':
			print('Changes applied to ' + os.path.basename(solutions) + ':')
			print(values_to_print)
		soltab_len   = '{:^' + str(max([ len(soltab_name) for soltab_name in soltab_names        ])) + '}'
		soltab_names = ' '.join([ soltab_len.format(soltab_name) for soltab_name in soltab_names ])
		print('Amount of flagged solutions per station and solution table:')
		print(antenna_len.format('Station') + ' ' + soltab_names)
		for antenna in antennas:
			values_to_print = []
			for solset_name in data.getSolsetNames():
				solset   = data.getSolset(solset_name)
				soltabs  = list([soltab.name for soltab in solset.getSoltabs()])
				for soltab_name in soltabs:
					try:
						values_to_print.append(soltab_len.format('{:6.2f}'.format(100 * flagged_solutions[soltab_name][antenna]) + '%'))
					except KeyError:
						values_to_print.append(soltab_len.format(' '))
			values_to_print = ' '.join(values_to_print)
			print(antenna_len.format(antenna) + ' ' + values_to_print)
		print('')
		data.close()

	state_len    = '{:^' + str(max([ len(state)       for state   in states                         ])) + '}'
	state_names  = ' '.join([ state_len.format(state) for state   in states                         ])
	print('Amount of flagged data per station at a given state:')
	print(antenna_len.format('Station') + '  ' + state_names)
	for antenna in antennas:
		values_to_print = []
		index = next(i for (i, item) in enumerate(station_statistics) if item['station'] == antenna)
		for state in states:
			try:
				values_to_print.append(state_len.format('{:6.2f}'.format(station_statistics[index]['percentage_flagged'][state]) + '%'))
			except KeyError:
				values_to_print.append(state_len.format(' '))
		values_to_print = ' '.join(values_to_print)
		print(antenna_len.format(antenna) + ' ' + values_to_print)

	## write bad antenna information into JSON file
	for bad_antenna in bad_antennas_list:
		station_statistics = json_output['metrics'][pipeline]['stations']
		try:
			index = [i for (i, item) in enumerate(station_statistics) if item['station'] == bad_antenna][0]
			json_output['metrics'][pipeline]['stations'][index] = {'station' : bad_antenna, 'removed' : 'yes'}
		except IndexError:
			json_output['metrics'][pipeline]['stations'].append(  {'station' : bad_antenna, 'removed' : 'yes'})

	## write JSON file
	if not output_fname:
		output_fname = (source + '_' + pipeline + '_'  + run_type + '_summary.json').lstrip('_')
	with open(output_fname, 'w') as fp:
		json.dump(json_output, fp)
	
	print('\n' + '*' * 10)
	print('Summary file is written to: ' + output_fname)
	print('Summary has been created.')
	
	return(0)


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description='Creates summary of a given prefactor3-CWL run.')
	parser.add_argument('flagFiles', nargs='+', help='List of flag files in JSON format')
	parser.add_argument('--pipeline', type=str, default='prefactor', help='Name of the pipeline.')
	parser.add_argument('--run_type', type=str, default='calibrator', help='Type of the pipeline')
	parser.add_argument('--filtered_antennas', type=str, default='[CR]S*&', help='Filter these antenna string from the processing.')
	parser.add_argument('--bad_antennas', type=str, default='[CR]S*&', help='Antenna string to be processed')
	parser.add_argument('--output_fname', type=str, default=None, help='Name of the output filename (default: summary.json)')
	parser.add_argument('--structure_file', type=str, default=None, help='Location of the structure function logfile')
	parser.add_argument('--Ateam_separation_file', type=str, default=None, help='Location of the Ateam_separation JSON file')
	parser.add_argument('--solutions', type=str, default=None, help='Location of the solutions h5parm file')
	parser.add_argument('--clip_sources', type=str, default='', help='Comma-separated list of sources that were clipped')
	parser.add_argument('--demix_sources', type=str, default='', help='Comma-separated list of sources that were demixed')
	parser.add_argument('--demix', type=bool, default=None, help='Tell the summary that demixing was enabled')
	parser.add_argument('--removed_bands', type=str, default='', help='Comma-separated list of bands that were removed from the data')
	parser.add_argument('--min_unflagged', type=float, default=0.5, help='Minimum fraction of unflagged data per band')
	parser.add_argument('--refant', type=str, default='', help='Reference antenna used')
	
	args = parser.parse_args()
	
	# start running script
	main(flagFiles = args.flagFiles, pipeline = args.pipeline, run_type = args.run_type, filtered_antennas = args.filtered_antennas, bad_antennas = args.bad_antennas, output_fname = args.output_fname, structure_file = args.structure_file, Ateam_separation_file = args.Ateam_separation_file, solutions = args.solutions, clip_sources = args.clip_sources, demix_sources = args.demix_sources, demix = args.demix, removed_bands = args.removed_bands, min_unflagged = args.min_unflagged, refant = args.refant)
	
	sys.exit(0)
