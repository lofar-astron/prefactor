#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import argparse
import subprocess, platform
import numpy, re
import multiprocessing

########################################################################
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
   return str_list

###############################################################################
def find_flagged_fraction(ms_file):
    
   outputs  = os.popen('DPPP msin=' + ms_file + ' msout=. steps=[count] count.type=counter count.warnperc=0.0000001 | grep NOTE').readlines()
   fraction_flagged = { output.split('(')[-1].rstrip(')\n'):output.split('%')[0][-5:].strip() for output in outputs if 'station' in output }
   return fraction_flagged

###############################################################################
def main(ID, observation_directory = '/data/share/pipeline/Observation', h5parmdb = 'solutions.h5', MSfile = '[]'):
	"""
	Creates statistics of a given ID (needs xml_parser.py and plot_statistics.py in the same directioy)
	
	Parameters
	----------
	ID : str
	ID to create statistical runtime plots
	
	"""
	summary = observation_directory + '/' + ID + '.log'
	f_summary = open(summary, 'w')
	
	# location of logfile
	print('Summary logfile is written to ' + summary)

	try:
		log_date = os.walk(observation_directory + ID + '/logs/').next()[1][0]
		pass
	except:
		f_summary.write(observation_directory + ID + '/logs/ does not exist. Please check your ID or you observation directory!')
		f_summary.close()
		return(1)
		pass      
	
	## get logfile/statistical information
	pipeline_log = observation_directory + ID + '/logs/' + log_date + '/pipeline.log'
	
	## check for the h5parm
	if os.path.exists(observation_directory  + ID + '/results/cal_values/' + h5parmdb):
		h5parmdb = observation_directory + ID + '/results/cal_values/' + h5parmdb
		pass
	elif os.path.exists(observation_directory + ID + '/results/cal_values/cal_' + h5parmdb):
		h5parmdb =  observation_directory + ID + '/results/cal_values/cal_' + h5parmdb
	else:
		f_summary.write('No h5parm solutions file found in ' + observation_directory + ID + '/results/cal_values/' + h5parmdb)
		f_summary.close()
		return(1)
		pass
	
	## convert stringlist to list
	mslist = input2strlist_nomapfile(MSfile)
	
	## get proper solset
	from losoto.h5parm import h5parm
	data = h5parm(h5parmdb, readonly = True)
	if 'target' in data.getSolsetNames():
		solset = 'target'
		pass
	elif 'calibrator' in data.getSolsetNames():
		solset = 'calibrator'
		pass
	else:
		f_summary.write('Neither calibrator or target solset has been found in ' + h5parmdb)
		f_summary.close()
		return(1)
		pass
	
	## get antenna information
	solset       = data.getSolset(solset)
	soltabs      = list([soltab.name for soltab in solset.getSoltabs()])
	source       = solset.obj._f_get_child('source')[0][0]
	antenna_len  = '{:<' + str(max([ len(antenna)     for antenna     in solset.getAnt()])) + '}'
	soltab_len   = '{:^' + str(max([ len(soltab_name) for soltab_name in soltabs        ])) + '}'
	
	## get software versions
	try:
		DPPP_version      =             subprocess.Popen(['DPPP', '-v'], stdout=subprocess.PIPE).stdout.read()
		losoto_version    = 'losoto ' + subprocess.Popen(['losoto', '--version'], stderr=subprocess.PIPE).stderr.read()
		lsmtool_version   =             subprocess.Popen(['lsmtool', '--version'], stdout=subprocess.PIPE).stdout.read()
		aoflagger_version =             subprocess.Popen(['aoflagger', '--version'], stdout=subprocess.PIPE).stdout.read()
		wsclean_version   =             subprocess.Popen(['wsclean', '--version'], stdout=subprocess.PIPE).stdout.read().split('\n')[1] + '\n'
		python_version    =             subprocess.Popen(['python', '--version'], stderr=subprocess.PIPE).stderr.read()
		
		import matplotlib, scipy, astropy
		modules_version   = 'matplotlib ' + str(matplotlib.__version__) + ', scipy ' + str(scipy.__version__) + ', astropy ' + str(astropy.__version__) + '\n'
		os_version        = ' '.join(platform.linux_distribution()) + ' ' + platform.release() + ' ' + platform.version() + '\n'
		f_summary.write('Software versions currently used:\n' + os_version        \
                                                           + DPPP_version      \
                                                           + aoflagger_version \
                                                           + losoto_version    \
                                                           + lsmtool_version   \
                                                           + wsclean_version   \
                                                           + python_version    \
                                                           + modules_version   )
	except OSError:
		f_summary.write('Could not find all LOFAR or related software packages. Could not determine the used software versions.\n')
		pass
	
	## get the list of removed stations:
	baseline_list = list(set([ str(line.split(";")[1:]) for line in open(pipeline_log) if 'baseline:' in line and ';' in line ]))
	bad_antennas  = list(set(re.sub("[\[\]\"\ '!*]", "", str(baseline_list)).replace('\\n','').replace('\\','').split(',')))
	if len(bad_antennas) == 1 and bad_antennas[0] == '':
		bad_antennas = 'NONE'
		pass
	
	## check for A-Team warning:
	Ateam_list = list(set([ str(line.split('source ')[-1]).split(' is')[0] for line in open(pipeline_log) if 'WARNING: The A-Team source' in line ]))
	if len(Ateam_list) == 0:
		Ateam_list = 'NONE'
		pass
	
	f_summary.write('\n')
	f_summary.write('Antennas removed from the data: ' + re.sub("[\[\]']", "", str(bad_antennas)) + '\n')
	f_summary.write('A-Team sources close to the phase reference center: ' + re.sub("[\[\]']", "", str(Ateam_list)) + '\n')
	f_summary.write('\n')
	
	## diffractive scale
	structure_file = observation_directory  + ID + '/results/inspection/' + source + '_structure.txt'
	if os.path.exists(structure_file):
		with open(structure_file, 'r') as infile:
			for line in infile:
				diffractive_scale_XX = float(line.split()[2])
				diffractive_scale_YY = float(line.split()[2])
		f_summary.write('XX diffractive scale: %3.1f km'%(diffractive_scale_XX/1000.) + '\n')
		f_summary.write('YY diffractive scale: %3.1f km'%(diffractive_scale_YY/1000.) + '\n') 
	
	## check whether reference has been used
	if 'bandpass' in soltabs:
		soltab  = solset.getSoltab('bandpass')
		history = soltab.getHistory()
		if not history.strip('\n') == '':
			f_summary.write(history + '\n')
			pass
	
	## get table of flagged solutions
	f_summary.write('\n')
	flagged_solutions = {}
	import losoto.lib_operations as losoto
	for soltab_name in soltabs:
		soltab   = solset.getSoltab(soltab_name)
		antennas = soltab.ant
		axes     = soltab.getAxesNames()
		axes.insert(0, axes.pop(axes.index('ant'))) ## put antenna table to the second place
		flagged_solutions[soltab_name] = {}
		for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axes, weight=True):
			vals    = losoto.reorderAxes(vals,    soltab.getAxesNames(), axes)
			weights = losoto.reorderAxes(weights, soltab.getAxesNames(), axes)
		for i, antenna in enumerate(antennas):
			flagged_fraction = 1. - float(numpy.mean(weights[:,i]))
			flagged_solutions[soltab_name][antenna] = soltab_len.format('{:.2f}'.format(100 * flagged_fraction) + '%')
	soltab_names = ' '.join([ soltab_len.format(soltab_name) for soltab_name in soltabs ])
	f_summary.write('Amount of flagged solutions per station and solution table:\n')
	f_summary.write(antenna_len.format('Station') + ' ' + soltab_names + '\n')
	for antenna in antennas:
		values_to_print = []
		for soltab_name in soltabs:
			try:
				values_to_print.append(flagged_solutions[soltab_name][antenna])
			except KeyError:
				values_to_print.append(soltab_len.format(' '))
		values_to_print = ' '.join(values_to_print)
		f_summary.write(antenna_len.format(antenna) + ' ' + values_to_print + '\n')
		
	
	## derive the fraction of flagged data of the entire observation
	pool = multiprocessing.Pool(processes = multiprocessing.cpu_count())
	flagged_fraction_list = pool.map(find_flagged_fraction, mslist)
	
	flagged_fraction_data = {}
	for entry in flagged_fraction_list:
		antennas = entry.keys()
		for antenna in antennas:
			try:
				flagged_fraction_data[antenna].append(float(entry[antenna]))
			except KeyError:
				flagged_fraction_data[antenna] = [float(entry[antenna])]
	
	f_summary.write('\n')
	f_summary.write('Overall amount of flagged data in the final data:\n')
	f_summary.write(antenna_len.format('Station') + '\n')
	for antenna in sorted(flagged_fraction_data.keys()):
		flagged_fraction_data[antenna] = sum(flagged_fraction_data[antenna]) / len(flagged_fraction_data[antenna])
		f_summary.write(antenna_len.format(antenna) + ' ' + '{:>10}'.format('{:.2f}'.format(flagged_fraction_data[antenna]) + '%') + '\n')
	
	print('Summary has been created.')
	return 0


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description='Creates statistics of a given ID (needs xml_parser.py and plot_statistics.py in the same directioy)') 
	parser.add_argument('ID', type=str, help='ID to evaluate statistics.')
	parser.add_argument('--obsdir', type=str, default='/data/share/pipeline/Observation', help='Directory (followed by ID) where to find the log subdirectory.')
	parser.add_argument('--h5parm', '--h5parm', type=str, default='solutions.h5', help='Name of the h5parm solutions file (default: solutions.h5)')
	parser.add_argument('--MSfile', '--MSfile', type=str, default='[]', help='List of MS to be analysed (default: [])')
	
	args = parser.parse_args()
	
	# start running script
	main(args.ID, args.obsdir, args.h5parm, args.MSfile)
	
	sys.exit(0)
	pass