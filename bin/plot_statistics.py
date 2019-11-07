#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys
import argparse
import logging
import numpy
if not 'matplotlib' in sys.modules:
	import matplotlib as mpl
	mpl.rcParams['xtick.labelsize'] = 20
	mpl.rcParams['font.size'] = 20
	mpl.use("Agg")
import matplotlib
import matplotlib.pyplot as plt

def add_coloring_to_emit_ansi(fn):

	def new(*args):
		levelno = args[0].levelno
		if(levelno>=50):
			color = '\x1b[31m' # red
			pass
		elif(levelno>=40):
			color = '\x1b[31m' # red
			pass
		elif(levelno>=30):
			color = '\x1b[33m' # yellow
			pass
		elif(levelno>=20):
			color = '\x1b[32m' # green
			pass
		elif(levelno>=10):
			color = '\x1b[35m' # pink
			pass
		else:
			color = '\x1b[0m' # normal
			pass
		args[0].msg = color + args[0].msg +  '\x1b[0m'  # normal
		return fn(*args)
		pass
	return new
	pass


def my_handler(type, value, tb):
	exception = logger.critical("{0}".format(str(value)))
	pass


########################################################################
def main(statistics_log, ID = None, logdate = None):
	"""
	Parse a statistics log file and plot their contents
	
	Parameters
	----------
	statistics_log : str
	UNIX-compatible location of the statistics log file
	
	ID : str
	provide ID for more verbose header
	
	lodate : str
	provide logdate for more verbose header
	"""
	
	## extract for overview plot
	try:
		steps        = [ line.split("STEP:")[-1].split(", duration:")[0].replace(' \x1b[35m','').replace('\x1b[32m','') for line in open(statistics_log) if 'STEP:' in line and 'PLUGIN' not in line ]
		duration     = [ float(line.split("STEP:")[-1].split(", duration:")[-1].replace(' \x1b[35m','').replace('s\x1b[0m\x1b[0m\n','')) for line in open(statistics_log) if 'STEP:' in line and 'PLUGIN' not in line ]
		runtime      = [ line.split("Full runtime is:")[-1].replace(' \x1b[35m','').replace('\x1b[0m\x1b[0m\n','') for line in open(statistics_log) if 'Full runtime is:' in line ][0]
		cluster_type = [ line.split("Type of cluster used for processing: ")[-1].replace('\x1b[35m','').replace('\x1b[0m\x1b[0m\n','') for line in open(statistics_log) if 'Type of cluster' in line ][0]
		xaxis        = numpy.arange(0,len(steps))
		pass
	except IndexError:
		logging.error('The provided logfile ' + statistics_log + ' does not provide all necessary information. Plotting will be skipped!')
		return 1
		pass
	
	## extract for node plots
	total_number     = [ int(line.split('step: ')[-1].replace('\x1b[35m','').replace('\x1b[0m\x1b[0m\n','')) for line in open(statistics_log) if 'Total number' in line ]
	total_nodes_used = [ line for line in open(statistics_log) if 'job(s) were run' in line ]
	multi_nodes      = [ number_nodes for number_nodes in total_number if number_nodes != 1 ]
	multi_steps      = [ step for step,number_nodes in zip(steps,total_number) if number_nodes != 1 ]
	multi_duration   = [ durations for durations,number_nodes in zip(duration,total_number) if number_nodes != 1 ]
	multi_nodes_used = [ nodes for nodes in total_nodes_used if int(nodes.split(' job(s) were run on: ')[0].replace('\x1b[1mINFO:\x1b[0m \x1b[32m\x1b[32m','')) != 1 ]
	
	jobs  = []
	nodes = []
	times = []
	index = 0
	
	for nodes_used in multi_nodes:
		jobs_list = [ int(multi_nodes_used[i].split(' job(s) were run on: ')[0].replace('\x1b[1mINFO:\x1b[0m \x1b[32m\x1b[32m','')) for i in numpy.arange(index,nodes_used + index) ] 
		node_list = [ multi_nodes_used[i].split(' job(s) were run on: ')[-1].split(' with an')[0].replace('\x1b[35m','').replace('\x1b[32m','') for i in numpy.arange(index,nodes_used + index) ] 
		time_list = [ float(multi_nodes_used[i].split('time of: ')[-1].split('s per job')[0].replace('\x1b[35m ','')) for i in numpy.arange(index,nodes_used + index) ] 
		jobs.append(jobs_list)
		nodes.append(node_list)
		times.append(time_list)
		index += nodes_used
		pass

	with open(statistics_log) as myfile:
		header = [next(myfile).replace('\x1b[1mINFO:\x1b[0m \x1b[32m\x1b[32m','').replace('\x1b[35m','').replace('\x1b[0m\x1b[0m','') for x in xrange(3)]
	
	legend = ''
	for i,step in enumerate(steps):
		legend += str(i) + ': ' + str(step) + '\n'
		pass
	
	runtime_int       = int(runtime.split('h')[0]) * 3600 + int(runtime.split('h')[-1].split('m')[0]) * 60 + int(runtime.split('m')[-1].rstrip('s'))
	relative_duration = numpy.array(duration) / float(runtime_int)
	
	title = ''
	if ID != None:
		title += 'ID: ' + ID + ' '
		pass
        if logdate != None:
		title += '(' + logdate + ')'
		pass
	if ID != None or logdate != None:
		title += '\n'
		pass
	fig = plt.figure(figsize=(20,18))
	ax  = plt.axes([0.1,0.1,0.68,0.85])
	
	plt.bar(xaxis, relative_duration, align='center', log=False)
	plt.style.use('classic')
	plt.title(title + 'Full runtime: ' + str(runtime) + ', cluster type: ' + cluster_type)
	plt.xlabel('Step ID'            , fontsize=14)
	plt.ylabel('fraction of runtime', fontsize=14)
	ax.set_xticks(xaxis)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	
	plt.gcf().text(0.80, 0.55 - len(steps)/1000., legend.rstrip('\n'), verticalalignment='center')#, bbox=dict(fc="none"))
	plt.gcf().text(0.05, 0.05, header[1] + header[2], verticalalignment='center')#, bbox=dict(fc="none"))
	
	output = os.path.dirname(statistics_log) + '/' + os.path.basename(statistics_log).split('.')[0] + '.png'
	logging.info('Plotting statistics to ' + output)
	plt.savefig(output)
	plt.close(fig)
	
	## plot the node distribution of the jobs
	for node, durations, step, time, job in zip(nodes, multi_duration, multi_steps, times, jobs):
		fig       = plt.figure(figsize=(20,18))
		ax        = plt.axes([0.1,0.1,0.68,0.85])
		xaxis     = numpy.arange(0,len(job))
		cmap      = matplotlib.cm.get_cmap('rainbow')
		normalize = matplotlib.colors.Normalize(vmin=numpy.mean(job)*0.5, vmax=numpy.mean(job)*1.5)
		colors    = [cmap(normalize(value)) for value in job]
		legend    = ''
		for i,ID in enumerate(node):
			legend += str(i) + ': ' + str(ID) + '\n'
			pass
		plt.bar(xaxis, time, align='center', log=False, color=colors)
		plt.style.use('classic')
		plt.title(title + 'Step: ' + step + ', duration: ' + str(durations) + 's')
		plt.xlabel('Node ID' , fontsize=14)
		plt.ylabel('average processing time per job [s]', fontsize=14)
		ax.set_xticks(xaxis)
		ax.tick_params(axis='x', labelsize=14)
		ax.tick_params(axis='y', labelsize=14)
		plt.gcf().text(0.80, 0.55 - len(nodes)/1000., legend.rstrip('\n'), verticalalignment='center')
		cax, _ = matplotlib.colorbar.make_axes(ax, fraction=0.15, shrink=1.0, aspect=30)
		cbar   = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize, label='#jobs')
		output = statistics_log.split('.')[:-1][0] + '_' + step + '.png'
		logging.info('Plotting statistics to ' + output)
		plt.savefig(output)
		plt.close(fig)
		pass


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description=' Plots the contents of a statistics file.') 
	parser.add_argument('statistics_file', type=str, help='A statistics file to plot.')
	parser.add_argument('--ID', type=str, default=None, help='Provide ID (optional).')
	parser.add_argument('--logdate', type=str, default=None, help='Provide log date (optional).')
	
	args = parser.parse_args()
	
	format_stream = logging.Formatter("\033[1m%(levelname)s:\033[0m %(message)s")
	format_file   = format_stream
	logging.root.setLevel(logging.INFO)
	
	log      = logging.StreamHandler()
	log.setFormatter(format_stream)
	log.emit = add_coloring_to_emit_ansi(log.emit)
	logging.root.addHandler(log)
	
	LOG_FILENAME = os.path.dirname(args.statistics_file) + '/' + os.path.basename(args.statistics_file).split('.')[0] + '.log'
	logfile = logging.FileHandler(LOG_FILENAME)
	logfile.setFormatter(format_file)
	logfile.emit = add_coloring_to_emit_ansi(logfile.emit)
	logging.root.addHandler(logfile)
	
	# install exception handler
	logger = logging.getLogger(LOG_FILENAME)
	sys.excepthook = my_handler
	
	# location of logfile
	logging.info('\033[0mLog file is written to ' + LOG_FILENAME)
	
	# start running script
	main(args.statistics_file, args.ID, args.logdate)
	
	sys.exit(0)
	pass
