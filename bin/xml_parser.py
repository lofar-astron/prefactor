#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import argparse
import logging
import xml.etree.ElementTree
import numpy
import re

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
def main(xml_input, pipeline_log):
	"""
	Parse an XML-file and create a logfile out of it
	
	Parameters
	----------
	xml_input : str
	UNIX-compatible location of the XML file to parse
	
        pipeline_log : str
	UNIX-compatible location of the genericpipeline logfile
	"""
	
	## get the tree and the root of it
	tree = xml.etree.ElementTree.parse(xml_input)
	root = tree.getroot()
	
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
            
	## count savestate steps:
	saved_state_steps = sum([ 1 for line in open(pipeline_log) if 'already exists in saved state; skipping' in line ])
	
	logging.info('Antennas removed from the data: \033[35m' + re.sub("[\[\]']", "", str(bad_antennas)))
	logging.info('A-Team sources close to the phase reference center: \033[35m' + re.sub("[\[\]']", "", str(Ateam_list)))
	logging.info('Overview on the time usage')
	full_time = 0.
	index = 0
	for child in root:   ## loop over all jobs
		try:
			full_time += float(child.attrib['duration'])
			if len(child) > 0:
				log_string = 'STEP: \033[35m' + child.tag + '\033[32m, duration: \033[35m %.2fs' % float(child.attrib['duration'])
				logging.info(log_string)
				if saved_state_steps == 0:
					index += 1
					pass
				else:
					logging.error('Logfile does not contain any node and time information about this step. Statistics may be wrong.')
					saved_state_steps -= 1
					raise KeyError
					pass
			else:
				log_string = 'PLUGIN STEP: \033[35m' + child.tag + '\033[32m, duration: \033[35m %.2fs' % float(child.attrib['duration'])
				logging.info(log_string)
				#continue
				pass
			for jobs in child:
				node_duration = {}
				number_jobs   = {}
				hosts         = []
				for job in jobs:
					job_id    = int(job.attrib['job_id'])
					errorcode = int(job.attrib['returncode'])
					if errorcode != 0:
						logging.error('Job \033[35m' + str(job_id) + '\033[31m returned errorcode: \033[35m' + str(errorcode))
						continue
						pass
					job_host = job.attrib['job_host']
					if job_host not in node_duration.keys():
						node_duration[job_host] = float(job.attrib['duration'])
						number_jobs[job_host]   = 1
						hosts.append(job_host)
						pass
					else:
						node_duration[job_host] += float(job.attrib['duration'])
						number_jobs[job_host]   += 1
						pass
				logging.info('Total number of nodes used for this step: \033[35m' + str(len(hosts)))
				for host in hosts:
					time_per_job = node_duration[host] / number_jobs[host]
					logging.info(str(number_jobs[host]) + ' job(s) were run on: \033[35m' + host + '\033[32m with an average processing time of: \033[35m %.2fs per job' % time_per_job)
					pass
				pass
			pass
		except KeyError:
			pass
		logging.info('---------------------------------------------')
		pass
	
	if 'cpu' in hosts[0]:
		cluster_type = 'CPU'
		pass
	elif 'gpu' in hosts[0]:
		cluster_type = 'GPU'
		pass
	else:
		cluster_type = 'UNKNOWN'
		pass
	logging.info('Type of cluster used for processing: \033[35m' + cluster_type)
	m, s = divmod(full_time, 60)
	h, m = divmod(m, 60)
	logging.info('Full runtime is: \033[35m%dh%02dm%02ds' % (h, m, int(s)))
	return 0
            

if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description=' Extracts statistics and important information out of the prefactor logfiles') 
	parser.add_argument('XMLfile', type=str, help='One XML-file to create a logfile out of it.')
	parser.add_argument('pipeline_log', type=str, help='Logfile from the genericpipeline run (pipeline.log).')
	
	args = parser.parse_args()
	
	format_stream = logging.Formatter("\033[1m%(levelname)s:\033[0m %(message)s")
	format_file   = format_stream
	logging.root.setLevel(logging.INFO)
	
	log      = logging.StreamHandler()
	log.setFormatter(format_stream)
	log.emit = add_coloring_to_emit_ansi(log.emit)
	logging.root.addHandler(log)
	
	LOG_FILENAME = os.path.dirname(args.XMLfile) + '/' + os.path.basename(args.XMLfile).split('.')[0] + '.log'
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
	main(args.XMLfile, args.pipeline_log)
	
	sys.exit(0)
	pass