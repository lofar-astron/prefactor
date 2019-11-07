#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import argparse
import logging

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


###############################################################################
def main(ID = 'Pre-Facet-Target', observation_directory = '/data/scratch/drabent/working_directory/'):
	"""
	Creates statistics of a given ID (needs xml_parser.py and plot_statistics.py in the same directioy)
	
	Parameters
	----------
	ID : str
	ID to create statistical runtime plots
	
	"""
	cwd = os.getcwd()
	logging.info('The current working directory is: ' + str(cwd))
	if not os.path.exists(cwd + '/xml_parser.py'):
		logging.error('xml_parser.py not found in current directory!')
		return 1
		pass
	elif not os.path.exists(cwd + '/plot_statistics.py'):
		logging.error('plot_statistics.py not found in current directory!')
		return 1
		pass
        
	try:
		log_date          = os.walk(observation_directory + ID + '/logs/').next()[1][0]
		pass
	except:
		logging.error(observation_directory + ID + '/logs/ does not exist. Please check your ID or you observation directory!')
		return 1
		pass

	logfile_directory = observation_directory + ID + '/logs/' + log_date
	
	if os.path.exists(logfile_directory + '/statistics.log'):
		logging.warning(logfile_directory + '/statistics.log will be overwritten!')
		os.remove(logfile_directory + '/statistics.log')
		pass
	os.system(cwd + '/xml_parser.py ' + logfile_directory + '/statistics.xml ' + logfile_directory + '/pipeline.log')
	os.system(cwd + '/plot_statistics.py ' + logfile_directory + '/statistics.log --ID ' + ID + ' --logdate ' + log_date)
        
        logging.info('Statistics have been created.')
        return 0


if __name__=='__main__':
    
	parser = argparse.ArgumentParser(description='Creates statistics of a given ID (needs xml_parser.py and plot_statistics.py in the same directioy)') 
	parser.add_argument('ID', type=str, default='Pre-Facet-Target', help='ID to evaluate statistics.')
	parser.add_argument('--obsdir', type=str, default='/data/scratch/drabent/working_directory/', help='Directory (followed by ID) where to find the log subdirectory.')
	
	args = parser.parse_args()
	
	format_stream = logging.Formatter("\033[1m%(levelname)s:\033[0m %(message)s")
	format_file   = format_stream
	logging.root.setLevel(logging.INFO)
	
	log      = logging.StreamHandler()
	log.setFormatter(format_stream)
	log.emit = add_coloring_to_emit_ansi(log.emit)
	logging.root.addHandler(log)
	
	LOG_FILENAME = args.ID + '.log'
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
	main(args.ID, args.obsdir)
	
	sys.exit(0)
	pass