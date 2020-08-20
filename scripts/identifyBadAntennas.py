#!/usr/bin/env python
"""
Identify fully flagged antennas
"""

import os
import multiprocessing
import logging

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

########################################################################
def find_flagged_antennas(ms_file):
    
   logging.info('Reading ' + str(ms_file))
   outputs = os.popen('DPPP msin=' + ms_file + ' msout=. steps=[count] count.type=counter count.warnperc=100 | grep NOTE').readlines()
   flaggedants = [ output.split('(')[-1].rstrip(')\n') for output in outputs if 'station' in output ]
   return(flaggedants)
########################################################################
def main(MSfile, filter = '*&'):

    mslist         = input2strlist_nomapfile(MSfile)

    pool = multiprocessing.Pool(processes = multiprocessing.cpu_count())
    flaggedants_list = pool.map(find_flagged_antennas, mslist)

    flagged_antenna_list = set.intersection(*map(set, flaggedants_list)) 

    for flagged_antenna in flagged_antenna_list:
        filter += ';!' + flagged_antenna + '*&&*'
        pass

    print('Identified bad antennas: ' + str(flagged_antenna_list))

    ## return results
    result = {'filter':str(filter)}
    logging.info('The following stations should be used for further processing: ' + str(result['filter']))
    return(result)

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Identify fully flagged antennas.')

    parser.add_argument('MSfiles', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the antenna list.')
    parser.add_argument('--filter', type=str, default='*&',
                        help='Stations to include in the check. Use LOFAR antenna selection syntax.')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)
    
    main(args.MSfile, filter = args.filter)
