#!/usr/bin/env python3
"""
Identify fully flagged antennas
"""

import os
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
   return(','.join(flaggedants))
########################################################################
def main(MSfile):

    ms_file         = input2strlist_nomapfile(MSfile)[0]

    flaggedants     = find_flagged_antennas(ms_file)

    ## return results
    result = {'flaggedants':str(flaggedants)}
    logging.info(str(flaggedants) + ' are flagged in ' + str(ms_file))
    return(result)

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Identify fully flagged antennas.')

    parser.add_argument('MSfiles', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the antenna list.')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)
    
    main(args.MSfile)
