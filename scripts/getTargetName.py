#!/usr/bin/env python
"""
Get name of the target field
"""

import pyrap.tables
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

def main(ms_input):
    
    msfile         = input2strlist_nomapfile(ms_input)[0]
    
    observationTable = pyrap.tables.table(msfile + '::OBSERVATION')
    targetName       = observationTable.getcol('LOFAR_TARGET')['array'][0]
        
    ## return results
    result = {'targetName':targetName}
    logging.info('The target name of ' + str(msfile) + ' is ' + str(result['targetName']))
    return(result)
    

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get target name of an MS file.')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the target name.')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)
    
    main(args.MSfile)
