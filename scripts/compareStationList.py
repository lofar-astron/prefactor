#!/usr/bin/env python
"""
Compare station list between MSfile and h5parm
"""

import pyrap.tables as pt
from losoto.h5parm import h5parm
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
def main(MSfile, h5parmdb, solset_name = 'calibrator', filter = '*&'):

    msfile = input2strlist_nomapfile(MSfile)[0]

    ## reading ANTENNA table of MS
    logging.info('Collecting information from the ANTENNA table of ' + str(msfile))
    antennaTable = pt.table(msfile + "::ANTENNA", ack = False)
    antennaNames = antennaTable.getcol('NAME')
    
    ## reading ANTENNA information of h5parm
    data   = h5parm(h5parmdb, readonly = True)
    solset = data.getSolset(solset_name)
    station_names = solset.getAnt().keys()
    
    #station_names = [station_name.decode('utf-8') for station_name in station_names]
    
    ## check whether there are more stations in the target than in the calibrator solutions
    missing_stations = list(set(antennaNames) - set(station_names))
    for missing_station in missing_stations:
        filter += ';!' + missing_station + '*'
        pass

    data.close()
    
    ## return results
    result = {'filter':str(filter)}
    logging.info('The following stations should be used for further processing: ' + str(result['filter']))
    return(result)
    

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Check whether station entries in MSfiles and h5parms match.')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the antenna list.')
    parser.add_argument('h5parmdb', type=str,
                        help='h5parmd for which you want to get the antenna list..')
    parser.add_argument('--solset_name', type=str, default='calibrator',
                        help='Name of the solset in the h5parmdb to check.')
    parser.add_argument('--filter', type=str, default='*&',
                        help='Stations to include in the check. Use LOFAR antenna selection syntax.')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)
    
    main(args.MSfile, args.h5parmdb, solset_name = args.solset_name, filter = args.filter)
