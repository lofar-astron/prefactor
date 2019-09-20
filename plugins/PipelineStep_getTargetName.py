#!/usr/bin/env python
"""
Identify fully flagged antennas
"""

from lofarpipe.support.data_map import DataMap, DataProduct
import pyrap.tables

def find_flagged_antennas(ms_file):
    
   print 'Reading ' + str(ms_file)
   outputs = os.popen('DPPP msin=' + ms_file + ' msout=. steps=[count] count.type=counter count.warnperc=100 | grep NOTE').readlines()
   flaggedants = [ output.split('(')[-1].rstrip(')\n') for output in outputs ]
   return flaggedants

def plugin_main(args, **kwargs):
    """
    Takes in list of targets and an h5parm solution set and returns a list of stations
    in the target data which mismatch the calibrator solutions antenna table
    
    Parameters
    ----------
    mapfile_in : str
        Mapfile for input measurement sets
    filter: str
        Default filter constrains for the ndppp_prep_target step (usually removing International Baselines)
    
    Returns
    -------
    result : dict
        Output station names to filter
    """
    mapfile_in     = kwargs['mapfile_in']
    data           = DataMap.load(mapfile_in)
    mslist         = [data[i].file for i in xrange(len(data))]
    msfile         = mslist[0]
    
    observationTable = pyrap.tables.table(msfile + '::OBSERVATION')
    targetName       = observationTable.getcol('LOFAR_TARGET')['array'][0]
        
    ## return results
    result = {'targetName':targetName}
    return result
    
    pass    
