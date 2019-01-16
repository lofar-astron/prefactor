#!/usr/bin/env python
import sys
import glob
import re
import pyrap.tables as pt
import numpy
import os
from lofarpipe.support.data_map import DataMap, DataProduct

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
def getsystemmemory():
   
   memory = int(os.popen('cat /proc/meminfo | grep MemAvailable').readlines()[0].split(':')[-1].split()[0])
   return memory
   pass

########################################################################
def getfilesize(MS):
   
   size = int(os.popen('du -cks ' + MS).readlines()[0].split()[0])
   return size

   pass

########################################################################
def main(ms_input, ms_output, min_length, filename=None, mapfile_dir=None):

    """
    Virtually concatenate subbands

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the calibrator MSs
    ms_output : str
        String from the outut concatenated MS
    filename: str
        Name of output mapfile
    mapfile_dir : str
        Directory for output mapfile

    """
    filelist      = input2strlist_nomapfile(ms_input)
    max_space     = int(getsystemmemory() / getfilesize(filelist[0]))
    max_length    = len(filelist) / ((len(filelist) / max_space) + 1)
    
    if max_length >= int(min_length):
        memory = '-memory-read'
        pass
    else:
        max_length = len(filelist)
        memory = '-indirect-read'
        pass
    
    print "The max_length value is: " + str(max_length)
    set_ranges    = list(numpy.arange(0, len(filelist), int(max_length)))
    set_ranges.append(len(filelist))

    map_out = DataMap([])
    for i in numpy.arange(len(set_ranges) - 1):
        f = ms_output + '_' + str(i)
        pt.msconcat(filelist[set_ranges[i]:set_ranges[i + 1]], f)
        map_out.data.append(DataProduct('localhost', f, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'concatmapfile': fileid, 'memory': memory}

    return result

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Virtually concat subbands')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) that we want to concatenate.')
    parser.add_argument('MSout', type=str,
                        help='Output MS file')
    parser.add_argument('-min_length', type=str,
                        help='Minimum amount of subbands to concatenate in frequency.')



    args = parser.parse_args()

    main(args.MSfile,args.MSout,args.min_length)
