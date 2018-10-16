#!/usr/bin/env python
import sys
import glob
import re
import pyrap.tables as pt
import numpy

#max_length = 147

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
def main(ms_input,ms_output,max_length):

    """
    Virtuall concatenate subbands 
  

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the calibrator MSs
    ms_output : str
        String from the outut concatenated MS

    """    
    filelist      = input2strlist_nomapfile(ms_input)
    set_ranges    = list(numpy.arange(0, len(filelist), int(max_length)))
    set_ranges.append(len(filelist))
    
    for i in numpy.arange(len(set_ranges) - 1):
        pt.msconcat(filelist[set_ranges[i]:set_ranges[i + 1]], ms_output + '_' + str(i))
        pass
                 
    return { 'pattern' : ms_output.split('/')[-1] + '*'}

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Virtually concat subbands')
    
    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) that we want to concatenate.')
    parser.add_argument('MSout', type=str, 
                        help='Output MS file')
    parser.add_argument('max_length', type=str, 
                        help='Max length of concatenation')
    
        
 
   
    args = parser.parse_args()
    
    main(args.MSfile,args.MSout,args.max_length)
