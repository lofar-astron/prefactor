#!/usr/bin/env python
import os,sys
import glob
import re
import pyrap.tables as pt



def plugin_main(args,**kwargs):
    """
    Extract for a map the full path of a file as a string
    
    Parameters
    ----------
    mapfile_dir : str
        Path to the mapfile containing the path to extract
        (NB: we suppose for now that the mapfile contains only one file)
              
    Returns
    -------
    {'FilePath': } : "dict"
        Path to the file
    """
    print kwargs['mapfile_dir']
    try :
        
        fi=open(kwargs['mapfile_dir'], "r")
        line=fi.readline() #we suppose the map contains only one line
        return {'FilePath':line.split('file\'')[-1].split("'")[1]}
    except IOError,e:
		print"***** PROBLEM WITH OPENING OF THE FILE"+kwargs['mapfile_dir']+"******* ",e
	

