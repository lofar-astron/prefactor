#!/usr/bin/env python
import os
import re
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    """
    Gets several mapfiles of length 1 as input, extracts names (pathes) from these
    mapfiles and returns those as mapfiles.
    Allows to access mapfiles written by a pythonplugin.

    Parameters
    ----------
    mapfile_* : str
      Takes one or more arguments that start with "mapfile_". The part after 
      the "_" will become the name of the mapfile in the output directory.
      (E.g. "mapfile_data_single" will become "data_single")

    Returns
    -------
    result : dict
        Output dictionary with the mapfile-names

    """
    result = {}
    mapfile_keys = []
    
    for keyname in kwargs.keys():
        if keyname[0:8] == "mapfile_":
            mapfile_keys.append(keyname)
            if len(keyname) < 10:
                raise ValueError("MapfilenamesFromMapfiles: Key: "+keyname+" is too short!")
        else:
            print "MapfilenamesFromMapfiles: input key:",keyname,"in unkown!"
    for keyname in mapfile_keys:
        inmap = DataMap.load(kwargs[keyname])
        if len(inmap) != 1:
            raise ValueError("MapfilenamesFromMapfiles: Length of mapfile "+kwargs[keyname]+" != 1. (From key: "+keyname+")")
        outname = keyname[8:]
        result[outname] = inmap[0].file
    return result

def string2bool(instring):
    if not isinstance(instring, basestring):
        raise ValueError('string2bool: Input is not a basic string!')
    if instring.upper() == 'TRUE' or instring == '1':
        return True
    elif instring.upper() == 'FALSE' or instring == '0': 
        return False
    else:
        raise ValueError('string2bool: Cannot convert string "'+instring+'" to boolean!')


