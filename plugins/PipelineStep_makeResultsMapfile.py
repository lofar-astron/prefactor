#!/usr/bin/env python
import os
import re
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# This functionality should go into the general "createMapfile"!

# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    result = {}
    
    if 'make_target_dir' in kwargs and string2bool(kwargs['make_target_dir']):
        if not os.path.isdir(kwargs['target_dir']):
            os.makedirs(kwargs['target_dir'])
    datamap = MiniMapfileManager.load(kwargs['mapfile_in'])
    datamap.change_directory(kwargs['target_dir'])
    if 'new_suffix' in kwargs:
        datamap.new_ext(kwargs['new_suffix'])

    if datamap:
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
        datamap.save(fileid)
        result['mapfile'] = fileid
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

class MiniMapfileManager(DataMap):
    # I didn't want to copy the full MapfileManager, but it should be possible to 
    # just copy&paste the methods in here into the "real" MapfileManager

    def change_directory(self, newDirectory):
	for i, item in enumerate(self._data):
            basename = os.path.basename(item.file)
            if not basename:
                raise ValueError('MiniMapfileManager.change_directory: basename of',item.file,'is empty!')
            self._data[i].file = os.path.join(newDirectory,basename)

    def new_ext(self, newextension):
        for i, item in enumerate(self._data):
            self._data[i].file = os.path.splitext(item.file)[0] + newextension
