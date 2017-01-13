import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by repeating max size in input mapfile items

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_out = DataMap([])

    # Find max size in x and y
    xsize_list = []
    ysize_list = []
    for item in map_in:
        xsize, ysize = [int(s) for s in item.file.split(' ')]
        xsize_list.append(xsize)
        ysize_list.append(ysize)
    maxsize = '{0} {1}'.format(max(xsize_list), max(ysize_list))

    for item in map_in:
        map_out.data.append(DataProduct(item.host, maxsize, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
