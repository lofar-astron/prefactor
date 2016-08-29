import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by expanding single input mapfile item into many items

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing single item
    mapfile_to_match : str
        Filename of datamap containing multiple items
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_to_match = kwargs['mapfile_to_match']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_match = DataMap.load(mapfile_to_match)
    map_out = DataMap([])

    map_match.iterator = DataMap.SkipIterator
    for item in map_match:
        map_out.data.append(DataProduct(item.host, map_in[0].file, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
