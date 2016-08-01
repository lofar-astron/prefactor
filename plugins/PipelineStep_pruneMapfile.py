import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Prunes entries from a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to trim
    prune_str : str
        Entries starting with this string will be removed.

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    prune_str = kwargs['prune_str'].lower()
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    prunelen = len(prune_str)

    map_out = DataMap([])
    map_in = DataMap.load(mapfile_in)

    for i, item in enumerate(map_in):
        if item.file[:prunelen].lower() != prune_str:
            map_out.data.append(DataProduct(item.host, item.file, item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
