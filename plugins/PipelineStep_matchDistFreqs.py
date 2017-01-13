import os
import casacore.tables as pt
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Matches a mapfile with one in which the MSs are distributed

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    mapfile_dist : str
        Filename of mapfile with distributed MS files
    mapfile_full : str
        Filename of mapfile with all MS files from which distributed one was
        made
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    num: int, optional
        Number of frequencies in output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_dist = kwargs['mapfile_dist']
    mapfile_full = kwargs['mapfile_full']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    fileid = os.path.join(mapfile_dir, filename)

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_dist = DataMap.load(mapfile_dist)
    map_dist.iterator = DataMap.SkipIterator
    map_full = DataMap.load(mapfile_full)
    map_full.iterator = DataMap.SkipIterator
    map_out = DataMap()
    map_out.data = []
    map_out._data = []

    # do not re-run if we already ran, and input files are deleted.
    if os.path.exists(fileid) and not os.path.exists(map_in[0].file):
        print 'PipelineStep_matchDistFreqs: Not re-running because output file exists, but input files don\'t!'
        return  {'mapfile': fileid}

    # find matches
    all_files_hosts = [(item.file, item.host) for item in map_full]
    dist_files = [item.file for item in map_dist]
    for i, (f, h) in enumerate(all_files_hosts):
        if f in dist_files:
            map_out.append(DataProduct(h, map_in[i].file, False))

    map_out.save(fileid)
    del(map_in)
    del(map_out)
    result = {'mapfile': fileid}

    return result
