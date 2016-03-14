import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Trims a string from filenames in a mapfile

    Note that everything from the last instance of the matching string to the
    end is trimmed.

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to trim
    trim_str : str
        String to remove
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    counter : int
        If counter is greater than 0, replace "image32" with "image42". This is
        a special argument for facetselfcal looping only

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']
    trim_str = kwargs['trim']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if 'counter' in kwargs:
        counter = int(kwargs['counter'])
    else:
        counter = 0

    map_out = DataMap([])
    map_in = DataMap.load(mapfile_in)

    for i, item in enumerate(map_in):
        index = item.file.rfind(trim_str)
        if index >= 0:
            item_trim = item.file[:index]
            if counter > 0:
                item_trim = item_trim.replace('image32', 'image42')
            map_out.data.append(DataProduct(item.host, item_trim,
                item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
