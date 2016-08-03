import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile for list of files

    Parameters
    ----------
    files : list or str
        List of files or mapfile with such a list as the only entry. May be
        given as a list of strings or as a string (e.g.,
        '[s1.skymodel, s2.skymodel]'
    hosts : list or str
        List of hosts/nodes. May be given as a list or as a string (e.g.,
        '[host1, host2]'
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    if type(kwargs['files']) is str:
        try:
            # Check if input is mapfile containing list as a string
            map_in = DataMap.load(kwargs['files'])
            in_files = [item.file for item in map_in]
            files = []
            for f in in_files:
                files += f.strip('[]').split(',')
        except:
            files = kwargs['files']
            files = files.strip('[]').split(',')
        files = [f.strip() for f in files]
    if type(kwargs['hosts']) is str:
        hosts = kwargs['hosts'].strip('[]').split(',')
        hosts = [h.strip() for h in hosts]
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    for i in range(len(files)-len(hosts)):
        hosts.append(hosts[i])

    map_out = DataMap([])
    for h, f in zip(hosts, files):
        map_out.data.append(DataProduct(h, f, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
