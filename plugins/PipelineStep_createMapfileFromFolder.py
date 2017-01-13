import os
import re
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile for files in a folder

    Parameters
    ----------
    folder : str
        Path to files to add to mapfile
    pattern : str
        Pattern to use to filter files to add
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    folder = kwargs['folder']
    pattern = kwargs['pattern']
    measurements = os.listdir(folder)
    measurements.sort()

    #convert pattern-string to RegExp format and make compiled RegExp
    rePattern = pattern.strip().replace('.','\.').replace('?','.').replace('*','.*')+'$'
    PatternReg = re.compile(rePattern)

    # Find files
    files = []
    for ms in measurements:
        if PatternReg.match(ms):
            files.append(folder + '/' + ms)

    # Determine hosts: use PBS_NODEFILE env variable if it exists; otherwise set
    # to localhost
    try:
        pbs_file = os.environ['PBS_NODEFILE']
        hosts = []
        with open(pbs_file, 'r') as f:
            for line in f:
                node_name = line.split()[0]
                if node_name not in hosts:
                    hosts.append(node_name)
    except KeyError:
        hosts = ['localhost']

    for i in range(len(files)-len(hosts)):
        hosts.append(hosts[i])

    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    map_out = DataMap([])
    for h, f in zip(hosts, files):
        map_out.data.append(DataProduct(h, f, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
