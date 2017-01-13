import os
import glob
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Updates the hosts in an input datamap

    Parameters
    ----------
    mapfile_in : str, optional
        Filename of datamap
    mapfile_dir: str, optional
        Directory containing mapfiles. All mapfiles in this directory will be
        updated

    """
    if 'mapfile_dir' in kwargs:
        mapfiles_in = glob.glob(os.path.join(kwargs['mapfile_dir'], '*'))
    else:
        mapfiles_in = [kwargs['mapfile_in']]

    if len(mapfiles_in) == 0:
        return

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

    for mapfile_in in mapfiles_in:
        try:
            map = DataMap.load(mapfile_in)
            for i in range(len(map)-len(hosts)):
                hosts.append(hosts[i])

            for item, host in zip(map, hosts):
                item.host = host

            map.save(mapfile_in)
        except:
            print('File {} does not appear to be a mapfile. Skipping it.'.format(mapfile_in))
