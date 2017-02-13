import os
import casacore.tables as pt
from lofarpipe.support.data_map import DataMap, DataProduct


def get_distributed_indices(start, end, n):
    """Returns list of evenly distributed indices for given range and number"""
    n = min(end+1, max(n, 2))
    step = (end-start)/float(n-1)
    return [int(round(start+x*step)) for x in xrange(n)]


def plugin_main(args, **kwargs):
    """
    Makes a mapfile with the MSs spread across the full bandwidth

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
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
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if 'num' in kwargs:
        num = int(kwargs['num'])
    else:
        num = 6
    fileid = os.path.join(mapfile_dir, filename)

    map_in = DataMap.load(mapfile_in)
    map_in.iterator = DataMap.SkipIterator
    map_out = DataMap()
    map_out.data = []
    map_out._data = []

    # do not re-run if we already ran, and input files are deleted.
    if os.path.exists(fileid) and not os.path.exists(map_in[0].file):
        print 'PipelineStep_selectDistFreqs: Not re-running because output file exists, but input files don\'t!'
        return  {'mapfile': fileid}

    #sort into frequency groups
    freq_groups = {}
    hosts = []
    for item in map_in:
        if '[' in item.file and ']' in item.file:
            files = item.file.strip('[]').split(',')
            files = [f.strip() for f in files]
            ms_file = files[0]
        else:
            ms_file = item.file

        # Get the frequency info from the MS file
        sw = pt.table(ms_file+'::SPECTRAL_WINDOW', ack=False)
        freq = int(sw.col('REF_FREQUENCY')[0])
        sw.close()
        if freq in freq_groups:
            freq_groups[freq].append(item.file)
        else:
            freq_groups[freq] = [item.file]
        if not item.host in hosts:
            hosts.append(item.host)

    # find maximum number of files per frequency-group
    maxfiles = max([len(group) for group in freq_groups.values()])

    # select frequencies
    freqs =  freq_groups.keys()
    freqs.sort()
    num_freqs = len(freqs)
    if num > num_freqs:
        print 'PipelineStep_selectDistFreqs: less than %d frequency groups found, contiunig with %d groups.'%(num, num_freqs)
        num = num_freqs
    dist_ind = get_distributed_indices(0, num_freqs-1, num)
    selfreqs = [freqs[ind] for ind in dist_ind]
    if len(selfreqs) < 1:
        print "PipelineStep_selectDistFreqs: Selected less than one frequency bands."
        raise ValueError("Selected less than one frequency bands.")

    all_files = []
    for selfreq in selfreqs:
        all_files.extend(freq_groups[selfreq])

    # extend the hosts-list
    for i in range(len(all_files)-len(hosts)):
        hosts.append(hosts[i])

    # fill the output-map
    for (host,fname) in zip(hosts,all_files):
        map_out.append(DataProduct(host, fname, False))

    map_out.save(fileid)
    del(map_in)
    del(map_out)
    result = {'mapfile': fileid}

    return result
