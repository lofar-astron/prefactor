import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a multi-mapfile by compressing input mapfile items

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing MS files
    nitems_to_compress: int
        Number of input items to compress into each output item. Set to zero or
        negative number to compress all input items to a single output item
        (default = -1). If greater than 0, the input map must have only one file
        per group (i.e., it must be a normal DataMap)
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
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']
    if 'nitems_to_compress' in kwargs:
        nitems_to_compress = int(float(kwargs['nitems_to_compress']))
    else:
        nitems_to_compress = -1

    map_in = MultiDataMap.load(mapfile_in)
    map_out = MultiDataMap([])
    map_in.iterator = DataMap.SkipIterator
    if nitems_to_compress > 0:
        all_files = []
        for item in map_in:
            if type(item.file) is list:
                all_files.extend(item.file)
            else:
                all_files.append(item.file)
        file_groups = [all_files[i:i+nitems_to_compress] for i  in range(0, len(all_files), nitems_to_compress)]
        all_hosts = [item.host for item in map_in]
        host_groups = [all_hosts[i:i+nitems_to_compress] for i  in range(0, len(all_hosts), nitems_to_compress)]
        for file_list, host_list in zip(file_groups, host_groups):
            map_out.data.append(MultiDataProduct(host_list[0], file_list, False))
    else:
        file_list = []
        for item in map_in:
            if type(item.file) is list:
                file_list.extend(item.file)
            else:
                file_list.append(item.file)
        host_list = [item.host for item in map_in]
        map_out.data.append(MultiDataProduct(host_list[0], file_list, False))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

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


class MultiDataProduct(DataProduct):
    """
    Class representing multiple files in a DataProduct.
    """
    def __init__(self, host=None, file=None, skip=True):
        super(MultiDataProduct, self).__init__(host, file, skip)
        if not file:
            self.file = list()
        else:
            self._set_file(file)

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': '{1}', 'skip': {2}}}".format(self.host,
                '[{}]'.format(','.join(self.file)), str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return ':'.join((self.host, str(self.file)))

    def _set_file(self, data):
        try:
            # Try parsing as a list
            if isinstance(data, list):
                self.file = data
            if isinstance(data, DataProduct):
                self._from_dataproduct(data)
            if isinstance(data, DataMap):
                self._from_datamap(data)

        except TypeError:
            raise DataProduct("No known method to set a filelist from %s" % str(file))

    def _from_dataproduct(self, prod):
        print 'setting filelist from DataProduct'
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print 'setting filelist from DataMap'
        filelist = {}
        for item in inmap:
            if not item.host in filelist:
                filelist[item.host] = []
            filelist[item.host].append(item.file)
        self.file = filelist['i am']

    def append(self, item):
        self.file.append(item)


class MultiDataMap(DataMap):
    """
    Class representing a specialization of data-map, a collection of data
    products located on the same node, skippable as a set and individually
    """
    def __init__(self, data=list(), iterator=iter):
        super(MultiDataMap, self).__init__(data, iterator)

    @classmethod
    def load(cls, filename):
        """Load a data map from file `filename`. Return a DataMap instance."""
        with open(filename) as f:
            datamap = eval(f.read())
            for i, d in enumerate(datamap):
                file_entry = d['file']
                if file_entry.startswith('[') and file_entry.endswith(']'):
                    file_list = [e.strip(' \'\"') for e in file_entry.strip('[]').split(',')]
                    datamap[i] = {'host': d['host'], 'file': file_list, 'skip': d['skip']}
            return cls(datamap)

    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            mdpdict = {}
            data.iterator = DataMap.SkipIterator
            for item in data:
                if not item.host in mdpdict:
                    mdpdict[item.host] = []
                mdpdict[item.host].append(item.file)
            mdplist = []
            for k, v in mdpdict.iteritems():
                mdplist.append(MultiDataProduct(k, v, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        elif isinstance(data, MultiDataProduct):
            self._set_data(data, dtype=MultiDataProduct)
        elif not data:
            pass
        else:
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)

