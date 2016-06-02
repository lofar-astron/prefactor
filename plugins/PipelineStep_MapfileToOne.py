#!/usr/bin/env python
import os
import re
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import argparse
from argparse import RawTextHelpFormatter


# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    result = {}
    datamap = None
    if kwargs['method'] == 'mapfile_all_to_one':
        datamap = _create_mapfile_ato(kwargs['mapfile_in'])
    if kwargs['method'] == 'dummy':
        datamap = MapfileManager()
        datamap.expand(kwargs['number'])

    if datamap:
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
        datamap.save(fileid)
        result['mapfile'] = fileid
    return result


# helper function
def _create_mapfile_ato(inmap):
    return MultiDataMap(DataMap.load(inmap))




class MapfileManager(DataMap):

    def __init__(self, folder=None, pattern=None, excludepattern=False):
        super(MapfileManager, self).__init__([])
        if folder:
            self.from_folder(folder, pattern, excludepattern)
        #self.map = DataMap([])

    def expand(self, number, hostlist=None, filelist=None):
        if hostlist:
            if len(hostlist) != number:
                print 'Error: length of hostlist should correspond to number of expansions'
                exit(1)
        else:
            print 'Info: no hostlist given. Will use "localhost" instead'
            hostlist = []
            for item in range(number):
                hostlist.append('localhost')

        if filelist:
            if len(filelist) != number:
                print 'Error: length of hostlist should correspond to number of expansions'
                exit(1)
        else:
            print 'Info: no filelist given. Will use "dummy" instead'
            filelist = []
            for item in range(number):
                filelist.append('dummy')

        prodlist = []
        for h, f in zip(hostlist, filelist):
            prodlist.append(DataProduct(h, f, False))

        self._set_data(prodlist)

    def insert(self, place, data):
        self._insert(place, data)

    def _insert(self, place, data, dtype=DataProduct):
        try:
            if isinstance(data, dtype):
                self._data.insert(place, data)
            elif isinstance(data, dict):
                self._data.insert(place, dtype(**data))
            elif all(isinstance(item, tuple) for item in data):
                self._data.insert(place, dtype(*data))
            else:
                raise TypeError
        except TypeError:
            raise DataMapError("Failed to validate data map: %s" % repr(data))

    def append(self, data):
        self._append(data)

    def _append(self, data, dtype=DataProduct):
        try:
            if isinstance(data, dtype):
                self._data.append(data)
            elif isinstance(data, dict):
                self._data.append(dtype(**data))
            elif all(isinstance(item, tuple) for item in data):
                self._data.append(dtype(*data))
            else:
                raise TypeError
        except TypeError:
            raise DataMapError("Failed to validate data map: %s" % repr(data))

    def delete(self, host=None, data=None, skip=None, pattern=None):
        if pattern:
            #convert pattern-strin to RegExp format and make compiled RegExp
            rePattern = pattern.strip().replace('.','\.').replace('?','.').replace('*','.*')+'$'
            PatternReg = re.compile(rePattern)        
        for i, item in enumerate(self._data):
            if item.host == host or item.file == data or item.skip == skip:
                del self._data[i]
            if pattern:
                if PatternReg.match(item.file):
                    del self._data[i]

    def from_folder(self, folder, pattern=None, exclude_pattern=False):
        measurements = os.listdir(folder)
        measurements.sort()
        if pattern:
            #convert pattern-strin to RegExp format and make compiled RegExp
            rePattern = pattern.strip().replace('.','\.').replace('?','.').replace('*','.*')+'$'
            PatternReg = re.compile(rePattern)
        for ms in measurements:
            if pattern:
                if not exclude_pattern and PatternReg.match(ms):
                    self._append(DataProduct('localhost', folder + '/' + ms, False))
                elif exclude_pattern and not PatternReg.match(ms):
                    self._append(DataProduct('localhost', folder + '/' + ms, False))
            else:
                self._append(DataProduct('localhost', folder + '/' + ms, False))

    def set_data_from_parset(self, parset, identifier):
        dps = parset.makeSubset(
            parset.fullModuleName('DataProducts') + '.'
        )
        self._set_data([
            tuple(os.path.join(location, filename).split(':')) + (skip,)
            for location, filename, skip in zip(
                dps.getStringVector(identifier + '.locations'),
                dps.getStringVector(identifier + '.filenames'),
                dps.getBoolVector(identifier + '.skip'))
        ])

    def from_parts(self, host='localhost', data='dummy', skip=False, ntimes=1):
        hostlist = self._input_to_list(host)
        datalist = self._input_to_list(data)
        skiplist = self._input_to_list(skip)
        if len(hostlist) is not len(datalist) or len(hostlist) is not len(skiplist) or len(hostlist) is not ntimes:
            print 'Length of parts is not equal. Will expand to max length given.'
            maxval = max(len(hostlist), len(datalist), len(skiplist), ntimes)
            lastval = hostlist[-1]
            if len(hostlist) is not maxval:
                for x in range(len(hostlist), maxval):
                    hostlist.append(lastval)

            lastval = datalist[-1]
            if len(datalist) is not maxval:
                for x in range(len(datalist), maxval):
                    datalist.append(lastval)

            lastval = skiplist[-1]
            if len(skiplist) is not maxval:
                for x in range(len(skiplist), maxval):
                    skiplist.append(lastval)
        prodlist = []
        for h, f, z in zip(hostlist, datalist, skiplist):
            prodlist.append(DataProduct(h, f, z))
        self._set_data(prodlist)

    def _input_to_list(self, data):
        datalist = []
        if isinstance(data, list):
            datalist = data
        else:
            datalist.append(data)
        return datalist


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
        print 'FILE: ', self.file

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
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
    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            datalist = []
            hostname = 'localhost'
            data.iterator = DataMap.SkipIterator
            for item in data:
                datalist.append(item.file)
                hostname = item.host
            mdplist = []
            mdplist.append(MultiDataProduct(hostname,datalist, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        else:
            print 'HELP: ', data
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)

if __name__ == '__main__':
    descriptiontext = "This script lets you create Mapfiles for the Lofar Pipeline Framework"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('name', help='Give the name of the file to output')
    parser.add_argument('-n', '--number',     help='Number of times to expand a dummy mapfile', type=int, default=0)
    parser.add_argument('-M', '--mapfile_in', help='Input mapfile to be changed.')
    parser.add_argument('-f', '--folder',     help='Path to the directory containing the data.', default=None)
    parser.add_argument('-p', '--pattern',    help='Pattern the files should contain (*.MS).')
    parser.add_argument('-e', '--exclude_pattern',
                        help='Files with the given pattern -p should be ignored', action='store_true')
    parser.add_argument('-s', '--suffix',     help='Suffix to add to the file names.')
    parser.add_argument('-m', '--method',     help='Method of creating the mapfile.', default='mapfile_from_folder')
    parser.add_argument('-d', '--directory',  help='Directory to store the mapfile in.', default='.')

    args = parser.parse_args()
    kwarg = {
        'filename': args.name,
        'number': args.number,
        'mapfile_in': args.mapfile_in,
        'folder': args.folder,
        'pattern': args.pattern,
        'exclude_pattern': args.exclude_pattern,
        'add_suffix_to_file': args.suffix,
        'method': args.method,
        'mapfile_dir': args.directory}

    if args.number != 0:
        kwarg['method'] = 'dummy'
    if args.suffix:
        kwarg['method'] = 'add_suffix_to_file'

    plugin_main([], **kwarg)
