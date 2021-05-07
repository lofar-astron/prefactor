#! /usr/bin/env python3
"""
Script to sort a list of MSs by into frequency groups by time-stamp
"""
import pyrap.tables as pt
import sys, os
import numpy as np
import logging
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

def _calc_edge_chans(inmap, numch, edgeFactor=32):
    """
    Generates a map with strings that can be used as input for NDPPP to flag the edges 
    of the input MSs during (or after) concatenation.
    
    inmap      - MultiDataMap (not mapfilename!) with the files to be concatenated.
    numch      - Number of channels per input file (All files are assumed to have the same number 
                 of channels.)
    edgeFactor - Divisor to compute how many channels are to be flagged at beginning and end. 
                 (numch=64 and edgeFactor=32 means "flag two channels at beginning and two at end")
    """
    outmap = DataMap([])
    for group in inmap:
        flaglist = []
        for i in range(len(group.file)):
            flaglist.extend(range(i*numch,i*numch+numch/edgeFactor))
            flaglist.extend(range((i+1)*numch-numch/edgeFactor,(i+1)*numch))
        outmap.append(DataProduct(group.host,str(flaglist).replace(' ',''),group.skip))
        print('_calc_edge_chans: flaglist:', str(flaglist).replace(' ',''))
    return(outmap)

def input2bool(invar):
    if invar == None:
        return(None)
    if isinstance(invar, bool):
        return(invar)
    elif isinstance(invar, str):
        if invar.upper() == 'TRUE' or invar == '1':
            return(True)
        elif invar.upper() == 'FALSE' or invar == '0':
            return(False)
        else:
            raise ValueError('input2bool: Cannot convert string "'+invar+'" to boolean!')
    elif isinstance(invar, int) or isinstance(invar, float):
        return(bool(invar))
    else:
        raise TypeError('input2bool: Unsupported data type:'+str(type(invar)))
  
def input2int(invar):
    if invar == None:
        return(None)
    if isinstance(invar, int):
        return(invar)
    elif isinstance(invar, float):
        return(int(invar))
    elif isinstance(invar, str):
        if invar.strip().upper() == 'NONE' or invar == 'FALSE':
            return(None)
        else:
            return(int(invar))
    else:
        raise TypeError('input2int: Unsupported data type:'+str(type(invar)))


    
def main(ms_input, filename=None, mapfile_dir=None, numSB=-1, hosts=None, NDPPPfill=True, target_path=None, stepname=None,
         mergeLastGroup=False, truncateLastSBs=True, firstSB=None):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    filename: str
        Name of output mapfile
    mapfile_dir : str
        Directory for output mapfile
    numSB : int, optional 
        How many files should go into one frequency group. Values <= 0 mean put 
        all files of the same time-step into one group.
        default = -1
    hosts : list or str
        List of hostnames or string with list of hostnames
    NDPPPfill : bool, optional
        Add dummy file-names for missing frequencies, so that NDPPP can
        fill the data with flagged dummy data.
        default = True
    target_path : str, optional
        Change the path of the "groups" files to this. (I.e. write output files 
        into this directory with the subsequent NDPPP call.)
        default = keep path of input files
    stepname : str, optional
        Add this step-name into the file-names of the output files.
    mergeLastGroup, truncateLastSBs : bool, optional
        mergeLastGroup = True, truncateLastSBs = True:
          not allowed
        mergeLastGroup = True, truncateLastSBs = False:
          put the files from the last group that doesn't have SBperGroup subbands 
          into the second last group (which will then have more than SBperGroup entries). 
        mergeLastGroup = False, truncateLastSBs = True:
          ignore last files, that don't make for a full group (not all files are used).
        mergeLastGroup = False, truncateLastSBs = False:
          keep inclomplete last group, or - with NDPPPfill=True - fill
          last group with dummies.      
    firstSB : int, optional
        If set, then reference the grouping of files to this station-subband. As if a file 
        with this station-subband would be included in the input files.
        (For HBA-low, i.e. 0 -> 100MHz, 55 -> 110.74MHz, 512 -> 200MHz)

    Returns
    -------
    result : dict
        Dict with the name of the generated mapfile

    """

    NDPPPfill = input2bool(NDPPPfill)
    mergeLastGroup = input2bool(mergeLastGroup)
    truncateLastSBs = input2bool(truncateLastSBs)
    firstSB = input2int(firstSB)
    numSB = int(numSB)

    if not filename or not mapfile_dir:
        raise ValueError('sort_times_into_freqGroups: filename and mapfile_dir are needed!')
    if mergeLastGroup and truncateLastSBs:
        raise ValueError('sort_times_into_freqGroups: Can either merge the last partial group or truncate at last full group, not both!')
#    if mergeLastGroup:
#        raise ValueError('sort_times_into_freqGroups: mergeLastGroup is not (yet) implemented!')
    if type(ms_input) is str:
        if ms_input.startswith('[') and ms_input.endswith(']'):
            ms_list = [f.strip(' \'\"') for f in ms_input.strip('[]').split(',')]
        else:
            map_in = DataMap.load(ms_input)
            map_in.iterator = DataMap.SkipIterator
            ms_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        ms_list.append(f.strip(' \'\"'))
                else:
                    ms_list.append(fname.strip(' \'\"'))  
    elif type(ms_input) is list:
        ms_list = [str(f).strip(' \'\"') for f in ms_input]
    else:
        raise TypeError('sort_times_into_freqGroups: type of "ms_input" unknown!')

    if type(hosts) is str:
        hosts = [h.strip(' \'\"') for h in hosts.strip('[]').split(',')]
    if not hosts:
        hosts = ['localhost']
    numhosts = len(hosts)
    print("sort_times_into_freqGroups: Working on",len(ms_list),"files (including flagged files).")

    time_groups = {}
    # sort by time
    for i, ms in enumerate(ms_list):
        # work only on files selected by a previous step
        if ms.lower() != 'none':
            # use the slower but more reliable way:
            obstable = pt.table(ms, ack=False)
            timestamp = int(round(np.min(obstable.getcol('TIME'))))
            #obstable = pt.table(ms+'::OBSERVATION', ack=False)
            #timestamp = int(round(obstable.col('TIME_RANGE')[0][0]))
            obstable.close()
            if timestamp in time_groups:
                time_groups[timestamp]['files'].append(ms)
            else:
                time_groups[timestamp] = {'files': [ ms ], 'basename' : os.path.splitext(ms)[0] }
    print("sort_times_into_freqGroups: found",len(time_groups),"time-groups")

    # sort time-groups by frequency
    timestamps = list(time_groups.keys())
    timestamps.sort()   # not needed now, but later
    first = True
    nchans = 0
    for time in timestamps:
        freqs = []
        for ms in time_groups[time]['files']:
            # Get the frequency info
            sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
            freq = sw.col('REF_FREQUENCY')[0]            
            if first:
                file_bandwidth = sw.col('TOTAL_BANDWIDTH')[0]
                nchans = sw.col('CHAN_WIDTH')[0].shape[0]
                chwidth = sw.col('CHAN_WIDTH')[0][0]
                freqset = set([freq])
                first = False
            else:
                assert file_bandwidth == sw.col('TOTAL_BANDWIDTH')[0]
                assert nchans == sw.col('CHAN_WIDTH')[0].shape[0]
                assert chwidth == sw.col('CHAN_WIDTH')[0][0]
                freqset.add(freq)
            freqs.append(freq)
            sw.close()
        time_groups[time]['freq_names'] = list(zip(freqs,time_groups[time]['files']))
        time_groups[time]['freq_names'].sort(key=lambda pair: pair[0])
        #time_groups[time]['files'] = [name for (freq,name) in freq_names]
        #time_groups[time]['freqs'] = [freq for (freq,name) in freq_names]
    print("sort_times_into_freqGroups: Collected the frequencies for the time-groups")

    freqliste = np.array(list(freqset))
    freqliste.sort()
    freq_width = np.min(freqliste[1:]-freqliste[:-1])
    if file_bandwidth > freq_width:
        raise ValueError("Bandwidth of files is larger than minimum frequency step between two files!")
    if file_bandwidth < (freq_width*0.51):
        #raise ValueError("Bandwidth of files is smaller than 51% of the minimum frequency step between two files! (More than about half the data is missing.)")
        logging.warning("Bandwidth of files is smaller than 51% of the minimum frequency step between two files! (More than about half the data is missing.)")
    #the new output map
    filemap = MultiDataMap()
    groupmap = DataMap()
    # add 1% of the SB badwidth in case maxfreq might be "exactly" on a group-border
    maxfreq = np.max(freqliste)+freq_width*0.51
    if firstSB != None:
        if freqliste[0] < 100e6:
            # LBA Data
            minfreq = (float(firstSB)/512.*100e6)-freq_width/2.
        elif freqliste[0] > 100e6 and freqliste[0] < 200e6:
            # HBA-Low
            minfreq = (float(firstSB)/512.*100e6)+100e6-freq_width/2.
        elif freqliste[0] > 200e6 and freqliste[0] < 300e6:
            # HBA-high
            minfreq = (float(firstSB)/512.*100e6)+200e6-freq_width/2.
        else:
            raise ValueError('sort_times_into_freqGroups: Frequency of lowest input data is higher than 300MHz!')        
        if np.min(freqliste) < minfreq:
            raise ValueError('sort_times_into_freqGroups: Frequency of lowest input data is lower than reference frequency!')
    else:
        minfreq = np.min(freqliste)-freq_width/2.
    groupBW = freq_width*numSB
    if groupBW < 1e6 and groupBW > 0:
        print('sort_times_into_freqGroups: ***WARNING***: Bandwidth of concatenated MS is lower than 1 MHz. This may cause conflicts with the concatenated file names!')
    if groupBW < 0:
    # this is the case for concatenating all subbands
       groupBW = maxfreq-minfreq
    truncateLastSBs = input2bool(False)
    NDPPPfill = input2bool(True)
    freqborders = np.arange(minfreq,maxfreq,groupBW)
    if mergeLastGroup:
        freqborders[-1] = maxfreq
    elif truncateLastSBs:
        pass #nothing to do! # left to make the logic more clear!
    elif not truncateLastSBs and NDPPPfill:
        freqborders = np.append(freqborders,(freqborders[-1]+groupBW))
    elif not truncateLastSBs and not NDPPPfill:
        freqborders = np.append(freqborders,maxfreq)

    freqborders = freqborders[freqborders>(np.min(freqliste)-groupBW)]
    ngroups = len(freqborders)-1
    if ngroups == 0:
        raise ValueError('sort_times_into_freqGroups: Not enough input subbands to create at least one full (frequency-)group!')
    
    print("sort_times_into_freqGroups: Will create",ngroups,"group(s) with",numSB,"file(s) each.")

    hostID = 0
    for time in timestamps:
        (freq,fname) = time_groups[time]['freq_names'].pop(0)
        for groupIdx in range(ngroups):
            files = []
            skip_this = True
            filefreqs_low = np.arange(freqborders[groupIdx],freqborders[groupIdx+1],freq_width)
            for lower_freq in filefreqs_low:
                if freq > lower_freq and freq < lower_freq+freq_width:
                    assert freq!=1e12
                    files.append(fname)
                    if len(time_groups[time]['freq_names'])>0:
                        (freq,fname) = time_groups[time]['freq_names'].pop(0)
                    else:
                        (freq,fname) = (1e12,'This_shouldn\'t_show_up')
                    skip_this = False
                elif NDPPPfill:
                    files.append('dummy.ms')
            if not skip_this:
                filemap.append(MultiDataProduct(hosts[hostID%numhosts], files, skip_this))
                freqID = int((freqborders[groupIdx]+freqborders[groupIdx+1])/2e6)
                groupname = time_groups[time]['basename']+'_%Xt_%dMHz.ms'%(time,freqID)
                if type(stepname) is str:
                    groupname += stepname
                if type(target_path) is str:
                    groupname = os.path.join(target_path,os.path.basename(groupname))
                groupmap.append(DataProduct(hosts[hostID%numhosts],groupname, skip_this))
        orphan_files = len(time_groups[time]['freq_names'])
        if freq < 1e12:
            orphan_files += 1
        if orphan_files > 0:
            print("sort_times_into_freqGroups: Had %d unassigned files in time-group %xt."%(orphan_files, time))
    filemapname = os.path.join(mapfile_dir, filename)
    filemap.save(filemapname)
    groupmapname = os.path.join(mapfile_dir, filename+'_groups')
    groupmap.save(groupmapname)
    # genertate map with edge-channels to flag
    flagmap = _calc_edge_chans(filemap, nchans)
    flagmapname = os.path.join(mapfile_dir, filename+'_flags')
    flagmap.save(flagmapname)
    result = {'mapfile': filemapname, 'groupmapfile': groupmapname, 'flagmapfile': flagmapname}
    return(result)

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
        return(
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return(':'.join((self.host, str(self.file))))

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
        print('setting filelist from DataProduct')
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print('setting filelist from DataMap')
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
            for i in range(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)



if __name__ == '__main__':
    import optparse
    import glob
    import random

    opt = optparse.OptionParser(usage='%prog [options] <MSPattern> \n')
    opt.add_option('-v', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-r', '--randomize', help='Randomize order of the input files. (default=False)', action='store_true', default=False)
    opt.add_option('-d', '--decimate', help='Remove every 10th file (after randomization if that is done). (default=False)', action='store_true', default=False)
    opt.add_option('-n', '--numbands', help='Number of how many files should be grouped together in frequency. (default=all files in one group)', type='int', default=-1)
    opt.add_option('-f', '--filename', help='Name for the mapfiles to write. (default=\"test.mapfile\")', type='string', default='test.mapfile')

    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 1:
        opt.print_help()
        sys.exit()

    # first argument: pattern for measurement-sets
    inMSs = glob.glob(args[0])
    if options.randomize:
        random.shuffle(inMSs)
    if options.decimate:
        for i in range((len(inMSs)-1),-1,-10):
            inMSs.pop(i)

    ergdict = main(inMSs, options.filename, '.', numSB=options.numbands, hosts=None, NDPPPfill=True)

    groupmap = DataMap.load(ergdict['groupmapfile'])
    filemap = MultiDataMap.load(ergdict['mapfile'])
    print("len(groupmap) : %d , len(filemap) : %d " % (len(groupmap),len(filemap)))
    if len(groupmap) != len(filemap):
        print("groupmap and filemap have different length!")
        sys.exit(1)
    for i in range(len(groupmap)):
        print("Group \"%s\" has %d entries."%(groupmap[i].file,len(filemap[i].file)))
