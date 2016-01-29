#!/usr/bin/env python
import os
import re
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# mandatory arguments:
# cmdline for type of mapfile creation
# options: mapfile-dir, filename, identifier(name in parsetparset)
def plugin_main(args, **kwargs):
    result = {}
    datamap = None
    
    mergeLastGroup = False
    startFromZero = False
    NDPPPfill = False
    truncateLastSBs = False
    if 'mergeLastGroup' in kwargs:
        mergeLastGroup = string2bool(kwargs['mergeLastGroup'])
    if 'startFromZero' in kwargs:
        startFromZero = string2bool(kwargs['startFromZero'])  
    if 'truncateLastSBs' in kwargs:
        truncateLastSBs = string2bool(kwargs['truncateLastSBs'])  
    if 'NDPPPfill' in kwargs:
        NDPPPfill = string2bool(kwargs['NDPPPfill'])
    (datamap , firstmap) = _sort_ObsSB(kwargs['mapfile_in'], int(kwargs['numSB']), filllNDPPPdummies=NDPPPfill, 
                                       mergeLastGroup=mergeLastGroup, startFromZero=startFromZero,
                                       truncateLastSBs=truncateLastSBs)
    if datamap:
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
        datamap.save(fileid)
        result['mapfile'] = fileid
    if firstmap:
        if 'target_path' in kwargs and 'target_ext' in kwargs:
            for item in firstmap:
                basename = os.path.splitext(os.path.basename(item.file))[0]
                item.file = os.path.join(kwargs['target_path'],basename+kwargs['target_ext'])
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename']+'_groups')
        firstmap.save(fileid)        
        result['groupmap'] = fileid
    if 'numCh_per_SB' in kwargs:
        flagmap = _calc_edge_chans(datamap,int(kwargs['numCh_per_SB']))
        fileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename']+'_flagch')
        flagmap.save(fileid)        
        result['flagmap'] = fileid
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
        for i in xrange(len(group.file)):
            flaglist.extend(range(i*numch,i*numch+numch/edgeFactor))
            flaglist.extend(range((i+1)*numch-numch/edgeFactor,(i+1)*numch))
        outmap.append(DataProduct(group.host,str(flaglist).replace(' ',''),group.skip))
        print str(flaglist).replace(' ','')
    return outmap

# helper function
def _sort_ObsSB(Mapfile, SBperGroup, filllNDPPPdummies=False, mergeLastGroup=False, startFromZero=False, 
                truncateLastSBs=False):
    """
    Sorts the MSs in "inmap" into groups with "SBperGroup" consecutive subband numbers for each 
    observation number. Expects files that contain a string like: "L123456_SAP007_SB890" 
    or: "L123456_SB890". The files in the groups will be sorted by subband number.
    The hostname is taken from the first file found for each group.

    Mapfile           - Name of the mapfile with the input data, should contain single MSs
    SBperGroup        - Group that many subbands into one group.
    filllNDPPPdummies - If True, add dummy file-names for missing subbands, so that NDPPP can
                        fill the data with flagged dummy data.
    mergeLastGroup    - True: put the files from the last group that doesn't have SBperGroup subbands 
                        into the second last group (which will then have more than SBperGroup entries). 
                        False: keep inclomplete last group, or - with filllNDPPPdummies=True - fill
                        last group with dummies.
    startFromZero     - Start grouping with 0, even if there is no SB000 in the map.
    """
    if mergeLastGroup and truncateLastSBs:
        raise ValueError('_sort_ObsSB: Can either merge the last partial group or truncate at last full group, not both!')
    inmap = MiniMapfileManager.load(Mapfile)
    # I'm going to need that:
    PatternReg = re.compile(r'(L\d+)(_SAP\d+)?_(SB\d+)')
    # get smallest and largest subband number: sort by subband and get SB-number of first and last entry
    inmap.sortMS_by_SBnum()
    if not startFromZero:
        minSB = int(PatternReg.search(inmap.data[0].file).group(3)[2:])
    else:
        minSB = 0
    maxSB = int(PatternReg.search(inmap.data[-1].file).group(3)[2:])
    # sort data into dictionary
    sortDict = {}
    for item in inmap.data:
        if not item.skip:
            regerg = PatternReg.search(item.file)
            Obs = regerg.group(1)
            SBnum = int(regerg.group(3)[2:])
            SBgroup = int((SBnum-minSB)/SBperGroup)
            if not Obs in sortDict:
                sortDict[Obs] = { }
            if not SBgroup in sortDict[Obs]:
                replacestring = Obs+'_SBgr%03d-%d' %(SBgroup,SBperGroup)
                reffile = PatternReg.sub(replacestring,item.file,1)
                sortDict[Obs][SBgroup] = { 'host' : item.host , 'files' : [], 'firstfile' : reffile }
            #the data is sorted by SBnum, so if files with lower SBnum are not already 
            #in the list, then they are missing!
            while filllNDPPPdummies and len(sortDict[Obs][SBgroup]['files']) < ((SBnum-minSB) % SBperGroup) :
                sortDict[Obs][SBgroup]['files'].append('dummy_entry')
            sortDict[Obs][SBgroup]['files'].append(item.file)
    # now go through the dictionary and put the data into the new map
    newmap = MultiDataMap()
    firstfileMap = DataMap()
    numGroups = (maxSB-minSB+1)/SBperGroup
    SBs_in_last = (maxSB-minSB+1)%SBperGroup
    obsNames = sortDict.keys()
    obsNames.sort()
    for obs in obsNames:
        obsDict = sortDict[obs]
        for SBgroup in xrange(numGroups-1):
            if SBgroup in obsDict:
                while filllNDPPPdummies and len(obsDict[SBgroup]['files']) < SBperGroup :
                    obsDict[SBgroup]['files'].append('dummy_entry')
                newmap.append(MultiDataProduct(obsDict[SBgroup]['host'], obsDict[SBgroup]['files'], False))
                firstfileMap.append(DataProduct(obsDict[SBgroup]['host'], obsDict[SBgroup]['firstfile'], False))
        #work on last full group
        SBgroup = numGroups-1
        if SBgroup in obsDict:
            while filllNDPPPdummies and len(obsDict[SBgroup]['files']) < SBperGroup :
                obsDict[SBgroup]['files'].append('dummy_entry')
        if mergeLastGroup and SBs_in_last != 0:
            lastList = []
            if numGroups in obsDict:
                lastList = obsDict[numGroups]['files']
            while filllNDPPPdummies and len(lastList) < SBs_in_last :
                lastList.append('dummy_entry')
            obsDict[SBgroup]['files'].extend(lastList)
        newmap.append(MultiDataProduct(obsDict[SBgroup]['host'], obsDict[SBgroup]['files'], False))
        firstfileMap.append(DataProduct(obsDict[SBgroup]['host'], obsDict[SBgroup]['firstfile'], False))
        #need to process incomplete group
        if SBs_in_last != 0 and not mergeLastGroup and not truncateLastSBs:
            if numGroups in obsDict:
                while filllNDPPPdummies and len(obsDict[numGroups]['files']) < SBs_in_last :
                    obsDict[numGroups]['files'].append('dummy_entry')
                newmap.append(MultiDataProduct(obsDict[numGroups]['host'], obsDict[numGroups]['files'], False))
                firstfileMap.append(DataProduct(obsDict[numGroups]['host'], obsDict[numGroups]['firstfile'], False))
    #That should be it!
    return (newmap,firstfileMap)


class MiniMapfileManager(DataMap):
    # I didn't want to copy the full MapfileManager, but it should be possible to 
    # just copy&paste the methods in here into the "real" MapfileManager

    def sort_by_string(self, reverse=False):
        getfile = lambda item : item.file
        self.data.sort(key=getfile,reverse=reverse)

    def sortMS_by_SBnum(self, reverse=False):
        self._sortMS(False,True,reverse=reverse)

    def sortMS_by_Obsnum(self, reverse=False):
        self._sortMS(True,False,reverse=reverse)
       
    def sortMS_by_Obs_SB(self, reverse=False):
        self._sortMS(True,True,reverse=reverse)
 
    def sortMS_by_SB_Obs(self, reverse=False):
        self._sortMS(True,True,obsfirst=False,reverse=reverse)

    def _sortMS(self, byObsnum, bySBnum, obsfirst=True, reverse=False):
        #regerg = re.search(r'(L\d+)(_SAP\d+)?_(SB\d+)',item.file)
        PatternReg = re.compile(r'(L\d+)(_SAP\d+)?_(SB\d+)')
        if byObsnum and bySBnum:
            if obsfirst:
                def keyfunc(item):
                    regerg = PatternReg.search(item.file)
                    return regerg.group(1)+regerg.group(3)
                self.data.sort(key=keyfunc,reverse=reverse)
            else:
                def keyfunc(item):
                    regerg = PatternReg.search(item.file)
                    return regerg.group(3)+regerg.group(1)
                self.data.sort(key=keyfunc,reverse=reverse)
        elif byObsnum: 
            def keyfunc(item):
                regerg = PatternReg.search(item.file)
                return regerg.group(1)
            self.data.sort(key=keyfunc,reverse=reverse)
        elif bySBnum:
            def keyfunc(item):
                regerg = PatternReg.search(item.file)
                return regerg.group(3)
            self.data.sort(key=keyfunc,reverse=reverse)
        else:
            raise ValueError("_sortMS: byObsnum and bySBnum both false!")


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

