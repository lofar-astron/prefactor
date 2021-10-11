#! /usr/bin/env python
"""
Script to sort a list of MSs by into frequency groups by time-stamp
"""
import pyrap.tables as pt
import sys, os
import numpy as np
import logging

########################################################################
def input2strlist_nomapfile(invar):
   """ 
   from bin/download_IONEX.py
   give the list of MSs from the list provided as a string
   """
   str_list = None
   if type(invar) is str:
       if invar.startswith('[') and invar.endswith(']'):
           str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
       else:
           str_list = [invar.strip(' \'\"')]
   elif type(invar) is list:
       str_list = [str(f).strip(' \'\"') for f in invar]
   else:
       raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
   return str_list
########################################################################
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
########################################################################
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


########################################################################
def main(MSfile, numSB=10, NDPPPfill=True, stepname=None, mergeLastGroup=False, truncateLastSBs=True, firstSB=None):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    MSfile : list or str
        List of MS filenames, or string with list, or path to a mapfile
    numSB : int, optional 
        How many files should go into one frequency group. Values <= 0 mean put 
        all files of the same time-step into one group.
        default = -1
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

    ms_list = input2strlist_nomapfile(MSfile)
    logging.info("Working on" + str(len(ms_list)) + "files (including flagged files).")


    ## here the real work starts
    time_groups = {}
    # sort by time
    for i, ms in enumerate(ms_list):
        # work only on files selected by a previous step
        if ms.lower() != 'none':
            # use the slower but more reliable way:
            obstable = pt.table(ms, ack=False)
            timestamp = int(round(np.min(obstable.getcol('TIME'))))
            obstable.close()
            if timestamp in time_groups:
                time_groups[timestamp]['files'].append(ms)
            else:
                time_groups[timestamp] = {'files': [ ms ], 'basename' : os.path.splitext(ms)[0] }
    logging.info("Found" + str(len(time_groups))  + "time-groups")

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
    logging.info("Collected the frequencies for the time-groups.")

    freqliste = np.array(list(freqset))
    freqliste.sort()
    freq_width = np.min(freqliste[1:]-freqliste[:-1])
    if file_bandwidth > freq_width:
        raise ValueError("Bandwidth of files is larger than minimum frequency step between two files!")
    if file_bandwidth < (freq_width*0.51):
        #raise ValueError("Bandwidth of files is smaller than 51% of the minimum frequency step between two files! (More than about half the data is missing.)")
        logging.warning("Bandwidth of files is smaller than 51% of the minimum frequency step between two files! (More than about half the data is missing.)")
    #the new output map

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
    
    logging.info("Will create " + str(ngroups) + " group(s) with " + str(numSB)  + "file(s) each.")

    groupnames = []
    filenames  = {}
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
                freqID = int((freqborders[groupIdx]+freqborders[groupIdx+1])/2e6)
                groupname = time_groups[time]['basename']+'_%Xt_%dMHz'%(time,freqID)
                if type(stepname) is str:
                    groupname += stepname
                groupname = os.path.basename(groupname)
                groupnames.append(groupname)
                filenames[groupname] = files

        orphan_files = len(time_groups[time]['freq_names'])
        if freq < 1e12:
            orphan_files += 1
        if orphan_files > 0:
            logging.info("Had %d unassigned files in time-group %xt."%(orphan_files, time))

    nr_of_groups = len(groupnames)
    total_bandwidth = nr_of_groups * groupBW
    results = {'filenames': filenames, 'groupnames': groupnames, 'total_bandwidth': total_bandwidth}

    return(results)
    
########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Identify fully flagged antennas.')

    parser.add_argument('MSfiles', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the antenna list.')
    parser.add_argument('--numbands', type=int, default=10,
                        help='Number of how many files should be grouped together in frequency. (default=10)')
    parser.add_argument('--NDPPPfill', type=bool, default=True,
                        help='Add dummy file-names for missing frequencies, so that NDPPP can fill the data with flagged dummy data.')
    parser.add_argument('--stepname', type=str, default=None,
                        help='Add this step-name into the file-names of the output files.')
    parser.add_argument('--mergeLastGroup', type=bool, default=False,
                        help='')
    parser.add_argument('--truncateLastSBs', type=bool, default=True,
                        help='')
    parser.add_argument('--firstSB', type=int, default=None,
                        help='If set, then reference the grouping of files to this station-subband. As if a file with this station-subband would be included in the input files.')



    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    # first argument: pattern for measurement-sets
    main(args.MSfiles, numSB=args.numbands, NDPPPfill=args.NDPPPfill, stepname=args.stepname, mergeLastGroup=args.mergeLastGroup, truncateLastSBs=args.truncateLastSBs, firstSB=args.firstSB)
