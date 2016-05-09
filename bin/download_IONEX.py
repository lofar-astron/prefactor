#! /usr/bin/env python
"""
Script to sort a list of MSs by into frequency groups by time-stamp
"""
import pyrap.tables as pt
import sys, os
import numpy as np
import RMextract.getIONEX as ionex;
import RMextract.PosTools as PosTools

def input2strlist_nomapfile(invar):
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

def download_IONEX(ms_input,server="ftp://ftp.unibe.ch/aiub/CODE/",prefix='CODG',ionexPath="IONEXdata/"):
    """
    Call getIONEX() once for every day that is covered by one of the input MSs

    
    Parameters
    ----------
    ms_input : list or str
        list of MS filenames, or string with list, or single MS filename
    server : str
        URL of the server where the IONEX files can be found
    prefix : str
        prefix of the IONEX files
    ionexPath : str
        directory where to store the IONEX files
    """
    mslist = input2strlist_nomapfile(ms_input)
    
    dates_done = []
    numfiles = len(mslist)
    if numfiles>100:
        printstat = True
    else:
        printstat = False
    for (i,MS) in enumerate(mslist):
        if printstat and (i+1)%100==0:
            print "Working on file %d out of %d"%(i+1,numfiles)
        obstable = pt.table(MS+'::OBSERVATION', ack=False)
        time_range = obstable.col('TIME_RANGE')[0]
        for timestamp in time_range:
            date_parms = PosTools.obtain_observation_year_month_day_fraction(timestamp)
            datestr = "%d-%0.2d-%0.2d" % date_parms[0:3]
            if datestr not in dates_done:
                print "Downloading IONEX file for:",datestr
                ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath)
                dates_done.append(datestr)

def main(msname, ionex_server="ftp://ftp.unibe.ch/aiub/CODE/", ionex_prefix='CODG', ionexPath="IONEXdata/"):
    # no download if ionex_server == None
    if ionex_server.strip(' []\'\"').lower() != 'none':
        download_IONEX(msname,server=ionex_server,prefix=ionex_prefix,ionexPath=ionexPath)
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Download IONEX files for given measurement sets.')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One or more MSs for which the IONEX data should be downloaded.')
    parser.add_argument('--server', type=str,                         
                        help='URL of the server to use. (default: \"ftp://ftp.unibe.ch/aiub/CODE/\")')
    parser.add_argument('--prefix', type=str,                         
                        help='Prefix of the IONEX files. (default: \"CODG\")')
    parser.add_argument('--destination', '--path', type=str,                         
                        help='Directory where to store the IONEX files. (default: \"IONEXdata/\"')

    args = parser.parse_args()
    
    server = 'ftp://ftp.unibe.ch/aiub/CODE/'
    if args.server:
        server = args.server
    prefix = 'CODG'
    if args.prefix:
        prefix = args.prefix
    path = 'IONEXdata/'
    if args.destination:
        path = args.destination

    download_IONEX(args.MSfile, server, prefix, path)
