#!/usr/bin/env python
# -* coding: utf-8 -*-

"""
Adds phases (=0) and amplitudes (=1) to any missing station if they appear in an h5parm, but not in a particular soltab. 

Created on Tue Mar 14 2019

@author: Alexander Drabent
"""

import argparse
import numpy as np
import os
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import logging

def main(h5parmfile, refh5 = None, solset='sol000', refsolset='sol000', soltab_in='phase000', soltab_out='GSMphase', filter='[CR]S*&', bad_antennas='[CR]S*&'):


    if refh5 == None:
        refh5 = h5parmfile
        pass
    
    if not os.path.exists(h5parmfile) or not os.path.exists(refh5):
        logging.error("H5parm file %s doesn't exist!" % h5parmfile)
        return(1)

    if soltab_in == soltab_out:
        logging.error("Output soltab has to be different from input soltab!")
        return(1)

    ### Open up the h5parm, get an example value
    data       = h5parm(h5parmfile, readonly = False)
    refdata    = h5parm(refh5, readonly = False)
    solset     = data.getSolset(solset)
    refsolset  = refdata.getSolset(refsolset)
    
    ### Get antenna information of the solset
    ref_station_names = sorted(refsolset.getAnt().keys())
    ref_station_names = [ref_station_name.decode() for ref_station_name in ref_station_names]
    bad_antennas_list = bad_antennas.lstrip(filter).replace('!','').replace('*','').replace('&','').split(';')
    new_station_names = [ ref_station_name for ref_station_name in ref_station_names if ref_station_name not in bad_antennas_list ]
    
    ### Load antenna list of input soltab
    logging.info('Processing solution table %s'%(soltab_in))
    soltab            = solset.getSoltab(soltab_in)
    old_station_names = list(soltab.ant[:])
    soltab_type       = soltab.getType()
    out_axes_vals     = []
    out_axes          = soltab.getAxesNames()
    out_axes.insert(1, out_axes.pop(out_axes.index('ant'))) ## put antenna table to the second place

    ### resort the val and weights array according to the new axes order
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=out_axes, weight=True):
        vals = reorderAxes( vals, soltab.getAxesNames(), out_axes)
        weights = reorderAxes( weights,  soltab.getAxesNames(), out_axes )

    ### setting the proper soltab_axis
    for axis in out_axes:
        if axis == 'time':
            out_axes_vals.append(soltab.time)
        elif axis == 'freq':
            out_axes_vals.append(soltab.freq)
        elif axis == 'ant':
            out_axes_vals.append(new_station_names)
        elif axis == 'pol':
            out_axes_vals.append(soltab.pol)
        elif axis == 'dir':
            out_axes_vals.append(soltab.dir)
        else:
            logging.error('Unknown axis in soltab: ' + str(axis))
            data.close()
            refdata.close()
            return 1           

    ### just copy if number of antennas is the same
    if len(new_station_names) == len(old_station_names):
        logging.warning('Station list in soltab ' + str(soltab_in) + ' matches the station list in the selected solset. Data will just be copied.')
        new_soltab = solset.makeSoltab(soltype=soltab_type, soltabName=soltab_out, axesNames=out_axes, axesVals=out_axes_vals, vals=vals, weights=weights)
    
    ### add missing stations with zero phase/uniform amplitudes and uniform weights
    elif len(new_station_names) > len(old_station_names):

        # define proper shape of array
        logging.info('Adding missing antennas to the soltab: ' + str(soltab_out))
        dimension      = list(np.shape(vals))
        dimension[1]   = len(new_station_names)
        added_stations = []
        
        if soltab_type == 'amplitude':
            new_vals    = np.ones(tuple(dimension))
        else:
            new_vals    = np.zeros(tuple(dimension))
        new_weights     = np.ones(tuple(dimension))
        
        for i, new_station in enumerate(new_station_names):
            if new_station in old_station_names:
                ant_index        = list(soltab.ant).index(new_station)
                new_vals[:,i]    = vals[:,ant_index]
                new_weights[:,i] = weights[:,ant_index]
            else:
                added_stations.append(new_station)
            
        new_soltab = solset.makeSoltab(soltype=soltab_type, soltabName=soltab_out, axesNames=out_axes, axesVals=out_axes_vals, vals=new_vals, weights=new_weights)
        new_soltab.addHistory('Added stations ' + str(added_stations).lstrip('[').rstrip(']') + ' with zero phases.')
        
    else:
        logging.error('There are fewer antennas in the solset than in the soltab ' + str(soltab_in))
        data.close()
        refdata.close()
        return(1)

    if len(bad_antennas_list) > 1:
        new_soltab.addHistory('Bad stations ' + str(bad_antennas_list[1:]).lstrip('[').rstrip(']') + ' have not been added back.')
    
    data.close()
    refdata.close()
    
    return(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Adds phases and amplitudes to any missing station if they appear in an h5parm, but not in a particular soltab.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('--refh5', type=str,
                        help='External H5parm from which the full list of antennas is used from.')
    parser.add_argument('--solset', type=str, default='sol000',
			help='Input calibration solutions')
    parser.add_argument('--refsolset', type=str, default='sol000',
			help='Input calibration solutions of the reference h5parm file')
    parser.add_argument('--soltab_in', type=str, default='phase000',
                        help='Input solution table')
    parser.add_argument('--soltab_out', type=str, default='GSMphase',
                        help='Output solution table (has to be different from input solution table)')


    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    main(h5parmfile=args.h5parm, refh5=args.refh5, solset=args.solset, refsolset=args.refsolset, soltab_in=args.soltab_in, soltab_out=args.soltab_out)

