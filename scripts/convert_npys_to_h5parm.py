#!/usr/bin/env python
import numpy as np
import pyrap.tables as pt
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import logging
from os.path import join

def makesolset(MS, data, solset_name):
    solset = data.makeSolset(solset_name)    

    antennaFile = MS + "/ANTENNA"
    logging.info('Collecting information from the ANTENNA table.')
    antennaTable = pt.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset.obj._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames,antennaPositions)))
    
    fieldFile = MS + "/FIELD"
    logging.info('Collecting information from the FIELD table.')
    fieldTable = pt.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    fieldTable.close()

    sourceTable = solset.obj._f_get_child('source')
    # add the field centre, that is also the direction for Gain and Common*
    sourceTable.append([('pointing',pointing)])

    return antennaNames
    pass

########################################################################
def main(MSfiles, h5parmfile, store_basename='caldata_transfer', store_directory='.', solset_out = 'calibrator'):

    # select on MS
    mslist = MSfiles.lstrip('[').rstrip(']').replace(' ','').replace("'","").split(',')
    
    if len(mslist) == 0:
        logging.error("Did not find any existing directory in input MS list!")
        return(1)
        pass
    else:
        MSfile = mslist[0]
        pass

    # Open up the h5parm, get an example value
    data   = h5parm(h5parmfile, readonly = False)
    
    # name (path) for parmdb to be written
    reference_station_names = makesolset(MSfile, data, solset_out)
    OutSolset = data.getSolset(solset_out)

    # load the numpy arrays written by the previous scripts
    # (filenames constructed in the same way as in these scripts)
    freqs_ampl    = np.load(join(store_directory, 'freqs_for_amplitude_array.npy'))
    amps_array    = np.load(join(store_directory, store_basename + '_amplitude_array.npy'))
    clock_array   = np.load(join(store_directory, 'fitted_data_dclock_' + store_basename + '_1st.npy'))
    tec_array     = np.load(join(store_directory, 'fitted_data_dTEC_' + store_basename + '_1st.npy'))
    freqs_phase   = np.load(join(store_directory, 'freqs_for_phase_array.npy'))
    phases_array  = np.load(join(store_directory, store_basename + '_phase_array.npy'))
    station_names = np.load(join(store_directory, store_basename + '_station_names.npy'))
    
    if len(reference_station_names) > len(station_names):
        missing_stations = set(reference_station_names) - set(station_names)
        logging.warning('CAUTION: For these stations now calibration solutions are available: ' + str(missing_stations))
        pass
    
    # create weights
    dimensions_amp   = np.shape(amps_array)
    dimensions_phase = np.shape(phases_array)
    dimensions_clock = np.shape(clock_array)
    dimensions_tec   = np.shape(tec_array)
    weights_amp      = np.ndarray(shape = dimensions_amp)
    weights_phase    = np.ndarray(shape = dimensions_phase)
    weights_clock    = np.ndarray(shape = dimensions_clock)
    weights_tec      = np.ndarray(shape = dimensions_tec)
    weights_amp.fill(1)
    weights_phase.fill(1)
    weights_clock.fill(1)
    weights_tec.fill(1)
    
    OutSolset.makeSoltab(soltype='clock', soltabName='clock',
                                        axesNames=['ant'],
                                        axesVals=[station_names],
                                        vals=clock_array[-1,:], weights=weights_clock[-1,:])
    OutSolset.makeSoltab(soltype='tec', soltabName='tec',
                                        axesNames=['ant'],
                                        axesVals=[station_names],
                                        vals=tec_array[-1,:], weights=weights_tec[-1,:])
    OutSolset.makeSoltab(soltype='amplitude', soltabName='amplitude',
                                        axesNames=['ant', 'freq', 'pol'],
                                        axesVals=[station_names, freqs_ampl, ['XX','YY']],
                                        vals=amps_array[:,-1,:,:], weights=weights_amp[:,-1,:,:])
    OutSolset.makeSoltab(soltype='phase', soltabName='phase',
                                        axesNames=['freq', 'ant'],
                                        axesVals=[freqs_phase, station_names],
                                        vals=phases_array, weights=weights_phase)

    return(0)


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create a parmDB with values from the calibrator.')

    parser.add_argument('MSfiles', type=str, 
                        help='One or more MSs to use as a reference.')
    parser.add_argument('h5parmfile', type=str,
                        help='Output H5parm file.')
    parser.add_argument('--basename', type=str, default='caldata_transfer',
                        help='Base-name of the numpy files with the calibrator values. (default: \"caldata_transfer\")')
    parser.add_argument('--storedir', type=str, default='.',
                        help='Directory of the numpy files with the calibrator values. (default: \".\")')
    parser.add_argument('--solset_out', type=str, default='calibrator',
                        help='Name of the solution set')

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    args = parser.parse_args()
    
    MSfiles = args.MSfiles
    logging.error(MSfiles)
    h5parmfile = args.h5parmfile

    main(MSfiles, h5parmfile, store_basename=args.basename, store_directory=args.storedir, solset_out=args.solset_out)
