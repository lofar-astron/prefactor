#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors:
# Francesco de Gasperin (H5parm_importer.py), 2014
# Andreas Horneffer (PipelineStep_losotoImporter.py), 2015

import os, sys
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

import numpy as np
import pyrap.tables as pt
import lofar.parmdb as pdb
from losoto.h5parm import h5parm, solWriter

import logging

# mandatory arguments:
# options: mapfile_in, mapfile-dir, filename 


# this is the code if you use it as a plugin
def plugin_main(args, **kwargs):
    fileid = kwargs['mapfile_in']
    datamap = DataMap.load(fileid)
    hdf5File = os.path.join(kwargs['hdf5_dir'],kwargs['hdf5file'])
    if kwargs.has_key('instrument'):
        instrument = kwargs['instrument']
    else:
        instrument = '/instrument'
    if kwargs.has_key('compression'):
        compression = int(kwargs['compression'])
    else:
        compression = 5
    if kwargs.has_key('solset'):
        solsetName = kwargs['solset']
    else:
        solsetName = None


    # Check is all the necessary files are available
    antennaFile = os.path.join(datamap[0].file,'ANTENNA')
    if not os.path.isdir(antennaFile):
        logging.critical('Missing ANTENNA table.')
        sys.exit(1)
    fieldFile = os.path.join(datamap[0].file,'FIELD')
    if not os.path.isdir(fieldFile):
        logging.critical('Missing FIELD table.')
        sys.exit(1)
    skydbFile = os.path.join(datamap[0].file,'sky')
    if not os.path.isdir(skydbFile):
        logging.warning('No sky table found. (Direction-dependent parameters will not work.)')
        skydbFile = None
        
    #generate list of parmDB-filenames
    parmDBnames = [ MS.file+instrument for MS in datamap ]

    #create and fill the hdf5-file:
    solset = parmDBs2h5parm(hdf5File, parmDBnames, antennaFile, fieldFile, skydbFile, compression=compression, solsetName=solsetName)

    # Add CREATE entry to history 
    h5parmDB = h5parm(hdf5File, readonly = False)
    soltabs = h5parmDB.getSoltabs(solset=solset)
    for st in soltabs:
        sw = solWriter(soltabs[st])
        sw.addHistory('CREATE (by PipelineStep_losotoImporter from %s / %s - %s)' % (os.path.abspath(''), 
                                   os.path.basename(parmDBnames[0]), os.path.basename(parmDBnames[-1]) ) )
    h5parmDB.close()

    #generate mapfile and wrap up
    mapfileentry = {}
    mapfileentry['host'] = 'localhost'
    mapfileentry['file'] = hdf5File
    mapfileentry['skip'] = False            
    outfileid = os.path.join(kwargs['mapfile_dir'], kwargs['filename'])
    outmap = open(outfileid, 'w')
    outmap.write(repr([mapfileentry]))
    outmap.close()
    result = {}
    result['mapfile'] = outfileid
    return result

# this is the code if you use it as a python-plugin step
def main(msfileslist, hdf5fileName, hdf5_dir='.', instrument='/instrument', solsetName=None, compression=5):
    tmp_msfiles = msfileslist.lstrip('[').rstrip(']').split(',')
    msfiles = [ MS.strip("\' \"") for MS in tmp_msfiles]
    hdf5File = os.path.join(hdf5_dir,hdf5fileName)
    compression = int(compression) #doesn't hurt if it is already an int
    instrument = instrument.strip()

    # Check is all the necessary files are available
    antennaFile = os.path.join(msfiles[0],'ANTENNA')
    if not os.path.isdir(antennaFile):
        logging.critical('Missing ANTENNA table.')
        sys.exit(1)
    fieldFile = os.path.join(msfiles[0],'FIELD')
    if not os.path.isdir(fieldFile):
        logging.critical('Missing FIELD table.')
        sys.exit(1)
    skydbFile = os.path.join(msfiles[0],'sky')
    if not os.path.isdir(skydbFile):
        logging.warning('No sky table found. (Direction-dependent parameters will not work.)')
        skydbFile = None
        
    #generate list of parmDB-filenames
    parmDBnames = [ MS+instrument for MS in msfiles ]

    #create and fill the hdf5-file:
    solset = parmDBs2h5parm(hdf5File, parmDBnames, antennaFile, fieldFile, skydbFile, compression=compression, solsetName=solsetName)

    # Add CREATE entry to history 
    h5parmDB = h5parm(hdf5File, readonly = False)
    soltabs = h5parmDB.getSoltabs(solset=solset)
    for st in soltabs:
        sw = solWriter(soltabs[st])
        sw.addHistory('CREATE (by losotoImporter from %s / %s - %s)' % (os.path.abspath(''), 
                                   os.path.basename(parmDBnames[0]), os.path.basename(parmDBnames[-1]) ) )
    h5parmDB.close()

    result = {}
    result['h5parm'] = hdf5File
    return result


def parmDBs2h5parm(h5parmName,parmDBs,antennaFile,fieldFile,skydbFile=None,compression=5,solsetName=None):
    """
    Write the contents of a list of parmDBs into a losoto-style hdf5-file
    h5parmName   - name (path) of the hdf5 file to generate
    parmDBs      - list with the names of the parmDBs
    antennaFile  - name (path) of an ANTENNA table of the observation
    fieldFile    - name (path) of a FIELD table of the observation
    skydbFile    - name (path) of a skydb table of the calibration run (Needed for direction dependent parameters)
    compresion   - compression level for the hdf5-file (0 - none ; 9 - highest)
    solsetName   - name of the solset to generate (default: "sol000")
    """

    # open/create the h5parm file and the solution-set
    h5parmDB = h5parm(h5parmName, readonly = False, complevel = compression)
    solset = h5parmDB.makeSolset(solsetName)

    #open the first instrument table, so that we know where to look for names and stuff
    firstInst = pdb.parmdb(parmDBs[0])

    # get unique list of solution types
    solTypes = list(set(x1.split(":")[0] for x1 in firstInst.getNames()))

    # rewrite solTypes in order to put together
    # Gain <-> DirectionalGain
    # CommonRotationAngle <-> RotationAngle
    # CommonScalarPhase <-> ScalarPhase
    # it also separate Real/Imag/Ampl/Phase into different solTypes
    if "Gain" in solTypes:
        solTypes.remove('Gain')
        solTypes.append('*Gain:*:Real')
        solTypes.append('*Gain:*:Imag')
        solTypes.append('*Gain:*:Ampl')
        solTypes.append('*Gain:*:Phase')
    if "DirectionalGain" in solTypes:
        solTypes.remove('DirectionalGain')
        solTypes.append('*Gain:*:Real')
        solTypes.append('*Gain:*:Imag')
        solTypes.append('*Gain:*:Ampl')
        solTypes.append('*Gain:*:Phase')
    if "RotationAngle" in solTypes:
        solTypes.remove('RotationAngle')
        solTypes.append('*RotationAngle')
    if "CommonRotationAngle" in solTypes:
        solTypes.remove('CommonRotationAngle')
        solTypes.append('*RotationAngle')
    if "RotationMeasure" in solTypes:
        solTypes.remove('RotationMeasure')
        solTypes.append('*RotationMeasure')
    if "ScalarPhase" in solTypes:
        solTypes.remove('ScalarPhase')
        solTypes.append('*ScalarPhase')
    if "CommonScalarPhase" in solTypes:
        solTypes.remove('CommonScalarPhase')
        solTypes.append('*ScalarPhase')
    if "CommonScalarAmplitude" in solTypes:
        solTypes.remove('CommonScalarAmplitude')
        solTypes.append('*ScalarAmplitude')
    # and remove duplicate entries
    solTypes = list(set(solTypes))

    for solType in solTypes:
        if len(firstInst.getNames(solType+':*')) == 0: continue
        pols = set(); dirs = set(); ants = set();
        freqs = set(); times = set(); ptype = set()

        for pDBname in parmDBs:
            instrumentdb = pdb.parmdb(pDBname)
            # create the axes grid, necessary if not all entries have the same axes lenght
            data = instrumentdb.getValuesGrid(solType+':*')
            for solEntry in data:
                pol, dir, ant, parm = parmdbToAxes(solEntry)
                if pol != None: pols |= set([pol])
                if dir != None: dirs |= set([dir])
                if ant != None: ants |= set([ant])
                freqs |= set(data[solEntry]['freqs'])
                times |= set(data[solEntry]['times'])
            #close the parmDB
            instrumentdb = 0 

        pols = np.sort(list(pols)); dirs = np.sort(list(dirs)); 
        ants = np.sort(list(ants)); freqs = np.sort(list(freqs)); 
        times = np.sort(list(times))
        shape = [i for i in (len(pols), len(dirs), len(ants), len(freqs), len(times)) if i != 0]
        vals = np.empty(shape)
        vals[:] = np.nan
        weights = np.zeros(shape)

        for pDBname in parmDBs:
            instrumentdb = pdb.parmdb(pDBname)
            # fill the values
            data = instrumentdb.getValuesGrid(solType+':*')
            if 'Real' in solType: dataIm = instrumentdb.getValuesGrid(solType.replace('Real','Imag')+':*')
            if 'Imag' in solType: dataRe = instrumentdb.getValuesGrid(solType.replace('Imag','Real')+':*')
            for solEntry in data:
                pol, dir, ant, parm = parmdbToAxes(solEntry)
                ptype |= set([solEntry.split(':')[0]]) # original parmdb solution type
                freq = data[solEntry]['freqs']
                time = data[solEntry]['times']
                val = data[solEntry]['values']
                # convert Real and Imag in Amp and Phase respectively
                if parm == 'Real':
                    solEntryIm = solEntry.replace('Real','Imag')
                    valI = dataIm[solEntryIm]['values']
                    val = np.sqrt((val**2)+(valI**2))
                if parm == 'Imag':
                    solEntryRe = solEntry.replace('Imag','Real')
                    valR = dataRe[solEntryRe]['values']
                    val = np.arctan2(val, valR)

                coords = []
                if pol != None:
                    polCoord = np.searchsorted(pols, pol)
                    coords.append(polCoord)
                if dir != None:
                    dirCoord = np.searchsorted(dirs, dir)
                    coords.append(dirCoord)
                if ant != None:
                    antCoord = np.searchsorted(ants, ant)
                    coords.append(antCoord)
                freqCoord = np.searchsorted(freqs, freq)
                timeCoord = np.searchsorted(times, time)
                vals[tuple(coords)][np.ix_(freqCoord,timeCoord)] = val.T
                weights[tuple(coords)][np.ix_(freqCoord,timeCoord)] = 1
            #close the parmDB
            instrumentdb = 0             

        vals = np.nan_to_num(vals) # replace nans with 0 (flagged later)

        if solType == '*RotationAngle':
            np.putmask(weights, vals == 0., 0) # flag where val=0
            h5parmDB.makeSoltab(solset, 'rotation', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*RotationMeasure':
            np.putmask(weights, vals == 0., 0) # flag where val=0
            h5parm.makeSoltab(solset, 'rotationmeasure', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*ScalarPhase':
            np.putmask(weights, vals == 0., 0)
            h5parmDB.makeSoltab(solset, 'scalarphase', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*ScalarAmplitude':
            np.putmask(weights, vals == 0., 0)
            h5parm.makeSoltab(solset, 'scalaramplitude', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == 'Clock':
            np.putmask(weights, vals == 0., 0)
            # clock may be diag or scalar
            if len(pols) == 0:
                h5parm.makeSoltab(solset, 'clock', axesNames=['ant','freq','time'], \
                    axesVals=[ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
            else:
                h5parm.makeSoltab(solset, 'clock', axesNames=['pol','ant','freq','time'], \
                    axesVals=[pol,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == 'TEC':
            np.putmask(weights, vals == 0., 0)
            # tec may be diag or scalar
            if len(pols) == 0:
                h5parm.makeSoltab(solset, 'tec', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
            else:
                h5parm.makeSoltab(solset, 'tec', axesNames=['pol','dir','ant','freq','time'], \
                    axesVals=[pols,dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*Gain:*:Real' or solType == '*Gain:*:Ampl':
            np.putmask(vals, vals == 0., 1) # nans end up into 1s (as BBS output, flagged next line)
            np.putmask(weights, vals == 1., 0) # flag where val=1
            h5parmDB.makeSoltab(solset, 'amplitude', axesNames=['pol','dir','ant','freq','time'], \
                    axesVals=[pols,dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*Gain:*:Imag' or solType == '*Gain:*:Phase':
            np.putmask(weights, vals == 0., 0) # falg where val=0
            h5parmDB.makeSoltab(solset, 'phase', axesNames=['pol','dir','ant','freq','time'], \
                    axesVals=[pols,dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))

    antennaTable = pt.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames,antennaPositions)))

    fieldTable = pt.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    fieldTable.close()

    sourceTable = solset._f_get_child('source')
    # add the field centre, that is also the direction for Gain and CommonRotationAngle
    sourceTable.append([('pointing',pointing)])

    dirs = []
    for tab in solset._v_children:
        c = solset._f_get_child(tab)
        if c._v_name != 'antenna' and c._v_name != 'source':
            dirs.extend(list(set(c.dir)))
    # remove duplicates
    dirs = list(set(dirs))
    # remove any pointing (already in the table)
    if 'pointing' in dirs:
        dirs.remove('pointing')

    if dirs != []:
        if skydbFile == None:
            logging.critical('No sky table given, but direction dependent parameters in parmDB. Exiting!')
            sys.exit(1)
        sourceFile = skydbFile + '/SOURCES'
        src_table = pt.table(sourceFile, ack=False)
        sub_tables = src_table.getsubtables()
        vals = []
        ra = dec = np.nan
        has_patches_subtable = False
        for sub_table in sub_tables:
            if 'PATCHES' in sub_table:
                has_patches_subtable = True
        if has_patches_subtable:
            # Read values from PATCHES subtable
            src_table.close()
            sourceFile = skydbFile + '/SOURCES/PATCHES'
            src_table = pt.table(sourceFile, ack=False)
            patch_names = src_table.getcol('PATCHNAME')
            patch_ras = src_table.getcol('RA')
            patch_decs = src_table.getcol('DEC')
            for source in dirs:
                try:
                    patch_indx = patch_names.index(source)
                    ra = patch_ras[patch_indx]
                    dec = patch_decs[patch_indx]
                except ValueError:
                    ra = np.nan
                    dec = np.nan
                vals.append([ra, dec])
            src_table.close()
        else:
            # Try to read default values from parmdb instead
            skydb = pdb.parmdb(skydbFile)
            vals = []
            ra = dec = np.nan

            for source in dirs:
                try:
                    ra = skydb.getDefValues('Ra:' + source)['Ra:' + source][0][0]
                    dec = skydb.getDefValues('Dec:' + source)['Dec:' + source][0][0]
                except KeyError:
                    # Source not found in skymodel parmdb, try to find components
                    ra = np.array(skydb.getDefValues('Ra:*' + source + '*').values())
                    dec = np.array(skydb.getDefValues('Dec:*' + source + '*').values())
                    if len(ra) == 0 or len(dec) == 0:
                        ra = np.nan
                        dec = np.nan
                    else:
                        ra = ra.mean()
                        dec = dec.mean()
                vals.append([ra, dec])
        sourceTable.append(zip(*(dirs,vals)))

    solsetname = solset._v_name
    # close the hdf5-file
    h5parmDB.close()
    return solsetname


def parmdbToAxes(solEntry):
    """
    Extract the information written as a string in the parmdb format
    """
    pol = None; pol1 = None; pol2 = None;
    dir = None; ant = None; parm = None

    thisSolType = solEntry.split(':')[0]

    # For CommonRotationAngle assuming [CommonRotationAngle:ant]
    if thisSolType == 'CommonRotationAngle':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For RotationAngle assuming [RotationAngle:ant:sou]
    elif thisSolType == 'RotationAngle':
        thisSolType, ant, dir = solEntry.split(':')

    # For RotationMeasure assuming [RotationMeasure:ant:sou]
    elif thisSolType == 'RotationMeasure':
        dir = 'pointing'
        try:
            thisSolType, ant = solEntry.split(':')
        except:
            thisSolType, ant, dir = solEntry.split(':')

    # For TEC assuming [TEC:ant or TEC:pol:ant]
    elif thisSolType == 'TEC':
        try:
            thisSolType, ant = solEntry.split(':')
        except:
            thisSolType, pol, ant = solEntry.split(':')
            pol1 = pol
            pol2 = pol
        dir = 'pointing'

    # For Clock assuming [Clock:ant or Clock:pol:ant]
    elif thisSolType == 'Clock':
        try:
            thisSolType, ant = solEntry.split(':')
        except:
            thisSolType, pol, ant = solEntry.split(':')
            pol1 = pol
            pol2 = pol
        dir = 'pointing'

    # For CommonScalarPhase assuming [CommonScalarPhase:ant]
    elif thisSolType == 'CommonScalarPhase':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For CommonScalarAmplitude assuming [CommonScalarAmplitude:ant]
    elif thisSolType == 'CommonScalarAmplitude':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For ScalarPhase assuming [ScalarPhase:ant:sou]
    elif thisSolType == 'ScalarPhase':
        thisSolType, ant, dir = solEntry.split(':')

    # For ScalarAmplitude assuming [ScalarAmplitude:ant:sou]
    elif thisSolType == 'ScalarAmplitude':
        thisSolType, ant, dir = solEntry.split(':')

    # For Gain assuming [Gain:pol1:pol2:parm:ant]
    elif thisSolType == 'Gain':
        thisSolType, pol1, pol2, parm, ant = solEntry.split(':')
        dir = 'pointing'

    # For DirectionalGain assuming [DirecitonalGain:pol1:pol2:parm:ant:sou]
    elif thisSolType == 'DirectionalGain':
        thisSolType, pol1, pol2, parm, ant, dir = solEntry.split(':')

    else:
        logging.error('Unknown solution type "'+thisSolType+'". Ignored.')

    if pol1 != None and pol2 != None:
        if pol1 == '0' and pol2 == '0': pol = 'XX'
        if pol1 == '1' and pol2 == '0': pol = 'YX'
        if pol1 == '0' and pol2 == '1': pol = 'XY'
        if pol1 == '1' and pol2 == '1': pol = 'YY'

    return pol, dir, ant, parm


if __name__=='__main__':
    # Options
    import optparse
    import glob

    opt = optparse.OptionParser(usage='%prog [-v] <H5parm> <MSPattern> \n'
                                '  <H5parm>    = (Path)name of the (new) H5parm file to be written.\n'
                                '  <MSPattern> = Search pattern for the measurement sets with instrument tables.\n'
                                '                (e.g. \"/data/scratch/MyObs/calibrator/L*.dppp\")\n'
                                '                Probably needs to be put in quotes when called from a shell!')
    opt.add_option('-i', '--instrument', dest="Instrument", type='string', default='/instrument',
                   help="Name of the instrument tables of the measurement sets. "
                   "If it starts with \"/\" -> instrument table is sub-directory within the MS directory. "
                   "(default=\"/instrument\")")
    opt.add_option('-v', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Solution-set name (default=sol###)', type='string', default=None)
    opt.add_option('-c', '--complevel', type='int', default='5',
                   help='Compression level from 0 (no compression, fast) to 9 (max compression, slow) (default=5)')

    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()

    # first argument: H5parm file name
    hdf5File = args[0]
    # second argument: pattern for measurement-sets
    inMSs = glob.glob(args[1])

    # options with default values
    instrument = options.Instrument
    compression = options.complevel
    solsetName = options.solset

    # Check is all the necessary files are available
    antennaFile = os.path.join(inMSs[0],'ANTENNA')
    if not os.path.isdir(antennaFile):
        logging.critical('Missing ANTENNA table.')
        sys.exit(1)
    fieldFile = os.path.join(inMSs[0],'FIELD')
    if not os.path.isdir(fieldFile):
        logging.critical('Missing FIELD table.')
        sys.exit(1)
    skydbFile = os.path.join(inMSs[0],'sky')
    if not os.path.isdir(skydbFile):
        logging.warning('No sky table found. (Direction-dependent parameters will not work.)')
        skydbFile = None
        
    #generate list of parmDB-filenames
    parmDBnames = [ MS.rstrip('/')+instrument for MS in inMSs ]

    #create and fill the hdf5-file:
    solset = parmDBs2h5parm(hdf5File, parmDBnames, antennaFile, fieldFile, skydbFile, compression=compression, solsetName=solsetName)

    # Add CREATE entry to history 
    h5parmDB = h5parm(hdf5File, readonly = False)
    soltabs = h5parmDB.getSoltabs(solset=solset)
    for st in soltabs:
        sw = solWriter(soltabs[st])
        sw.addHistory('CREATE (by PipelineStep_losotoImporter from %s / %s - %s)' % (os.path.abspath(''), 
                                   os.path.basename(parmDBnames[0]), os.path.basename(parmDBnames[-1]) ) )
    h5parmDB.close()

