#!/usr/bin/env python
"""
Script to merge selfcal parmdbs
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import numpy as np
import lofar.parmdb as pdb

# Make an empty parmDB with only the defaults and return the parmdb object
def make_empty_parmdb(outname):
    myParmdb=pdb.parmdb(outname,create=True)
    myParmdb.addDefValues("Gain:0:0:Ampl",1.)
    myParmdb.addDefValues("Gain:1:1:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Ampl",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Ampl",1.)
    myParmdb.addDefValues("Gain:0:0:Real",1.)
    myParmdb.addDefValues("Gain:1:1:Real",1.)
    myParmdb.addDefValues("DirectionalGain:0:0:Real",1.)
    myParmdb.addDefValues("DirectionalGain:1:1:Real",1.)
    myParmdb.addDefValues("AntennaOrientation",5.497787144)
    myParmdb.addDefValues("RotationMeasure",1e-6)    
    return myParmdb


def main(parmdb_in, parmdb_out ,debug=False):
    """
    Merges facet selfcal parmdbs into a parmdb for a single band

    Parameters
    ----------
    parmdb_in : str
        Name of (path to) the parmdb with Gain in polar coordinates.
    parmdb_out : str
        Name of the new parmdb that will be written
    """

    pdb_in =  pdb.parmdb(parmdb_in,create=False)
    pdb_out = make_empty_parmdb(parmdb_out)
    inparms = pdb_in.getValuesGrid('*')
    names_in = inparms.keys()
    for parmname in names_in:
        nameparts = parmname.split(':')
        if (len(nameparts) == 5 
            and (nameparts[0] == 'Gain' or nameparts[0] == 'DirectionalGain')): 
            if nameparts[3] == 'Phase' :
                amplname = nameparts[0]+':'+nameparts[1]+':'+nameparts[2]+':Ampl:'+nameparts[4]
                if amplname in names_in:
                    amps = inparms[amplname]['values']
                else:
                    amps = np.ones_like(inparms[parmname]['values'])
                realname = nameparts[0]+':'+nameparts[1]+':'+nameparts[2]+':Real:'+nameparts[4]
                imagname = nameparts[0]+':'+nameparts[1]+':'+nameparts[2]+':Imag:'+nameparts[4]
                realvals = amps * np.cos(inparms[parmname]['values'])
                imagvals = amps * np.sin(inparms[parmname]['values'])
                if debug: print 'Adding parameters:',realname,imagname
                ValueHolder = pdb_out.makeValue(values=realvals,
                                                sfreq=inparms[parmname]['freqs'],
                                                efreq=inparms[parmname]['freqwidths'],
                                                stime=inparms[parmname]['times'],
                                                etime=inparms[parmname]['timewidths'],
                                                asStartEnd=False)
                pdb_out.addValues(realname,ValueHolder)  
                ValueHolder = pdb_out.makeValue(values=imagvals,
                                                sfreq=inparms[parmname]['freqs'],
                                                efreq=inparms[parmname]['freqwidths'],
                                                stime=inparms[parmname]['times'],
                                                etime=inparms[parmname]['timewidths'],
                                                asStartEnd=False)
                pdb_out.addValues(imagname,ValueHolder)  
            elif nameparts[3] == 'Real' or nameparts[3] == 'Imag' :
                print 'Gain value:',parmname,'not in polar coordinates!'
        else:
            ValueHolder = pdb_out.makeValue(values=inparms[parmname]['values'],
                                            sfreq=inparms[parmname]['freqs'],
                                            efreq=inparms[parmname]['freqwidths'],
                                            stime=inparms[parmname]['times'],
                                            etime=inparms[parmname]['timewidths'],
                                            asStartEnd=False)
            pdb_out.addValues(parmname,ValueHolder)  

    pdb_out.flush()

if __name__ == '__main__':
    descriptiontext = "Merge parmdbs in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('parmdb_in', help='Name of (path to) the parmdb with Gain in polar coordinates.')
    parser.add_argument('parmdb_out', help='Name of the new parmdb that will be written')
    parser.add_argument('-d', '--debug', help='Generate dubug output', action='store_true', default=False)

    args = parser.parse_args()
    main(args.parmdb_in, args.parmdb_out,args.debug)
