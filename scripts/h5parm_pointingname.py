#!/usr/bin/env python

from losoto.h5parm import h5parm
from losoto.lib_operations import *

import numpy
import logging
########################################################################
def main(h5parmdb, solsetName='sol000', pointing='POINTING'):


    logging.info("Setting " + str(pointing) + " as pointing direction for the solset " + str(solsetName) + " in " + str(h5parmdb))
    
    data          = h5parm(h5parmdb, readonly = False)       
    solset        = data.getSolset(solsetName)
    
    sources = solset.obj._f_get_child('source')
    direction = list(sources[0][-1])
    
    for i in numpy.arange(len(sources)):
        sources[i] = (pointing, direction)
    
    data.close()
    return(0)

    
########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Define the pointing direction of an h5parm.')

    parser.add_argument('h5parmdb', type=str,
                        help='H5parm of which pointing will be set.')
    parser.add_argument('--solsetName', '--solset', type=str, default='sol000',
                        help='Name of the h5parm solution set (default: sol000)')
    parser.add_argument('--pointing', type=str, default='POINTING',
                        help='Name of the h5parm solution set (default: POINTING)')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)
    
    h5parmdb = args.h5parmdb
    main(h5parmdb, solsetName=args.solsetName, pointing=args.pointing)
