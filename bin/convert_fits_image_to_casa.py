#!/usr/bin/env python
"""
Convert a FITS image to a CASA image

For WSClean 1.7, use force_stokes_I = True to overide incorrect Stokes
keyword in model image headers. The restfreq is also set to avoid
problems with casapy2bbs.py
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.images as pim
import pyrap.tables as pt
import numpy as np


def main(fitsimage, outfilename, force_stokes_i=False):
    """
    Convert a fits image to a CASA image

    Parameters
    ----------
    fitsimage : str
        Name of FITS image
    outfilename : str
        Name of output CASA image
    force_stokes_i : bool, optional
        If True, force Stokes axis to be 'I'

    """
    casaimage = pim.image(fitsimage)
    casaimage.saveas(outfilename, overwrite=True)

    if type(force_stokes_i) is str:
        if force_stokes_i.lower() == 'true':
            force_stokes_i = True
        else:
            force_stokes_i = False

    if force_stokes_i:
        coords = casaimage.coordinates().dict()
        coords['stokes1']['stokes'] = ['I']
        freq = coords['spectral2']['wcs']['crval']
        coords['spectral2']['restfreqs'] = np.array([freq])
        outtable = pt.table(outfilename, readonly=False, ack=False)
        outtable.putkeywords({'coords': coords})
        outtable.done()


if __name__ == '__main__':
    descriptiontext = "Convert a fits image to a CASA image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fitsimage', help='name of the full skymodel')
    parser.add_argument('outfilename', help='name for the output')
    parser.add_argument('-f', '--force_stokes_i', help='Force Stokes I (for WSClean v1.7)',
        type=bool, default=False)
    args = parser.parse_args()

    main(args.fitsimage, args.outfilename, args.force_stokes_i)
