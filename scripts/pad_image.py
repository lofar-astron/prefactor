#! /usr/bin/env python
"""
Script to pad a FITS image
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits as pyfits
import numpy as np
import sys
import os


def main(infile, xypadsize):
    """Pad a fits image with zeros to padsize"""
    pad_xsize, pad_ysize = [int(s) for s in xypadsize.split(' ')]

    hdu = pyfits.open(infile)
    imdata = hdu[0].data[0, 0]
    (ysize, xsize) = imdata.shape
    if xsize > pad_xsize or ysize > pad_ysize:
        raise ValueError('pad_image: padded size is smaller than current size!')
    if xsize == pad_xsize and ysize == pad_ysize:
        return()

    xoffset = (pad_xsize - xsize) / 2
    yoffset = (pad_ysize - ysize) / 2
    newdata=np.zeros((1, 1, pad_ysize, pad_xsize))

    newdata[0, 0, yoffset:yoffset+ysize, xoffset:xoffset+xsize] = imdata
    hdu[0].data = newdata
    hdu[0].header['CRPIX1'] += xoffset
    hdu[0].header['CRPIX2'] += yoffset
    hdu.writeto(infile, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Pad FITS image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('infile', help='filename of input image')
    parser.add_argument('xypadsize', help='padded size as "xsize ysize"')
    args = parser.parse_args()

    main(args.infile, args.xypadsize)
