#! /usr/bin/env python
"""
Script to make a sky model from fits model images
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits
from astropy import wcs
import numpy as np
import scipy.interpolate
import casacore.tables as pt
import sys
import os
import glob


def ra2hhmmss(deg):
    """Convert RA coordinate (in degrees) to HH MM SS"""

    from math import modf
    if deg < 0:
        deg += 360.0
    x, hh = modf(deg/15.)
    x, mm = modf(x*60)
    ss = x*60

    return (int(hh), int(mm), ss)


def dec2ddmmss(deg):
    """Convert DEC coordinate (in degrees) to DD MM SS"""

    from math import modf
    sign = (-1 if deg < 0 else 1)
    x, dd = modf(abs(deg))
    x, ma = modf(x*60)
    sa = x*60

    return (int(dd), int(ma), sa, sign)


def convert_radec_str(ra, dec):
    """Takes ra, dec in degrees and returns makesourcedb strings"""
    ra = ra2hhmmss(ra)
    sra = str(ra[0]).zfill(2)+':'+str(ra[1]).zfill(2)+':'+str("%.3f" % (ra[2])).zfill(6)
    dec = dec2ddmmss(dec)
    decsign = ('-' if dec[3] < 0 else '+')
    sdec = decsign+str(dec[0]).zfill(2)+'.'+str(dec[1]).zfill(2)+'.'+str("%.3f" % (dec[2])).zfill(6)

    return sra, sdec


def main(fits_models, ms_file, skymodel, fits_masks, min_flux_jy=0.005, interp='linear'):
    """
    Make a makesourcedb sky model for input MS from WSClean fits model images

    Parameters
    ----------
    fits_models : str
        Filename of FITS model images. Can be a list of files (e.g.,
        '[mod1.fits,mod2.fits,...]')
    ms_file : str
        Filename of MS for which sky model is to be made. Can be a list of files
        (e.g., '[ms1,ms2,...]', in which case they should all have the same
        frequency
    skymodel : str
        Filename of the output makesourcedb sky model
    fits_masks : str
        Filename of FITS mask images. Can be a list of files (e.g.,
        '[msk1.fits,msk2.fits,...]')
    min_flux_jy : float, optional
        Minimum value of flux in Jy of a source to include in output model
    interp : str, optional
        Interpolation method. Can be any supported by scipy.interpolate.interp1d:
            'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'

    """
    min_flux_jy = float(min_flux_jy)

    # Get filenames of model images and masks
    if '[' in fits_models and ']' in fits_models:
        fits_models = fits_models.strip('[] ').split(',')
        fits_models = [f.strip('\'\" ') for f in fits_models]
    else:
        fits_models = [fits_models.strip('\'\" ')]
    if '[' in fits_masks and ']' in fits_masks:
        fits_masks = fits_masks.strip('[]').split(',')
        fits_masks = [f.strip('\'\" ') for f in fits_masks]
    else:
        fits_masks = [fits_masks.strip('\'\" ')]

    # Read (first) MS file and get the frequency info
    if '[' in ms_file and ']' in ms_file:
        files = ms_file.strip('[] ').split(',')
        #files = [f.strip() for f in files]
        ms_file = files[0].strip('\'\" ')
    sw = pt.table(ms_file+'::SPECTRAL_WINDOW', ack=False)
    ms_freq = sw.col('REF_FREQUENCY')[0]
    ms_freq_low = sw.col('CHAN_FREQ')[0][0]
    ms_freq_high = sw.col('CHAN_FREQ')[0][-1]
    sw.close()

    # Get frequencies and data of model images and masks
    freqs = []
    model_images = []
    mask_images = []
    for f in fits_models:
        hdr = fits.getheader(f, 0)
        freqs.append(hdr['CRVAL3']) # Hz
        model_images.append(fits.getdata(f, 0))
    for f in fits_masks:
        mask_images.append(fits.getdata(f, 0))

    # Sort by freq
    sorted_ind = np.argsort(freqs)
    freqs = np.array(freqs)[sorted_ind]
    fits_models = np.array(fits_models)[sorted_ind]
    model_images = np.array(model_images)[sorted_ind]
    mask_images = np.array(mask_images)[sorted_ind]

    # Check if there is a model at the ms frequency. If so, just use that one
    ind = np.where( np.logical_and(freqs >= ms_freq_low, freqs <= ms_freq_high) )
    if len(ind[0]) == 1:
        freqs = freqs[ind]
        fits_models = fits_models[ind]
        model_images = model_images[ind]
        mask_images = mask_images[ind]

    # Set the WCS reference
    hdr = fits.getheader(fits_models[0], 0)
    w = wcs.WCS(hdr)

    # Find nonzero pixels in stacked image
    stacked_model = np.zeros(model_images[0].shape)
    stacked_mask = np.zeros(mask_images[0].shape)
    for im, ma in zip(model_images, mask_images):
        stacked_model += im
        stacked_mask += ma
    nonzero_ind = np.where((stacked_model != 0.0) & (stacked_mask > 0))

    # Interpolate the fluxes to the frequency of the MS
    nsources = len(nonzero_ind[0])
    fluxes = []
    names = []
    ras = []
    decs = []
    for i in range(nsources):
        flux_array = []
        freq_array = []
        for ind, (im, ma) in enumerate(zip(model_images, mask_images)):
            index = [nonzero_ind[j][i] for j in range(4)]
            if ma[tuple(index)] > 0:
                flux_array.append(im[tuple(index)])
                freq_array.append(freqs[ind])

        if len(flux_array) > 0:
            # Only create a source entry if MS frequency is within +/- 4 MHz of
            # the sampled frequency range to prevent large extrapolations when
            # the primary beam may cause rapid spectral changes
            if ms_freq > freq_array[0]-4e6 and ms_freq < freq_array[-1]+4e6:
                # If MS frequency lies outside range, just use nearest freq
                if ms_freq < freq_array[0]:
                    flux = flux_array[0]
                elif ms_freq > freq_array[-1]:
                    flux = flux_array[-1]
                else:
                    # Otherwise interpolate
                    flux = scipy.interpolate.interp1d(freq_array, flux_array, kind=interp)(ms_freq)
                if flux > min_flux_jy:
                    index.reverse() # change to WCS coords
                    ras.append(w.wcs_pix2world(np.array([index]), 0, ra_dec_order=True)[0][0])
                    decs.append(w.wcs_pix2world(np.array([index]), 0, ra_dec_order=True)[0][1])
                    names.append('cc{}'.format(i))
                    fluxes.append(flux)

    # Write sky model
    with open(skymodel, 'w') as outfile:
        outfile.write('FORMAT = Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency\n')
        for name, ra, dec, flux in zip(names, ras, decs, fluxes):
            ra_str, dec_str = convert_radec_str(ra, dec)
            outfile.write('{0}, POINT, {1}, {2}, {3}, 0.0, 0.0, 0.0, {4}\n'
                .format(name, ra_str, dec_str, flux, ms_freq))


if __name__ == '__main__':
    descriptiontext = "Make a makesourcedb sky model from WSClean fits model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fits_models', help='Model images')
    parser.add_argument('msfile', help='MS file')
    parser.add_argument('skymodel', help='Filename of output sky model')
    parser.add_argument('fits_masks', help='Mask images')

    args = parser.parse_args()
    main(args.fits_models, args.msfile, args.skymodel, args.fits_masks)
