#!/usr/bin/env python
"""
plot image overlaid with mask contours
Original code from Reinout
"""
import argparse
from argparse import RawTextHelpFormatter
#import pyrap.images as pim
#import pyrap.tables as pt
#import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import aplpy
import numpy
import astropy.io.fits

def meanclip(indata, clipsig=4.0, maxiter=10, converge_num=0.001, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.

   .. note:: MYMEANCLIP routine from ACS library.

   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.

   Examples
   --------
   >>> mean, sigma = meanclip(indata)

   Parameters
   ----------
   indata: array_like
       Input data.

   clipsig: float
       Number of sigma at which to clip.

   maxiter: int
       Ceiling on number of clipping iterations.

   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.

   verbose: {0, 1}
       Print messages to screen?

   Returns
   -------
   mean: float
       N-sigma clipped mean.

   sigma: float
       Standard deviation of remaining pixels.

   """
   # Flatten array
   skpix = indata.reshape( indata.size, )
 
   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0
 
   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = numpy.median(skpix)
       sig = numpy.std(skpix)
       wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop
 
   mean  = numpy.mean( skpix )
   sigma = robust_sigma( skpix )
 
   if verbose:
       prf = 'MEANCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)
 
   return mean, sigma


def find_imagenoise(data):
  mean, rms =  meanclip(data)
  return rms

def robust_sigma(in_y, zero=0):
    """
    Calculate a resistant estimate of the dispersion of
    a distribution. For an uncontaminated distribution,
    this is identical to the standard deviation.

    Use the median absolute deviation as the initial
    estimate, then weight points using Tukey Biweight.
    See, for example, Understanding Robust and
    Exploratory Data Analysis, by Hoaglin, Mosteller
    and Tukey, John Wiley and Sons, 1983.

    .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.

    Examples
    --------
    >>> result = robust_sigma(in_y, zero=1)

    Parameters
    ----------
    in_y : array_like
        Vector of quantity for which the dispersion is
        to be calculated

    zero : int
        If set, the dispersion is calculated w.r.t. 0.0
        rather than the central value of the vector. If
        Y is a vector of residuals, this should be set.

    Returns
    -------
    out_val : float
        Dispersion value. If failed, returns -1.

    """
    # Flatten array
    y = in_y.ravel()

    eps = 1.0E-20
    c1 = 0.6745
    c2 = 0.80
    c3 = 6.0
    c4 = 5.0
    c_err = -1.0
    min_points = 3

    if zero:
        y0 = 0.0
    else:
        y0 = numpy.median(y)

    dy    = y - y0
    del_y = abs( dy )

    # First, the median absolute deviation MAD about the median:

    mad = numpy.median( del_y ) / c1

    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps:
        mad = del_y.mean() / c2
    if mad < eps:
        return 0.0

    # Now the biweighted value:
    u  = dy / (c3 * mad)
    uu = u * u
    q  = numpy.where(uu <= 1.0)
    count = len(q[0])
    if count < min_points:
        module_logger.warn('ROBUST_SIGMA: This distribution is TOO WEIRD! '
                           'Returning {}'.format(c_err))
        return c_err

    numerator = numpy.sum( (y[q] - y0)**2.0 * (1.0 - uu[q])**4.0 )
    n    = y.size
    den1 = numpy.sum( (1.0 - uu[q]) * (1.0 - c4 * uu[q]) )
    siggma = n * numerator / ( den1 * (den1 - 1.0) )

    if siggma > 0:
        out_val = numpy.sqrt( siggma )
    else:
        out_val = 0.0

    return out_val




def main(fitsimage, maskfits, outfilename):

    #imagenoise = 5E-3
    
    #f = aplpy.FITSFigure(fitsimage,slices=[0,0])
    #f.show_colorscale(vmax=10*imagenoise,vmin=-3*imagenoise)
    ##outplotname = fitsimage.replace('.fits','.png')
    #f.add_colorbar()
    #f.save('%s.png'%(outfilename),dpi=750)
    #f.show_contour(maskfits, levels=[0.5], colors='m')
    #f.save('%s_mask.png'%(outfilename),dpi=750)
    
    ##f = aplpy.FITSFigure(maskfits,slices=[0,0])
    ##f.show_colorscale(vmax=0,vmin=1)
    ###f.show_contour(maskfits, levels=[0.5], colors='w')
    ###outplotname = fitsimage.replace('.fits','.png')
    ##f.add_colorbar()
    ##f.save('%s-mask.png'%(outfilename),dpi=500)



    #outplotname = fitsimage.replace('.fits','.png')

    # find image noise
    hdulist = astropy.io.fits.open(fitsimage)
    data    = hdulist[0].data
    imagenoise = find_imagenoise(data)
    hdulist.close()
    print 'Image noise is: ', imagenoise*1e3, ' mJy'


    f = aplpy.FITSFigure(fitsimage,slices=[0,0])
    f.show_colorscale(vmax=16*imagenoise,vmin=-4*imagenoise, cmap='bone')


    #f.add_beam()
    #f.beam.set_color('white')
    #f.beam.set_hatch('.')


    f.add_grid()
    f.grid.set_color('white')
    f.grid.set_alpha(0.5)
    f.grid.set_linewidth(0.2)

    f.add_colorbar()


    f.show_contour(maskfits, colors='red', levels=[0.0],filled=False, smooth=1,alpha=0.6, linewidths=0.25)
    f.colorbar.set_axis_label_text('Flux (Jy beam$^{-1}$)')

    f.save('%s.png'%(outfilename),dpi=400)




if __name__ == '__main__':
    descriptiontext = "Convert a fits image to a CASA image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fitsimage', help='name of the image')
    parser.add_argument('maskfits', help='name of the mask')
    parser.add_argument('outname', help='name for the output')
    args = parser.parse_args()

    main(args.fitsimage, args.maskfits, args.outname)
