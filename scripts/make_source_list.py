#! /usr/bin/env python
"""
Script to make a source catalog for an image
"""
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
try:
    import matplotlib
    matplotlib.use('Agg')
except (RuntimeError, ImportError):
    pass
try:
    import bdsf
except ImportError:
    from lofar import bdsm as bdsf
import lsmtool


def main(image_name, catalog_name, atrous_do=False, threshisl=3.0, threshpix=5.0, rmsbox=(60, 20),
         rmsbox_bright=(35, 7), adaptive_rmsbox=False, atrous_jmax=6, adaptive_thresh=150.0,
         compare_dir=None, format='png'):
    """
    Make a source catalog for an image

    Parameters
    ----------
    image_name : str
        Filename of input image from which catalog will be made
    catalog_name : str
        Filename of output source catalog in FITS format
    atrous_do : bool, optional
        Use wavelet module of PyBDSF?
    threshisl : float, optional
        Value of thresh_isl PyBDSF parameter
    threshpix : float, optional
        Value of thresh_pix PyBDSF parameter
    rmsbox : tuple of floats, optional
        Value of rms_box PyBDSF parameter
    rmsbox_bright : tuple of floats, optional
        Value of rms_box_bright PyBDSF parameter
    adaptive_rmsbox : tuple of floats, optional
        Value of adaptive_rms_box PyBDSF parameter
    atrous_jmax : int, optional
        Value of atrous_jmax PyBDSF parameter
    adaptive_thresh : float, optional
        If adaptive_rmsbox is True, this value sets the threshold above
        which a source will use the small rms box
    compare_dir : str, optional
        Directory for output plots/info of comparison of source properties to surveys
    format : str, optional
        Format of output plots: 'png' or 'pdf'
    """
    if type(rmsbox) is str:
        rmsbox = eval(rmsbox)

    if type(rmsbox_bright) is str:
        rmsbox_bright = eval(rmsbox_bright)

    if type(atrous_do) is str:
        if atrous_do.lower() == 'true':
            atrous_do = True
            threshisl = 4.0  # override user setting to ensure proper source fitting
        else:
            atrous_do = False

    if type(adaptive_rmsbox) is str:
        if adaptive_rmsbox.lower() == 'true':
            adaptive_rmsbox = True
        else:
            adaptive_rmsbox = False

    if not os.path.exists(image_name):
        print('ERROR: input image not found')
        sys.exit(1)

    atrous_jmax = int(atrous_jmax)
    threshpix = float(threshpix)
    threshisl = float(threshisl)
    adaptive_thresh = float(adaptive_thresh)

    img = bdsf.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                             thresh_pix=threshpix, thresh_isl=threshisl,
                             atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                             adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=adaptive_thresh,
                             rms_box_bright=rmsbox_bright, rms_map=True, quiet=True,
                             atrous_jmax=atrous_jmax)
    srcroot = 'prefactor'
    img.write_catalog(outfile=catalog_name, bbs_patches='source', srcroot=srcroot, clobber=True)

    # Use LSMTool to make some basic comparisons to surveys
    if compare_dir is not None:
        s = lsmtool.load(catalog_name)
        s.setPatchPositions(method='wmean')
        _, _, refRA, refDec = s._getXY()
        def_dict = s.getDefaultValues()
        if 'ReferenceFrequency' in def_dict:
            ref_freq_hz = def_dict['ReferenceFrequency']
        else:
            ref_freq_hz = 0.0
        if ref_freq_hz > 140e6 and ref_freq_hz < 160e6:
            # use TGSS
            vocat = 'TGSS'
            ignspec = None
        else:
            # too far in freq from TGSS -> use GSM, ignoring sources that lack spectral
            # information (spec_indx = -0.7)
            vocat = 'GSM'
            ignspec = -0.7
        s_vo = lsmtool.load(vocat, VOPosition=[refRA, refDec], VORadius='5 deg')
        s_vo.group('threshold', FWHM='30 arcsec', threshold=0.001)
        s.compare(s_vo, radius='30 arcsec', excludeMultiple=False, outDir=compare_dir,
                  name1='LOFAR', name2=vocat, format=format, ignoreSpec=ignspec)


if __name__ == '__main__':
    descriptiontext = "Make a source catalog from an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_name', help='Image name')
    parser.add_argument('catalog_name', help='Output catalog name')
    parser.add_argument('-a', '--atrous_do', help='use wavelet fitting', type=bool, default=False)
    parser.add_argument('-i', '--threshisl', help='', type=float, default=3.0)
    parser.add_argument('-p', '--threshpix', help='', type=float, default=5.0)
    parser.add_argument('-r', '--rmsbox', help='rms box width and step (e.g., "(60, 20)")',
                        type=str, default='(60, 20)')
    parser.add_argument('--rmsbox_bright', help='rms box for bright sources(?) width and step (e.g., "(60, 20)")',
                        type=str, default='(60, 20)')
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)
    parser.add_argument('-j', '--atrous_jmax', help='Max wavelet scale', type=int, default=3)
    parser.add_argument('-c', '--compare_dir', help='', type=str, default=None)
    parser.add_argument('-f', '--format', help='', type=str, default='pdf')

    args = parser.parse_args()
    main(args.image_name, args.catalog_name, atrous_do=args.atrous_do,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
         rmsbox_bright=args.rmsbox_bright, adaptive_rmsbox=args.adaptive_rmsbox,
         atrous_jmax=args.atrous_jmax, compare_dir=args.compare_dir, format=args.format)
