#! /usr/bin/env python
"""
Script to make a clean mask
"""
import argparse
from argparse import RawTextHelpFormatter
from lofar import bdsm
import casacore.images as pim
from astropy.io import fits as pyfits
from astropy.coordinates import Angle
import pickle
import numpy as np
import sys
import os
#from factor.lib.polygon import Polygon

class Polygon:
    """
    Generic polygon class

    Polygons are used to define the facet boundaries.

    Parameters
    ----------
    x : array
        A sequence of nodal x-coords.
    y : array
        A sequence of nodal y-coords.

    """
    def __init__(self, x, y):
        if len(x) != len(y):
            raise IndexError('x and y must be equally sized.')
        self.x = np.asfarray(x)
        self.y = np.asfarray(y)

        # Closes the polygon if were open
        x1, y1 = x[0], y[0]
        xn, yn = x[-1], y[-1]
        if x1 != xn or y1 != yn:
            self.x = np.concatenate((self.x, [x1]))
            self.y = np.concatenate((self.y, [y1]))

        # Anti-clockwise coordinates
        if _det(self.x, self.y) < 0:
            self.x = self.x[::-1]
            self.y = self.y[::-1]


    def is_inside(self, xpoint, ypoint, smalld=1e-12):
        """
        Check if point is inside a general polygon.

        An improved version of the algorithm of Nordbeck and Rydstedt.

        REF: SLOAN, S.W. (1985): A point-in-polygon program. Adv. Eng.
        Software, Vol 7, No. 1, pp 45-47.

        Parameters
        ----------
        xpoint : array or float
            The x-coord of the point to be tested.
        ypoint : array or float
            The y-coords of the point to be tested.
        smalld : float
            Tolerance within which point is considered to be on a side.

        Returns
        -------
        mindst : array or float
            The distance from the point to the nearest point of the polygon:

            If mindst < 0 then point is outside the polygon.
            If mindst = 0 then point in on a side of the polygon.
            If mindst > 0 then point is inside the polygon.

        """
        xpoint = np.asfarray(xpoint)
        ypoint = np.asfarray(ypoint)

        # Scalar to array
        if xpoint.shape is tuple():
            xpoint = np.array([xpoint], dtype=float)
            ypoint = np.array([ypoint], dtype=float)
            scalar = True
        else:
            scalar = False
        # Check consistency
        if xpoint.shape != ypoint.shape:
            raise IndexError('x and y  must be equally sized.')

        # If snear = True: Dist to nearest side < nearest vertex
        # If snear = False: Dist to nearest vertex < nearest side
        snear = np.ma.masked_all(xpoint.shape, dtype=bool)

        # Initialize arrays
        mindst = np.ones_like(xpoint, dtype=float) * np.inf
        j = np.ma.masked_all(xpoint.shape, dtype=int)
        x = self.x
        y = self.y
        n = len(x) - 1  # Number of sides/vertices defining the polygon

        # Loop over each side defining polygon
        for i in range(n):
            d = np.ones_like(xpoint, dtype=float) * np.inf

            # Start of side has coords (x1, y1)
            # End of side has coords (x2, y2)
            # Point has coords (xpoint, ypoint)
            x1 = x[i]
            y1 = y[i]
            x21 = x[i + 1] - x1
            y21 = y[i + 1] - y1
            x1p = x1 - xpoint
            y1p = y1 - ypoint

            # Points on infinite line defined by
            #     x = x1 + t * (x1 - x2)
            #     y = y1 + t * (y1 - y2)
            # where
            #     t = 0    at (x1, y1)
            #     t = 1    at (x2, y2)
            # Find where normal passing through (xpoint, ypoint) intersects
            # infinite line
            t = -(x1p * x21 + y1p * y21) / (x21 ** 2 + y21 ** 2)
            tlt0 = t < 0
            tle1 = (0 <= t) & (t <= 1)

            # Normal intersects side
            d[tle1] = ((x1p[tle1] + t[tle1] * x21) ** 2 +
                       (y1p[tle1] + t[tle1] * y21) ** 2)

            # Normal does not intersects side
            # Point is closest to vertex (x1, y1)
            # Compute square of distance to this vertex
            d[tlt0] = x1p[tlt0] ** 2 + y1p[tlt0] ** 2

            # Store distances
            mask = d < mindst
            mindst[mask] = d[mask]
            j[mask] = i

            # Point is closer to (x1, y1) than any other vertex or side
            snear[mask & tlt0] = False

            # Point is closer to this side than to any other side or vertex
            snear[mask & tle1] = True

        if np.ma.count(snear) != snear.size:
            raise IndexError('Error computing distances')
        mindst **= 0.5

        # Point is closer to its nearest vertex than its nearest side, check if
        # nearest vertex is concave.
        # If the nearest vertex is concave then point is inside the polygon,
        # else the point is outside the polygon.
        jo = j.copy()
        jo[j == 0] -= 1
        area = _det([x[j + 1], x[j], x[jo - 1]], [y[j + 1], y[j], y[jo - 1]])
        mindst[~snear] = np.copysign(mindst, area)[~snear]

        # Point is closer to its nearest side than to its nearest vertex, check
        # if point is to left or right of this side.
        # If point is to left of side it is inside polygon, else point is
        # outside polygon.
        area = _det([x[j], x[j + 1], xpoint], [y[j], y[j + 1], ypoint])
        mindst[snear] = np.copysign(mindst, area)[snear]

        # Point is on side of polygon
        mindst[np.fabs(mindst) < smalld] = 0

        # If input values were scalar then the output should be too
        if scalar:
            mindst = float(mindst)
        return mindst



def read_vertices(filename):
    """
    Returns facet vertices stored in input file
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def read_casa_polys(filename, image):
    """
    Reads casa region file and returns polys

    Note: only regions of type "poly" are supported
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    polys = []
    for line in lines:
        if line.startswith('poly'):
            poly_str_temp = line.split('[[')[1]
            poly_str = poly_str_temp.split(']]')[0]
            poly_str_list = poly_str.split('], [')
            ra = []
            dec = []
            for pos in poly_str_list:
                RAstr, Decstr = pos.split(',')
                ra.append(Angle(RAstr, unit='hourangle').to('deg').value)
                dec.append(Angle(Decstr.replace('.', ':', 2), unit='deg').to('deg').value)
            poly_vertices = [np.array(ra), np.array(dec)]

            # Convert to image-plane polygon
            xvert = []
            yvert = []
            for RAvert, Decvert in zip(np.array(ra), np.array(dec)):
                try:
                    pixels = image.topixel([0, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                except:
                    pixels = image.topixel([1, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                xvert.append(pixels[2]) # x -> Dec
                yvert.append(pixels[3]) # y -> RA
            polys.append(Polygon(xvert, yvert))

        elif line.startswith('ellipse'):
            ell_str_temp = line.split('[[')[1]
            if '], 0.0' not in ell_str_temp and '], 90.0' not in ell_str_temp:
                print('Only position angles of 0.0 and 90.0 are support for CASA '
                    'regions of type "ellipse"')
                sys.exit(1)
            if '], 0.0' in ell_str_temp:
                ell_str = ell_str_temp.split('], 0.0')[0]
                pa = 0
            else:
                ell_str = ell_str_temp.split('], 90.0')[0]
                pa = 90
            ell_str_list = ell_str.split('], [')

            # Ellipse center
            RAstr, Decstr = ell_str_list[0].split(',')
            ra_center = Angle(RAstr, unit='hourangle').to('deg').value
            dec_center = Angle(Decstr.replace('.', ':', 2), unit='deg').to('deg').value
            pixels = image.topixel([0, 1, dec_center*np.pi/180.0,
                ra_center*np.pi/180.0])
            x_center = pixels[2] # x -> Dec
            y_center = pixels[3] # y -> RA

            # Ellipse semimajor and semiminor axes
            a_str, b_str = ell_str_list[1].split(',')
            a_deg = float(a_str.split('arcsec')[0])/3600.0
            b_deg = float(b_str.split('arcsec')[0])/3600.0
            pixels1 = image.topixel([0, 1, (dec_center-a_deg/2.0)*np.pi/180.0,
                ra_center*np.pi/180.0])
            a_pix1 = pixels1[2]
            pixels2 = image.topixel([0, 1, (dec_center+a_deg/2.0)*np.pi/180.0,
                ra_center*np.pi/180.0])
            a_pix2 = pixels2[2]
            a_pix = abs(a_pix2 - a_pix1)
            ex = []
            ey = []
            for th in range(0, 360, 1):
                if pa == 0:
                    # semimajor axis is along x-axis
                    ex.append(a_pix * np.cos(th * np.pi / 180.0)
                        + x_center) # x -> Dec
                    ey.append(a_pix * b_deg / a_deg * np.sin(th * np.pi / 180.0) + y_center) # y -> RA
                elif pa == 90:
                    # semimajor axis is along y-axis
                    ex.append(a_pix * b_deg / a_deg * np.cos(th * np.pi / 180.0)
                        + x_center) # x -> Dec
                    ey.append(a_pix * np.sin(th * np.pi / 180.0) + y_center) # y -> RA
            polys.append(Polygon(ex, ey))

        elif line.startswith('box'):
            poly_str_temp = line.split('[[')[1]
            poly_str = poly_str_temp.split(']]')[0]
            poly_str_list = poly_str.split('], [')
            ra = []
            dec = []
            for pos in poly_str_list:
                RAstr, Decstr = pos.split(',')
                ra.append(Angle(RAstr, unit='hourangle').to('deg').value)
                dec.append(Angle(Decstr.replace('.', ':', 2), unit='deg').to('deg').value)
            ra.insert(1, ra[0])
            dec.insert(1, dec[1])
            ra.append(ra[2])
            dec.append(dec[0])
            poly_vertices = [np.array(ra), np.array(dec)]

            # Convert to image-plane polygon
            xvert = []
            yvert = []
            for RAvert, Decvert in zip(np.array(ra), np.array(dec)):
                try:
                    pixels = image.topixel([0, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                except:
                    pixels = image.topixel([1, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                xvert.append(pixels[2]) # x -> Dec
                yvert.append(pixels[3]) # y -> RA
            polys.append(Polygon(xvert, yvert))

        elif line.startswith('#'):
            pass

        else:
            print('Only CASA regions of type "poly", "box", or "ellipse" are supported')
            sys.exit(1)

    return polys


def make_template_image(image_name, reference_ra_deg, reference_dec_deg,
    imsize=512, cellsize_deg=0.000417):
    """
    Make a blank image and save it to disk

    Parameters
    ----------
    image_name : str
        Filename of output image
    reference_ra_deg : float, optional
        RA for center of output mask image
    reference_dec_deg : float, optional
        Dec for center of output mask image
    imsize : int, optional
        Size of output image
    cellsize_deg : float, optional
        Size of a pixel in degrees

    """
    shape_out = [1, 1, imsize, imsize]
    hdu = pyfits.PrimaryHDU(np.zeros(shape_out, dtype=np.float32))
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header

    # Add WCS info
    header['CRVAL1'] = reference_ra_deg
    header['CDELT1'] = -cellsize_deg
    header['CRPIX1'] = imsize/2.0
    header['CUNIT1'] = 'deg'
    header['CTYPE1'] = 'RA---SIN'
    header['CRVAL2'] = reference_dec_deg
    header['CDELT2'] = cellsize_deg
    header['CRPIX2'] = imsize/2.0
    header['CUNIT2'] = 'deg'
    header['CTYPE2'] = 'DEC--SIN'

    # Add STOKES info
    header['CRVAL3'] = 1.0
    header['CDELT3'] = 1.0
    header['CRPIX3'] = 1.0
    header['CUNIT3'] = ''
    header['CTYPE3'] = 'STOKES'

    # Add frequency info
    header['RESTFRQ'] = 15036
    header['CRVAL4'] = 150e6
    header['CDELT4'] = 3e8
    header['CRPIX4'] = 1.0
    header['CUNIT4'] = 'HZ'
    header['CTYPE4'] = 'FREQ'
    header['SPECSYS'] = 'TOPOCENT'

    # Add equinox
    header['EQUINOX'] = 2000.0

    # Add telescope
    header['TELESCOP'] = 'LOFAR'

    hdulist[0].header = header

    hdulist.writeto(image_name, clobber=True)
    hdulist.close()



def main(image_name, mask_name, atrous_do=False, threshisl=0.0, threshpix=0.0, rmsbox=None,
         rmsbox_bright=(35, 7), iterate_threshold=False, adaptive_rmsbox=False, img_format='fits',
         threshold_format='float', trim_by=0.0, vertices_file=None, atrous_jmax=6,
         pad_to_size=None, skip_source_detection=False, region_file=None, nsig=1.0,
         reference_ra_deg=None, reference_dec_deg=None, cellsize_deg=0.000417,
         use_adaptive_threshold=False, adaptive_thresh=150.0):
    """
    Make a clean mask and return clean threshold

    Parameters
    ----------
    image_name : str
        Filename of input image from which mask will be made. If the image does
        not exist, a template image with center at (reference_ra_deg,
        reference_dec_deg) will be made internally
    mask_name : str
        Filename of output mask image
    atrous_do : bool, optional
        Use wavelet module of PyBDSM?
    threshisl : float, optional
        Value of thresh_isl PyBDSM parameter
    threshpix : float, optional
        Value of thresh_pix PyBDSM parameter
    rmsbox : tuple of floats, optional
        Value of rms_box PyBDSM parameter
    rmsbox_bright : tuple of floats, optional
        Value of rms_box_bright PyBDSM parameter
    iterate_threshold : bool, optional
        If True, threshold will be lower in 20% steps until
        at least one island is found
    adaptive_rmsbox : tuple of floats, optional
        Value of adaptive_rms_box PyBDSM parameter
    img_format : str, optional
        Format of output mask image (one of 'fits' or 'casa')
    threshold_format : str, optional
        Format of output threshold (one of 'float' or 'str_with_units')
    trim_by : float, optional
        Fraction by which the perimeter of the output mask will be
        trimmed (zeroed)
    vertices_file : str, optional
        Filename of file with vertices (must be a pickle file containing
        a dictionary with the vertices in the 'vertices' entry)
    atrous_jmax : int, optional
        Value of atrous_jmax PyBDSM parameter
    pad_to_size : int, optional
        Pad output mask image to a size of pad_to_size x pad_to_size
    skip_source_detection : bool, optional
        If True, source detection is not run on the input image
    region_file : str, optional
        Filename of region file in CASA format. If given, no mask image
        is made (the region file is used as the clean mask)
    nsig : float, optional
        Number of sigma of returned threshold value
    reference_ra_deg : float, optional
        RA for center of output mask image
    reference_dec_deg : float, optional
        Dec for center of output mask image
    cellsize_deg : float, optional
        Size of a pixel in degrees
    use_adaptive_threshold : bool, optional
        If True, use an adaptive threshold estimated from the negative values in
        the image
    adaptive_thresh : float, optional
        If adaptive_rmsbox is True, this value sets the threshold above
        which a source will use the small rms box

    Returns
    -------
    result : dict
        Dict with nsig-sigma rms threshold

    """
    if rmsbox is not None and type(rmsbox) is str:
        rmsbox = eval(rmsbox)

    if type(rmsbox_bright) is str:
        rmsbox_bright = eval(rmsbox_bright)

    if pad_to_size is not None and type(pad_to_size) is str:
        pad_to_size = int(pad_to_size)

    if type(atrous_do) is str:
        if atrous_do.lower() == 'true':
            atrous_do = True
            threshisl = 4.0 # override user setting to ensure proper source fitting
        else:
            atrous_do = False

    if type(iterate_threshold) is str:
        if iterate_threshold.lower() == 'true':
            iterate_threshold = True
        else:
            iterate_threshold = False

    if type(adaptive_rmsbox) is str:
        if adaptive_rmsbox.lower() == 'true':
            adaptive_rmsbox = True
        else:
            adaptive_rmsbox = False

    if type(skip_source_detection) is str:
        if skip_source_detection.lower() == 'true':
            skip_source_detection = True
        else:
            skip_source_detection = False

    if type(use_adaptive_threshold) is str:
        if use_adaptive_threshold.lower() == 'true':
            use_adaptive_threshold = True
        else:
            use_adaptive_threshold = False

    if reference_ra_deg is not None and reference_dec_deg is not None:
        reference_ra_deg = float(reference_ra_deg)
        reference_dec_deg = float(reference_dec_deg)

    if not os.path.exists(image_name):
        print('Input image not found. Making empty image...')
        if not skip_source_detection:
            print('ERROR: Source detection cannot be done on an empty image')
            sys.exit(1)
        if reference_ra_deg is not None and reference_dec_deg is not None:
            image_name = mask_name + '.tmp'
            make_template_image(image_name, reference_ra_deg, reference_dec_deg,
                cellsize_deg=float(cellsize_deg))
        else:
            print('ERROR: if image not found, a refernce position must be given')
            sys.exit(1)

    trim_by = float(trim_by)
    atrous_jmax = int(atrous_jmax)
    threshpix = float(threshpix)
    threshisl = float(threshisl)
    nsig = float(nsig)
    adaptive_thresh = float(adaptive_thresh)
    threshold = 0.0

    if not skip_source_detection:
        if vertices_file is not None:
            # Modify the input image to blank the regions outside of the polygon
            temp_img = pim.image(image_name)
            image_name += '.blanked'
            temp_img.saveas(image_name, overwrite=True)
            input_img = pim.image(image_name)
            data = input_img.getdata()

            vertices = read_vertices(vertices_file)
            RAverts = vertices[0]
            Decverts = vertices[1]
            xvert = []
            yvert = []
            for RAvert, Decvert in zip(RAverts, Decverts):
                pixels = input_img.topixel([1, 1, Decvert*np.pi/180.0,
                    RAvert*np.pi/180.0])
                xvert.append(pixels[2]) # x -> Dec
                yvert.append(pixels[3]) # y -> RA
            poly = Polygon(xvert, yvert)

            # Find masked regions
            masked_ind = np.where(data[0, 0])

            # Find distance to nearest poly edge and set to NaN those that
            # are outside the facet (dist < 0)
            dist = poly.is_inside(masked_ind[0], masked_ind[1])
            outside_ind = np.where(dist < 0.0)
            if len(outside_ind[0]) > 0:
                data[0, 0, masked_ind[0][outside_ind], masked_ind[1][outside_ind]] = np.nan

            # Save changes
            input_img.putdata(data)

        if use_adaptive_threshold:
            # Get an estimate of the rms
            img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                     thresh_pix=threshpix, thresh_isl=threshisl,
                                     atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=adaptive_thresh,
                                     rms_box_bright=rmsbox_bright, rms_map=True, quiet=True,
                                     atrous_jmax=atrous_jmax, stop_at='isl')

            # Find min and max pixels
            max_neg_val = abs(np.min(img.ch0_arr))
            max_neg_pos = np.where(img.ch0_arr == np.min(img.ch0_arr))
            max_pos_val = abs(np.max(img.ch0_arr))
            max_pos_pos = np.where(img.ch0_arr == np.max(img.ch0_arr))

            # Estimate new thresh_isl from min pixel value's sigma, but don't let
            # it get higher than 1/2 of the peak's sigma
            threshisl_neg = 2.0 * max_neg_val / img.rms_arr[max_neg_pos][0]
            max_sigma = max_pos_val / img.rms_arr[max_pos_pos][0]
            if threshisl_neg > max_sigma / 2.0:
                threshisl_neg = max_sigma / 2.0

            # Use the new threshold only if it is larger than the user-specified one
            if threshisl_neg > threshisl:
                threshisl = threshisl_neg

        if iterate_threshold:
            # Start with given threshold and lower it until we get at least one island
            nisl = 0
            while nisl == 0:
                img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                         thresh_pix=threshpix, thresh_isl=threshisl,
                                         atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                         adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=adaptive_thresh,
                                         rms_box_bright=rmsbox_bright, rms_map=True, quiet=True,
                                         atrous_jmax=atrous_jmax)
                nisl = img.nisl
                threshpix /= 1.2
                threshisl /= 1.2
                if threshpix < 5.0:
                    break
        else:
            img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                     thresh_pix=threshpix, thresh_isl=threshisl,
                                     atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=adaptive_thresh,
                                     rms_box_bright=rmsbox_bright, rms_map=True, quiet=True,
                                     atrous_jmax=atrous_jmax)

        if img.nisl == 0:
            if region_file is None or region_file == '[]':
                print('No islands found. Clean mask cannot be made.')
                sys.exit(1)
            else:
                # Continue on and use user-supplied region file
                skip_source_detection = True
                threshold = nsig * img.clipped_rms

        # Check if there are large islands preset (indicating that multi-scale
        # clean is needed)
        has_large_isl = False
        for isl in img.islands:
            if isl.size_active > 100:
                # Assuming normal sampling, a size of 100 pixels would imply
                # a source of ~ 10 beams
                has_large_isl = True

    if (region_file is not None and region_file != '[]' and skip_source_detection):
        # Copy region file and return if source detection was not done
        os.system('cp {0} {1}'.format(region_file.strip('[]"'), mask_name))
        if threshold_format == 'float':
            return {'threshold_5sig': threshold}
        elif threshold_format == 'str_with_units':
            # This is done to get around the need for quotes around strings in casapy scripts
            # 'casastr/' is removed by the generic pipeline
            return {'threshold_5sig': 'casastr/{0}Jy'.format(threshold)}
    elif not skip_source_detection:
        img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
                         img_format=img_format, clobber=True)

    if (vertices_file is not None or trim_by > 0 or pad_to_size is not None
        or (region_file is not None and region_file != '[]')
        or skip_source_detection):
        # Alter the mask in various ways
        if skip_source_detection:
            # Read the image
            mask_im = pim.image(image_name)
        else:
            # Read the PyBDSM mask
            mask_im = pim.image(mask_name)
        data = mask_im.getdata()
        coordsys = mask_im.coordinates()
        if reference_ra_deg is not None and reference_dec_deg is not None:
            values = coordsys.get_referencevalue()
            values[2][0] = reference_dec_deg/180.0*np.pi
            values[2][1] = reference_ra_deg/180.0*np.pi
            coordsys.set_referencevalue(values)
        imshape = mask_im.shape()
        del(mask_im)

        if pad_to_size is not None:
            imsize = pad_to_size
            coordsys['direction'].set_referencepixel([imsize/2, imsize/2])
            pixmin = (imsize - imshape[2]) / 2
            if pixmin < 0:
                print("The padded size must be larger than the original size.")
                sys.exit(1)
            pixmax = pixmin + imshape[2]
            data_pad = np.zeros((1, 1, imsize, imsize), dtype=np.float32)
            data_pad[0, 0, pixmin:pixmax, pixmin:pixmax] = data[0, 0]
            new_mask = pim.image('', shape=(1, 1, imsize, imsize), coordsys=coordsys)
            new_mask.putdata(data_pad)
        else:
            new_mask = pim.image('', shape=imshape, coordsys=coordsys)
            new_mask.putdata(data)

        data = new_mask.getdata()

        if skip_source_detection:
            # Mask all pixels
            data[:] = 1

        if vertices_file is not None:
            # Modify the clean mask to exclude regions outside of the polygon
            vertices = read_vertices(vertices_file)
            RAverts = vertices[0]
            Decverts = vertices[1]
            xvert = []
            yvert = []
            for RAvert, Decvert in zip(RAverts, Decverts):
                try:
                    pixels = new_mask.topixel([0, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                except:
                    pixels = new_mask.topixel([1, 1, Decvert*np.pi/180.0,
                                               RAvert*np.pi/180.0])
                xvert.append(pixels[2]) # x -> Dec
                yvert.append(pixels[3]) # y -> RA
            poly = Polygon(xvert, yvert)

            # Find masked regions
            masked_ind = np.where(data[0, 0])

            # Find distance to nearest poly edge and unmask those that
            # are outside the facet (dist < 0)
            dist = poly.is_inside(masked_ind[0], masked_ind[1])
            outside_ind = np.where(dist < 0.0)
            if len(outside_ind[0]) > 0:
                data[0, 0, masked_ind[0][outside_ind], masked_ind[1][outside_ind]] = 0

        if trim_by > 0.0:
            sh = np.shape(data)
            margin = int(sh[2] * trim_by / 2.0 )
            data[0, 0, 0:sh[2], 0:margin] = 0
            data[0, 0, 0:margin, 0:sh[3]] = 0
            data[0, 0, 0:sh[2], sh[3]-margin:sh[3]] = 0
            data[0, 0, sh[2]-margin:sh[2], 0:sh[3]] = 0

        if region_file is not None and region_file != '[]':
            # Merge the CASA regions with the mask
            casa_polys = read_casa_polys(region_file.strip('[]"'), new_mask)
            for poly in casa_polys:
                # Find unmasked regions
                unmasked_ind = np.where(data[0, 0] == 0)

                # Find distance to nearest poly edge and mask those that
                # are inside the casa region (dist > 0)
                dist = poly.is_inside(unmasked_ind[0], unmasked_ind[1])
                inside_ind = np.where(dist > 0.0)
                if len(inside_ind[0]) > 0:
                    data[0, 0, unmasked_ind[0][inside_ind], unmasked_ind[1][inside_ind]] = 1

        # Save changes
        new_mask.putdata(data)
        if img_format == 'fits':
            new_mask.tofits(mask_name, overwrite=True)
        elif img_format == 'casa':
            new_mask.saveas(mask_name, overwrite=True)
        else:
            print('Output image format "{}" not understood.'.format(img_format))
            sys.exit(1)

    if not skip_source_detection:
        if threshold_format == 'float':
            return {'threshold_5sig': nsig * img.clipped_rms, 'multiscale': has_large_isl}
        elif threshold_format == 'str_with_units':
            # This is done to get around the need for quotes around strings in casapy scripts
            # 'casastr/' is removed by the generic pipeline
            return {'threshold_5sig': 'casastr/{0}Jy'.format(nsig * img.clipped_rms),
                'multiscale': has_large_isl}
    else:
        return {'threshold_5sig': '0.0'}


if __name__ == '__main__':
    descriptiontext = "Make a clean mask.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_name', help='Image name')
    parser.add_argument('mask_name', help='Mask name')
    parser.add_argument('-a', '--atrous_do', help='use wavelet fitting', type=bool, default=False)
    parser.add_argument('-i', '--threshisl', help='', type=float, default=3.0)
    parser.add_argument('-p', '--threshpix', help='', type=float, default=5.0)
    parser.add_argument('-r', '--rmsbox', help='rms box width and step (e.g., "(60, 20)")',
        type=str, default='(60, 20)')
    parser.add_argument('--rmsbox_bright', help='rms box for bright sources(?) width and step (e.g., "(60, 20)")',
        type=str, default='(60, 20)')
    parser.add_argument('-t', '--iterate_threshold', help='iteratively decrease threshold until at least '
        'one island is found', type=bool, default=False)
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)
    parser.add_argument('-f', '--img_format', help='format of output mask', type=str, default='casa')
    parser.add_argument('-d', '--threshold_format', help='format of return value', type=str, default='float')
    parser.add_argument('-b', '--trim_by', help='Trim masked region by this number of pixels', type=float, default=0.0)
    parser.add_argument('-v', '--vertices_file', help='file containing facet polygon vertices', type=str, default=None)
    parser.add_argument('--region_file', help='File containing casa regions to be merged with the detected mask', type=str, default=None)
    parser.add_argument('-j', '--atrous_jmax', help='Max wavelet scale', type=int, default=3)
    parser.add_argument('-z', '--pad_to_size', help='pad mask to this size', type=int, default=None)
    parser.add_argument('-s', '--skip_source_detection', help='skip source detection', type=bool, default=False)

    args = parser.parse_args()
    erg = main(args.image_name, args.mask_name, atrous_do=args.atrous_do,
               threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
               rmsbox_bright=args.rmsbox_bright,
               iterate_threshold=args.iterate_threshold,
               adaptive_rmsbox=args.adaptive_rmsbox, img_format=args.img_format,
               threshold_format=args.threshold_format, trim_by=args.trim_by,
               vertices_file=args.vertices_file, atrous_jmax=args.atrous_jmax,
               pad_to_size=args.pad_to_size, skip_source_detection=args.skip_source_detection,
               region_file=args.region_file)
    print erg
