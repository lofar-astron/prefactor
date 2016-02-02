#! /usr/bin/env python
"""
Script to make a clean mask from an image
"""
import argparse
from argparse import RawTextHelpFormatter
try:
    # Temp version of PyBDSM that uses faster fftconvolve
    from lofar2 import bdsm
except ImportError:
    from lofar import bdsm
import pyrap.images as pim
import pickle
import numpy as np
import sys
import os

class Polygon:
    """
    Polygon object.

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
        '''Check if point is inside a general polygon.

        Input parameters:

        xpoint -- The x-coord of the point to be tested.
        ypoint -- The y-coords of the point to be tested.
        smalld -- A small float number.

        xpoint and ypoint could be scalars or array-like sequences.

        Output parameters:

        mindst -- The distance from the point to the nearest point of the
                  polygon.
                  If mindst < 0 then point is outside the polygon.
                  If mindst = 0 then point in on a side of the polygon.
                  If mindst > 0 then point is inside the polygon.

        Notes:

        An improved version of the algorithm of Nordbeck and Rydstedt.

        REF: SLOAN, S.W. (1985): A point-in-polygon program. Adv. Eng.
             Software, Vol 7, No. 1, pp 45-47.

        '''
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
            raise IndexError('x and y has different shapes')
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
    Returns facet vertices
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def main(image_name, mask_name, atrous_do=False, threshisl=0.0, threshpix=0.0, rmsbox=None,
         iterate_threshold=False, adaptive_rmsbox=False, img_format='fits',
         threshold_format='float', trim_by=0.0, vertices_file=None, atrous_jmax=6,
         pad_to_size=None, skip_source_detection=False, region_file=None):
    """
    Run PyBDSM to make an island clean mask

    Parameters
    ----------
    TODO

    Returns
    -------
    result : dict
        Dict with 5-sigma rms threshold

    """
    if region_file is not None:
        if region_file != '[]':
            # Copy the CASA region file (stripped of brackets, etc.) and return
            os.system('cp {0} {1}'.format(region_file.strip('[]"'), mask_name))
            return {'threshold_5sig': '0.0'}

    if atrous_do:
        threshisl = 4.0

    if rmsbox is not None and type(rmsbox) is str:
        rmsbox = eval(rmsbox)

    if pad_to_size is not None and type(pad_to_size) is str:
        pad_to_size = int(pad_to_size)

    if type(atrous_do) is str:
        if atrous_do.lower() == 'true':
            atrous_do = True
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

    trim_by = float(trim_by)
    atrous_jmax = int(atrous_jmax)
    threshpix = float(threshpix)
    threshisl = float(threshisl)

    if not skip_source_detection:
        if iterate_threshold:
            # Start with given threshold and lower it until we get at least one island
            nisl = 0
            while nisl == 0:
                img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                         thresh_pix=threshpix, thresh_isl=threshisl,
                                         atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                         adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=150,
                                         rms_box_bright=(35,7), rms_map=True, quiet=True,
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
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=150,
                                     rms_box_bright=(35,7), rms_map=True, quiet=True,
                                     atrous_jmax=atrous_jmax)

        if img.nisl == 0:
            print('No islands found. Clean mask cannot be made.')
            sys.exit(1)

        img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
                         img_format=img_format, clobber=True)

    if vertices_file is not None or trim_by > 0 or pad_to_size is not None or skip_source_detection:
        # Alter the mask in various ways
        if skip_source_detection:
            # Read the image
            mask_im = pim.image(image_name)
        else:
            # Read the PyBDSM mask
            mask_im = pim.image(mask_name)
        data = mask_im.getdata()
        coordsys = mask_im.coordinates()
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
                pixels = new_mask.topixel([0, 1, Decvert*np.pi/180.0,
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

        # Save changes
        new_mask.putdata(data)
        if img_format == 'fits':
            new_mask.tofits(mask_name, overwrite=True)
        else:
            new_mask.saveas(mask_name, overwrite=True)

    if not skip_source_detection:
        if threshold_format == 'float':
            return {'threshold_5sig': 5.0 * img.clipped_rms}
        elif threshold_format == 'str_with_units':
            # This is done to get around the need for quotes around strings in casapy scripts
            # 'casastr/' is removed by the generic pipeline
            return {'threshold_5sig': 'casastr/{0}Jy'.format(5.0 * img.clipped_rms)}
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
    parser.add_argument('-t', '--iterate_threshold', help='iteratively decrease threshold until at least '
        'one island is found', type=bool, default=False)
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)
    parser.add_argument('-f', '--img_format', help='format of output mask', type=str, default='casa')
    parser.add_argument('-d', '--threshold_format', help='format of return value', type=str, default='float')
    parser.add_argument('-b', '--trim_by', help='Trim masked region by this number of pixels', type=float, default=0.0)
    parser.add_argument('-v', '--vertices_file', help='file containing facet polygon vertices', type=str, default=None)
    parser.add_argument('-j', '--atrous_jmax', help='Max wavelet scale', type=int, default=3)
    parser.add_argument('-z', '--pad_to_size', help='pad mask to this size', type=int, default=None)
    parser.add_argument('-s', '--skip_source_detection', help='skip source detection', type=bool, default=False)

    args = parser.parse_args()
    main(args.image_name, args.mask_name, atrous_do=args.atrous_do,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
         iterate_threshold=args.iterate_threshold,
         adaptive_rmsbox=args.adaptive_rmsbox, img_format=args.img_format,
         threshold_format=args.threshold_format, trim_by=args.trim_by,
         vertices_file=args.vertices_file, atrous_jmax=args.atrous_jmax,
         pad_to_size=args.pad_to_size, skip_source_detection=args.skip_source_detection)
