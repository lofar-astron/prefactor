#!/usr/bin/env python
"""
Script to merge two sky models
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import lsmtool
import sys


def main(inmodel1, inmodel2, outmodel, match_by='name', radius=0.0, keep='all'):
    """
    Creates mosaic

    Parameters
    ----------
    inmodel1 : str
        Filename of input model 1
    inmodel2 : str
        Filename of input model 2
    outmodel : str
        Filename of output model
    match_by : str, optional
        The match_by parameter of LSMTool
    radius : float, optional
        Radius in degrees for cross match
    keep : str, optional
        The keep parameter of LSMTool

    """
    s1 = lsmtool.load(inmodel1)
    s2 = lsmtool.load(inmodel2)

    s1.concatenate(s2, matchBy=match_by, radius=float(radius), keep=keep,
        inheritPatches=True)
    s1.group('every')
    s1.write(fileName=outmodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Merge two sky models.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('inmodel1', help='input model 1 filename')
    parser.add_argument('inmodel2', help='input model 2 filename')
    parser.add_argument('outmodel', help='output model filename')
    parser.add_argument('-m', '--match_by', help='match type', type=str, default='name')
    parser.add_argument('-r', '--radius', help='radius in degrees for matching', type=float, default=0.0)
    parser.add_argument('-k', '--keep', help='', type=str, default='all')

    args = parser.parse_args()
    main(args.inmodel1, args.inmodel2, args.outmodel, matchBy=args.match_by,
         radius=args.radius, keep=args.keep)
