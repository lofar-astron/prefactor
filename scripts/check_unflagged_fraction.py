#!/usr/bin/env python
"""
Flag MS if unflagged fraction is too low.
"""
import argparse
import os
import numpy as np

def find_unflagged_fraction(ms_file):
    """
    Finds the fraction of data that is unflagged
    Parameters
    ----------
    ms_file : str
        Filename of input MS
    Returns
    -------
    unflagged_fraction : float
        Fraction of unflagged data
    """
    import subprocess

    # Call taql. Note that we do not use pt.taql(), as pt.taql() can cause
    # hanging/locking issues on some systems
    p = subprocess.Popen("taql 'CALC sum([select nfalse(FLAG) from {0}]) / "
        "sum([select nelements(FLAG) from {0}])'".format(ms_file),
        shell=True, stdout=subprocess.PIPE)
    r = p.communicate()
    unflagged_fraction = float(r[0])

    return unflagged_fraction

def main(ms_file, min_fraction=0.01, print_fraction=True):
    """
    Flag MS if unflagged fraction is too low.

    Parameters
    ----------
    ms_file : str
        Name (path) of input MS
    min_fraction : float , optional
        minimum fraction of unflagged data needed to keep this MS
    print_fraction : bool, optional
        print the actual fration of unflagged data

    Returns
    -------
    result : dict
        Dict with the name of the input MS or "None"

    """
    min_fraction = float(min_fraction)

    unflagged_fraction = find_unflagged_fraction(ms_file)
    if print_fraction:
        print("File %s has %.2f%% unflagged data."%(os.path.basename(ms_file.rstrip('/')),unflagged_fraction * 100.))
    if unflagged_fraction < min_fraction:
        print('Unflagged fraction of {0} is: {1}%, '
              'removing file.'.format(os.path.basename(ms_file.rstrip('/')), str(unflagged_fraction * 100)))
        return {'flagged': 'None',  'unflagged_fraction':  unflagged_fraction, 'filename': ms_file}
    else:
        return {'flagged': ms_file, 'unflagged_fraction':  unflagged_fraction, 'filename': ms_file}

if __name__ == '__main__':
    descriptiontext = "Check a MS for a minimum fraction of unflagged data.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inputms', help='name of the input MS')
    parser.add_argument('-f', '--min_fraction', help='Minimum fraction of unflagged data needed to keep this MS '
                        '(default 0.01 = "keep if at least 1%% is not flagged")',  type=float, default=0.01)
    args = parser.parse_args()

    erg = main(args.inputms, args.min_fraction,print_fraction=True)
