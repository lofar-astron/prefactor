#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import logging
import casacore.tables as ct
import numpy as np


def find_flagged_fraction(ms_file):
    """
    Finds flagged fraction for all stations

    Parameters
    ----------
    ms_file : str
        Name (path) of input MS

    Returns
    -------
    result : dict
        Dict of station_name:flagged_fraction for all stations
    """
    # Use taql to find the flagged fraction per station
    ant = ct.table(ms_file+'::ANTENNA', ack=False)
    station_names = ant.col('NAME')[:]
    ant.close()

    t = ct.table(ms_file, ack=False)
    fraction_flagged = {}
    for i, station_name in enumerate(station_names):
        t1 = t.query('ANTENNA1 == {0} or ANTENNA2 == {0}'.format(i))
        flags_per_element = t1.calc('ntrue(FLAG)')
        nelements = t1.calc('nelements(FLAG)')[0]  # = number of channels * number of pols
        fraction_flagged[station_name] = float(np.sum(flags_per_element)) / nelements / len(t1)
        t1.close()
    t.close()

    return fraction_flagged


def main(ms_file):
    # derive the fraction of flagged data of the entire observation
    logging.info('Reading data.')
    flagged_fraction_dict = find_flagged_fraction(ms_file)

    return flagged_fraction_dict


if __name__ == '__main__':
    descriptiontext = "Check the flagged fraction of stations in a MS.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inputms', help='name of the input MS')
    args = parser.parse_args()

    erg = main(ms_file=args.inputms)
