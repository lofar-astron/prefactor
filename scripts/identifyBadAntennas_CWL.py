#!/usr/bin/env python
"""
Identify fully flagged antennas
"""
import casacore.tables as ct
import numpy as np
import logging


def input2strlist_nomapfile(invar):
    """
    from bin/download_IONEX.py
    give the list of MSs from the list provided as a string
    """
    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list


def find_flagged_antennas(ms_file, min_fraction=1.0):
    """
    Finds antennas with flagged fractions equal to or above the specified min

    Parameters
    ----------
    ms_file : str
        Name (path) of input MS
    min_fraction : float, optional
        Minimum fraction of flagged data: stations with flagged fractions equal to
        or above this value are considered to be flagged

    Returns
    -------
    result : str
        Comma-separated string of the names of the flagged stations. If no stations
        are flagged, an empty string is returned
    """
    # Use taql to find the flagged fraction per station
    logging.info('Reading ' + str(ms_file))
    ant = ct.table(ms_file+'::ANTENNA', ack=False)
    station_names = ant.col('NAME')[:]
    ant.close()

    t = ct.table(ms_file, ack=False)
    flaggedants = []
    for i, station_name in enumerate(station_names):
        t1 = t.query('ANTENNA1 == {0} or ANTENNA2 == {0}'.format(i))
        flags_per_element = t1.calc('ntrue(FLAG)')
        nelements = t1.calc('nelements(FLAG)')[0]  # = number of channels * number of pols
        flagged_fraction = float(np.sum(flags_per_element)) / nelements / len(t1)
        if flagged_fraction >= min_fraction:
            flaggedants.append(station_name)
        t1.close()
    t.close()

    return ','.join(flaggedants)


def main(MSfile):
    ms_file = input2strlist_nomapfile(MSfile)[0]

    flaggedants = find_flagged_antennas(ms_file)

    # return results
    result = {'flaggedants': str(flaggedants)}
    logging.info(str(flaggedants) + ' are flagged in ' + str(ms_file))
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Identify fully flagged antennas.')

    parser.add_argument('MSfiles', type=str, nargs='+',
                        help='One (or more MSs) for which you want to get the antenna list.')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s", "%Y-%m-%d %H:%M:%S")
    format_file = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    main(args.MSfile)
