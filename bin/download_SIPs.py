#! /usr/bin/env python

import get_MOM_data as MOM

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Download SIPs for given measurement sets from MOM.')

    parser.add_argument('MSfiles', type=str, nargs='+',
                        help='One or more MSs for which the SIPs should be downloaded.')
    parser.add_argument('--cache_file', type=str,
                        default='/media/scratch/test/horneff/Pipeline-Test/feedback_test/sip_cache.pkl',
                        help='File to store the downloaded SIPs. (Will be updated if it exists.)')
    parser.add_argument('-v','--verbose', action="store_true", default=False,
                        help="Generate more verbose output.")

    args = parser.parse_args()
    
    MOM.init_cache(args.cache_file)

    for msfile in args.MSfiles:
        MOM.get_SIP_from_MSfile(msfile, verbose=args.verbose)
