#!/usr/bin/env python
import os,sys
import glob
import pyrap.tables as pt
import numpy as np
import time
import lsmtool

########################################################################
def grab_coord_MS(MS):
    """
    Read the coordinates of a field from one MS corresponding to the selection given in the parameters

    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field

    Returns
    -------
    RA, Dec : "tuple"
        coordinates of the field (RA, Dec in deg , J2000)
    """

    # reading the coordinates ("position") from the MS
    # NB: they are given in rad,rad (J2000)
    [[[ra,dec]]] = pt.table(MS+'::FIELD', readonly=True, ack=False).getcol('PHASE_DIR')

    # RA is stocked in the MS in [-pi;pi]
    # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
    if ra<0:
        ra=ra+2*np.pi

    # convert radians to degrees
    ra_deg =  ra/np.pi*180.
    dec_deg = dec/np.pi*180.

    # and sending the coordinates in deg
    return ra_deg,dec_deg


########################################################################
def input2strlist_nomapfile(invar):
    """ from bin/download_IONEX.py
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


########################################################################
def main(ms_input, SkymodelPath, Radius="5.", DoDownload="True", Source="TGSS"):
    """
    Download the skymodel for the target field

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    SkymodelPath : str
        Full name (with path) to the skymodel; if YES is true, the skymodel will be downloaded here
    Radius : string with float (default = "5.")
        Radius for the TGSS/GSM cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the TGSS skymodel or GSM.
        "Force": download skymodel from TGSS or GSM, delete existing skymodel if needed.
        "True" or "Yes": use existing skymodel file if it exists, download skymodel from
                         TGSS or GSM if it does not.
        "False" or "No": Do not download skymodel, raise an exception if skymodel
                         file does not exist.

    """

    FileExists = os.path.isfile(SkymodelPath)
    if (not FileExists and os.path.exists(SkymodelPath)):
        raise ValueError("download_tgss_skymodel_target: Path: \"%s\" exists but is not a file!"%(SkymodelPath))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(SkymodelPath)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising skymodel in "+ SkymodelPath
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising skymodel in "+ SkymodelPath
            return
         else:
            raise ValueError("download_tgss_skymodel_target: Path: \"%s\" does not exist and skymodel download is disabled!"%(SkymodelPath))

    # If we got here, then we are supposed to download the skymodel.
    assert download_flag is True # Jaja, belts and suspenders...
    print "DOWNLOADING skymodel for the target into "+ SkymodelPath

    # Reading a MS to find the coordinate (pyrap)
    [RATar,DECTar]=grab_coord_MS(input2strlist_nomapfile(ms_input)[0])

    # Downloading the skymodel, skip after five tries
    errorcode = 1
    tries     = 0
    while errorcode != 0 and tries < 5:
        if Source == 'TGSS':
            errorcode = os.system("wget -O "+SkymodelPath+ " \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv3.cgi?coord="+str(RATar)+","+str(DECTar)+"&radius="+str(Radius)+"&unit=deg&deconv=y\' ")
            pass
        elif Source == 'GSM':
            errorcode = os.system("wget -O "+SkymodelPath+ " \'http://172.104.228.177/cgi-bin/gsmv1.cgi?coord="+str(RATar)+","+str(DECTar)+"&radius="+str(Radius)+"&unit=deg&deconv=y\' ")
            pass
        time.sleep(5)
        tries += 1
        pass

    if not os.path.isfile(SkymodelPath):
        raise IOError("download_tgss_skymodel_target: Path: \"%s\" does not exist after trying to download the skymodel."%(SkymodelPath))

    # Treat all sources as one group (direction)
    skymodel = lsmtool.load(SkymodelPath)
    skymodel.group('single')
    skymodel.write(clobber=True)
    
    return


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=' Download the TGSS or GSM skymodel for the target field')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which a TGSS/GSM skymodel will be download.')
    parser.add_argument('SkyTar', type=str,
                        help='Full name (with path) to the skymodel; the TGSS/GSM skymodel will be downloaded here')
    parser.add_argument('--Radius', type=float, default=5.,
                        help='Radius for the TGSS/GSM cone search in degrees')
    parser.add_argument('--Source', type=str, default='TGSS',
                        help='Choose source for skymodel: TGSS or GSM')
    parser.add_argument('--DoDownload', type=str, default="True",
                        help='Download or not the TGSS skymodel or GSM ("Force" or "True" or "False").')

    args = parser.parse_args()
    radius=5
    if args.Radius:
        radius=args.Radius

    main(args.MSfile, args.SkyTar, str(radius), args.DoDownload, args.Source)
