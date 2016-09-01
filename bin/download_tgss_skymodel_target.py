#!/usr/bin/env python
import os,sys
import glob
import re
import pyrap.tables as pt
#import pylab
from astropy import units as u
from astropy.coordinates import SkyCoord as sc


########################################################################
def grab_coo_MS(MS):
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
    [[[ra,dec]]] = pt.table(MS+'/FIELD', readonly=True, ack=False).getcol('PHASE_DIR')
    
    # for a small conversion to degrees for the tgss query
    coo_tar=sc(ra,dec,unit=(u.rad, u.rad),frame="icrs")
    
    # and sending the coordinates in deg
    return coo_tar.ra.deg,coo_tar.dec.deg


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
def main(ms_input,SkymodelTar,Radius,Yes="True"):
    """
    Download the TGSS skymodel for the target field

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    SkymodelTar : str
        Full name (with path) to the skymodel; if YES is true, the TGSS skymodel will be downloaded here
    Radius : float
        Radius for the TGSS cone search in degrees
    YES : str ("True" or "False")
        Download or not the TGSS skymodel.
        (If False, one supposes the pipeline will use the skymodel provided in SkymodelTar)
    
    """
    
    if Yes=="True":# [TODO: (How to) send directly a bolean)
        print "DOWNLOADING TGSS Skymodel for the target into "+ SkymodelTar
        
        # Reading a MS to find the coordinate (pyrap)
        [RATar,DECTar]=grab_coo_MS(input2strlist_nomapfile(ms_input)[0])
        
        # Downloading the skymodel
        os.system("wget -O "+SkymodelTar+ " \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv2.cgi?coord="+str(RATar)+","+str(DECTar)+"&radius="+Radius+"&unit=deg&deconv=y\' ")
    else:
        print "USING the skymodel provided in "+ SkymodelTar
        
    
########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Find automatically between skymodels the one to use (for the Calibrator)')
    
    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which we search matching skymodel.')
    parser.add_argument('SkyTar', type=str, 
                        help='Full name (with path) to the skymodel; the TGSS skymodel will be downloaded here')
    parser.add_argument('--Radius', type=float, 
                        help='Radius for the TGSS cone search in degrees')
        
 
   
    args = parser.parse_args()
    radius=5
    if args.Radius:
        radius=args.Radius
    
    main(args.MSfile,args.SkyTar, str(radius))
    

