#!/usr/bin/env python
import os
import glob
import pyrap.tables as pt
import math
import lsmtool
import numpy
    
def grab_pointing(MS):
    """
    Read the name of the calibrator from one MS corresponding to the selection given in the parameters

    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field


    Returns
    -------
    nameCal : str
        Name of the calibrator as provided in the MS
        NB: we suppose that all the calibrators' observations have this field filled in (MS/Observation, column LOFAR_TARGET)
    """
    
    [ra, dec] = pt.table(MS+'/FIELD', readonly=True, ack=False).getcol('PHASE_DIR')[0][0] * 180 / math.pi
    return ra, dec

    

########################################################################

def check_skymodel(skymodel, ra, dec, max_separation_arcmin = 1.0):
    """
    Searches for an appropriate sky model
    """
    s = lsmtool.load(skymodel)
    dist_deg = s.getDistance(ra, dec)
    if any(dist_deg * 60.0 < max_separation_arcmin):
        patch_position = int(numpy.where(dist_deg * 60 < max_separation_arcmin)[0][0])
        patch_name = s.getPatchNames()[patch_position]
        return (True, patch_name)
        pass
    else:
        return (False, '')
        pass
    pass

########################################################################
def find_skymodel(ra, dec, PathSkyMod, extensionSky = ".skymodel", max_separation_arcmin = 1.0):
    """
    Find in the provided folder the correponding skymodel for the given source

    Parameters
    ----------
    NameSouce : str
        Name of the source for which we want to find the skymodel
    PathSkyMod : str
        Path to the folder containing the skymodels in use
        (with Pattern, will enable to read and get the info from a MS)
    [extensionSky] :  str
        Default: ".skymodel"
        extension of the skymodel files
 
    NB: we suppose that the name of the calibrator is included in the name of the skymodel

    Returns
    -------
    list_skymodel[0] : str
        Full name (with path) to the matching skymodel 
    """
    

    skymodels = glob.glob(PathSkyMod + "/*" + extensionSky)
    
    # remove any Ateam models from the listing (only in filenames)
    #skymodels = [s for s in skymodels if 'Ateam' not in s]
    #skymodels = [s for s in skymodels if 'A-Team' not in s]
    
    for skymodel in skymodels:
        check = check_skymodel(skymodel, ra, dec, max_separation_arcmin)
        if check[0]:
            print "The following skymodel will be used for the calibrator: " + skymodel.split("/")[-1] + " (in " + PathSkyMod + ")"
            return (skymodel, check[-1])
            pass
        else:
            pass
        pass
    
    raise TypeError('find_skymodel: SKYMODEL FOR THE CALIBRATOR NOT FOUND IN ' + PathSkyMod)
    pass
            
########################################################################
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

########################################################################
def main(ms_input, DirSkymodelCal, extensionSky=".skymodel", max_separation_arcmin = 1.0):

    """
    Find automatically the skymodel to use for the Calibrator 
  

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the calibrator MSs
    DirSkymodelCal : str
        Path to the skymodel file, or to the folder where the skymodels are stored    
    [extensionSky] :  str
        Default: ".skymodel"
        extension of the skymodel files
 
        
    Returns
    -------
    {'SkymodelCal':skymodelCal} : "dict"
        Path to the skymodel of the calibrator
    """    

    if os.path.isfile(DirSkymodelCal):
        print "Using the skymodel provided: " + DirSkymodelCal
        return { 'SkymodelCal' : DirSkymodelCal }
        
    elif os.path.isdir(DirSkymodelCal):
        ra, dec = grab_pointing(input2strlist_nomapfile(ms_input)[0])
        skymodelCal, skymodelName  = find_skymodel(ra, dec, DirSkymodelCal, extensionSky, max_separation_arcmin)
        return { 'SkymodelCal' : skymodelCal, 'SkymodelName': skymodelName}
    else:
        raise ValueError("find_skymodel_cal: The path \"%s\" is neither a file nor a directory!"%(DirSkymodelCal))


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Find automatically between skymodels the one to use (for the Calibrator)')
    
    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which we search the matching skymodel.')
    parser.add_argument('DirSky', type=str, 
                        help='Path to the skymodel file, or to the folder where the skymodels are stored.')
    parser.add_argument('--extsky', type=str, 
                        help='extension of the skymodel files. (default: \".skymodel\")')
        
    args = parser.parse_args()
    extensionSky='.skymodel'
    if args.extsky:
        extensionSky=args.extsky
    
    main(args.MSfile,args.DirSky, extensionSky)
    
    pass
