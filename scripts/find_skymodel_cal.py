#!/usr/bin/env python
import os,sys
import glob
import re
import pyrap.tables as pt

    
def grab_name_MS(MS):
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
    
    [nameCal] =pt.table(MS+'/OBSERVATION', readonly=True, ack=False).getcol('LOFAR_TARGET')['array']
    return nameCal

    
########################################################################        
def find_skymodel(NameSource,PathSkyMod,extensionSky=".skymodel"):
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
    

    list_skymodel=[] #"list" of corresponding skymodels for the calibrator 
    
    for skymodel in glob.glob(PathSkyMod+"/*"+extensionSky):
        if skymodel.split("/")[-1].lower().find(NameSource.lower())!=-1:
            list_skymodel.append(skymodel)
        
    # checking that one and only one skymodel is found
    nbr_skymodel=len(list_skymodel)
    if nbr_skymodel==1:
        print "The following skymodel will be used for the calibrator: "+list_skymodel[0].split("/")[-1]+ " (in "+PathSkyMod+")"
        return list_skymodel[0]
    elif nbr_skymodel==0:
        raise TypeError('find_skymodel: /!\ SKYMODEL FOR THE CALIBRATOR NOT FOUND IN '+PathSkyMod)
    elif nbr_skymodel>1:
        raise TypeError('find_skymodel: Several skymodels match. \n TODO: Please implement here which one is generally prefered \n Calibrator Skymodels are :\n'+'\n\t'.join(list_skymodel))

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
def main(ms_input,DirSkymodelCal,extensionSky=".skymodel"):

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
        print "Using the skymodel provided :" +DirSkymodelCal
        return { 'SkymodelCal' : DirSkymodelCal }
        
    elif os.path.isdir(DirSkymodelCal):
        # Getting the name of the Calibrator from the information stored in a MS
        # NB: we suppose that all the calibrators observations have this field filled in (MS/Observation, column LOFAR_TARGET)
        nameCal=grab_name_MS(input2strlist_nomapfile(ms_input)[0])
        
        # Looking in the folder DirSkymodelCal to find the corresponding skymodel
        # NB: we suppose that the name of the calibrator is included in the name of the skymodel
        skymodelCal=find_skymodel(nameCal,DirSkymodelCal,extensionSky)
        return { 'SkymodelCal' : skymodelCal }
    else:
        raise ValueError("find_skymodel_cal: The path \"%s\" is neither a file nor a directory!"%(DirSkymodelCal) )


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
    
