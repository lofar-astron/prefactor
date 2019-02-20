#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Written by Bas van der Tol (vdtol@strw.leidenuniv.nl), March 2011.
# Adapted for prefactor by Alexander Drabent (alex@tls-tautenburg.de), February 2019.

#from pylab import *
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import pylab
import pyrap.quanta as qa
import pyrap.tables as pt
import pyrap.measures as pm
import sys
import numpy
import os

targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'Sun'},
            {'name' : 'Jupiter'},
            {'name' : 'Moon'}]

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
def main(ms_input, min_separation = 30, outputimage = None):

    """
    Print seperation of the phase reference center of an input MS 
  

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the calibrator MSs
        
    Returns
    -------
    0 --just for printing
    """    

    msname = input2strlist_nomapfile(ms_input)[0]
  

    # Create a measures object
    me = pm.measures()

    # Open the measurement set and the antenna and pointing table
    ms = pt.table(msname)  

    # Get the position of the first antenna and set it as reference frame
    ant_table = pt.table(msname + '/ANTENNA')  
    ant_no = 0
    pos = ant_table.getcol('POSITION')
    x = qa.quantity( pos[ant_no,0], 'm' )
    y = qa.quantity( pos[ant_no,1], 'm' )
    z = qa.quantity( pos[ant_no,2], 'm' )
    position =  me.position( 'wgs84', x, y, z )
    me.doframe( position )
    ant_table.close()

    # Get the first pointing of the first antenna
    field_table = pt.table(msname + '/FIELD')
    field_no = 0
    direction = field_table.getcol('PHASE_DIR')
    ra = direction[ ant_no, field_no, 0 ]
    dec = direction[ ant_no, field_no, 1 ]
    targets.insert(0, {'name' : 'Pointing', 'ra' : ra, 'dec' : dec})
    field_table.close()

    # Get a ordered list of unique time stamps from the measurement set
    time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
    time = time_table.getcol('TIME')
    time1 = time/3600.0
    time1 = time1 - pylab.floor(time1[0]/24)*24

    ra_qa  = qa.quantity( targets[0]['ra'], 'rad' )
    dec_qa = qa.quantity( targets[0]['dec'], 'rad' )
    pointing =  me.direction('j2000', ra_qa, dec_qa)

    separations = []

    print 'SEPERATION from A-Team sources'
    print '------------------------------'
    print 'The minimal expected distance to an A-Team source is: ' + str(min_separation) + ' deg.'
    for target in targets:
   
        t = qa.quantity(time[0], 's')
        t1 = me.epoch('utc', t)
        me.doframe(t1)

        if 'ra' in target.keys():
            ra_qa  = qa.quantity( target['ra'], 'rad' )
            dec_qa = qa.quantity( target['dec'], 'rad' )
            direction =  me.direction('j2000', ra_qa, dec_qa)
            pass
        else :
            direction =  me.direction(target['name'])
            pass
      
        separations.append(me.separation(pointing, direction))

        # Loop through all time stamps and calculate the elevation of the pointing
        el = []
        for t in time:
            t_qa = qa.quantity(t, 's')
            t1 = me.epoch('utc', t_qa)
            me.doframe(t1)
            a = me.measure(direction, 'azel')
            elevation = a['m1']
            el.append(elevation['value']/pylab.pi*180)
            pass
        
        el = numpy.array(el)
        pylab.plot(time1, el)
        
        if target['name'] != 'Pointing':
            print target['name'] + ': ' + str(me.separation(pointing, direction))
            if int(float(min_separation)) > int(float(str(me.separation(pointing, direction)).split(' deg')[0])):
                print 'WARNING: The A-Team source ' + target['name'] + ' is closer than ' + str(min_separation) + ' deg to the phase reference center. Calibration might not perform as expected.'
                pass
            pass
        
        pass
    print '------------------------------'
    pylab.title('Pointing Elevation')
    pylab.title('Elevation')
    pylab.ylabel('Elevation (deg)');
    pylab.xlabel('Time (h)');
    pylab.legend( [ target['name'] + '(' + separation.to_string() + ')' for target, separation in zip(targets, separations) ])

    if outputimage != None:
        inspection_directory = os.path.dirname(outputimage)
        if not os.path.exists(inspection_directory):
            os.makedirs(inspection_directory)
            print "Directory" , inspection_directory ,  "created." 
            pass
        else:
            print("Directory", inspection_directory,  "already exists.")
            pass
        print 'Plotting A-Team elevation to: ' + outputimage
        pylab.savefig(outputimage)
        pass
    return 0
    pass

   
########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Print seperation of the phase reference center of an input MS to an A-team source')
    
    parser.add_argument('MSfile', type=str, nargs='+', help='One (or more MSs).')
    parser.add_argument('--min_separation', type=int, default=30, help='minimal accepted distance to an A-team source on the sky in degrees (will raise a WARNING). Default: 30')
    parser.add_argument('--outputimage', type=str, default=None, help='location of the elevation plot of the A-Team sources.')
        
    args = parser.parse_args()
    
    main(args.MSfile, args.min_separation, args.outputimage)
    
    pass