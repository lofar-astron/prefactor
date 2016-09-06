#!/usr/bin/env python

# ToDo: make this backward-compatible to pyrap
from casacore.measures import measures
import casacore.tables as pt
import casacore.quanta as qa
import numpy as np
import os, sys


def main(msnames, sourcelist="CasA CygA VirA TauA", elevlimit='10deg'):
    def get_mean_time(ms):
        obstable = pt.table(ms+'::OBSERVATION', ack=False)
        meantime = np.mean(obstable.col('TIME_RANGE')[0])
        return meantime

    if msnames[0] == '[':
        mslist = msnames.lstrip('[').rstrip(']').split(',')
        for ms in mslist:
            ms = ms.strip("\'\"")
            if os.path.isdir(ms):
                msname = ms
                break
    else:
        msname = msnames
    print 'gen_demixing_sources.py Using MS:',msname
    reftime = int(round(get_mean_time(msname)))
    dm = measures()
    #position of CS004
    dm.doframe(dm.position('ITRF','3.82659e+06m', '460866.m', '5.0649e+06m'))
    dm.doframe(dm.epoch('utc', qa.quantity(reftime, 's') ))
    ele_rad = qa.quantity(elevlimit.strip(' []\'\"')).get_value('rad')
    outstring = "["
    sources = [source.strip(' []\'\"') for source in sourcelist.split()]
    for source in sources:
        if source == 'CasA':
            CasA_dir = dm.direction('J2000', '23h23m27.9s','+58d48m42s')
            azeldir =  dm.measure(CasA_dir,'AZEL')
            if qa.quantity(azeldir['m1']).get_value('rad') > ele_rad :
                if len(outstring) > 1:
                    outstring +=','
                outstring += 'CasA'
        elif source == 'CygA':
            CygA_dir = dm.direction('J2000', '19h59m28.3s','+40d44m02s')
            azeldir =  dm.measure(CygA_dir,'AZEL')
            if qa.quantity(azeldir['m1']).get_value('rad') > ele_rad :
                if len(outstring) > 1:
                    outstring +=','
                outstring += 'CygA'
        elif source == 'VirA':
            VirA_dir = dm.direction('J2000', '12h30m49.4233s', '+12d23m28.043s')
            azeldir =  dm.measure(VirA_dir,'AZEL')
            if qa.quantity(azeldir['m1']).get_value('rad') > ele_rad :
                if len(outstring) > 1:
                    outstring +=','
                outstring += 'VirA'
        elif source == 'TauA':
            TauA_dir = dm.direction('J2000', '05h34m32.0s','+22d00m52s')
            azeldir =  dm.measure(TauA_dir,'AZEL')
            if qa.quantity(azeldir['m1']).get_value('rad') > ele_rad :
                if len(outstring) > 1:
                    outstring +=','
                outstring += 'TauA'
        elif source == 'HerA':
            HerA_dir = dm.direction('J2000', '16h51m08.1s','+04d59m33s')
            azeldir =  dm.measure(HerA_dir,'AZEL')
            if qa.quantity(azeldir['m1']).get_value('rad') > ele_rad :
                if len(outstring) > 1:
                    outstring +=','
                outstring += 'HerA'
        else:
            raise ValueError('gen_demixing_sources.py: Unknown Source '+source+'\n'
                             'Supported sources: CasA, CygA, VirA, TauA, HerA')
    outstring += "]"
    outrec = { 'demixsources' : outstring }
    #outrec['debugout'] = outstring+' ; '+qa.quantity(reftime, 's').formatted('YMD')
    return outrec
