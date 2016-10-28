#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, time
import siplib
import query
import re

import ssl

ssl._create_default_https_context = ssl._create_unverified_context

sip_cache = {}

def get_obsID_from_filename(path):
    filename = os.path.basename(path)
    obsReg = re.compile(r'L\d+')
    obsmatch = obsReg.match(filename)
    if (not obsmatch):
        print "get_obsID_from_filename: Could not find obs-ID in filename %s!"%(filename)
        raise ValueError("Could not find obs-ID in filename!")
    obsID = obsmatch.group(0)[1:]
    return obsID

def get_SIPs_from_obsID(obsID, projectID, verbose=False):
    if verbose:
        print "Downloading all SIPs for \"observation\" %s in project %s."%(obsID, projectID)
    dpids = query.getltadataproductids_fromlta_byprojectandsasid(projectID, obsID)
    if not dpids or len(dpids) < 1:
        print "get_SIPs_from_obsID: failed to get dataproduct-IDs for obs %s in project %s!"%(obsID, projectID)
    starttime = time.time()
    for (num,input_dpid) in enumerate(dpids):
        xml = query.getsip_fromlta_byprojectandltadataproductid('LC2_009',input_dpid)
        new_sip = siplib.Sip.from_xml(xml)
        filename = new_sip.sip.dataProduct.fileName
        sip_cache[filename] = new_sip        
        if verbose and ((num+1)%10)==0:
            duration = time.time()-starttime
            ETA = duration/(num+1)*(len(dpids)-num-1)
            print "read %d of %d SIPs, elapsed: %.2f seconds, ETA: %.2f seconds"%((num+1),len(dpids),duration,ETA)

def get_projectID_from_MSfile(path):
    import pyrap.tables as pt
    t = pt.table(path+"::OBSERVATION", readonly=True, ack=False)
    project = t.getcell('PROJECT',0)
    return project
            
def get_SIP_from_MSfile(path, verbose=False):
    filename = os.path.basename(path)
    nameparts = filename.split('.MS')
    if verbose and len(nameparts) < 2:
        print "get_SIP_from_MSfile: Filename \"%s\" does not look like a standard LOFAR MS file name. This may cause problems later on."%(filename)
    fileID = nameparts[0]
    strlen = len(fileID)
    for keyname in sip_cache.keys():
        if keyname[:strlen] == fileID:
            return sip_cache[keyname]
    if verbose:
        print "Cannot fine SIP for %s in cache."%(filename)        
    obsID = get_obsID_from_filename(filename)
    projectID = get_projectID_from_MSfile(path)
    get_SIPs_from_obsID(obsID, projectID, verbose)
    for keyname in sip_cache.keys():
        if keyname[:strlen] == fileID:
            return sip_cache[keyname]
    print "Cannot fine SIP for %s in cache after downloading SIPs for obs %s!"%(filename, obsID)
    raise ValueError("Failed to download or identify SIP!")

    
