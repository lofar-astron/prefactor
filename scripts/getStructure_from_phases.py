#!/usr/bin/env python

##########################################################################
#
# Written my Maaijke Mevius -- see Radio Science publication DOI: 10.1002/2016RS006028
#
##########################################################################
# History
# 20/07/2016 Outlier rejection and argparse added -- Tim Shimwell

import matplotlib as mpl
mpl.use("Agg")
import casacore.tables as pt
import lofar.parmdb as pd
from pylab import *
import scipy
from scipy.optimize import leastsq
import argparse

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def model(t, coeffs):
    coeffs[0] = abs(coeffs[0])
    return coeffs[0] + (t**coeffs[2])/coeffs[1]
def residuals(coeffs, y, t):
    coeffs[0] = abs(coeffs[0])
    return y - model(t, coeffs)

def getStructure(pdname,antennas,nr_grid=1,doplot=True,outbasename='ionosphere',dofilter=True):
    outfile = open(outbasename +'_structure.txt','w')
    myt=pd.parmdb(pdname)
    val=myt.getValuesGrid("*")
    stat=pt.table(antennas)
    names=stat.getcol('NAME')
    pos=stat.getcol('POSITION')
    stations={}
    for i,j in zip(names,pos):
        stations[i]=j
    phx=[]
    phy=[]
    allposx=[]
    allposy=[]
    nrtimes=val[val.keys()[0]]['times'].shape[0]
    timestep=int(nrtimes/nr_grid)
    for i in sorted(val.keys()):
        if not "CS" in i:
            continue
        #if 'CS201' in i:
        #    continue
        #if 'CS030' in i:
        #    continue
        if "0:0" in i:
            phx.append(val[i]['values'].flatten())
            allposx.append(stations[i.split(":")[-1]])
        if "1:1" in i:
            phy.append(val[i]['values'].flatten())
            allposy.append(stations[i.split(":")[-1]])

    phx=np.array(phx)
    phy=np.array(phy)
    allposx=np.array(allposx)
    D=allposx[:,np.newaxis,:]-allposx[np.newaxis,:]
    D2=np.sqrt(np.sum(D**2,axis=-1))
    allposy=np.array(allposy)
    Dy=allposy[:,np.newaxis,:]-allposy[np.newaxis,:]
    D2y=np.sqrt(np.sum(Dy**2,axis=-1))
    S0s=[]
    betas=[]
    for itime in xrange(nr_grid+1):
        tm=[0,1e9]
        if itime<nr_grid:
            tm[0]=itime*timestep
            tm[1]=tm[0]+timestep
        dphx=phx[:,np.newaxis,tm[0]:tm[1]]-phx[np.newaxis,:,tm[0]:tm[1]]
        dphy=phy[:,np.newaxis,tm[0]:tm[1]]-phy[np.newaxis,:,tm[0]:tm[1]]
        dphx=np.unwrap(np.remainder(dphx-dphx[:,:,0][:,:,np.newaxis]+np.pi,2*np.pi))
        dvarx=np.var(dphx,axis=-1)
        dphy=np.unwrap(np.remainder(dphy-dphy[:,:,0][:,:,np.newaxis]+np.pi,2*np.pi))
        dvary=np.var(dphy,axis=-1)
        myselect=np.logical_and(D2>0,np.logical_and(np.any(np.logical_and(dvarx>1e-7,dvarx<.1),axis=0)[np.newaxis],np.any(np.logical_and(dvarx>1e-7,dvarx<.1),axis=1)[:,np.newaxis]))
        if dofilter:
            x=D2[myselect]
            y=dvarx[myselect]
            # Seems like real values dont occur above 1.0
            flagselect = np.where(y > 1.0)
            xplot = np.delete(x,flagselect) 
            yplot = np.delete(y,flagselect)
            
            bins = np.logspace(np.log10(np.min(xplot)),np.log10(np.max(xplot)),10)
            binys = []
            binxs = []
            for i in range(0,len(bins)-1):
                binmin = bins[i]
                binmax = bins[i+1]
                inbin = np.intersect1d(np.where(xplot > binmin),np.where(xplot < binmax))
                binxs.append((binmin+binmax)/2.0)
                #binys.append(np.percentile(y[inbin],90))
                binys.append(np.median(yplot[inbin])+10*mad(yplot[inbin]))
            x0 = [0.1,1.0,1.0]
            xfit, flag = scipy.optimize.leastsq(residuals, x0, args=(binys,binxs))
            flagselect = np.where(y > model(x,xfit))
            print xfit,'fitting'
        if doplot:
            x=D2[myselect]
            y=dvarx[myselect]
            subplot(2,1,1)
            scatter(x,y,color='b')
            if dofilter:
                y2 = y[flagselect]
                x2 = x[flagselect]
                scatter(x2,y2,color='r')
                scatter(x,model(x,xfit),color='gray',linestyle='-')#,'gray',linestyle='-',linewidth=3)
        x=D2[np.logical_and(D2>1e3,myselect)]
        y=dvarx[np.logical_and(D2>1e3,myselect)]
        if dofilter:
            flagselect = np.where(y > model(x,xfit))
            x = np.delete(x,flagselect)
            y = np.delete(y,flagselect)
        A=np.ones((2,x.shape[0]),dtype=float)
        A[1,:]=np.log10(x)
        par=np.dot(np.linalg.inv(np.dot(A,A.T)),np.dot(A,np.log10(y)))
        S0=10**(-1*np.array(par)[0]/np.array(par)[1])
        if doplot:
            plot(x,10**(par[0]+par[1]*np.log10(x)),'r-')
            xscale("log")
            yscale("log")
            #ylim(1e-3,2)
            xlim(30,4000)
        myselect=np.logical_and(D2y>0,np.logical_and(np.any(np.logical_and(dvary>1e-7,dvary<.1),axis=0)[np.newaxis],np.any(np.logical_and(dvary>1e-7,dvary<.1),axis=1)[:,np.newaxis]))
        if dofilter:
            x=D2y[myselect]
            y=dvary[myselect]
            # seems like real values dont occur above 1.0
            flagselect = np.where(y > 1.0)
            xplot = np.delete(x,flagselect) 
            yplot = np.delete(y,flagselect)
            bins = np.logspace(np.log10(np.min(xplot)),np.log10(np.max(xplot)),20)
            binys = []
            binxs = []
            for i in range(0,len(bins)-1):
                binmin = bins[i]
                binmax = bins[i+1]
                inbin = np.intersect1d(np.where(xplot > binmin),np.where(xplot < binmax))
                binxs.append((binmin+binmax)/2.0)
                #binys.append(np.percentile(y[inbin],90))
                binys.append(np.median(yplot[inbin])+10*mad(yplot[inbin]))
            x0 = [0.1,1.0,1.0]
            xfit, flag = scipy.optimize.leastsq(residuals, x0, args=(binys,binxs))
            flagselect = np.where(y > model(x,xfit))
            print xfit,'fitting'
        if doplot:
            x=D2y[myselect]
            y=dvary[myselect]
            subplot(2,1,2)
            scatter(x,y,color='b')
            if dofilter:
                y2 = y[flagselect]
                x2 = x[flagselect]
                scatter(x2,y2,color='r')
                scatter(x,model(x,xfit),color='gray',linestyle='-')#,'gray',linestyle='-',linewidth=3)
        x=D2y[np.logical_and(D2y>1e3,myselect)]
        y=dvary[np.logical_and(D2y>1e3,myselect)]
        if dofilter:
            flagselect = np.where(y > model(x,xfit))
            x = np.delete(x,flagselect)
            y = np.delete(y,flagselect)
        A=np.ones((2,x.shape[0]),dtype=float)
        A[1,:]=np.log10(x)
        pary=np.dot(np.linalg.inv(np.dot(A,A.T)),np.dot(A,np.log10(y)))
        S0y=10**(-1*np.array(pary)[0]/np.array(pary)[1])
        if doplot:
            plot(x,10**(pary[0]+pary[1]*np.log10(x)),'r-')
            xscale("log")
            yscale("log")
            #ylim(1e-3,2)
            xlim(30,4000)
            savefig('%s_structure.png'%outbasename)
            close()
            cla()
        S0s.append([S0,S0y])
        betas.append([par[1],pary[1]])
        outfile.write('S0s ****%s**** %s beta %s %s\n'%(S0,S0y,par[1],pary[1]))
    return S0s,betas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('parmdb',help='Name of solution parmdb')
    parser.add_argument('anttab',help='Name of antenna table')
    parser.add_argument('outbasename',help='Outfile base name')
    parser.add_argument('-p',dest='doplot',default=True,help='Do plotting')
    parser.add_argument('-f',dest='dofilter',default=True,help='Attempt to filter out odd solutions')
    parser.add_argument('-nr_grid',dest='nr_grid',default=1,help='Time step')

    args = parser.parse_args()
    parmdb = args.parmdb
    anttab = args.anttab
    dofilter = args.dofilter
    doplot = args.doplot
    outbasename = args.outbasename
    nr_grid = args.nr_grid
    print parmdb,anttab,dofilter,doplot,outbasename
    
    getStructure(parmdb,anttab,nr_grid,doplot,outbasename,dofilter)
