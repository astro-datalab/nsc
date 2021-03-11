#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
from astropy.time import Time
from astropy.wcs import WCS
#import healpy as hp
from dlnpyutils import utils as dln, coords
import subprocess
import time
from argparse import ArgumentParser
#import socket
#from dustmaps.sfd import SFDQuery
#from astropy.coordinates import SkyCoord
from dl import queryClient as qc
import matplotlib
import matplotlib.pyplot as plt 
from glob import glob

#def cutout(exposure,ccdnum,ra=None,dec=None,fov=None):
#    """ Get an NSC cutout."""
#
#    expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposures.fits.gz',1)
#    #expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposure.fits.gz',1)
#    ind, = np.where(expstr['base']==exposure)
#    nind = len(ind)
#    if nind==0:
#        print(exposure+' NOT FOUND')
#        return
#    # Load the data
#    fluxfile = expstr['file'][ind]
#    ccdnum = expstr['ccdnum'][ind]
#    im,head = fits.getdata(fluxfile,ccdnum,header=True)

def cutout(im,xcen,ycen,size=101,missing=0.0):
    """ Make a cutout from an image."""

    # Make the plot
    ny,nx = im.shape
    cutim = np.zeros((size,size),float) + missing
    half = size//2
    # ximlo/hi are the *exact* indexes we will use for the input image
    # xcutlo/hi are the *exact* corresponding indexes for the cutout image
    ximlo = int(np.round(xcen)-half)
    xcutlo = 0
    if ximlo<0:
        xcutlo = int(np.abs(ximlo))
        ximlo = 0
    ximhi = int(np.round(xcen)+half)
    xcuthi = size-1
    if ximhi>(nx-1):
        xcuthi = int(size-1-(ximhi-(nx-1)))
        ximhi = nx-1
    yimlo = int(np.round(ycen)-half)
    ycutlo = 0
    if yimlo<0:
        ycutlo = int(np.abs(yimlo))
        yimlo = 0
    yimhi = int(np.round(ycen)+half)
    ycuthi = size-1
    if yimhi>(ny-1):
        ycuthi = int(size-1-(yimhi-(ny-1)))
        yimhi = ny-1
    cutim[ycutlo:ycuthi+1,xcutlo:xcuthi+1] = im[yimlo:yimhi+1,ximlo:ximhi+1]
    xr = [int(np.round(xcen)-half), int(np.round(xcen)+half)]
    yr = [int(np.round(ycen)-half), int(np.round(ycen)+half)]

    return cutim,xr,yr

def cutoutfig(im,meas,figfile):
    """ Make a figure of an image."""

    matplotlib.use('Agg')
    if os.path.exists(figfile): os.remove(figfile)
    fig = plt.gcf()   # get current graphics window                                                                  
    fig.clf()         # clear  

    ax = fig.subplots()
    fig.set_figheight(figsize*0.5)
    fig.set_figwidth(figsize)
    plt.imshow(im,origin='lower')
    leg = ax.legend(loc='upper left', frameon=False)
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.xlim(xr)
    #plt.ylim(yr)
    ax.annotate(r'S/N=%5.1f',xy=(np.mean(xr), yr[0]+dln.valrange(yr)*0.05),ha='center')
    plt.save(figfile)



def meascutout(meas,obj,size=101):
    """ Input the measurements and create cutouts. """

    expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposures.fits.gz',1)
    #expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposure.fits.gz',1)                                                              
    decam = Table.read('/home/dnidever/projects/delvered/data/decam.txt',format='ascii')

    objid = obj['objectid'][0]

    # Sort by MJD
    si = np.argsort(meas['mjd'])
    meas = meas[si]

    ind1,ind2 = dln.match(expstr['base'],meas['exposure'])
    nind = len(ind1)
    if nind==0:
        print('No matches')
        return
    # Sort by input meas catalog
    si = np.argsort(ind2)
    ind1 = ind1[si]
    ind2 = ind2[si]
    # Load the data
    instrument = expstr['instrument'][ind1]
    fluxfile = expstr['file'][ind1]
    fluxfile = fluxfile.replace('/net/mss1/','/mss1/')  # for thing/hulk
    ccdnum = meas['ccdnum'][ind2]
    figfiles = []
    for i in range(nind):
        if instrument[i]=='c4d':
            dind, = np.where(decam['CCDNUM']==ccdnum[i])
            extname = decam['NAME'][dind[0]]
            im,head = fits.getdata(fluxfile[i],header=True,extname=extname)
        else:
            im,head = fits.getdata(fluxfile[i],ccdnum[i],header=True)
        w = WCS(head)
        # Object X/Y position
        xobj,yobj = w.all_world2pix(obj['ra'],obj['dec'],0)
        # Get the cutout
        xcen = meas['x'][ind2[i]]-1   # convert to 0-indexes
        ycen = meas['y'][ind2[i]]-1
        smim = dln.gsmooth(im,2)
        # use the object coords for centering
        cutim,xr,yr = cutout(smim,xobj,yobj,51)
        #cutim,xr,yr = cutout(smim,xcen,ycen,51)

        # exposure_ccdnum, filter, MJD, delta_MJD, mag
        print(str(i+1)+' '+meas['exposure'][ind2[i]]+' '+str(ccdnum[i])+' '+str(meas['x'][ind2[i]])+' '+str(meas['y'][ind2[i]])+' '+str(meas['mag_auto'][ind2[i]]))

        figdir = '/net/dl2/dnidever/nsc/instcal/v3/hpm2/'
        figfile = figdir
        figfile += '%s_%04d_%s_%02d.jpg' % (str(obj['objectid'][0]),i+1,meas['exposure'][ind2[i]],ccdnum[i])
        figfiles.append(figfile)
        matplotlib.use('Agg')
        plt.rcParams.update({'font.size': 11})
        if os.path.exists(figfile): os.remove(figfile)
        fig = plt.gcf()   # get current graphics window                                                                  
        fig.clf()         # clear  

        figsize = 6.0
        ax = fig.subplots()
        fig.set_figheight(figsize*0.8)
        fig.set_figwidth(figsize)
        med = np.nanmedian(smim)
        sig = dln.mad(smim)
        bigim,xr2,yr2 = cutout(smim,xcen,ycen,151,missing=med)
        lmed = np.nanmedian(bigim)
        vmin = lmed-3*sig
        vmax = lmed+5*sig
        plt.imshow(cutim,origin='lower',aspect='auto',interpolation='none',extent=(xr[0],xr[1],yr[0],yr[1]),
                   vmin=vmin,vmax=vmax,cmap='Greys')   # viridis, Greys, jet
        plt.colorbar()
        # Meas X/Y position
        plt.scatter([xcen],[ycen],c='r',marker='+',s=100)
        # Object X/Y position
        #xobj,yobj = w.all_world2pix(obj['ra'],obj['dec'],0)
        #plt.scatter(xobj,yobj,marker='o',s=200,facecolors='none',edgecolors='y',linewidth=3)
        plt.scatter(xobj,yobj,c='y',marker='+',s=100)
        #leg = ax.legend(loc='upper left', frameon=False)
        plt.xlabel('X')
        plt.ylabel('Y')
        #plt.xlim(xr)
        #plt.ylim(yr)
        #ax.annotate(r'S/N=%5.1f',xy=(np.mean(xr), yr[0]+dln.valrange(yr)*0.05),ha='center')
        ax.annotate('%s  %02d  %s  %6.1f  ' % (meas['exposure'][ind2[i]],ccdnum[i],meas['filter'][ind2[i]],expstr['exptime'][ind1[i]]),
                    xy=(np.mean(xr), yr[0]+dln.valrange(yr)*0.05),ha='center',color='blue')
        ax.annotate('%10.5f  %10.5f  ' % (meas['mjd'][ind2[i]],meas['mjd'][ind2[i]]-np.min(meas['mjd'])),
                    xy=(xr[0]+dln.valrange(xr)*0.05, yr[1]-dln.valrange(yr)*0.05),ha='left',color='blue')
        ax.annotate('%5.2f +/- %4.2f' % (meas['mag_auto'][ind2[i]], meas['magerr_auto'][ind2[i]]),
                    xy=(xr[1]-dln.valrange(xr)*0.05, yr[1]-dln.valrange(yr)*0.05),ha='right',color='blue')
        plt.savefig(figfile)
        print('Cutout written to '+figfile)

        #import pdb; pdb.set_trace()

    # Make the animated gif
    animfile = figdir+str(objid)+'_cutouts.gif'
    print('Creating animated gif '+animfile)
    if os.path.exists(animfile): os.remove(animfile)
    ret = subprocess.run('convert -delay 100 '+figdir+str(objid)+'_*.jpg '+animfile,shell=True)
    #import pdb; pdb.set_trace()



def objcutouts(objid):
    """ Make cutouts for all the measurements of one object."""
    
    obj = qc.query(sql="select * from nsc_dr2.object where objectid='%s'" % objid,fmt='table',profile='db01')
    meas = qc.query(sql="select * from nsc_dr2.meas where objectid='%s'" % objid,fmt='table',profile='db01')
    nmeas = len(meas)
    print(str(nmeas)+' measurements for '+objid)
    meascutout(meas,obj)


if __name__ == "__main__":
    parser = ArgumentParser('Make cutout animated gif for one object.')
    parser.add_argument('objid', type=str, nargs='*', help='Exposure name')
    #parser.add_argument('fov', type=str, nargs=1, default='', help='Field of View (deg)')
    args = parser.parse_args()
    objid = args.objid

    if len(args.objid)==1:
        if os.path.exists(objid):
            objid = dln.readlines(objid)
            nobj = dln.size(objid)
    else:        
        nobj = len(args.objid)

    if type(objid) is not list: objid=[objid]

    for i in range(nobj):
        print(str(i+1)+' '+objid[i])
        try:
            objcutouts(objid[i])
        except Exception as e:
            print('Failed on '+objid[i]+' '+str(e))
            traceback.print_exc()
