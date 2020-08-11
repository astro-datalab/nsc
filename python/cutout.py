#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
from astropy.time import Time
import healpy as hp
from dlnpyutils import utils, coords
import subprocess
import time
from argparse import ArgumentParser
import socket
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord

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

def cutout(im,xcen,ycen,size=51,missing=0.0):
    """ Make a cutout from an image."""

    # Make the plot
    nx,ny = im.shape
    cutout = np.zeros((size,size),float) + missing
    half = size//2
    # ximlo/hi are the *exact* indexes we will use for the input image
    # xcutlo/hi are the *exact* corresponding indexes for the cutout image
    ximlo = np.round(xcen)-half
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
    cutout[xcutlo:xcuthi+1,ycutlo:ycuthi+1] = im[ximlo:ximhi+1,yimlo:yimhi+1]

    return cutout

def meascutout(meas,size=51):
    """ Input the measurements and create cutouts. """

    expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposures.fits.gz',1)
    #expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposure.fits.gz',1)                                                              
    decam = Table.read('/home/dnidever/projects/delvered/data/decam.txt',format='ascii')


    ind1,ind2 = dln.match(expstr['base'],meas['exposure'])
    nind = len(ind)
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
    ccdnum = meas['ccdnum'][ind2]
    for i in range(nind):
        if instrument[i]=='c4d':
            dind, = np.where(decam['CCDNUM']==ccdnum[i])
            extname = decam['NAME'][dind[0]]
            im,head = fits.getdata(fluxfile[i],header=True,extname=extname)
        else:
            im,head = fits.getdata(fluxfile[i],ccdnum[i],header=True)
        # Get the cutout
        xcen = meas['x'][ind2[i]]-1   # convert to 0-indexes
        ycen = meas['y'][ind2[i]]-1
        cutim = cutout(im,xcen,ycen,51)

        # exposure_ccdnum, filter, MJD, delta_MJD, mag


if __name__ == "__main__":
    parser = ArgumentParser(description='Fix pms in healpix object catalogs.')
    parser.add_argument('exposure', type=str, nargs=1, help='Exposure name')
    parser.add_argument('ra', type=str, nargs=1, default='', help='Right Ascension (deg)')
    parser.add_argument('dec', type=str, nargs=1, default='', help='Declination (deg)')
    parser.add_argument('fov', type=str, nargs=1, default='', help='Field of View (deg)')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    exposure = args.exposure[0]
    ra = args.ra[0]
    if ra=='': ra=None
    dec = args.dec[0]
    if dec=='': dec=None
    fov = args.fov[0]
    if fov=='': fov=None

    
