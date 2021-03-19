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
import traceback
from reproject import reproject_interp
import shutil

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

def getfitsext(filename,extname,header=True):
    """ Get FITS extension by name."""

    try:
        out = fits.getdata(filename,header=header,extname=extname)        

    # Sometimes getdata throws an error and says it can't find the extension but it is there
    except:
        hdu = fits.open(filename)
        extnames = []
        nhdu = len(hdu)
        for i in range(nhdu):
            extnames.append(hdu[i].header.get('extname'))
        gd, = np.where(np.char.array(extnames).astype(str)==extname)
        ngd = len(gd)
        if ngd==0:
            raise ValueError('Extension '+extname+' not found in '+filename)
        ext = gd[0]
        out = fits.getdata(filename,ext,header=header)

    return out

def meascutout(meas,obj,size=10,outdir='./'):
    """ Input the measurements and create cutouts. """

    expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposures.fits.gz',1)
    #expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposure.fits.gz',1)                                                              
    decam = Table.read('/home/dnidever/projects/delvered/data/decam.txt',format='ascii')

    objid = obj['id'][0]

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

    # Create the reference WCS
    wref = WCS(naxis=2)
    pixscale = 0.26   # DECam, "/pix
    npix = round(size/pixscale)
    if npix % 2 ==0:  # must be odd
        npix += 1
    hpix = npix//2  # center of image
    wref.wcs.ctype = ['RA---TAN','DEC--TAN']
    wref.wcs.crval = [obj['ra'][0],obj['dec'][0]]
    wref.wcs.crpix = [npix//2,npix//2]
    wref.wcs.cd = np.array([[pixscale/3600.0, 0.0],[0.0, pixscale/3600]])
    wref.array_shape = (npix,npix)
    refheader = wref.to_header()
    refheader['NAXIS'] = 2
    refheader['NAXIS1'] = npix
    refheader['NAXIS2'] = npix

    # Load the data
    instrument = expstr['instrument'][ind1]
    fluxfile = expstr['file'][ind1]
    fluxfile = fluxfile.replace('/net/mss1/','/mss1/')  # for thing/hulk
    maskfile = expstr['maskfile'][ind1]
    maskfile = maskfile.replace('/net/mss1/','/mss1/')  # for thing/hulk
    ccdnum = meas['ccdnum'][ind2]
    figfiles = []
    xmeas = []
    ymeas = []
    cutimarr = np.zeros((npix,npix,nind),float)
    for i in range(nind):
    #for i in range(3):
        try:
            if instrument[i]=='c4d':
                dind, = np.where(decam['CCDNUM']==ccdnum[i])
                extname = decam['NAME'][dind[0]]
                im,head = getfitsext(fluxfile[i],extname,header=True)
                mim,mhead = getfitsext(maskfile[i],extname,header=True)
                #im,head = fits.getdata(fluxfile[i],header=True,extname=extname)
                #mim,mhead = fits.getdata(maskfile[i],header=True,extname=extname)
            else:
                im,head = fits.getdata(fluxfile[i],ccdnum[i],header=True)
                mim,mhead = fits.getdata(maskfile[i],ccdnum[i],header=True)
        except:
            print('error')
            import pdb; pdb.set_trace()

        # Get chip-level information
        exposure = os.path.basename(fluxfile[i])[0:-8]  # remove fits.fz
        chres = qc.query(sql="select * from nsc_dr2.chip where exposure='"+exposure+"' and ccdnum="+str(ccdnum[i]),fmt='table')

        w = WCS(head)
        # RA/DEC correction for the object
        lon = obj['ra'][0]-chres['ra'][0]
        lat = obj['dec'][0]-chres['dec'][0]
        racorr = chres['ra_coef1'][0] + chres['ra_coef2'][0]*lon + chres['ra_coef3'][0]*lon*lat + chres['ra_coef4'][0]*lat
        deccorr = chres['dec_coef1'][0] + chres['dec_coef2'][0]*lon + chres['dec_coef3'][0]*lon*lat + chres['dec_coef4'][0]*lat
        # apply these offsets to the header WCS CRVAL
        #w.wcs.crval += [racorr,deccorr]
        #head['CRVAL1'] += racorr
        #head['CRVAL2'] += deccorr
        print(racorr,deccorr)

        # Object X/Y position
        xobj,yobj = w.all_world2pix(obj['ra'],obj['dec'],0)
        # Get the cutout
        xcen = meas['x'][ind2[i]]-1   # convert to 0-indexes
        ycen = meas['y'][ind2[i]]-1
        smim = dln.gsmooth(im,2)
        # use the object coords for centering
        #cutim,xr,yr = cutout(smim,xobj,yobj,size)

        # Mask the bad pixels
        badmask = (mim>0)
        im[badmask] = np.nanmedian(im[~badmask])

        # Create a common TAN WCS that each image gets interpoled onto!!!
        #hdu1 = fits.open(fluxfile[i],extname=extname)
        smim1 = dln.gsmooth(im,1.5)
        hdu = fits.PrimaryHDU(smim1,head)
        cutim, footprint = reproject_interp(hdu, refheader, order='bicubic')  # biquadratic
        cutim[footprint==0] = np.nanmedian(im[~badmask])   # set out-of-bounds to background
        #xr = [0,npix-1]
        #yr = [0,npix-1]
        xr = [-hpix*pixscale,hpix*pixscale]
        yr = [-hpix*pixscale,hpix*pixscale]

        # exposure_ccdnum, filter, MJD, delta_MJD, mag
        print(str(i+1)+' '+meas['exposure'][ind2[i]]+' '+str(ccdnum[i])+' '+str(meas['x'][ind2[i]])+' '+str(meas['y'][ind2[i]])+' '+str(meas['mag_auto'][ind2[i]]))

        #figdir = '/net/dl2/dnidever/nsc/instcal/v3/hpm2/cutouts/'
        figfile = outdir
        figfile += '%s_%04d_%s_%02d.jpg' % (str(obj['id'][0]),i+1,meas['exposure'][ind2[i]],ccdnum[i])
        figfiles.append(figfile)
        matplotlib.use('Agg')
        plt.rc('font',size=15)
        plt.rc('axes',titlesize=20)
        plt.rc('axes',labelsize=20)
        plt.rc('xtick',labelsize=20)
        plt.rc('ytick',labelsize=20)
        #plt.rcParams.update({'font.size': 15})
        #plt.rcParams.update({'axes.size': 20})
        #plt.rcParams.update({'xtick.size': 20})
        #plt.rcParams.update({'ytick.size': 20})
        if os.path.exists(figfile): os.remove(figfile)
        fig = plt.gcf()   # get current graphics window                                                                  
        fig.clf()         # clear  

        gskw = dict(width_ratios=[30,1])
        fig, ax = plt.subplots(ncols=2, nrows=1, gridspec_kw=gskw)

        figsize = 8.0 #6.0
        figheight = 8.0
        figwidth = 9.0
        #ax = fig.subplots()  # projection=wcs
        #fig.set_figheight(figsize*0.8)
        fig.set_figheight(figheight)
        fig.set_figwidth(figwidth)
        med = np.nanmedian(smim)
        sig = dln.mad(smim)
        bigim,xr2,yr2 = cutout(smim,xcen,ycen,151,missing=med)
        lmed = np.nanmedian(bigim)

        # Get the flux of the object and scale each image to the same height
        #meas.mag_aper1 = cat1.mag_aper[0] + 2.5*alog10(exptime) + chstr[i].zpterm
        #cmag = mag_auto + 2.5*alog10(exptime) + zpterm
        instmag = meas['mag_auto'][ind2[i]] - 2.5*np.log10(chres['exptime'][0]) - chres['zpterm'][0]
        #mag = -2.5*log(flux)+25.0
        instflux = 10**((25.0-instmag)/2.5)
        print('flux = '+str(instflux))
        # Get height of object
        #  flux of 2D Gaussian is ~2*pi*height*sigma^2
        pixscale1 = np.max(np.abs(w.wcs.cd))*3600
        fwhm = chres['fwhm'][0]/pixscale1
        instheight = instflux/(2*3.14*(fwhm/2.35)**2)
        print('height = '+str(instheight))
        # Scale the images to the flux level of the first image
        cutim -= lmed
        if i==0:
            instflux0 = instflux.copy()
            instheight0 = instheight.copy()
        else:
            scale = instflux0/instflux
            #scale = instheight0/instheight
            cutim *= scale
            print('scaling image by '+str(scale))


        #vmin = lmed-8*sig  # 3*sig
        #vmax = lmed+12*sig  # 5*sig
        if i==0:
            vmin = -8*sig  # 3*sig   
            #vmax = 12*sig  # 5*sig   
            vmax = 0.5*instheight # 0.5
            vmin0 = vmin
            vmax0 = vmax
        else:
            vmin = vmin0
            vmax = vmax0

        print('vmin = '+str(vmin))
        print('vmax = '+str(vmax))

        cutimarr[:,:,i] = cutim.copy()

        ax[0].imshow(cutim,origin='lower',aspect='auto',interpolation='none',extent=(xr[0],xr[1],yr[0],yr[1]),
                     vmin=vmin,vmax=vmax,cmap='viridis')   # viridis, Greys, jet
        #plt.imshow(cutim,origin='lower',aspect='auto',interpolation='none',
        #           vmin=vmin,vmax=vmax,cmap='viridis')   # viridis, Greys, jet
        #plt.colorbar()

        # show one vertical, one horizontal line pointing to the center but offset
        # then a small dot on the meas position
        # 13, 8
        ax[0].plot(np.array([0,0]),np.array([-0.066*npix,0.066*npix])*pixscale,c='white',alpha=0.7)
        ax[0].plot(np.array([-0.066*npix,0.066*npix])*pixscale,np.array([0,0]),c='white',alpha=0.7)

        # Meas X/Y position
        xmeas1,ymeas1 = wref.all_world2pix(meas['ra'][ind2[i]],meas['dec'][ind2[i]],0)
        xmeas.append(xmeas1)
        ymeas.append(ymeas1)
        ax[0].scatter([(xmeas1-hpix)*pixscale],[(ymeas1-hpix)*pixscale],c='r',marker='+',s=20)
        #plt.scatter([xmeas],[ymeas],c='r',marker='+',s=100)
        #plt.scatter([xcen],[ycen],c='r',marker='+',s=100)
        # Object X/Y position
        #xobj,yobj = w.all_world2pix(obj['ra'],obj['dec'],0)
        xobj,yobj = wref.all_world2pix(obj['ra'],obj['dec'],0)
        #plt.scatter(xobj,yobj,marker='o',s=200,facecolors='none',edgecolors='y',linewidth=3)
        #plt.scatter(xobj,yobj,c='y',marker='+',s=100)
        #leg = ax.legend(loc='upper left', frameon=False)
        ax[0].set_xlabel(r'$\Delta$ RA (arcsec)')
        ax[0].set_ylabel(r'$\Delta$ DEC (arcsec)')
        ax[0].set_xlim((xr[1],xr[0]))  # sky right
        ax[0].set_ylim(yr)
        #plt.xlabel('X')
        #plt.ylabel('Y')
        #plt.xlim(xr)
        #plt.ylim(yr)
        #ax.annotate(r'S/N=%5.1f',xy=(np.mean(xr), yr[0]+dln.valrange(yr)*0.05),ha='center')
        co = 'white' #'lightgray' # blue
        ax[0].annotate('%s  %02d  %s  %6.1f  ' % (meas['exposure'][ind2[i]],ccdnum[i],meas['filter'][ind2[i]],expstr['exptime'][ind1[i]]),
                       xy=(np.mean(xr), yr[0]+dln.valrange(yr)*0.05),ha='center',color=co)
        ax[0].annotate('%10.2f    $\Delta$t=%7.2f  ' % (meas['mjd'][ind2[i]],meas['mjd'][ind2[i]]-np.min(meas['mjd'])),
                       xy=(xr[0]+dln.valrange(xr)*0.05, yr[1]-dln.valrange(yr)*0.05),ha='left',color=co)
        ax[0].annotate('%s = %5.2f +/- %4.2f' % (meas['filter'][ind2[i]], meas['mag_auto'][ind2[i]], meas['magerr_auto'][ind2[i]]),
                       xy=(xr[1]-dln.valrange(xr)*0.05, yr[1]-dln.valrange(yr)*0.05),ha='right',color=co)

        # Progress bar
        frameratio = (i+1)/float(nind)
        timeratio = (meas['mjd'][ind2[i]]-np.min(meas['mjd']))/dln.valrange(meas['mjd'])
        #ratio = frameratio
        ratio = timeratio
        print('ratio = '+str(100*ratio))
        barim = np.zeros((100,100),int)
        ind = dln.limit(int(round(ratio*100)),1,99)
        barim[:,0:ind] = 1
        ax[1].imshow(barim.T,origin='lower',aspect='auto',cmap='Greys')
        ax[1].set_xlabel('%7.1f \n days' % (meas['mjd'][ind2[i]]-np.min(meas['mjd'])))
        #ax[1].set_xlabel('%d/%d' % (i+1,nind))
        ax[1].set_title('%d/%d' % (i+1,nind))
        ax[1].axes.xaxis.set_ticks([])
        #ax[1].axes.xaxis.set_visible(False)
        ax[1].axes.yaxis.set_visible(False)
        #ax[1].axis('off')
        right_side = ax[1].spines['right']
        right_side.set_visible(False)
        left_side = ax[1].spines['left']
        left_side.set_visible(False)
        top_side = ax[1].spines['top']
        top_side.set_visible(False)

        plt.savefig(figfile)
        print('Cutout written to '+figfile)

        #import pdb; pdb.set_trace()

    avgim = np.sum(cutimarr,axis=2)/nind
    avgim *= instheight0/np.max(avgim)
    medim = np.median(cutimarr,axis=2)

    # Make a single blank file at the end so you know it looped
    figfile = outdir
    figfile += '%s_%04d_%s.jpg' % (str(obj['id'][0]),i+2,'path')
    figfiles.append(figfile)
    matplotlib.use('Agg')
    if os.path.exists(figfile): os.remove(figfile)
    fig = plt.gcf()   # get current graphics window                                                                  
    fig.clf()         # clear  
    gskw = dict(width_ratios=[30,1])
    fig, ax = plt.subplots(ncols=2, nrows=1, gridspec_kw=gskw)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    ax[0].imshow(avgim,origin='lower',aspect='auto',interpolation='none',extent=(xr[0],xr[1],yr[0],yr[1]),
                 vmin=vmin,vmax=vmax,cmap='viridis')   # viridis, Greys, jet
    ax[0].plot(np.array([0,0]),np.array([-0.066*npix,0.066*npix])*pixscale,c='white',alpha=0.7,zorder=1)
    ax[0].plot(np.array([-0.066*npix,0.066*npix])*pixscale,np.array([0,0]),c='white',alpha=0.7,zorder=1)
    xmeas = np.array(xmeas)
    ymeas = np.array(ymeas)
    ax[0].plot((xmeas-hpix)*pixscale,(ymeas-hpix)*pixscale,c='r')
    #plt.scatter((xmeas-hpix)*pixscale,(ymeas-hpix)*pixscale,c='r',marker='+',s=30)
    ax[0].set_xlabel(r'$\Delta$ RA (arcsec)')
    ax[0].set_ylabel(r'$\Delta$ DEC (arcsec)')
    ax[0].set_xlim((xr[1],xr[0]))   # sky -right
    ax[0].set_ylim(yr)
    ax[1].axis('off')
    plt.savefig(figfile)
    # Make four copies
    for j in np.arange(2,11):
        #pathfile = figfile.replace('path1','path'+str(j))
        pathfile = figfile.replace('%04d'%(i+2),'%04d'%(i+1+j))
        if os.path.exists(pathfile): os.remove(pathfile)
        shutil.copyfile(figfile,pathfile)
        figfiles.append(pathfile)

    # Make the animated gif
    animfile = outdir+str(objid)+'_cutouts.gif'
    print('Creating animated gif '+animfile)
    if os.path.exists(animfile): os.remove(animfile)
    #ret = subprocess.run('convert -delay 100 '+figdir+str(objid)+'_*.jpg '+animfile,shell=True)
    ret = subprocess.run('convert -delay 20 '+' '.join(figfiles)+' '+animfile,shell=True)
    #import pdb; pdb.set_trace()
    dln.remove(figfiles)


def objcutouts(objid,size=40,outdir='./'):
    """ Make cutouts for all the measurements of one object."""
    
    obj = qc.query(sql="select * from nsc_dr2.object where id='%s'" % objid,fmt='table',profile='db01')
    meas = qc.query(sql="select * from nsc_dr2.meas where objectid='%s'" % objid,fmt='table',profile='db01')
    nmeas = len(meas)
    print(str(nmeas)+' measurements for '+objid)
    meascutout(meas,obj,size=size,outdir=outdir)


if __name__ == "__main__":
    parser = ArgumentParser('Make cutout animated gif for one object.')
    parser.add_argument('objid', type=str, nargs='*', help='Exposure name')
    parser.add_argument('--fov', type=str, nargs=1, default=40, help='Field of View (arcsec)')
    parser.add_argument('--outdir', type=str, nargs=1, default='', help='Output directory')
    args = parser.parse_args()
    objid = args.objid
    nobj = len(args.objid)
    fov = float(dln.first_el(args.fov))
    outdir = dln.first_el(args.outdir)
    if outdir=='':
        outdir = '/net/dl2/dnidever/nsc/instcal/v3/hpm2/cutouts2/'

    if len(args.objid)==1:
        if os.path.exists(objid[0]):
            objid = dln.readlines(objid)
            nobj = dln.size(objid)
    else:        
        nobj = len(args.objid)

    if type(objid) is not list: objid=[objid]

    for i in range(nobj):
        print(str(i+1)+' '+objid[i])
        try:
            objcutouts(objid[i],size=fov,outdir=outdir)
        except Exception as e:
            print('Failed on '+objid[i]+' '+str(e))
            traceback.print_exc()
