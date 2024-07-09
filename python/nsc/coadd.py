#!/usr/bin/env python
#
# COADD.PY   This is used to create coadds.
#

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20170911'  # yyyymmdd


"""
    Software to create an NSC coadd.

"""


import os
import sys
import numpy as np
import shutil
#import scipy
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy.interpolate import RectBivariateSpline
#from skimage import measure, morphology
#from scipy.signal import argrelmin
#import scipy.ndimage.filters as filters
import time
import glob
import logging
import sep
import healpy as hp
from reproject import reproject_interp
import tempfile
import subprocess
import glob
from dlnpyutils import utils as dln, coords

    
def brickwcs(ra,dec,npix=3600,step=0.262):
    """ Create the WCS and header for a brick."""

    # This creates a brick WCS given the brick structure

    # Make the tiling file
    #---------------------
    # Lines with the tiling scheme first
    nx = npix
    ny = npix
    step = step / npix
    xref = nx//2
    yref = ny//2

    w = WCS()
    w.wcs.crpix = [xref+1,yref+1]
    w.wcs.cdelt = np.array([step,step])
    w.wcs.crval = [ra,dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.array_shape = (npix,npix)

    #  Make the header as well
    hdu = fits.PrimaryHDU()
    head = hdu.header

    head['NAXIS'] = 2
    head['NAXIS1'] = nx
    head['CDELT1'] = step
    head['CRPIX1'] = xref+1
    head['CRVAL1'] = ra
    head['CTYPE1'] = 'RA---TAN'
    head['NAXIS2'] = ny
    head['CDELT2'] = step
    head['CRPIX2'] = yref+1
    head['CRVAL2'] = dec
    head['CTYPE2'] = 'DEC--TAN'
    #head['BRAMIN'] = brickstr.ra1,'RA min of unique brick area'
    #head['BRAMAX'] = brickstr.ra2,'RA max of unique brick area'
    #head['BDECMIN'] = brickstr.dec1,'DEC min of unique brick area'
    #head['BDECMAX'] = brickstr.dec2,'DEC max of unique brick area'

    return w,head

def inNativeByteOrder(im):
    """ Put image in native byte order."""
    if ((im.dtype.byteorder=='<') & (sys.byteorder=='big')) | ((im.dtype.byteorder=='>') & (sys.byteorder=='little')):
        return im.byteswap(inplace=True).newbyteorder()
    else:
        return im

def doImagesOverlap(head1,head2):
    """ Do the two images overlap."""

    if isinstance(head1,WCS):
        wcs1 = head1
    else:
        wcs1 = WCS(head1)
    ny1,nx1 = wcs1.array_shape
    ra1,dec1 = wcs1.wcs_pix2world(nx1/2,ny1/2,0)
    vra1,vdec1 = wcs1.wcs_pix2world([0,nx1-1,nx1-1,0],[0,0,ny1-1,ny1-1],0)
    vlon1,vlat1 = coords.rotsphcen(vra1,vdec1,ra1,dec1,gnomic=True)

    if isinstance(head2,WCS):
        wcs2 = head2
    else:
        wcs2 = WCS(head2)
    ny2,nx2 = wcs2.array_shape
    vra2,vdec2 = wcs2.wcs_pix2world([0,nx2-1,nx2-1,0],[0,0,ny2-1,ny2-1],0)
    vlon2,vlat2 = coords.rotsphcen(vra2,vdec2,ra1,dec1,gnomic=True)

    olap = coords.doPolygonsOverlap(vlon1,vlat1,vlon2,vlat2)

    return olap

 
def adxyinterp(head,rr,dd,nstep=10,xyad=False):
    """
    Instead of transforming the entire large 2D RA/DEC 
    arrays to X/Y do a sparse grid and perform linear 
    interpolation to the full size. 
 
    Parameters
    ----------
    head : header
       The FITS header with the WCS. 
    rr : numpy array
       The 2D RA array. 
    dd  : numpy array
       The 2D DEC array. 
    nstep : int, optional
       The undersampling to use, default nstep=10. 
    xyad : bool, optional
       Perform X/Y->RA/DEC conversion instead. 
          In this case the meaning of the coordinate 
          arrays are flipped, i.e. rr/dd<->xx/yy 
 
    Returns
    -------
    xx : numpy array
       The 2D X array. 
    yy : numpy array
       The 2D Y array. 
 
    Example
    -------

    xx,yy = adxyinterp(head,ra,dec,nstep=10)
            
    By D. Nidever  Oct. 2016 
    Translated to Python by D. Nidever, April 2022
    """

    nx,ny = rr.shape
    nxs = (nx-1)//nstep + 1 
    nys = (ny-1)//nstep + 1 

    wcs = WCS(head)
     
    # Small dimensions, just transform the full grid 
    if nx <= nstep or ny <= nstep: 
        if xyad==False:
            xx,yy = wcs.world_to_pixel(SkyCoord(ra=rr,dec=dd,unit='deg'))
        else:
            coo = wcs.pixel_to_world(rr,dd)
            xx = coo.ra.deg
            yy = coo.dec.deg

    # Subsample the RA/DEC arrays 
    rrs = rr[0:nx:nstep,0:ny:nstep] 
    dds = dd[0:nx:nstep,0:ny:nstep] 
    if xyad==False:
        xxs,yys = wcs.world_to_pixel(SkyCoord(ra=rrs,dec=dds,unit='deg'))
    else:
        coos = wcs.pixel_to_world(rrs,dds)
        xxs = coos.ra.deg
        yys = coos.dec.deg
     
    # Start final arrays 
    xx = np.zeros((nx,ny),float)
    yy = np.zeros((nx,ny),float)
     
    # Use CONGRID to perform the linear interpolation 
    #   congrid normally does nx_out/nx_in scaling 
    #   if /minus_one set it does (nx_out-1)/(nx_in-1) 

    xx0,yy0 = np.arange(xxs.shape[0]),np.arange(xxs.shape[1])
    xx1 = np.arange((nxs-1)*nstep+1)/((nxs-1)*nstep)*(xxs.shape[0]-1)
    yy1 = np.arange((nys-1)*nstep+1)/((nys-1)*nstep)*(xxs.shape[1]-1)
    ixx = RectBivariateSpline(xx0,yy0,xxs,kx=1,ky=1)(xx1,yy1)
    iyy = RectBivariateSpline(xx0,yy0,yys,kx=1,ky=1)(xx1,yy1)
    xx[0:(nxs-1)*nstep+1,0:(nys-1)*nstep+1] = ixx 
    yy[0:(nxs-1)*nstep+1,0:(nys-1)*nstep+1] = iyy 
     
    # Deal with right edge 
    if (nxs-1)*nstep+1 < nx: 
        # Make a short grid in X at the right edge 
        rrs_rt = rr[nx-nstep-1:nx:nstep,0:ny:nstep] 
        dds_rt = dd[nx-nstep-1:nx:nstep,0:ny:nstep] 
        if xyad==False:
            xxs_rt,yys_rt = wcs.world_to_pixel(SkyCoord(ra=rrs_rt,dec=dds_rt,unit='deg'))
        else: 
            coo_rt = wcs.pixel_to_world(rrs_rt,dds_rt)
            xxs_rt = coo_rt.ra.deg
            yys_rt = coo_rt.dec.deg

        xx0_rt,yy0_rt = np.arange(xxs_rt.shape[0]),np.arange(xxs_rt.shape[1])
        xx1_rt = np.arange(nstep+1)/nstep*(xxs_rt.shape[0]-1)
        yy1_rt = np.arange((nys-1)*nstep+1)/((nys-1)*nstep)*(xxs_rt.shape[1]-1)
        ixx_rt = RectBivariateSpline(xx0_rt,yy0_rt,xxs_rt,kx=1,ky=1)(xx1_rt,yy1_rt)
        iyy_rt = RectBivariateSpline(xx0_rt,yy0_rt,yys_rt,kx=1,ky=1)(xx1_rt,yy1_rt)
        xx[nx-nstep-1:nx,0:(nys-1)*nstep+1] = ixx_rt
        yy[nx-nstep-1:nx,0:(nys-1)*nstep+1] = iyy_rt 

    # Deal with top edge 
    if (nys-1)*nstep+1 < ny: 
        # Make a short grid in Y at the top edge 
        rrs_tp = rr[0:nx:nstep,ny-nstep-1:ny:nstep] 
        dds_tp = dd[0:nx:nstep,ny-nstep-1:ny:nstep] 
        if xyad==False:
            xxs_tp, yys_tp = wcs.world_to_pixel(SkyCoord(ra=rrs_tp,dec=dds_tp,unit='deg'))
        else: 
            coo_tp = wcs.pixel_to_world(rrs_tp,dds_tp)
            xxs_tp = coo_tp.ra.deg
            yys_tp = coo_tp.dec.deg

        xx0_tp,yy0_tp = np.arange(xxs_tp.shape[0]),np.arange(xxs_tp.shape[1])
        xx1_tp = np.arange((nxs-1)*nstep+1)/((nxs-1)*nstep)*(xxs_tp.shape[0]-1)
        yy1_tp = np.arange(nstep+1)/nstep*(xxs_tp.shape[1]-1)
        ixx_tp = RectBivariateSpline(xx0_tp,yy0_tp,xxs_tp,kx=1,ky=1)(xx1_tp,yy1_tp)
        iyy_tp = RectBivariateSpline(xx0_tp,yy0_tp,yys_tp,kx=1,ky=1)(xx1_tp,yy1_tp)
        xx[0:(nxs-1)*nstep+1,ny-nstep-1:ny] = ixx_tp 
        yy[0:(nxs-1)*nstep+1,ny-nstep-1:ny] = iyy_tp 

    # Deal with top/right corner 
    if (nxs-1)*nstep+1 < nx and (nys-1)*nstep+1 < ny: 
        # Make a short grid in X and Y at the top-right corner 
        rrs_tr = rr[nx-nstep-1:nx:nstep,ny-nstep-1:ny:nstep] 
        dds_tr = dd[nx-nstep-1:nx:nstep,ny-nstep-1:ny:nstep] 
        if xyad==False:
            xxs_tr, yys_tr = wcs.world_to_pixel(SkyCoord(ra=rrs_tr,dec=dds_tr,unit='deg'))
        else: 
            coo_tr = wcs.pixel_to_world(rrs_tr,dds_tr)
            xxs_tr = coo_tr.ra.deg
            yys_tr = coo_tr.dec.deg

        xx0_tr,yy0_tr = np.arange(xxs_tr.shape[0]),np.arange(xxs_tr.shape[1])
        xx1_tr = np.arange(nstep+1)/nstep*(xxs_tr.shape[0]-1)
        yy1_tr = np.arange(nstep+1)/nstep*(xxs_tr.shape[1]-1)
        ixx_tr = RectBivariateSpline(xx0_tr,yy0_tr,xxs_tr,kx=1,ky=1)(yy1_tr,xx1_tr)
        iyy_tr = RectBivariateSpline(xx0_tr,yy0_tr,yys_tr,kx=1,ky=1)(yy1_tr,xx1_tr)
        xx[nx-nstep-1:nx,ny-nstep-1:ny] = ixx_tr
        yy[nx-nstep-1:nx,ny-nstep-1:ny] = iyy_tr

    return xx,yy

def image_reproject_bilinear(im,head,outhead,wtim=None,verbose=False):
    """
    Resample image using python bilinear RectBivariateSpline.

    Parameters
    ----------
    im : numpy array
       Input image.
    head : Header
       Header for the input image that containts the WCS.
    outhead : Header
       Header for the output image that contains the WCS.
    wtim : numpy array, optional
       Weight image.
    verbose : boolean, optional
       Verbose output to the screen.

    Returns
    -------
    oim : numpy array
       Resampled image.
    ohead : Header
       Header for resampled image.
    owtim : numpy array
       Resampled weight image.  Only if wtim is input.

    Example
    -------

    out = image_reproject_bilinear(im,head,outhead)

    """

    t0 = time.time()

    if isinstance(head,WCS):
        wcs = head
    else:
        wcs = WCS(head)
    if isinstance(outhead,WCS):
        outwcs = outhead
    else:
        outwcs = WCS(outhead)

    ny,nx = im.shape
    fnx,fny = outwcs.array_shape

    # Make the wtim if not input
    #   1-good, 0-bad
    if wtim is None:   
        wtim = np.ones(ny,nx,float)
        saturate = head.get('saturate')
        if satuate is not None:
            gdpix = (im < saturate)
            bdpix = (im >= saturate)
            if np.sum(bdpix) > 0:
                background = np.nanmedian(im[gdpix])
                wtim[bdpix] = 0.0
                im[bdpix] = background

    # Get image ra/dec vertices
    vcoo = wcs.pixel_to_world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
    vertices_ra = vcoo.ra.deg
    vertices_dec = vcoo.dec.deg

    # 2D RA/DEC arrays for final image
    xb = np.zeros(fny,float).reshape(-1,1) + np.arange(fnx).reshape(1,-1)
    yb = np.arange(fny).reshape(-1,1) + np.zeros(fnx).reshape(1,-1)
    bcoo = outwcs.pixel_to_world(xb,yb)
    rab = bcoo.ra.deg
    decb = bcoo.dec.deg

    # Get X/Y range for this image in the final coordinate system
    vx,vy = outwcs.world_to_pixel(SkyCoord(ra=vertices_ra,dec=vertices_dec,unit='deg'))
    xout = [np.maximum(np.floor(np.min(vx))-2, 0),
            np.minimum(np.ceil(np.max(vx))+2, fnx-2)+1]
    xout = np.array(xout).astype(int)
    nxout = int(xout[1]-xout[0])
    yout = [np.maximum(np.floor(np.min(vy))-2, 0),
            np.minimum(np.ceil(np.max(vy))+2, fny-2)+1]
    yout = np.array(yout).astype(int)
    nyout = int(yout[1]-yout[0])
    rr = rab[yout[0]:yout[1],xout[0]:xout[1]]
    dd = decb[yout[0]:yout[1],xout[0]:xout[1]]
    xx,yy = adxyinterp(head,rr,dd,nstep=10)

    # The x/y position to bilinear need to be in the original system, ~1sec
    rim = np.zeros(xx.shape,float)
    rwtim = np.zeros(xx.shape,bool)
    good = ((xx>=0) & (xx<=im.shape[1]-1) & (yy>=0) & (yy<=im.shape[0]-1))
    if np.sum(good)>0:
        rim[good] = RectBivariateSpline(np.arange(im.shape[0]),np.arange(im.shape[1]),im,kx=1,ky=1).ev(yy[good],xx[good])
        rwtim[good] = RectBivariateSpline(np.arange(im.shape[0]),np.arange(im.shape[1]),wtim,kx=1,ky=1).ev(yy[good],xx[good])
 
    # Contruct final image
    oim = np.zeros((fny,fnx),float)
    oim[yout[0]:yout[1],xout[0]:xout[1]] = rim
    owtim = np.zeros((fny,fnx),bool)
    owtim[yout[0]:yout[1],xout[0]:xout[1]] = rwtim
 

    # Contruct the final header
    ohead = head.copy()
    # Delete any previous WCS keywords
    for n in ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CTYPE1','CTYPE2','CD1_1','CD1_2','CD2_1','CD2_2']:
        if n in ohead:
            del ohead[n]
    cards = [f[0] for f in ohead.cards]
    pvnames = dln.grep(cards,'PV[0-9]_+[0-9]')
    for p in pvnames:
        del ohead[p]
    # Add the new WCS
    whead = outwcs.to_header()
    ohead.extend(whead)
    ohead['NAXIS1'] = fnx
    ohead['NAXIS2'] = fny

    dt = time.time()-t0
    if verbose:
        print('dt = %8.2f sec' % dt)

    return oim,ohead,owtim
        

def image_reproject_swarp(im,head,outhead,wtim=None,tmproot='.',verbose=False):
    """
    Reproject image onto new projection with interpolation with Swarp.

    Parameters
    ----------
    im : numpy array
       Input image.
    head : Header
       Header for the input image that containts the WCS.
    outhead : Header
       Header for the output image that contains the WCS.
    wtim : numpy array, optional
       Weight image.
    tmproot : str
       Temporary directory to use for the work.  Default is the current directory.
    verbose : boolean, optional
       Verbose output to the screen.

    Returns
    -------
    oim : numpy array
       Resampled image.
    ohead : Header
       Header for resampled image.
    owtim : numpy array
       Resampled weight image.  Only if wtim is input.

    Example
    -------

    out = image_reproject_swarp(im,head,outhead)

    """

    t0 = time.time()

    if isinstance(head,WCS):
        wcs = head
    else:
        wcs = WCS(head)
    if isinstance(outhead,WCS):
        outwcs = outhead
    else:
        outwcs = WCS(outhead)

    # Make temporary directory and CD to it
    tmpdir = os.path.abspath(tempfile.mkdtemp(prefix="swrp",dir=tmproot))
    origdir = os.path.abspath(os.curdir)
    os.chdir(tmpdir)

    # Write image to temporary file
    imfile = 'image.fits'
    fits.PrimaryHDU(im,head).writeto(imfile,overwrite=True)
    if wtim is not None:
        wtfile = 'wt.fits'
        fits.PrimaryHDU(wtim,head).writeto(wtfile,overwrite=True)

    # Use swarp to do the resampling

    # Create configuration file
    # fields to modify: IMAGEOUT_NAME, WEIGHTOUT_NAME, WEIGHT_IMAGE, CENTER, PIXEL_SCALE, IMAGE_SIZE, GAIN?
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(os.path.dirname(os.path.dirname(fil)))
    paramdir = codedir+'/params/'
    shutil.copyfile(paramdir+"swarp.config",tmpdir+"/swarp.config")
    configfile = "swarp.config"
    clines = dln.readlines(configfile)

    imoutfile = 'image.out.fits'
    if wtim is not None:
        wtoutfile = 'wt.out.fits'

    # IMAGEOUT_NAME
    ind = dln.grep(clines,'^IMAGEOUT_NAME ',index=True)[0]
    clines[ind] = 'IMAGEOUT_NAME        '+imoutfile+'     # Output filename'
    # Weight parameteres
    if wtim is not None:
        # WEIGHTOUT_NAME
        ind = dln.grep(clines,'^WEIGHTOUT_NAME ',index=True)[0]
        clines[ind] = 'WEIGHTOUT_NAME    '+wtoutfile+' # Output weight-map filename'
        # WEIGHT_TYPE
        ind = dln.grep(clines,'^WEIGHT_TYPE ',index=True)[0]
        clines[ind] = 'WEIGHT_TYPE         MAP_WEIGHT            # BACKGROUND,MAP_RMS,MAP_VARIANCE'
        # WEIGHT_IMAGE
        ind = dln.grep(clines,'^WEIGHT_IMAGE ',index=True)[0]
        clines[ind] = 'WEIGHT_IMAGE      '+wtfile+'    # Weightmap filename if suffix not used'
    else:
        # Remove WEIGHT keywords
        ind = dln.grep(clines,'^WEIGHT_TYPE ',index=True)[0]
        clines[ind] = 'WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE'
        #ind = dln.grep(clines,'^WEIGHT',index=True)
        #clines = dln.remove_indices(clines,ind)
    # CENTER
    ind = dln.grep(clines,'^CENTER ',index=True)[0]
    clines[ind] = 'CENTER         %9.6f, %9.6f  # Coordinates of the image center' % (outwcs.wcs.crval[0],outwcs.wcs.crval[1])
    # PIXEL_SCALE
    onx,ony = outwcs.array_shape
    r,d = outwcs.all_pix2world([onx//2,onx//2],[ony//2,ony//2+1],0)
    pixscale = (d[1]-d[0])*3600
    ind = dln.grep(clines,'^PIXEL_SCALE ',index=True)[0]
    clines[ind] = 'PIXEL_SCALE           %6.3f             # Pixel scale' % pixscale
    # IMAGE_SIZE
    ind = dln.grep(clines,'^IMAGE_SIZE ',index=True)[0]
    clines[ind] = 'IMAGE_SIZE             %d,%d               # Image size (0 = AUTOMATIC)' % outwcs.array_shape
    # GAIN??
    #ind = dln.grep(clines,'^GAIN_KEYWORD ',index=True)[0]
    #clines[ind] = 'GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)'

    # Write the updated configuration file
    dln.writelines(configfile,clines,overwrite=True)

    # Run swarp
    if verbose is False:
        #slogfile = "swarp.log"
        #sf = open(slogfile,'w')
        #retcode = subprocess.call(["swarp",imfile,"-c",configfile],stdout=sf,stderr=subprocess.STDOUT,shell=False)
        #sf.close()
        #slines = dln.readlines(slogfile)
        retcode = subprocess.check_output(["swarp",imfile,"-c",configfile],shell=False,stderr=subprocess.STDOUT)        
    else:
        retcode = subprocess.run(["swarp",imfile,"-c",configfile],shell=False,stderr=subprocess.STDOUT)        
                
    # Load the output file
    oim,shead = fits.getdata(imoutfile,header=True)
    # By default swarp makes it sky-right, flip
    if shead['CD1_1'] < 0:
        if verbose:
            print('Flipping swarp image in RA axis')
        oim = oim[:,::-1]
        shead['CD1_1'] = -shead['CD1_1']

    # Contruct the final header
    ohead = head.copy()
    # Delete any previous WCS keywords
    for n in ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CTYPE1','CTYPE2','CD1_1','CD1_2','CD2_1','CD2_2']:
        if n in ohead:
            del ohead[n]
    cards = [f[0] for f in ohead.cards]
    pvnames = dln.grep(cards,'PV[0-9]_+[0-9]')
    for p in pvnames:
        del ohead[p]
    # Add the new Swarp header
    ohead.extend(shead)

    out = (oim,ohead)
    if wtim is not None:
        owtim,owthead = fits.getdata(wtoutfile,header=True)
        out = (oim,ohead,owtim)

    # Delete temporary directory and files??
    tmpfiles = glob.glob('*')
    for f in tmpfiles:
        os.remove(f)
    os.rmdir(tmpdir)

    # Go back to the original direcotry
    os.chdir(origdir)

    dt = time.time()-t0
    if verbose:
        print('dt = %8.2f sec' % dt)

    return out

def image_reproject(im,head,outhead,wtim=None,kind='swarp',tmproot='.',verbose=False):
    """
    Reproject image onto new projection with interpolation.
    Wrapper for the _swarp and _bilinear functions
    """

    if kind == 'swarp':
        return image_reproject_swarp(im,head,outhead,wtim=wtim,tmproot=tmproot,verbose=verbose)
    elif kind == 'bilinear':
        return image_reproject_swarp(im,head,outhead,wtim=wtim,verbose=verbose)
    else:
        raise ValueError(kind,' not supported')
    

def image_interp(imagefile,outhead,weightfile=None,masknan=False,verbose=False):
    """
    Interpolate a single image (can be multi-extension) to the output WCS.

    Parameters
    ----------
    imagefile : str
       Filename for the flux file.  Can be multi-extension.
    outhead : Header
       Header for the output image that contains the WCs.
    weightfile : str, optional
       Filename for the weight file.  Optional.
    masknan : boolean, optional
       Flag to mask NaNs in image.  Default is False.
    verbose : boolean, optional
       Verbose output to the screen.  Default is False.

    Returns
    -------
    oim : numpy array
       Resampled flux image.
    ohead : Header
       Header for final resampled image.
    obg : numpy array
       The resampled background image.
    owtim : numpy array
       Resample weight image if the weight file was input.

    Example
    -------

    oim,ohead,obg,owtim = image_interp('image.fits',outhead)

    """

    t0 = time.time()

    if os.path.exists(imagefile) is False:
        raise ValueError(imagefile+" NOT FOUND")
    if weightfile is not None:
        if os.path.exists(weightfile) is False:
            raise ValueError(weightfile+" NOT FOUND")    

    if isinstance(outhead,WCS):
        brickwcs = outhead
    else:
        brickwcs = WCS(outhead)

    # Output vertices
    bricknx,brickny = brickwcs.array_shape
    brickra,brickdec = brickwcs.wcs_pix2world(bricknx/2,brickny/2,0)
    brickvra,brickvdec = brickwcs.wcs_pix2world([0,bricknx-1,bricknx-1,0],[0,0,brickny-1,brickny-1],0)
    brickvlon,brickvlat = coords.rotsphcen(brickvra,brickvdec,brickra,brickdec,gnomic=True)

    # How many extensions
    hdulist = fits.open(imagefile)
    nimhdu = len(hdulist)
    hdulist.close()
    if weightfile is not None:
        hdulist = fits.open(weightfile)
        nwthdu = len(hdulist)
        hdulist.close()    
        if nimhdu != nwthdu:
            raise ValueError(imagefile+' and '+weightfile+' do NOT have the same number of extensions.')

    if verbose:
        print('Flux file = '+imagefile)
        if weightfile is not None:
            print('Weight file = '+weightfile)
        print(str(nimhdu)+' extensions')

    # Open the files
    imhdulist = fits.open(imagefile)
    if weightfile is not None:
        wthdulist = fits.open(weightfile)

    # Initialize final images
    oim = np.zeros((brickny,bricknx),float)
    owtim = np.zeros((brickny,bricknx),float)
    obg = np.zeros((brickny,bricknx),float)
        
    # Loop over the HDUs
    noverlap = 0
    for i in range(nimhdu):           
        # Just get the header
        head = imhdulist[i].header
        if head['NAXIS']==0:    # no image
            continue
        wcs = WCS(head)
        nx1 = head['NAXIS1']
        ny1 = head['NAXIS2']

        # Check that it overlaps the final area
        if doImagesOverlap(brickwcs,wcs) is False:
            continue
        print('Extension '+str(i)+' overlaps')
            
        mask = None
        # Flux image
        im = imhdulist[i].data
        head = imhdulist[i].header
        im = inNativeByteOrder(im)      # for sep need native byte order  
        nx1,ny1 = im.shape
        # Weight image
        if weightfile is not None:
            wt = wthdulist[i].data
            whead = wthdulist[i].header
            wt = inNativeByteOrder(wt)
            mask = (wt<=0)

        # Mask NaNs/Infs
        if masknan is True:
            if mask is not None:
                mask = (mask==True) | ~np.isfinite(im)
            else:
                mask = ~np.isfinite(im)
            im[mask] = np.median(im[~mask])


        # Step 1. Background subtract the image
        bkg = sep.Background(im, mask=mask, bw=64, bh=64, fw=3, fh=3)
        bkg_image = bkg.back()
        im -= bkg_image
        im[mask] = 0

        # Step 2. Reproject the image
        if weightfile is not None:
            newim,newhead,newwt = image_reproject(im,head,outhead,wtim=wt)
        else:
            newim,newhead = image_reproject(im,head,outhead)
        ohead = newhead.copy()   # initial the output header
        newbg, bghead = image_reproject(bkg_image,head,outhead)

        # Step 3. Add to final images
        oim += newim
        if weightfile is not None:
            owtim += newwt
        obg += newbg

        noverlap += 1

    ohead['NAXIS1'] = bricknx
    ohead['NAXIS2'] = brickny

    dt = time.time()-t0
    if verbose:
        print('dt = %8.2f sec' % dt)

    if noverlap==0:
        print('No overlap')
        return None,None,None

    ohead['NOVERLAP'] = (noverlap,'number of chips overlap')
    if weightfile is not None:
        return oim,ohead,obg,owtim
    else:
        return oim,ohead,obg

    
def meancube(imcube,wtcube,weights=None,crreject=False,statistic='mean'):
    """ This does the actual stack of an image cube.  The images must already be background-subtracted and scaled."""
    # Weights should be normalized, e.g., sum(weights)=1
    # pixels to be masked in an image should have wtcube = 0

    ny,nx,nimages = imcube.shape

    # Unweighted
    if weights is None:
        weights = np.ones(nimages,float)/nimages
    
    # Do the weighted average
    finaltot = np.zeros((ny,nx),float)
    finaltotwt = np.zeros((ny,nx),float)        
    totvarim = np.zeros((ny,nx),float)
    for i in range(nimages):
        mask = (wtcube[:,:,i] > 0)
        var = np.zeros((ny,nx),float)  # sig^2
        var[mask] = 1/wtcube[:,:,i][mask]
        finaltot[mask] += imcube[:,:,i][mask]*weights[i]
        finaltotwt[mask] += weights[i]
        # Variance in each pixel for noise images and the scalar weights
        totvarim[mask] += weights[i]*var[mask]
    # Create the weighted average image
    finaltotwt[finaltotwt<=0] = 1
    final = finaltot/finaltotwt
    # Create final error image
    error = np.sqrt(totvarim)

    # Sum
    if statistic == 'sum':
        final *= nimages
        error *= nimages
    
    # CR rejection
    if crreject is True:
        pass
    
    return final,error


def reassemble(filename):
    """ Reassemble an image that has been split into subregions."""

    if os.path.exists(filename)==False:
        raise ValueError(filename+' not found')
    
    hdu = fits.open(filename)
    head0 = hdu[0].header
    fnx = head0['ONAXIS1']
    fny = head0['ONAXIS2']    
    image = np.zeros((fny,fnx),float)
    
    for i in range(len(hdu)-1):
        head1 = hdu[i+1].header
        im1 = hdu[i+1].data        
        x0 = head1['SUBX0']
        x1 = head1['SUBX1']
        nx = head1['SUBNX']
        y0 = head1['SUBY0']
        y1 = head1['SUBY1']
        ny = head1['SUBNY']
        # Stuff into final image
        image[y0:y1+1,x0:x1+1] = im1

    head = head0.copy()
    head['NAXIS1'] = fnx
    head['NAXIS2'] = fny    

    # Delete temporary header keys
    for c in ['SUBX0','SUBX1','SUBNX','SUBY0','SUBY1','SUBNY','ONAXIS1','NAXIS2']:
        del head[c]

    hdu.close()
        
    return image,head


def stack(meta,statistic='mean'):
    """
    Actually do the stacking/averaging of multiple images already reprojected.

    Parameters
    ----------
    meta : table
       Table that contains all of the information to perform the stacking.
         Required columns are : "flxfile", "wtfile", "bgfile", "weight"
    statistic : str, optional
       The statistic to use when combining the images: 'mean' or 'sum'.  Default
         is 'mean'.

    Returns
    -------
    final : numpy array
       Final coadded flux image.
    error : numpy array
       Final coadded error image.

    Example
    -------

    final,error = stack(meta)

    """

    nimages = len(meta)
    imagefiles = meta['flxfile']
    weightfiles = meta['wtfile']
    bgfiles = meta['bgfile']
    weights = meta['weight']
    scales = meta['scale']

    # DO NOT use the error maps for the weighted average.  Use scalar weights for each exposure.
    #  otherwise you'll get screwy images
    # Only use the error maps to generate the final error image
    
    # imagefiles/weightfiles, list of filenames

    # The images should already be background subtracted and scaled
    
    # How many subimages
    file1 = imagefiles[0]
    hdu = fits.open(file1)
    nbin = len(hdu)
    hdu.close()

    # Original image size
    head = fits.getheader(imagefiles[0],0)
    fnx = head['ONAXIS1']
    fny = head['ONAXIS2']
    
    # Get the sizes and positions of the subimages
    dtype = np.dtype([('X0',int),('X1',int),('Y0',int),('Y1',int),('NX',int),('NY',int)])
    bintab = np.zeros(nbin,dtype=dtype)
    for b in range(nbin):
        head = fits.getheader(imagefiles[0],b)
        bintab['X0'][b] = head['SUBX0']
        bintab['X1'][b] = head['SUBX1']
        bintab['Y0'][b] = head['SUBY0']
        bintab['Y1'][b] = head['SUBY1']        
        bintab['NX'][b] = head['SUBNX']
        bintab['NY'][b] = head['SUBNY']

    # Final image
    final = np.zeros((fny,fnx),float)
    error = np.zeros((fny,fnx),float)    
        
    # Loop over bins
    for b in range(nbin):
        imcube = np.zeros((bintab['NY'][b],bintab['NX'][b],nimages),float)
        wtcube = np.zeros((bintab['NY'][b],bintab['NX'][b],nimages),float)        
        # Loop over images
        for f in range(nimages):
            im,head = fits.getdata(imagefiles[f],b,header=True)
            wt,whead = fits.getdata(weightfiles[f],b,header=True)
            # Deal with NaNs
            wt[~np.isfinite(im)] = 0

            # Scale the image
            #  divide image by "scale"
            im /= scales[f]
            #  wt = 1/err^2, need to perform same operation on err as on image
            wt *= scales[f]**2

            # Stuff into the cube
            imcube[:,:,f] = im
            wtcube[:,:,f] = wt
        # Do the weighted combination
        avgim,errim = meancube(imcube,wtcube,weights=weights,statistic=statistic)
        # Stuff into final image
        final[bintab['Y0'][b]:bintab['Y1'][b]+1,bintab['X0'][b]:bintab['X1'][b]+1] = avgim
        error[bintab['Y0'][b]:bintab['Y1'][b]+1,bintab['X0'][b]:bintab['X1'][b]+1] = errim

    return final,error


def mktempfile(im,head,bg,wt,outhead,nbin=2):
    """ Break up into bins and save to temporary file """

    ny,nx = im.shape

    # Set up temporary file n ames
    tid,tfile = tempfile.mkstemp(prefix="timage",dir="/tmp")
    os.close(tid)  # close open file
    tbase = os.path.basename(tfile)
    timfile = "/tmp/"+tbase+"_flx.fits"
    twtfile = "/tmp/"+tbase+"_wt.fits"
    tbgfile = "/tmp/"+tbase+"_bg.fits"

    timhdu = fits.HDUList()
    twthdu = fits.HDUList()
    tbghdu = fits.HDUList()

    xbin = ybin = nbin
    dx = nx // xbin
    dy = ny // ybin

    # Put information in header
    # ONAXIS1, ONAXIS2, SUBX0, SUBX1, SUBY0, SUBY1, SUBNX, SUBNY
    for i in range(xbin):
        x0 = i*dx
        x1 = x0 + dx
        if i==(xbin-1): x1=nx
        for j in range(ybin):
            y0 = j*dy
            y1 = y0 + dy
            if j==(ybin-1): y1=ny
            newhead = head.copy()
            newhead['ONAXIS1'] = newhead['NAXIS1']
            newhead['ONAXIS2'] = newhead['NAXIS2']
            newhead['SUBX0'] = x0
            newhead['SUBX1'] = x1-1
            newhead['SUBNX'] = x1-x0               
            newhead['SUBY0'] = y0
            newhead['SUBY1'] = y1-1
            newhead['SUBNY'] = y1-y0
            # Flux
            subim = im[y0:y1,x0:x1].copy()
            hdu1 = fits.PrimaryHDU(subim,newhead.copy())
            timhdu.append(hdu1)
            # Background
            subbg = bg[y0:y1,x0:x1].copy()
            hdu2 = fits.PrimaryHDU(subbg,newhead.copy())
            tbghdu.append(hdu2)                
            # Weight
            subwt = wt[y0:y1,x0:x1].copy()
            hdu3 = fits.PrimaryHDU(subwt,newhead.copy())
            twthdu.append(hdu3)

    timhdu.writeto(timfile,overwrite=True)
    timhdu.close()
    tbghdu.writeto(tbgfile,overwrite=True)
    tbghdu.close()        
    twthdu.writeto(twtfile,overwrite=True)
    twthdu.close()
    if os.path.exists(tfile): os.remove(tfile)

    return timfile,tbgfile,twtfile

    
def coadd(imagefiles,weightfiles,meta,outhead,statistic='mean',
          nbin=2,outfile=None,verbose=False):
    """
    Create a coadd given a list of images.

    Parameters
    ----------
    imagefiles : list
       List of flux image filenames.
    weightfiles : list
       List of weight image filenames.
    meta : table
       Table of information for each image.  Must have "zperm",
         "exptime", and "fwhm".
    outhead : header or WCS
       Header with projection for the output image.
    statistic : str, optional
       Statistic to use for coaddition: 'mean' or 'sum'.  Default is 'mean'.
    nbin : int, optional
       Number of bins to use (in X and Y) when splitting up the
         image for the temporary files.  Default is 2.
    outfile : str, optional
       Name of output FITS filename.
    verbose : boolean, optional
       Verbose output to the screen.  Default is False.

    Returns
    -------
    final : numpy array
       Final coadded flux image.
    error : numpy array
       Final coadded error image.

    Example
    -------

    final, error = coadd(imagefiles,weightfiles,meta,outhead)

    """

    t0 = time.time()

    # meta should have zpterm, exptime, fwhm
    
    nimages = dln.size(imagefiles)

    # Figure out scales and weights
    # F_trans = 10^(-0.8*(delta_mag-0.2))    
    scales = meta['exptime'] * 10**(-0.8*(meta['zpterm']-0.2))
    meta['scale'] = scales

    # Use weight~S/N
    # S/N goes as sqrt(exptime) and in background-dominated regime S/N ~ 1/FWHM
    # so maybe something like weight ~ sqrt(scaling)/FWHM
    weights = np.sqrt(scales)/meta['fwhm']
    weights /= np.sum(weights)    # normalize
    meta['weight'] = weights

    # Loop over the images
    meta['flxfile'] = 100*' '  # add columns for temporary file names
    meta['wtfile'] = 100*' '
    meta['bgfile'] = 100*' '
    for f in range(nimages):
        if verbose:
            print(str(f+1)+' '+imagefiles[f]+' '+weightfiles[f])

        # Interpolate image
        fim, fhead, fbg, fwt = image_interp(imagefiles[f],outhead,weightfile=weightfiles[f],verbose=verbose)
        ny,nx = fim.shape

        # Break up image and save to temporary files
        tflxfile,tbgfile,twtfile = mktempfile(fim,fhead,fbg,fwt,outhead,nbin=nbin)
        meta['flxfile'][f] = tflxfile
        meta['bgfile'][f] = tbgfile
        meta['wtfile'][f] = twtfile
        
    # Stack the images
    #   this does the scaling and weighting
    final,error = stack(meta,statistic=statistic)


    
    # Delete temporary files
    for i in range(len(meta)):
        if os.path.exists(meta['flxfile'][i]): os.remove(meta['flxfile'][i])
        if os.path.exists(meta['bgfile'][i]): os.remove(meta['bgfile'][i])
        if os.path.exists(meta['wtfile'][i]): os.remove(meta['wtfile'][i])

    # Final header
    #  scales, weights, image names, mean backgrounds
    fhead = outhead.copy()
    fhead['NIMAGES'] = len(imagefiles)
    for i in range(len(meta)):
        fhead['FILE'+str(i+1)] = imagefiles[i]
        fhead['WTFL'+str(i+1)] = weightfiles[i]
        fhead['WEIT'+str(i+1)] = meta['weight'][i]
        fhead['SCAL'+str(i+1)] = meta['scale'][i]
        fhead['EXPT'+str(i+1)] = meta['exptime'][i]
        fhead['FWHM'+str(i+1)] = meta['fwhm'][i]

    # Output final file
    if outfile is not None:
        hdu = fits.HDUList()
        hdu.append(fits.PrimaryHDU(final,fhead))
        hdu[0].header['COMMENT'] = 'Flux image'
        # units, number of images, names, weights, scales, 
        hdu.append(fits.ImageHDU(error))
        hdu[1].header['COMMENT'] = 'Weight/error image'
        hdu.writeto(outfile,overwrite=True)
        hdu.close()

    dt = time.time()-t0
    if verbose:
        print('dt = %8.2f sec' % dt)

    return final,error


    # OLD NOTES
    #-give list of FITS files (possible with extensions) and mask/noise file info, WCS/output header and coadd method
    #-loop through each file and resample onto final projection (maybe write to temporary file)
    #  do same for error and mask images
    #-perform the coadd method, simplest is just sum/average or median
    #-write output file
    # weighted sum/average/median
    
    # Need to deal with background level and scale
    # calculate background value myself and maybe use photometric zero-point for the scale.
    # does that take the exposure time into account as well?  I think so.
    
    # Use the astropy reproject package


    
    # NEW NOTES
    # When creating the NSC stacks, also create a deep/multi-band stack.

    # coadd steps:
    # -loop over each exposure that overlaps the image
    #   -homogenize masks
    #   -fix "bad" pixels in wt and flux arrays
    #   -fit and subtract 2D background
    #   -resample each overlapping chip onto the final brick WCS (flux preserving)
    #      be careful when interpolating near bad pixels
    #   -save flux/wt/mask/sky to a temporary directory to use later
    # -once all of the resampling is done, figure out the relative weights for all the exposures
    #   and the final scaling
    # -use NSC zeropoints to figure out relative weighting between the exposures
    # -need to be sure that the images are SCALED properly as well (takes into exptime/zeropoint)
    # -depending on the number of overlapping exposures, break up the break into subregions which
    #   will be combined separately.  need to make sure to use consistent weighting/scaling/zero
    #   otherwise there'll be jumps at these boundaries
    # -weighted/scaled combine in each subregion taking the mask and wt image into account properly.
    # -how about outlier rejection?
    # -fpack compress the final stacked images (flux, mask, weight) similar to how the CP files are done
    # -PSF-matching???
    # -with the exposure zpterm and exptime it should be possible to calculate the SCALING of the exposure
    # -should I use that also for the weight?  should take FWHM/SEEING into account.
    #   in allframe_calcweights/getweights I use weight~S/N
    #   S/N goes as sqrt(exptime) and in background-dominated regime S/N ~ 1/FWHM
    #   so maybe something like weight ~ sqrt(scaling)/FWHM
    #   how about teff?  see eqn. 4, pg.10 of Morganson+2018
    #   teff = (FWHM_fid/FwHM)^2 * (B_fid/B) * F_trans
    #   B = background
    #   F_trans = atmospheric transmission relative to a nearly clear night, basically zeropoint
    #   F_trans = 10^(-0.8*(delta_mag-0.2))
    #
    #   teff is the ratio between the actual exposure time and the exposure time necessary to achieve
    #   the same signal-to-noise for point sources observed in nominal conditions. An exposure taken
    #   under “fiducial” conditions has teff = 1.
    #
    #   see section 6.3 for how DES performs the coadds/stacks, they create "chi-mean" coadds of r/i/z
    #   also see Drlica-Wagner+2018, DES Y1 Cosmoogy results
    #
    # -ccdproc.combine
    # https://ccdproc.readthedocs.io/en/latest/image_combination.html#id1
    # -reproject.reproject_interp
    # https://reproject.readthedocs.io/en/stable/
    #   need to check the reprojection algorithms, doesn't support lanczos
    # reproject also has a coadd/mosaicking sofware and can take weights, and can do background matching
    # from reproject.mosaicking import reproject_and_coadd
    # could also use swarp
    # https://www.astromatic.net/pubsvn/software/swarp/trunk/doc/swarp.pdf
