#!/usr/bin/env python
#
# NSC_COADD.PY   This is used to create coadds for the NSC.
#

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20170911'  # yyyymmdd


"""
    Software to create an NSC coadd.

"""


import os
import sys
import numpy as np
#import scipy
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
#import photutils
#from skimage import measure, morphology
#from scipy.cluster import vq
#import gaps
#import matplotlib.pyplot as plt
#import pylab
#from scipy.signal import argrelmin
#import scipy.ndimage.filters as filters
import time
import glob
import logging
import sep
from reproject import reproject_interp
import tempfile
from dlnpyutils import utils as dln

def rootdirs():
    # Returns dldir, mssdir, localdir
    host = socket.gethostname()
    shost = host.split('.')[0]

    if shost == 'thing' or shost == 'hulk':
        return ('/dl1','/mss1/','/d0/')

    if shost == 'gp07' or shost == 'gp08' or shost == 'gp09':
        return ('/dl1','/net/mss1/','/data0/')

def meancube(imcube,wtcube,weights=None,crreject=False):
    """ This does the actual stack of an image cube.  The images must already be background-subtracted and scaled."""
    # Weights should be normalized, e.g., sum(weights)=1
    # pixels to be masked in an image should have wtcube = 0

    nx,ny,nimages = imcube.shape

    # Unweighted
    if weights is None:
        weights = np.ones(nimages,float)/nimages
    
    # Do the weighted average
    finaltot = np.zeros((nx,ny),float)
    finaltotwt = np.zeros((nx,ny),float)        
    totvarim = np.zeros((nx,ny),float)
    for i in range(nimages):
        mask = (wtcube[:,:,i] > )
        var = np.zeros((nx,ny),float)  # sig^2
        var[mask] = 1/wtcube[:,:,mask]
        finaltot[mask] += imcube[:,:,mask]*weights[i]
        finaltotwt[mask] += weights[i]
        # Variance in each pixel for noise images and the scalar weights
        totvarim[mask] += weights[i]*var
    # Create the weighted average image
    finaltotwt[finaltotwt<=0] = 1
    final = finaltot/finaltotwt
    # Create final error image
    error = np.sqrt(finalvarim)
    
    # CR rejection
    if crreject is True:
        pass
    
    return final,error


def stack(imagefiles,errorfiles,bgfiles,weights=None):
    """ Actually do the stacking/averaging of multiple images already reprojected."""

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
    head = fits.getheader(imagegfiles[0],0)
    fnx = head['ONAXIS1']
    fny = head['ONAXIS2']
    
    # Get the sizes and positions of the subimages
    dtype = np.dtype([('SUBX0',int),('SUBX1',int),('SUBY0',int),('SUBY1',int),('NX',int),('NY',int)])
    binstr = np.zeros(nbin,dtype=dtype)
    for b in range(nbin):
        head = fits.getheader(imagefiles[0],b)
        substr['X0'][b] = head['SUBX0']
        substr['X1'][b] = head['SUBX1']
        substr['Y0'][b] = head['SUBY0']
        substr['Y1'][b] = head['SUBY1']        
        binstr['NX'][b] = head['SUBNX']
        binstr['NY'][b] = head['SUBNY']

    # Final image
    final = np.zeros((fnx,fny),float)
    error = np.zeros((fnx,fny),float)    
        
    # Loop over bins
    for b in range(nbin):
        imcube = np.zeros((substr['NX'][b],substr['NY'][b],nimages),np.float64)
        wtcube = np.zeros((substr['NX'][b],substr['NY'][b],nimages),np.float64)        
        # Loop over images
        for f in range(nimages):
            im,head = fits.getheader(imagefiles[f],b,header=True)
            wt,whead = fits.getheader(weightfiles[f],b,header=True)
            # Background subtraction here???
            # Deal with NaNs
            wt[~np.isfinite(im)] = 0
            # Stuff into the cube
            imcube[:,:,f] = im
            wtcube[:,:,f] = wt
        # Do the weighted combination
        avgim,errim = meancube(imcube,wtcube,weights=weights)
        # Stuff into final image
        final[substr['X0'][b]:substr['X1'][b]+1,substr['Y0'][b]:substr['Y1'][b]+1] = avgim
        error[substr['X0'][b]:substr['X1'][b]+1,substr['Y0'][b]:substr['Y1'][b]+1] = errim

    return final,error

    
def coadd(imagefiles,weightfiles,meta,outhead,coaddtype='average'):
    """ Create a coadd given a list of images. """

    # meta should have zpterm, exptime, fwhm
    
    nimages = dln.size(imagefiles)

    # Figure out scales and weights
    # F_trans = 10^(-0.8*(delta_mag-0.2))    
    scales = cat['exptime'] * 10**(-0.8*(cat['zpterm']-0.2))

    # Use weight~S/N
    # S/N goes as sqrt(exptime) and in background-dominated regime S/N ~ 1/FWHM
    # so maybe something like weight ~ sqrt(scaling)/FWHM
    weights = np.sqrt(scales)/cat['FWHM']
    weights /= np.sum(weights)    # normalize
    
    # Loop over the images
    tempfiles = []
    tempwtfiles = []
    tempbgfiles = []
    for f in range(nimages):

        im,head = fits.getheader(imagefiles[f],header=True)
        im = im.byteswap(inplace=True).newbyteorder()    # for sep need native byte order
        wt,whead = fits.getheader(weightfiles[f],header=True)
        wt = wt.byteswap(inplace=True).newbyteorder()
        nx1,ny1 = im.shape
        
        # Step 1. Background subtract the image
        bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
        bkg_image = bkg.back()
        im -= bkg_image
        
        # Step 2. Reproject the image
        newim, footprint = reproject_interp((im,head), outhead)
        newwt, wfootprint = reproject_interp((wt,whead), outhead)
        newbg, bfootprint = reproject_interp((bkg_image,head), outhead)        
        fnx,fny = newim.shape

        # Step 3. Scale the image
        #  divide image by "scales"
        newim /= scales[i]
        #  wt = 1/err^2, need to perform same operation on err as on image
        newwt *= scales[i]**2
        
        # Step 4. Break up into bins and save to temporary file
        tid,tfile = tempfile.mkstemp(prefix="timage",dir="/tmp")
        tbase = os.path.basename(tfile)
        timfile = "/tmp/"+tbase+"_flx.fits"
        twtfile = "/tmp/"+tbase+"_wt.fits"
        tbgfile = "/tmp/"+tbase+"_bg.fits"

        timhdu = fits.open(timfile)
        twthdu = fits.open(twtfile)
        tbghdu = fits.open(tbgfile)        

        xbin = ybin = 2
        dx = fnx // xbin
        dy = fny // ybin
        
        # Put information in header
        # ONAXIS1, ONAXIS2, SUBX0, SUBX1, SUBY0, SUBY1, SUBNX, SUBNY
        for i in range(xbin):
            x0 = i*dx
            x1 = x0 + dx-1
            if i==(xbin-1): x1=(fnx-1)
            for j in range(ybin):
                y0 = j*dy
                y1 = y0 + dy-1
                if j==(ybin-1): y1=(fny-1)
                newhead = outhead.copy()
                newhead['SCALE'] = scales[i]
                newhead['WEIGHT'] = weights[i]
                newhead['ONAXIS1'] = newhead['NAXIS1']
                newhead['ONAXIS2'] = newhead['NAXIS2']
                newhead['SUBX0'] = x0
                newhead['SUBX1'] = x1
                newhead['SUBNX'] = x1-x0+1               
                newhead['SUBY0'] = y0
                newhead['SUBY1'] = y1
                newhead['SUBNY'] = y1-y0+1
                # Flux
                subim = newim[x0:x1+1,y0:y1+1].copy()
                hdu1 = fits.PrimaryHDU(subim,newhead.copy())
                timhdu.append(hdu1)
                # Weight
                subwt = newwt[x0:x1+1,y0:y1+1].copy()
                hdu1 = fits.PrimaryHDU(subwt,newhead.copy())
                twthdu.append(hdu1)
                # Background
                subbg = newbg[x0:x1+1,y0:y1+1].copy()
                hdu1 = fits.PrimaryHDU(subbg,newhead.copy())
                timhdu.append(hdu1)                
        timhdu.writeto(timfile,overwrite=True)
        timhdu.close()
        twthdu.writeto(twtfile,overwrite=True)
        twthdu.close()
        tbghdu.writeto(tbgfile,overwrite=True)
        tbghdu.close()        
                
        tempfiles.append(timfile)
        tempwtfiles.append(twtfile)
        tempbgfiles.append(tbgfile)

        
    # Step 4. Stack the images
    final = stack(tempfiles,tempwtfiles,tempbgfiles,weights)

    # Delete temporary files

    
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
