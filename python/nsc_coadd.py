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
from astropy.wcs import WCS
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
from dlnpyutils import utils as dln, coords


def rootdirs():
    # Returns dldir, mssdir, localdir
    host = socket.gethostname()
    shost = host.split('.')[0]

    if shost == 'thing' or shost == 'hulk':
        return ('/net/dl2','/mss1','/data0')

    if shost == 'gp06' or shost == 'gp07' or shost == 'gp08' or shost == 'gp09':
        return ('/net/dl2','/net/mss1','/data0')

    
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

def getbrickexposures(brick,band=None,version='v3'):
    """ Get exposures information that overlap a brick."""

    dldir, mssdir, localdir = rootdirs()
    brick_dbfile = dldir+'/dnidever/nsc/instcal/'+version+'/lists/nsc_bricks.db'
    brickdata = db.query(brick_dbfile, table='bricks', cols='*', where='brickname="'+brick+'"')

    # Healpix information
    pix128 = hp.ang2pix(128,brickdata['ra'],brickdata['dec'],lonlat=True)
    # neighbors
    neipix = hp.get_all_neighbours(128,pix128)
    
    # Get all of the exposures overlapping this region
    meta_dbfile = dldir+'/dnidever/nsc/instcal/'+version+'/lists/nsc_meta.db'
    allpix = np.hstack((neipix.flatten(),pix128))
    whr = ' or '.join(['ring128=='+h for h in allpix.astype(str)])
    chipdata = db.query(meta_dbfile, table='chip', cols='*', where=whr)

    # Do more overlap checking
    brick_vra = np.hstack((brickdata['ra1'],brickdata['ra2'],brickdata['ra2'],brickdata['ra1']))
    brick_vdec = np.hstack((brickdata['dec1'],brickdata['dec1'],brickdata['dec2'],brickdata['dec2']))    
    olap = np.zeros(len(chipdata),bool)
    for i in range(len(chipdata)):
        vra = np.hstack((chipdata['vra1'][i],chipdata['vra2'][i],chipdata['vra2'][i],chipdata['vra1'][i]))
        vdec = np.hstack((chipdata['vdec1'][i],chipdata['vdec1'][i],chipdata['vdec2'][i],chipdata['vdec2'][i]))
        olap[i] = coords.doPolygonsOverlap(vra,vdec,brick_vra,brick_vdec)
    ngdch = np.sum(olap)
    if ngdch==0:
        print('No exposures overlap brick '+brick)
        return None
    chipdata = chipdata[olap]
    # Check band
    if band is not None:
        gband, = np.where(chipdata['filter']==band)
        if len(gband)==0:
            print('No '+band+' exposures overlap brick '+brick)
            return None
        chipdata = chipdata[gband]
    exposure = np.unique(chipdata['exposure'])
    nexp = len(exposure)
    print(str(nexp)+' exposures overlap brick '+brick)

    # Get the exosure data
    whr = ' or '.join(['exposure=="'+e+'"' for e in exposure.astype(str)])
    expdata = db.query(meta_dbfile, table='exposure', cols='*', where=whr)

    return expdata
    
    
def nsc_coadd(brick,band='g',version='v3'):
    pass
    # This creates a coadd for one NSC brick

    # Make sure to fix the WCS using the coefficients I fit with Gaia DR2
    #  that are in the meta files.

    expdata = getbrickexposures(brick,band=band,version=version)
    
    dldir, mssdir, localdir = rootdirs()
    brick_dbfile = dldir+'/dnidever/nsc/instcal/'+version+'/lists/nsc_bricks.db'
    brickdata = db.query(brick_dbfile, table='bricks', cols='*', where='brickname="'+brick+'"')

    # Healpix information
    pix128 = hp.ang2pix(128,brickdata['ra'],brickdata['dec'],lonlat=True)
