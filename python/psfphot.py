#!/usr/bin/env python
#
# PSFPHOT.PY - PSF photometry routines
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20200323'  # yyyymmdd

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
#from astropy.table import Table, Column
import time
#import shutil
#import subprocess
#import logging
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve
#import astropy.stats
import sep
from dlnpyutils import utils as dln

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


# Run basic Source Extraction
def sextractor(im,err=None,mask=None,nsig=5.0,gain=1.0):

    # Check byte order, SEP needs little endian
    if im.dtype.byteorder == '>':
        data = im.byteswap().newbyteorder()
    else:
        data = im

    # Background estimation and subtraction
    bkg = sep.Background(data, mask, bw=256, bh=256, fw=3, fh=3)
    bkg_image = bkg.back()
    data_sub = data-bkg
    #data_sub[data>50000]=0.0
    # Detect and extract objects
    if err is None:
        objects = sep.extract(data_sub, nsig, err=bkg.globalrms, mask=mask)
    else:
        objects = sep.extract(data_sub, nsig, err=err, mask=mask)
        
    # Get mag_auto in 2 steps
    kronrad, krflag = sep.kron_radius(data_sub, objects['x'], objects['y'], objects['a'], objects['b'],
                                      objects['theta'], 6.0, mask=mask)
    flux, fluxerr, flag = sep.sum_ellipse(data_sub, objects['x'], objects['y'], objects['a'], objects['b'],
                                          objects['theta'], 2.5*kronrad, subpix=1, err=err, mask=mask, gain=gain)
    flag |= krflag  # combine flags into 'flag'

    # Use circular aperture if Kron radius is too small
    r_min = 1.75  # minimum diameter = 3.5
    use_circle = kronrad * np.sqrt(objects['a'] * objects['b']) < r_min
    if np.sum(use_circle)>0:
        cflux, cfluxerr, cflag = sep.sum_circle(data_sub, objects['x'][use_circle], objects['y'][use_circle],
                                                r_min, subpix=1, err=err, mask=mask, gain=gain)
        flux[use_circle] = cflux
        fluxerr[use_circle] = cfluxerr
        flag[use_circle] = cflag
    mag_auto = -2.5*np.log10(flux)+25.0

    # Make the final catalog
    newdt = np.dtype([('kronrad',float),('flux_auto',float),('mag_auto',float)])
    cat = dln.addcatcols(objects,newdt)
    cat['flag'] |= flag
    cat['kronrad'] = kronrad
    cat['flux_auto'] = flux
    cat['mag_auto'] = mag_auto
    
    return cat


