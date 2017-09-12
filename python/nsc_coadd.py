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
import shutil
import re
import subprocess
import glob
import logging
import socket

def coadd_healpix(pix,coaddtype='average'):
    ''' Create a coadd for a single healpix pixel.
    '''

    #-give it healpix pixel number and any flags/constraints to use
    #-it will figure out which exposures/chips overlap that healpix coadd region
    #   (large rectangular box with buffer)
    #-impose the flags/constraints
    #-set up output WCS/header
    #-call nsc_coadd.py with list of files, WCS header, coadd method


def coadd(images,outhead,coaddtype='average'):
    ''' Create a coadd given a list of images.
    '''

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
