#!/usr/bin/env python

import os
import sys
import numpy as np
import scipy
from astropy.io import fits
import photutils
from skimage import measure, morphology
from scipy.cluster import vq
#import gaps
import matplotlib.pyplot as plt
import pylab
from scipy.signal import argrelmin
import scipy.ndimage.filters as filters
import time

def loadpsf(filename):
    """ Load a DAOPHOT .psf file"""
    # Check if the file exists
    if os.path.exists(filename) is False:
        raise ValueError(filename+" NOT FOUND")

class Psf:
    """ DAOPHOT PSF class."""
    
    def __init__(self,filename):
        # Check if the file exists
        if os.path.exists(filename) is False:
            raise ValueError(filename+" NOT FOUND")
        # Load the file

        # Initalize the psf object


    def call(self,x,y,mag):
        """ Create a PSF image."""

    def __str__(self):
        pass

    def write(self,outfile):
        """ Write the PSF to a file."""
        pass
