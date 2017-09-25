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

def rootdirs():
    # Returns dldir, mssdir, localdir
    host = socket.gethostname()
    shost = host.split('.')[0]

    if shost == 'thing' or shost == 'hulk':
        return ('/dl1','/mss1/','/d0/')

    if shost == 'gp07' or shost == 'gp08' or shost == 'gp09':
        return ('/dl1','/net/mss1/','/data0/')

def coadd_healpix(pix,coaddtype='average'):
    ''' Create a coadd for a single healpix pixel.
    '''

    #-give it healpix pixel number and any flags/constraints to use
    #-it will figure out which exposures/chips overlap that healpix coadd region
    #   (large rectangular box with buffer)
    #-impose the flags/constraints
    #-set up output WCS/header
    #-call nsc_coadd.py with list of files, WCS header, coadd method

    # Load the list    
    dldir, mssdir, localdir = rootdirs()
    listfile = localdir+'dnidever/nsc/instcal/nsc_healpix_list.fits'
    if not os.path.exists(self.home):
        print(listfile,' NOT FOUND')
        return
    healstr = fits.getdata(listfile,1)
    index = fits.getdata(listfile,2)
    # Find our pixel
    ind, = np.where(index['pix'] == pix)
    nind = len(ind)
    if nind eq 0:
        print('No entries for Healpix pixel "',%s,'" in the list' % pix)
        return
    ind = ind[0]
    list = healstr[index[ind]['lo'][0]:index[ind]['hi'][0]+1]
    nlist = len(list)

    # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
    #  so we can deal with the edge cases
    #NEIGHBOURS_RING,nside,pix,neipix,nneipix
    #for i=0,nneipix-1 do begin
    #  ind = where(index.pix eq neipix[i],nind)
    #  if nind gt 0 then begin
    #    ind = ind[0]
    #    list1 = healstr[index[ind].lo:index[ind].hi]
    #    push,list,list1
    #  endif
    #endfor
    ## Get unique values
    #ui = uniq(list.file,sort(list.file))
    #list = list[ui]
    #nlist = n_elements(list)
    #print,strtrim(nlist,2),' exposures that overlap this pixel and neighbors'

    # Get the boundary coordinates
    #   healpy.boundaries but not sure how to do it in IDL
    #   pix2vec_ring/nest can optionally return vertices but only 4
    #     maybe subsample myself between the vectors
    # Expand the boundary to include a "buffer" zone
    #  to deal with edge cases
    #PIX2VEC_RING,nside,pix,vec,vertex

    # Use python code to get the boundary
    #  this takes ~2s mostly from import statements
    #tempfile = MKTEMP('bnd')
    #file_delete,tempfile+'.fits',/allow
    #step = 100
    #pylines = 'python -c "from healpy import boundaries; from astropy.io import fits;'+$
    #          ' v=boundaries('+strtrim(nside,2)+','+strtrim(pix,2)+',step='+strtrim(step,2)+');'+$
    #          " fits.writeto('"+tempfile+".fits'"+',v)"'
    #spawn,pylines,out,errout
    #vecbound = MRDFITS(tempfile+'.fits',0,/silent)
    #file_delete,[tempfile,tempfile+'.fits'],/allow
    #VEC2ANG,vecbound,theta,phi
    #rabound = phi*radeg
    #decbound = 90-theta*radeg

    # Expand the boundary by the buffer size
    #PIX2ANG_RING,nside,pix,centheta,cenphi
    #cenra = cenphi*radeg
    #cendec = 90-centheta*radeg
    # reproject onto tangent plane
    #ROTSPHCEN,rabound,decbound,cenra,cendec,lonbound,latbound,/gnomic
    # expand by a fraction, it's not an extact boundary but good enough
    #buffsize = 10.0/3600. ; in deg
    #radbound = sqrt(lonbound^2+latbound^2)
    #frac = 1.0 + 1.5*max(buffsize/radbound)
    #lonbuff = lonbound*frac
    #latbuff = latbound*frac
    #buffer = {cenra:cenra,cendec:cendec,lon:lonbuff,lat:latbuff}




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
