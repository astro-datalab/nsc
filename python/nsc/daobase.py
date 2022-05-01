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

class Exposure:
    """ An astronomical image or exposure which contains
       -flux image
       -noise/variance image
       -mask image
       -background image
       -FITS header
       -psf
       -wcs?
    """
    def __init__(self,flux,noise=None,mask=None,background=None,header=None):
        # Flux
        self.flux = flux
        ny, nx = flux.shape
        # Noise
        if noise is None:
            self.noise = np.zeros([ny,nx],'f')
            self.noise[:] = mad(flux)
        else:
            self.noise = noise
        # Mask, True means BAD
        if mask is None:
            self.mask = np.zeros([ny,nx],'b')+False
        else:
            self.mask = np.zeros([ny,nx],'b')+False
            self.mask[:] = mask
        # Background
        if background is None:
            self.background = np.zeros([ny,nx],'f')
            self.background[:] = np.median(flux)
        else:
            self.background = background
        # Header
        if header is None:
            self.header = ''
        else:
            self.header = header

    def __str__(self):
        txt = "Flux: ["+','.join(map(str,self.flux.shape))+"]\n"
        #txt += "  med="+np.median(self.flux).astype('a')+' sigma='+mad(self.flux).astype('a')+"\n"
        txt += "Noise: ["+','.join(map(str,self.noise.shape))+"]\n"
        #txt += "  med="+np.median(self.noise).astype('a')+' sigma='+mad(self.noise).astype('a')+"\n"
        txt += "Mask : ["+','.join(map(str,self.mask.shape))+"]\n"
        txt += "Backg : ["+','.join(map(str,self.background.shape))+"]\n"
        txt += "Header : ["+str(len(self.header))+"]"
        return txt
        
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
        It needs to be scaled by 1.4826x to be on the same "scale" as std. dev.
    """
    med = np.median(arr)
    return 1.4826*np.median(np.abs(arr - med))

def doPolygonsOverlap(xPolygon1, yPolygon1, xPolygon2, yPolygon2):
    """Returns True if two polygons are overlapping."""

    # How to determine if two polygons overlap.
    # If a vertex of one of the polygons is inside the other polygon
    # then they overlap.
    
    n1 = len(xPolygon1)
    n2 = len(xPolygon2)
    isin = False

    # Loop through all vertices of second polygon
    for i in range(n2):
        # perform iterative boolean OR
        # if any point is inside the polygon then they overlap   
        isin = isin or isPointInPolygon(xPolygon1, yPolygon1, xPolygon2[i], yPolygon2[i])

    # Need to do the reverse as well, not the same
    for i in range(n1):
        isin = isin or isPointInPolygon(xPolygon2, yPolygon2, xPolygon1[i], yPolygon1[i])

    return isin

def isPointInPolygon(xPolygon, yPolygon, xPt, yPt):
    """Returns boolean if a point is inside a polygon of vertices."""
    
    # How to tell if a point is inside a polygon:
    # Determine the change in angle made by the point and the vertices
    # of the polygon.  Add up the delta(angle)'s from the first (include
    # the first point again at the end).  If the point is inside the
    # polygon, then the total angle will be +/-360 deg.  If the point is
    # outside, then the total angle will be 0 deg.  Points on the edge will
    # outside.
    # This is called the Winding Algorithm
    # http://geomalgorithms.com/a03-_inclusion.html

    n = len(xPolygon)
    # Array for the angles
    angle = np.zeros(n)

    # add first vertex to the end
    xPolygon1 = np.append( xPolygon, xPolygon[0] )
    yPolygon1 = np.append( yPolygon, yPolygon[0] )

    wn = 0   # winding number counter

    # Loop through the edges of the polygon
    for i in range(n):
        # if edge crosses upward (includes its starting endpoint, and excludes its final endpoint)
        if yPolygon1[i] <= yPt and yPolygon1[i+1] > yPt:
            # if (P is  strictly left of E[i])    // Rule #4
            if isLeft(xPolygon1[i], yPolygon1[i], xPolygon1[i+1], yPolygon1[i+1], xPt, yPt) > 0: 
                 wn += 1   # a valid up intersect right of P.x

        # if edge crosses downward (excludes its starting endpoint, and includes its final endpoint)
        if yPolygon1[i] > yPt and yPolygon1[i+1] <= yPt:
            # if (P is  strictly right of E[i])    // Rule #4
            if isLeft(xPolygon1[i], yPolygon1[i], xPolygon1[i+1], yPolygon1[i+1], xPt, yPt) < 0: 
                 wn -= 1   # a valid up intersect right of P.x

    # wn = 0 only when P is outside the polygon
    if wn == 0:
        return False
    else:
        return True

def isLeft(x1, y1, x2, y2, x3, y3):
    # isLeft(): test if a point is Left|On|Right of an infinite 2D line.
    #   From http://geomalgorithms.com/a01-_area.html
    # Input:  three points P1, P2, and P3
    # Return: >0 for P3 left of the line through P1 to P2
    # =0 for P3 on the line
    # <0 for P3 right of the line
    return ( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) )

def create_mask(im,head):
    """ Create the mask using the image and the
        saturation level set in the header
    """
    # Create the mask
    #  mask=1   bad
    #  mask=0   good
    try:
        saturate = head['saturate']
    except:
        saturate = 60000.0
    mask = im > 0.9*saturate
    return mask
    
def create_noise(im,head,mask=None):
    """ This creates the noise image using the image
        and GAIN and RDNOISE from the header
        If the mask is input then it also sets
        noise to 999999. for masked pixels.
    """

    # Create the noise mask
    try:
        # electrons
        rdnoise = head['rdnoise']
    except:
        rdnoise = 0.0
    try:
        # electrons/ADU
        gain = head['gain']
    except:
        # Sky sigma should be sqrt(im*gain)/gain = sqrt(im/gain)
        # so we can roughly (ignoring rdnoise) back out gain
        # gain ~ im/skysig^2
        sig = mad(im)
        skysig = sig
        skymed = np.median(im)
        gain = skymed/skysig**2
        # Now try to remove rdnoise, normally in electrons
        if rdnoise > 0.0:
            skysig = np.sqrt(sig**2 - (rdnoise/gain)**2)
            gain = skymed/skysig**2
    # Make the simple noise model
    noise = np.sqrt(im/gain+(rdnoise/gain)**2)
    # Set noise to high value for masked pixels
    if mask is not None:
        noise[mask > 0] = 999999.

    return noise

def get_subim(im,xind,yind,hwidth,mask=None,noise=None,sky=0.0):
    """ This function returns a subimage
        xind, yind   the indices for the im array
        hwidth is the half-width.  Total width = 2*hwidth+1
        sky  the sky level to use, default=0.0
        It turns an exposure object
    """

    # Getting the image size
    ny, nx = im.shape

    # Initializing the subimage
    subim = np.zeros([2*hwidth+1,2*hwidth+1],'f') + sky
    submask = np.zeros([2*hwidth+1,2*hwidth+1],'b')
    subnoise = np.zeros([2*hwidth+1,2*hwidth+1],'f')+1.0
    # set noise=1.0 in out of bounds region so weights
    # won't get messed up
    
    # Checking that we're getting any part of the image
    # The center can be off the edge
    if ( (xind+hwidth) >= 0 ) & ( (xind-hwidth) <= (nx-1) ) \
       & ( (yind+hwidth) >= 0 ) & ( (yind-hwidth) <= (ny-1) ):

       # Indices for the big image
       xbg0 = (xind-hwidth) if (xind-hwidth) > 0 else 0
       xbg1 = (xind+hwidth) if (xind+hwidth) < (nx-1) else (nx-1)
       ybg0 = (yind-hwidth) if (yind-hwidth) > 0 else 0
       ybg1 = (yind+hwidth) if (yind+hwidth) < (ny-1) else (ny-1)

       # Indices for the subimage
       xsm0 = hwidth+xbg0-xind
       xsm1 = hwidth+xbg1-xind
       ysm0 = hwidth+ybg0-yind
       ysm1 = hwidth+ybg1-yind

       # Getting part of the image
       # must add+1 to end indices
       subim[ysm0:ysm1+1,xsm0:xsm1+1] = im[ybg0:ybg1+1,xbg0:xbg1+1]

       # Fill in the mask
       if mask is not None:
           submask[ysm0:ysm1+1,xsm0:xsm1+1] = mask[ybg0:ybg1+1,xbg0:xbg1+1]

       # Fill in the noise
       if noise is not None:
           subnoise[ysm0:ysm1+1,xsm0:xsm1+1] = noise[ybg0:ybg1+1,xbg0:xbg1+1]
           
    # Create the sub exposure
    subexp = Exposure(subim,mask=submask,noise=subnoise)
       
    return subexp

def get_subim_rec(im,coord,boundary,mask=None,noise=None,sky=0.0):
    """ This function returns a rectangular subimage
        im        The flux image
        coord     [xcenter, ycenter]
        boundary  Boundary values of the subimage, [ymin,xmin,ymax,xmax].
                     They can go beyond the boundary of IM.
        mask      Mask for IM, default=None
        noise     Noise for IM, default=None
        sky       Sky level to use, default=0.0
        It returns an exposure object
    """

    # Getting the image size
    ny, nx = im.shape

    # Get boundary limits
    ymin = boundary[0]
    xmin = boundary[1]
    ymax = boundary[2]
    xmax = boundary[3]
    nysm = ymax-ymin+1
    nxsm = xmax-xmin+1
    
    # Initializing the subimage
    subim = np.zeros([nysm,nxsm],'f') + sky
    submask = np.zeros([nysm,nxsm],'b')
    subnoise = np.zeros([nysm,nxsm],'f')+1.0
    # set noise=1.0 in out of bounds region so weights
    # won't get messed up
    
    # Checking that we're getting any part of the image
    # The center can be off the edge
    if ( xmax >= 0 ) & ( xmin <= (nx-1) ) & ( ymax >= 0 ) & ( ymin <= (ny-1) ):

       # Indices for the big image
       xbg0 = xmin if xmin > 0 else 0
       xbg1 = xmax if xmax < (nx-1) else (nx-1)
       ybg0 = ymin if ymin > 0 else 0
       ybg1 = ymax if ymax < (ny-1) else (ny-1)

       # Indices for the subimage
       xsm0 = xbg0-xmin
       xsm1 = xbg1-xmin
       ysm0 = ybg0-ymin
       ysm1 = ybg1-ymin

       # Getting part of the image
       # must add+1 to end indices
       subim[ysm0:ysm1+1,xsm0:xsm1+1] = im[ybg0:ybg1+1,xbg0:xbg1+1]

       # Fill in the mask
       if mask is not None:
           submask[ysm0:ysm1+1,xsm0:xsm1+1] = mask[ybg0:ybg1+1,xbg0:xbg1+1]

       # Fill in the noise
       if noise is not None:
           subnoise[ysm0:ysm1+1,xsm0:xsm1+1] = noise[ybg0:ybg1+1,xbg0:xbg1+1]
           
    # Create the sub exposure
    subexp = Exposure(subim,mask=submask,noise=subnoise)
       
    return subexp

def get_fluxcenter(subim,mask=None,noise=None):
    """ Calculate the flux-weighted center     
        Returns (xcen,ycen,xcenerr,ycenerr) tuple.
        If "noise" is input then errors in the centroids
        are computed, other xcenerr/ycenerr are NAN
    """
    
    # Center-of-Mass like centroid.  Weight by flux
    # Add up all position along columns since they will have the same
    # X value. Similar for Y
    ny, nx = subim.shape
    # Use input mask and mask out negative pixels
    if mask is None:
        gmask = subim >= 0.0
    else:
        gmask = (subim >= 0.0) & (mask == False)
    totflux = np.sum(subim*gmask)
    xcen = np.sum( np.sum(subim*gmask,axis=0)*np.arange(nx) )/totflux
    ycen = np.sum( np.sum(subim*gmask,axis=1)*np.arange(ny) )/totflux

    # Calculate uncertainties in flux-weighted center
    #  if noise is input
    if noise is not None:
        # xc = sum_i(sum_j(flux[i,j]*x[j])) / sum_i(sum_j(flux[i,j]))
        # Using eqns. 32+33 from SExtractor manual
        xcenerr = np.sqrt( np.sum( np.sum((noise**2)*gmask,axis=0)*np.arange(nx)**2 )/totflux**2 )
        ycenerr = np.sqrt( np.sum( np.sum((noise**2)*gmask,axis=1)*np.arange(ny)**2 )/totflux**2 )
    else:
        xcenerr = np.nan
        ycenerr = np.nan
        
    return (xcen,ycen,xcenerr,ycenerr)

def detect_segment(exp,nsig=5.0):
    """ Detect with image segmentation
    """
    
    # Parameter checking
    if nsig <= 0:
        print "Nsig must be >0"
        return

    # Subtract background from image
    im = exp.flux - exp.background
    sigma = exp.noise
    ny, nx = im.shape
    
    # image segmentation
    #skimage.measure.label
    #http://www.scipy-lectures.org/packages/scikit-image/
    #t0 = time.time()
    blobs = (im > nsig*sigma) & (exp.mask == False)
    all_labels = measure.label(blobs,background=0)
    # remove regions that are too small
    labels = morphology.remove_small_objects(all_labels,5)
    all_props = measure.regionprops(labels,im)
    #  area, bbox, centroid, coords, eccentricity, image, intensity_image
    #  major_axis_length, minor_axis_length, moments_central
    #  orientation, perimeter, weighted_centroid, weighted_moments_central
    #
    # there's also the scipy.ndimage.label and measurements function
    # figure out which one's faster, phoutils uses scipy.ndimage
    # I think skimage extends the capabilities of scipy.ndimage
    #http://scikit-image.org/docs/dev/overview.html
    #print time.time()-t0
    # 0.28 sec
    # I could use this to just compute what I want

    # Create peak structure
    nreg = len(all_props)
    dt = np.dtype([('id',int),('x',float),('y',float),('area',int),('bbox_x0',int),('bbox_x1',int),
                   ('bbox_y0',int),('bbox_y1',int),('major_axis',float),('minor_axis',float),
                   ('theta',float),('eccentricity',float)])
    cat = np.zeros(nreg,dtype=dt)
    for i in xrange(nreg):
        cat['id'][i] = all_props[i].label
        centroid = all_props[i].weighted_centroid
        cat['x'][i] = centroid[1]
        cat['y'][i] = centroid[0]
        cat['area'][i] = all_props[i].area
        bbox = all_props[i].bbox
        cat['bbox_y0'][i] = bbox[0]
        cat['bbox_x0'][i] = bbox[1]
        cat['bbox_y1'][i] = bbox[2]
        cat['bbox_x1'][i] = bbox[3]
        cat['major_axis'][i] = all_props[i].major_axis_length
        cat['minor_axis'][i] = all_props[i].minor_axis_length
        cat['theta'][i] = np.rad2deg(all_props[i].orientation)
        cat['eccentricity'][i] = all_props[i].eccentricity

    return cat
        
def detect_delta(exp,nsig=5.0,fluxfrac=0.5):
    """ This detects peaks in an image
    """

    # Parameter checking
    if nsig <= 0:
        print "Nsig must be >0"
        return
    if fluxfrac < 0:
        print "Fluxfrac must be >0"
        return

    # Subtract background from image
    im = exp.flux - exp.background
    sigma = exp.noise
    ny, nx = im.shape
    # Smooth with a small Gaussian
    #   mask the image
    sigma0 = 1.0
    fmask = 1-exp.mask.astype('f')
    smim_masked = filters.gaussian_filter(im*fmask,sigma=sigma0,mode='mirror',truncate=2.0)
    smwt = filters.gaussian_filter(fmask,sigma=sigma0,mode='mirror',truncate=2.0)
    nogoodpix = smwt <= 0.0
    if np.sum(nogoodpix) > 0:
        smwt[nogoodpix] = 1.0
        smim = smim_masked / smwt
        smim[nogoodpix] = 0.0
    else:
        smim = smim_masked / smwt

    # Getting maxima points
    #  fix edges, assume the perimeter is 0.5*fluxfrac of last column
    #  this will allow it to pass through the neighbor fluxfraction
    #  criteria for detection below.
    diffx1 = smim - np.roll(smim,1,axis=1)
    diffx1[:,0] = 0.5*fluxfrac*smim[:,0] * (smim[:,0] > 0.0)
    diffx2 = smim - np.roll(smim,-1,axis=1)
    diffx2[:,nx-1] = 0.5*fluxfrac*smim[:,nx-1] * (smim[:,nx-1] > 0.0)
    diffy1 = smim - np.roll(smim,1,axis=0)
    diffy1[0,:] = 0.5*fluxfrac*smim[0,:] * (smim[0,:] > 0.0)
    diffy2 = smim - np.roll(smim,-1,axis=0)
    diffy2[ny-1,:] = 0.5*fluxfrac*smim[ny-1,:] * (smim[ny-1,:] > 0.0)

    # Try the maximum_filter
    #neighborhood_size = 5
    #smim_max = filters.maximum_filter(smim, neighborhood_size)
    #smim_maxima = (smim == smim_max)
    # This works but finds MANY peaks, needs to apply sigma threshold
    # problems is that it's harder to easily check the flux
    # ratio/fraction of the neighbors
    
    # Set the threshold for the flux difference for neighbors
    diffth = smim*fluxfrac * (smim > 0.0)
    
    # Detect
    #  first line - it's a peak
    #  second line - neighbor flux fraction threshold
    #       they must be >=fluxfrac*peak
    #  third line - Nsigma threshold in original image
    #  fourth line - mask not set for pixel
    detect = (diffx1 >= 0) & (diffx2 >= 0) & (diffy1 >= 0) & (diffy2 >= 0) \
             & (diffx1 < diffth) & (diffx2 < diffth) & (diffy1 < diffth) & (diffy2 < diffth) \
             & (im >= nsig*sigma) \
             & (exp.mask == False)
    src = np.where(detect == True)
    ndetect = len(src[0])
    xdetect = np.zeros(ndetect,'i')
    xdetect[:] = src[1]
    ydetect = np.zeros(ndetect,'i')
    ydetect[:] = src[0]

    # Create peak structure
    dt = np.dtype([('xcen',int),('ycen',int),('nsig',float)])
    peaks = np.zeros(ndetect,dtype=dt)
    peaks['nsig'] = 0.0
    
    # Stuff in some information
    #peaks['xcen'] = xdetect
    #peaks['ycen'] = ydetect
    #peaks['nsig'] = im[ydetect,xdetect] / sigma[ydetect,xdetect]

    #pylab.ion()
    #plt.imshow(smim,interpolation='nearest',origin='lower')
    #plt.plot(xdetect,ydetect,'r+')
    #import pdb; pdb.set_trace()

    # I'M GETTING "HOT PIXELS" COMING THROUGH
    
    # Loop through the peaks and make sure
    #   that we have the brightest peak within +/-2 pixels
    #   this is to make sure that we don't get any duplicates.
    #   Can't make this too large or we won't detect fainter
    #   neighbor peaks that are close by
    detmask = np.zeros(im.shape,'b')
    detbuff = 2
    for i in xrange(ndetect):
        # Buffer x/y ranges
        xlo = (xdetect[i]-detbuff) if (xdetect[i]-detbuff) > 0 else 0
        xhi = (xdetect[i]+detbuff) if (xdetect[i]+detbuff) < (nx-1) else (nx-1)
        ylo = (ydetect[i]-detbuff) if (ydetect[i]-detbuff) > 0 else 0
        yhi = (ydetect[i]+detbuff) if (ydetect[i]+detbuff) < (ny-1) else (ny-1)
        maxsubsmim = np.max(smim[ylo:yhi+1,xlo:xhi+1])
        ndetmask = np.sum(detmask[ylo:yhi+1,xlo:xhi+1])
        # Brightest value within +/-2 pixels (use smoothed image)
        #  and no peak previously detected nearby
        if (smim[ydetect[i],xdetect[i]] >= maxsubsmim) & (ndetmask == 0):
            # Since we find peaks in the smoothed image they might
            #  not correspond to peaks in the original image
            #  there might be a +/-1 pixel offset
            peaks['xcen'][i] = xdetect[i]
            peaks['ycen'][i] = ydetect[i]
            peaks['nsig'][i] = im[ydetect[i],xdetect[i]] / sigma[ydetect[i],xdetect[i]]
            detmask[ylo:yhi+1,xlo:xhi+1] = 1
            
    # Only keep good peaks
    gdpeaks, = np.where(peaks['nsig'] >= nsig)
    ngdpeaks = len(gdpeaks)
    if ngdpeaks > 0:
        peaks = peaks[gdpeaks]
    else:
        peaks = None

    print ngdpeaks, "sources detected in image"
        
    return peaks

def detect_peaksegment(exp,nsig=5.0):
    """ Detect with peaks and use image segmentation
        to find their footprint
    """
    
    # Parameter checking
    if nsig <= 0:
        print "Nsig must be >0"
        return

    # Subtract background from image
    im = exp.flux - exp.background
    sigma = exp.noise
    ny, nx = im.shape

    # Detect peaks
    peaks1 = detect_delta(exp,nsig=nsig)
    
    # Image segmentation
    #skimage.measure.label
    #http://www.scipy-lectures.org/packages/scikit-image/
    #t0 = time.time()
    blobs = (im > nsig*sigma) & (exp.mask == False)
    labels = measure.label(blobs,background=0)
    # remove regions that are too small
    #labels = morphology.remove_small_objects(all_labels,5)
    props = measure.regionprops(labels,im)
    nreg = len(props)
    #  area, bbox, centroid, coords, eccentricity, image, intensity_image
    #  major_axis_length, minor_axis_length, moments_central
    #  orientation, perimeter, weighted_centroid, weighted_moments_central
    #
    # there's also the scipy.ndimage.label and measurements function
    # figure out which one's faster, photutils uses scipy.ndimage
    # I think skimage extends the capabilities of scipy.ndimage
    #http://scikit-image.org/docs/dev/overview.html
    #print time.time()-t0
    # 0.28 sec
    # I could use this to just compute what I want

    # Create peak structure
    npeaks = len(peaks1)
    dt = np.dtype([('xcen',int),('ycen',int),('nsig',float),('maxflux',float),
                   ('area',int),('bbox_x0',int),('bbox_x1',int),('bbox_y0',int),('bbox_y1',int)])
    peaks = np.zeros(npeaks,dtype=dt)
    peaks['xcen'] = peaks1['xcen']
    peaks['ycen'] = peaks1['ycen']
    peaks['nsig'] = peaks1['nsig']
    peaks['maxflux'] = im[peaks1['ycen'],peaks1['xcen']]
    
    # Get region labels for all the peaks
    peak_labels = np.zeros(npeaks,'i')
    peak_labels[:] = labels[peaks1['ycen'],peaks1['xcen']]
    # Indices for all the peaks
    allpeaksind = np.arange(npeaks)
    
    # Now loop through the regions and resegment
    # with halfmax
    for i in xrange(nreg):
        # Get bbox footprint for this region
        bbox = props[i].bbox
        subim = im[bbox[0]:bbox[2]+1,bbox[1]:bbox[3]+1]
        submask = exp.mask[bbox[0]:bbox[2]+1,bbox[1]:bbox[3]+1]
        nys, nxs = subim.shape

        # Figure out which peaks have this label
        labelbool = np.in1d(peak_labels,props[i].label)
        if np.sum(labelbool) > 0:
            peakind = allpeaksind[labelbool]
            npeakind = len(peakind)
            # Loop over the peaks in this region
            for j in xrange(npeakind):
                peakind1 = peakind[j]
                # Make the halfmax mask image
                #  make the threshold level slightly smaller
                #  so that later on we can get the whole
                #  half-max region with the contour
                halfmax_subim = (subim >= 0.4*peaks['maxflux'][peakind1]) & (submask == False)
                # Now segment the halfmax_subim
                sublabels = measure.label(halfmax_subim,background=0)
                subprops = measure.regionprops(sublabels,subim)
                nsubreg = len(subprops)
                # Find the sublabel for this peak
                xsub = peaks['xcen'][peakind1]-bbox[1]
                ysub = peaks['ycen'][peakind1]-bbox[0]
                peak_sublabel = sublabels[ysub,xsub]
                # NOT SURE IF THIS INDEXING WILL BE RIGHT
                peak_subprop = subprops[peak_sublabel-1]
                bbox_sub = peak_subprop.bbox
                peaks[peakind1]['area'] = peak_subprop.area
                peaks[peakind1]['bbox_y0'] = bbox_sub[0]+bbox[0]
                peaks[peakind1]['bbox_x0'] = bbox_sub[1]+bbox[1]
                peaks[peakind1]['bbox_y1'] = bbox_sub[2]+bbox[0]
                peaks[peakind1]['bbox_x1'] = bbox_sub[3]+bbox[1]

                #pylab.ion()
                #plt.imshow(smim,interpolation='nearest',origin='lower')
                #plt.plot(xdetect,ydetect,'r+')
                #import pdb; pdb.set_trace()
                
                #cat['id'][i] = all_props[i].label
                #centroid = all_props[i].weighted_centroid
                #cat['x'][i] = centroid[1]
                #cat['y'][i] = centroid[0]
                #cat['area'][i] = all_props[i].area
                #bbox = all_props[i].bbox
                #cat['bbox_y0'][i] = bbox[0]
                #cat['bbox_x0'][i] = bbox[1]
                #cat['bbox_y1'][i] = bbox[2]
                #cat['bbox_x1'][i] = bbox[3]
                #cat['major_axis'][i] = all_props[i].major_axis_length
                #cat['minor_axis'][i] = all_props[i].minor_axis_length
                #cat['theta'][i] = np.rad2deg(all_props[i].orientation)
                #cat['eccentricity'][i] = all_props[i].eccentricity

    # Deal with peaks not associated with any label
    # label=0
        
    return peaks
        

def get_rec_footprint(im,coord):
    """ This returns the rectangular coordinates of the approximate
        1/2 maximum footprint for the peak input.
        im     The image
        coord  The [y,x] coordinates of the peak to find the footprint
    """

    ny,nx = im.shape
    ycen = coord[0]
    xcen = coord[1]
    hwidth = 30
    maxim = im[ycen,xcen]
    halfmax = 0.5*maxim
   
    # Check the marginal sums, do +/-1 and +/-30 for now
    # positive vertical stripe, thin in x, long in y
    if (ycen < (ny-1)):
        pvx0 = (xcen-1) if (xcen-1)>=0 else 0
        pvx1 = (xcen+1) if (xcen+1)<=(nx-1) else (nx-1)
        pvy0 = (ycen+1) if (ycen+1)>=0 else 0
        pvy1 = (ycen+hwidth) if (ycen+hwidth)<=(ny-1) else (ny-1)
        pvert = im[pvy0:pvy1+1,pvx0:pvx1+1]
        mpvert = np.sum(pvert,axis=1)/(pvx1-pvx0+1)
        # Find where it drops below 1/2 max
        pvloind, = np.where(mpvert < halfmax)
        if len(pvloind) > 0:
            ymax = pvy0+np.min(pvloind)
        else:
            # no half crossing found, but all elements
            #  are lower than max, so slowly declining
            pvminima, = argrelmin(mpvert,mode='wrap')
            # take the first minimum
            ymax = pvy0 + pvminima[0]
    else:
        ymax = ycen
    # negative vertical stripe, thin in x, long in y
    if (ycen > 0):
        nvx0 = (xcen-1) if (xcen-1)>=0 else 0
        nvx1 = (xcen+1) if (xcen+1)<=(nx-1) else (nx-1)
        nvy0 = (ycen-hwidth) if (ycen-hwidth)>=0 else 0
        nvy1 = (ycen-1) if (ycen-1)<=(ny-1) else (ny-1)
        nvert = im[nvy0:nvy1+1,nvx0:nvx1+1]
        mnvert = np.sum(nvert,axis=1)/(nvx1-nvx0+1)
        # Find where it drops below 1/2 max
        nvloind, = np.where(mnvert < halfmax)
        if len(nvloind) > 0:
            ymin = nvy0+np.max(nvloind)
        else:
            nvminima, = argrelmin(mnvert,mode='wrap')
            # take the last minimum
            ymin = nvy0 + np.max(nvminima)
    else:
        ymin = ycen
    # positive horizontal stripe, thin in y, long in x
    if (xcen < (nx-1)):
        phx0 = (xcen+1) if (xcen+1)>=0 else 0
        phx1 = (xcen+hwidth) if (xcen+hwidth)<=(nx-1) else (nx-1)
        phy0 = (ycen-1) if (ycen-1)>=0 else 0
        phy1 = (ycen+1) if (ycen+1)<=(ny-1) else (ny-1)
        phori = im[phy0:phy1+1,phx0:phx1+1]
        mphori = np.sum(phori,axis=0)/(phy1-phy0+1)
        # Find where it drops below 1/2 max
        phloind, = np.where(mphori < halfmax)
        if len(phloind) > 0:
            xmax = phx0+np.min(phloind)
        else:
            phminima, = argrelmin(mphori,mode='wrap')
            # take the first minimum
            xmax = phx0 + np.min(phminima)
    else:
        xmax = xcen
    # negative horizontal stripe, thin in y, long in x
    if (xcen > 0):
        nhx0 = (xcen-hwidth) if (xcen-hwidth)>=0 else 0
        nhx1 = (xcen-1) if (xcen-1)<=(nx-1) else (nx-1)
        nhy0 = (ycen-1) if (ycen-1)>=0 else 0
        nhy1 = (ycen+1) if (ycen+1)<=(ny-1) else (ny-1)
        nhori = im[nhy0:nhy1+1,nhx0:nhx1+1]
        mnhori = np.sum(nhori,axis=0)/(nhy1-nhy0+1)
        # Find where it drops below 1/2 max
        nhloind, = np.where(mnhori < halfmax)
        if len(nhloind) > 0:
            xmin = nhx0+np.max(nhloind)
        else:
            nhminima, = argrelmin(mnhori,mode='wrap')
            # take the last minimum
            xmin = nhx0 + np.max(nhminima)
    else:
        xmin = xcen
            
    #pylab.ion()
    #plt.imshow(smim,interpolation='nearest',origin='lower')
    #plt.plot(xdetect,ydetect,'r+')
    #import pdb; pdb.set_trace()

    # Now do similar procedure for marginal sums over 1.5-2.0x larger area
    # Only do this if we NEED to.  Check the perimeter pixels (xmin,xmax,ymin,ymax)
    # if they are all BELOW halfmax then we don't need to do the next step
    # of marginal sums
    
    boundary = [ymin,xmin,ymax,xmax]
    
    return boundary
    

def get_contour(im,level,coord):
    """ This returns the contour in IM at LEVEL around the COORD coordinates
        im     The image to get the contour for.
        level  The flux level for the contour.
        coord  The [y,x] coordinates for the center about which
                 the contour is desired.
    """
    ycen = coord[0]
    xcen = coord[1]
    # Getting the contours
    allcontours = measure.find_contours(im, level,positive_orientation='high')
    # list of (n,2) contours
    # need to find the one that encloses the center
    # check my vertex overlap functions in printVisitSkyMap.py to see
    #  which contour encloses the center or flux center
    isinpoly = np.zeros(len(allcontours))
    for f in xrange(len(allcontours)):
        isinpoly[f] = isPointInPolygon(allcontours[f][:,0],allcontours[f][:,1],ycen,xcen)
    gdcont = np.where(isinpoly == 1)[0]
    ngdcont = len(gdcont)
    if ngdcont > 0:
        contour = allcontours[gdcont[0]]
    else:
        # return empty list
        contour = []

    return contour

def get_morph_single(im,exp,morph,hwidth=10,hwidthS=3,noerrors=False):
    """ Measure morphological parameters of a source in a small image
    """
    
    # Get rectangular footprint
    #boundary = get_rec_footprint(im,[morph['y0'],morph['x0']])

    # Get small sub-image just to compute the flux center
    #  just with the pixels near the peak
    if noerrors is False:
        subexpS = get_subim(im,morph['x0'],morph['y0'],hwidthS,mask=exp.mask,noise=exp.noise)
    else:
        subexpS = get_subim(im,morph['x0'],morph['y0'],hwidthS,mask=exp.mask)

    # Getting flux center from small subimage
    if noerrors is False:
        xcenS, ycenS, xcenerrS, ycenerrS = get_fluxcenter(subexpS.flux,mask=subexpS.mask,noise=subexpS.noise)
    else:
        xcenS, ycenS, xcenerrS, ycenerrS = get_fluxcenter(subexpS.flux,mask=subexpS.mask)
    morph['x'] = xcenS + morph['x0'] - hwidthS
    morph['y'] = ycenS + morph['y0'] - hwidthS
    morph['xerr'] = xcenerrS
    morph['yerr'] = ycenerrS
    
    # Getting maximum from small subimage
    maxim = np.max(subexpS.flux*(1-subexpS.mask))
    morph['max'] = maxim
    
    # Get larger sub-image
    # If bbox exists then use that to get the width/subimage
    if 'bbox_x0' in morph.dtype.names:
        hwidth1 = np.max([  np.abs([morph['bbox_x0'],morph['bbox_x1']]-morph['x0']),
                            np.abs([morph['bbox_y0'],morph['bbox_y1']]-morph['y0']) ])
        hwidth1 = hwidth1 if (hwidth1>hwidth) else 5
    # Using preset subimage size
    else:
        hwidth1 = hwidth
    # Get the large subexposure
    if noerrors is False:
        subexp = get_subim(im,morph['x0'],morph['y0'],hwidth1,mask=exp.mask,noise=exp.noise)
    else:
        subexp = get_subim(im,morph['x0'],morph['y0'],hwidth1,mask=exp.mask)
        
    # Get the mask of the "good" pixels
    goodmask = (subexp.flux > 0.0) & (subexp.mask == False)

    # Make xcen/ycen for the larger subimage
    xcen = xcenS + (hwidth1-hwidthS)
    ycen = ycenS + (hwidth1-hwidthS)
    
    # Construct X- and Y-arrays
    #xx = np.zeros([2*hwidth+1,2*hwidth+1],'f')
    #for j in xrange(2*hwidth+1):
    #    xx[:,j] = j
    #yy = np.zeros([2*hwidth+1,2*hwidth+1],'f')
    #for j in xrange(2*hwidth+1):
    #    yy[j,:] = j   
    yy, xx = np.indices([2*hwidth1+1,2*hwidth1+1],'f')
    #print "need to test this np.indices code"
    
    # Calculate the flux in the subimage
    morph['flux'] = np.sum(subexp.flux * goodmask)

    # Getting the contour at 1/2 maximum flux
    #  add a perimeter of zeros to ensure that
    #  we get back a contour even if it hits the edge
    ny1, nx1 = subexp.flux.shape
    tflux = np.zeros([ny1+2,nx1+2],'f')
    tflux[1:ny1+1,1:nx1+1] = subexp.flux
    contour = get_contour(tflux,maxim*0.5,[ycen+1,xcen+1])
    # Offset the coordinates to the original image
    if len(contour) > 0:
        contour -= 1
            
    # Good contour, make the measurements
    if len(contour) > 0:
        # Getting the path
        xpath = contour[:,1]
        ypath = contour[:,0]
        xmnpath = np.mean(xpath)
        ymnpath = np.mean(ypath)
    
        # Calculating the FWHM
        dist = np.sqrt((xpath-xmnpath)**2.0 + (ypath-ymnpath)**2.0)  
        fwhm1 = 2.0 * np.mean(dist)
        morph['contour_fwhm'] = fwhm1
            
        # Measuring "ellipticity", (1-a/b)
        elip = 2*np.std(dist-fwhm1)/fwhm1
        morph['contour_elip'] = elip
    
        # Calculate the position angle
        #  angle for point where dist is maximum
        #  angle from positive x-axis
        maxind = np.where(dist == dist.max())[0]
        theta = np.rad2deg( np.arctan2(ypath[maxind]-ymnpath,xpath[maxind]-xmnpath) )
        theta = theta[0] % 360
        # want values between -90 and +90
        if theta > 180:
            theta -= 180
        if theta > 90:
            theta -= 180
        morph['contour_theta'] = theta
            
    else:
        morph['contour_fwhm'] = np.nan
        morph['contour_elip'] = np.nan
        morph['contour_theta'] = np.nan
        # THIS CAN HAPPEN IF THE CONTOUR HITS THE EDGE

    # Computing the "round" factor
    # round = difference of the heights of the two 1D Gaussians
    #         -------------------------------------------------
    #                 average of the two 1D Gaussians
    #
    # Where the 1D Gaussians are of the marginal sums, i.e. sum
    # along either the x or y dimensions
    # round~0 is good
    # round<0 object elongated in x-direction
    # round>0 object elongated in y-direction
    htx = np.max(np.sum(subexp.flux*(1-subexp.mask),axis=0))
    hty = np.max(np.sum(subexp.flux*(1-subexp.mask),axis=1))
    morph['round'] = (hty-htx)/np.mean([htx,hty])

    # 2D Gaussian fitting???

    # Create a "window" mask, only including pixels within 1/2 maximum contour
    if len(contour) > 0:
        contmask = measure.grid_points_in_poly(subexp.flux.shape,contour)
        cgoodmask = goodmask*contmask
    else:
        cgoodmask = np.copy(goodmask)
            
    # Second MOMENTS of the windowed image
    #  The factor of 3.33 corrects for the fact that we are missing
    #  a decent chunk of the 2D Gaussian.  I derived this empirically
    #  using simulations.
    posflux = np.sum( subexp.flux*cgoodmask )
    ixx = np.sum( subexp.flux*cgoodmask * (xx-xcen)**2 ) / posflux * 3.33
    iyy = np.sum( subexp.flux*cgoodmask * (yy-ycen)**2 ) / posflux * 3.33
    ixy = np.sum( subexp.flux*cgoodmask * (xx-xcen) * (yy-ycen) ) / posflux * 3.33
    morph['ixx'] = ixx
    morph['iyy'] = iyy
    morph['ixy'] = ixy
        
    # Computing semi-major, semi-minor and theta from the moments
    #  The SExtractor manual has the solutions to this on pg.31
    # semi-major axis = (ixx+iyy)/2 + sqrt( ((ixx-iyy)/2)^2 + ixy^2 ) = sigx^2
    # semi-minor axis = (ixx+iyy)/2 - sqrt( ((ixx-iyy)/2)^2 + ixy^2 ) = sigy^2
    # tan(2*theta) = 2*ixy/(ixx-iyy)
    # two solutions between -90 and +90 with opposite signs
    # by definition, theta is the position angle for which ixx_theta is maximized
    # (in rotated frame), counter-clockwise from positive x-axis.
    # so theta is the solution the tan equation that has the same sign as ixy
    siga = np.sqrt( (ixx+iyy)/2 + np.sqrt( ((ixx-iyy)/2)**2 + ixy**2 ) )
    sigb = np.sqrt( (ixx+iyy)/2 - np.sqrt( ((ixx-iyy)/2)**2 + ixy**2 ) ) \
           if np.sqrt( ((ixx-iyy)/2)**2 + ixy**2 ) <= (ixx+iyy)/2 else 0.1
    if ixx != iyy:
        theta = np.rad2deg( np.arctan2(2*ixy,ixx-iyy) / 2 )
        theta = np.abs(theta)*np.sign(ixy)
    else:
        theta = 0.0
    morph['siga'] = siga
    morph['sigb'] = sigb
    morph['theta'] = theta
    # THETA is more accurate than CONTOUR_THETA
    
    # Construct Gaussian model for Gaussian weighted photometry
    #  https://en.wikipedia.org/wiki/Gaussian_function
    #  theta in the wikipedia equation is CLOCKWISE so add - before theta
    thetarad = np.deg2rad(theta)
    a = ((np.cos(-thetarad)**2) / (2*siga**2)) + ((np.sin(-thetarad)**2) / (2*sigb**2))
    b = -((np.sin(-2*thetarad)) / (4*siga**2)) + ((np.sin(-2*thetarad)) / (4*sigb**2))
    c = ((np.sin(-thetarad)**2) / (2*siga**2)) + ((np.cos(-thetarad)**2) / (2*sigb**2))
    g = np.exp(-(a*(xx-xcen)**2 + 2*b*(xx-xcen)*(yy-ycen) + c*(yy-ycen)**2))
    g /= np.sum(g)

    # Gaussian weighted flux
    #  from Valdes (2007)
    gausswtflux = np.sum( g*subexp.flux*cgoodmask/subexp.noise**2 )
    gausswtflux /= np.sum( g**2 * cgoodmask/subexp.noise**2 )
    morph['gausswtflux'] = gausswtflux
    
    # Compute gaussian scaling factor using a 
    # weighted mean of the fraction flux/gaussian
    #  weight by (S/N)^2
    wt = np.zeros(subexp.flux.shape,'f')
    wt[:,:] = (subexp.flux/subexp.noise)**2
    # only use values where the gaussian is large enough
    # and the image isn't masked out
    gmask = (g > np.max(g)*0.05) & (cgoodmask == True)
    wt *= gmask
    wt /= np.sum(wt)
    # set values of the gaussian in the denominator
    #   to 1.0 where they are too low
    gdenom = np.copy(g)
    gdenom[gmask == False] = 1
    gauss_scale = np.sum( subexp.flux*wt / gdenom )
    morph['gaussflux'] = gauss_scale
        
    # compute chi-squared
    chisq = np.sum( ((subexp.flux-g*gausswtflux)*gmask / subexp.noise)**2 ) / np.sum(gmask)
    morph['chisq'] = chisq
        
    return morph

def get_morph(exp,peaks,noerrors=False):
    """ Measure morphological parameters of a source in a small image
    """
    subim = exp.flux - exp.background
    npeaks = len(peaks)
    hwidthS = 3
    hwidth = 10
    # Create morph structure
    # peak segment
    if 'bbox_x0' in peaks.dtype.names:
        dt = np.dtype([('x0',int),('y0',int),('nsig',float),('maxflux',float),('area',int),
                       ('bbox_x0',int),('bbox_y0',int),('bbox_x1',int),('bbox_y1',int),
                       ('x',float),('xerr',float),('y',float),
                       ('yerr',float),('flux',float),('round',float),('max',float),('contour_fwhm',float),
                       ('contour_elip',float),('contour_theta',float),('ixx',float),('iyy',float),
                       ('ixy',float),('siga',float),('sigb',float),('theta',float),
                       ('gausswtflux',float),('gaussflux',float),('chisq',float)])
        morph = np.zeros(npeaks,dtype=dt)
        morph['x0'] = peaks['xcen']
        morph['y0'] = peaks['ycen']
        morph['nsig'] = peaks['nsig']
        morph['maxflux'] = peaks['maxflux']
        morph['area'] = peaks['area']
        morph['bbox_x0'] = peaks['bbox_x0']
        morph['bbox_y0'] = peaks['bbox_y0']
        morph['bbox_x1'] = peaks['bbox_x1']
        morph['bbox_y1'] = peaks['bbox_y1']
    # Delta peaks
    else:
        dt = np.dtype([('x0',int),('y0',int),('nsig',float),('x',float),('xerr',float),('y',float),
                       ('yerr',float),('flux',float),('round',float),('max',float),('contour_fwhm',float),
                       ('contour_elip',float),('contour_theta',float),('ixx',float),('iyy',float),
                       ('ixy',float),('siga',float),('sigb',float),('theta',float),
                       ('gausswtflux',float),('gaussflux',float),('chisq',float)])
        morph = np.zeros(npeaks,dtype=dt)
        morph['x0'] = peaks['xcen']
        morph['y0'] = peaks['ycen']
        morph['nsig'] = peaks['nsig']
    
    # Loop through the peaks
    for i in range(npeaks):
        morph1in = morph[i]

        #pylab.ion()
        #plt.imshow(smim,interpolation='nearest',origin='lower')
        #plt.plot(xdetect,ydetect,'r+')
        #import pdb; pdb.set_trace()
    
        #print morph1in.contour_fwhm
        morph1out = get_morph_single(subim,exp,morph1in,hwidth=hwidth,noerrors=noerrors)
        #print morph1out.contour_fwhm
        # The box is too small, try twice the size
        #if (np.isfinite(morph1out['contour_fwhm']) is False) | (morph1out['siga'] > (hwidth-1)) \
        #   | (morph1out['sigb'] > (hwidth-1)):
        #    morph1out = get_morph_single(subim,exp,morph1in,hwidth=2*hwidth,noerrors=noerrors)
        morph[i] = morph1out

    return morph

def get_morphclusters(morph,ngroups=3):
    """ Find clusters in the shapes of the sources
    """

    # Get sources with "good" values for everything
    gd, = np.where(np.isfinite(morph['fwhm']) & np.isfinite(morph['round']) &
                   np.isfinite(morph['elip']) & np.isfinite(morph['theta']) &
                   np.isfinite(morph['ixx']) & np.isfinite(morph['iyy']) &
                   np.isfinite(morph['ixy']))           
    ngd = len(gd)
    
    # Use K-means
    #  Try three groups: CRs, stars, and galaxies
    features = np.zeros([ngd,7],'f')
    features[:,0] = morph['fwhm'][gd]
    features[:,1] = morph['round'][gd]
    features[:,2] = morph['elip'][gd]
    features[:,3] = morph['theta'][gd]
    features[:,4] = morph['ixx'][gd]
    features[:,5] = morph['iyy'][gd]
    features[:,6] = morph['ixy'][gd]
    #features[np.isnan(features)] = 999999.
    # maybe only use sources with non-NAN values
    whitened = vq.whiten(features)
    #gp = gaps.gap(whitened)
    #print gp
    res, idx = vq.kmeans2(whitened,ngroups,missing='warn')
    # Create cluster structure
    nclusters = np.max(idx)+1
    dt = np.dtype([('nsources',int),('med_fwhm',float),('sig_fwhm',float),('med_round',float),
                   ('sig_round',float),('med_elip',float),('sig_elip',float),
                   ('med_theta',float),('sig_theta',float),('med_ixx',float),
                   ('sig_ixx',float),('med_iyy',float),('sig_iyy',float),
                   ('med_ixy',float),('sig_ixy',float)])
    clusters = np.zeros(nclusters,dtype=dt)
    for i in xrange(nclusters):
        grp, = np.where(idx == i)
        clusters['nsources'][i] = len(grp)
        if len(grp) > 1:
            clusters['med_fwhm'][i] = np.median(features[grp,0])
            clusters['sig_fwhm'][i] = mad(features[grp,0])
            clusters['med_round'][i] = np.median(features[grp,1])
            clusters['sig_round'][i] = mad(features[grp,1])
            clusters['med_elip'][i] = np.median(features[grp,2])
            clusters['sig_elip'][i] = mad(features[grp,2])
            clusters['med_theta'][i] = np.median(features[grp,3])
            clusters['sig_theta'][i] = mad(features[grp,3])
            clusters['med_ixx'][i] = np.median(features[grp,4])
            clusters['sig_ixx'][i] = mad(features[grp,4])
            clusters['med_iyy'][i] = np.median(features[grp,5])
            clusters['sig_iyy'][i] = mad(features[grp,5])
            clusters['med_ixy'][i] = np.median(features[grp,6])
            clusters['sig_ixy'][i] = mad(features[grp,6])
    # Get the non-empty clusters
    gdclusters, = np.where(clusters['nsources'] > 1)
    clusters = clusters[gdclusters]
    
    return clusters
        
def imfwhm(exp):
    """ Measure the image PSF FWHM
        bkg has the image (data), background (background) and noise in the
        background (background_rms)
    """
    # Step 1. Detect with delta function
    peaks = detect_delta(exp)
    # Step 2. Aperture photometry
    #  not sure this is needed
    # Step 3. Measure morphology/moments
    morph = get_morph(exp,peaks)
    # Step 4. Get cluster of stars and measure FWHM
    clusters = get_morphclusters(morph)
    # Step 5. Now decide which cluster to use for "stars" and measure FWHM
    # three types of groups
    # 1.) CRs, hot/bad pixels
    #   -large variety of sizes, but lots will be very "peaky"
    #     and have small widths
    # 2.) stars
    #   -wider than CRs but thinner than galaxies
    #   -should be very uniform in all parameters
    # 3.) galaxies
    #   -should be the widest
    #   -should a larger dispersion of shape, size, ellipticity and angle
    # But how do we know if all of these three groups are presented?
    
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run DAOPHOT PSF photometry on FITS image.")

    parser.add_argument('--file', '-f', action="store", help="The FITS file to process", default=None)
    #parser.add_argument('--datarepo', '-d', action="store", help="The data repository directory", default="/data/lsst/decam/redux/cp/cosmos/")
    #parser.add_argument('--outfile', '-o', action="store", help="The output filename for the metrics.", default="qametrics.csv")
    parser.add_argument('--verbose', '-v', action="store_true", help="Print out the data as it is gathered.", default=False)
    parser.add_argument('--clobber', '-c', action="store_true", help="Overwrite the output file if it already exists.", default=False)

    args = parser.parse_args()
    file = args.file
    
    print "Running DAOPHOT PSF photometry on ", file

    # Figure out the output file
    dir = os.path.dirname(file)
    if dir == '':
        dir = '.'
    base, ext = os.path.splitext(os.path.basename(file))
    outfile = dir+'/'+base+'_cat.fits'
    if os.path.exists(outfile) and not args.clobber:
        print outfile," EXISTS and --clobber not set"
        sys.exit()
    
    # Load the file
    im, head = fits.getdata(file,0,header=True)
    nx, ny = im.shape
    print "Dimensions", im.shape

    # Make new "image" or "exposure" class that has:
    # -flux image
    # -noise/variance image
    # -mask image
    # -background image
    # -FITS header
    # -psf
    # -wcs?
    
    # Get the background
    print "Computing background image"
    # THIS TAKES WAY TOO LONG  
    bkg = photutils.Background2D(im, (nx/10, ny/10), filter_size=1,method='median')
    subim = im-bkg.background

    # Create the mask
    #  mask=1   bad
    #  mask=0   good
    #try:
    #    saturate = head['saturate']
    #except:
    #    saturate = 60000.0
    #mask = im > 0.9*saturate
    mask = create_mask(im,head)
    
    # Create the noise mask
    noise = create_noise(im,head,mask)
    #try:
    #    # electrons
    #    rdnoise = head['rdnoise']
    #except:
    #    rdnoise = 0.0
    #try:
    #    # electrons/ADU
    #    gain = head['gain']
    #except:
    #    # Sky sigma should be sqrt(im*gain)/gain = sqrt(im/gain)
    #    # so we can roughly (ignoring rdnoise) back out gain
    #    # gain ~ im/skysig^2
    #    sig = mad(im)
    #    skysig = sig
    #    skymed = np.median(im)
    #    gain = skymed/skysig**2
    #    # Now try to remove rdnoise, normally in electrons
    #    if rdnoise > 0.0:
    #        skysig = np.sqrt(sig**2 - (rdnoise/gain)**2)
    #        gain = skymed/skysig**2
    ## Make the simple noise model
    #noise = np.sqrt(im/gain+(rdnoise/gain)**2)
    #noise[mask > 0] = 999999.
    
    # Create exposure object
    #  background_rms is way too high
    #  NEED A BETTER "NOISE" IMAGE!!
    #exp = Exposure(im,bkg.background_rms,mask,bkg.background,head)
    exp = Exposure(im,noise,mask,bkg.background,head)
    
    # Measure PSF FWHM
    #print "Measuring PSF FWHM"
    #fwhm = imfwhm(exp)

    # Step 1. Detect with delta function
    peaks = detect_delta(exp)
    # Step 2. Aperture photometry

    # Step 3. Measure morphology/moments
    morph = get_morph(exp,peaks)
    
    # output to csv file
    print "Writing output catalog to ", outfile
    # Delete output file if it exists and clobber set
    if os.path.exists(outfile) and args.clobber:
        os.remove(outfile)
    fits.writeto(outfile,morph)
    #print "Writing outputs to", args.outfile
    #data.tofile(args.outfile,sep='\n')

    # Add header line
    #f = open(args.outfile,'r')  # read it all back in first
    #temp = f.read()
    #f.close()
    ##  now write out with header line
    #f = open(args.outfile, 'w')
    #f.write(morph(dt.names)+'\n')
    #f.write(temp)
    #f.close()
