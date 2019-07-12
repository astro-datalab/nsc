#!/usr/bin/env python
#
# UTILS.PY - Utility functions.
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180823'  # yyyymmdd

import re
import logging
import os
#import sys
import numpy as np
#import warnings
#from astropy.io import fits
#from astropy.utils.exceptions import AstropyWarning
#from astropy.table import Table, Column
#import time
#import shutil
#import glob
#import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve
#import astropy.stats
#import struct
#import tempfile

# Ignore these warnings, it's a bug
#warnings.filterwarnings("ignore", message="numpy.dtype size changed")
#warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def lt(x,limit):
    """Takes the lesser of x or limit"""
    if np.array(x).size>1:
        out = [i if (i<limit) else limit for i in x]
    else:
        out = x if (x<limit) else limit
    if type(x) is np.ndarray: return np.array(out)
    return out
    
def gt(x,limit):
    """Takes the greater of x or limit"""
    if np.array(x).size>1:
        out = [i if (i>limit) else limit for i in x]
    else:
        out = x if (x>limit) else limit
    if type(x) is np.ndarray: return np.array(out)
    return out        

def limit(x,llimit,ulimit):
    """Require x to be within upper and lower limits"""
    return lt(gt(x,llimit),ulimit)

# Standard grep function that works on string list
def grep(lines,expr,index=False):
    '''
    Similar to the standard unit "grep" but run on a list of strings.
    Returns a list of the matching lines unless index=True is set,
    then it returns the indices.

    Parameters
    ----------
    lines : list
          The list of string lines to check.
    expr  : str
          Scalar string expression to search for.
    index : bool, optional
          If this is ``True`` then the indices of matching lines will be
          returned instead of the actual lines.  index is ``False`` by default.

    Returns
    -------
    out : list
        The list of matching lines or indices.

    Example
    -------

    Search for a string and return the matching lines:

    .. code-block:: python

        mlines = grep(lines,"hello")

    Search for a string and return the indices of the matching lines:

    .. code-block:: python

        index = grep(lines,"hello",index=True)

    '''
    out = []
    cnt = 0L
    for l in lines:
        m = re.search(expr,l)
        if m != None:
            if index is False:
                out.append(l)
            else:
                out.append(cnt)
        cnt = cnt+1
    return out


# Read in all lines of files
def readlines(fil=None):
    '''
    Read in all lines of a file.
    
    Parameters
    ----------
    file : str
         The name of the file to load.
   
    Returns
    -------
    lines : list
          The list of lines from the file

    Example
    -------

    .. code-block:: python

       lines = readlines("file.txt")

    '''
    if fil is None:
        print("File not input")
        return
    f = open(fil,'r')
    lines = f.readlines()
    f.close()
    return lines


# Write all lines to file
def writelines(filename=None,lines=None,overwrite=True):
    '''
    Write a list of lines to a file.
    
    Parameters
    ----------
    filename : str
        The filename to write the lines to.
    lines : list
         The list of lines to write to a file.
    overwrite : bool, optional, default is True
        If the output file already exists, then overwrite it.

    Returns
    -------
    Nothing is returned.  The lines are written to `fil`.

    Example
    -------

    .. code-block:: python

       writelines("file.txt",lines)

    '''
    # Not enough inputs
    if lines is None:
        print("No lines input")
        return
    if filename is None:
        print("No file name input")
        return
    # Check if the file exists already
    if os.path.exists(filename):
        if overwrite is True:
            os.remove(filename)
        else:
            print(filename+" already exists and overwrite=False")
            return
    # Write the file
    f = open(filename,'w')
    f.writelines(lines)
    f.close()


# Remove indices from a list
def remove_indices(lst,index):
    '''
    This will remove elements from a list given their indices.

    Parameters
    ----------
    lst : list
          The list from which to remove elements.
    index : list or array
          The list or array of indices to remove.

    Returns
    -------
    newlst : list
           The new list with indices removed.

    Example
    -------

    Remove indices 1 and 5 from array `arr`.

    .. code-block:: python

        index = [1,5]
        arr  = range(10)
        arr2 = remove_indices(arr,index)
        print(arr)
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    '''
    newlst = []
    for i in range(len(lst)):
       if i not in index: newlst.append(lst[i])
    return newlst


# Little function used by numlines
def blocks(files, size=65536):
    '''
    This is a small utility function used by numlines()
    '''
    while True:
        b = files.read(size)
        if not b: break
        yield b


# Read number of lines in a file
def numlines(fil):
    '''
    This function quickly counts the number of lines in a file.

    Parameters
    ----------
    fil : str
          The filename to check the number of lines.

    Returns
    -------
    nlines : int
           The number of lines in `fil`.

    Example
    -------

    .. code-block:: python

        n = numlines("file.txt")

    '''
    with open(fil, "r") as f:
        return (sum(bl.count("\n") for bl in blocks(f)))

    # Could also use this
    #count=0
    #for line in open(fil): count += 1


# Set up basic logging to screen
def basiclogger(name=None):
    '''
    This sets up a basic logger that writes just to the screen.
    '''
    if name is None: name = "log"
    logger = logging.getLogger(name)
    # Only add a handler if none exists
    #  the logger might already have been created
    if len(logger.handlers)==0:
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(levelname)-2s %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger

def rotsph(lon,lat,clon,clat,anode=None,reverse=False,original=False):
    '''
    This rotates a spherical coordinate system to a new pole

    I got the equations for this from the paper
    Calabretta et al. 2002, A&A, 395, 1077
    Equation 5.

    Also, see rotate_lb.pro that uses a matrix method
    and for which you can specify the equator you'd like.
    rotsph.pro is faster than rotate_lb.pro
    By default, the origin is the point where the two equators
    cross (unless =ANODE is set).
    This should give you the same result (within ~1E-10")
    rotate_lb,lon,lat,[clon,clat],[clon+90,0.],nlon,nlat

    Parameters
    ----------
    lon       Array of longitudes to be rotated
    lat       Array of latitudes to be rotated
    clon      Longitude of the new NORTH POLE in the old coordinate system
    clat      Latitude of the new NORTH POLE in the old coordinate system
    =anode    The "Ascending Node" which is the longitude (in the new
             system) of the first point where the old equator cross
             the new one.  This sets the zero-point of the new
             longitude.  By default the zero-point of the new
             coordinate system is the point where the two equators
             cross.
    /original Set the new longitude zero-point to be clon (if clat>0)
             and clon+180 (if clat<0).  This is the way it was
             originally done.  DON'T USE WITH "ANODE"
    /stp      Stop at the end of the program
    /reverse  The reverse operation.  In that case (nlon,nlat) should be input
           as (lon,lat). E.g.

           rotsph,ra,dec,cra,cdec,nlon,nlat
           rotsph,nlon,nlat,cra,cdec,nra,ndec,/reverse
           
           (ra,dec) and (nra,ndec) should be identical to 1E-10.

    Returns
    -------
    nlon  Array of rotated longitudes
    nlat  Array of rotated latitudes

    '''

    radeg = 180.0/np.pi

    alphap = np.array(clon/radeg)
    deltap = np.array(clat/radeg)
    phip = np.array(90.0/radeg)
    if original: phip = np.array(180.0/radeg)   # original way
    thetap = np.array(90.0/radeg)

    # By default the origin of the new coordinate system is the point
    # where the two equators cross for the first time
    #  Unless /original is set.  Then the zero-point is at clon
    #   (if clat>0) and clon+180 (if clat<0)

    # NORMAL
    if reverse is False:
        alpha = np.array(lon/radeg)
        delta = np.array(lat/radeg)

        # arg(x,y) but atan(y,x)
        phi = phip + np.arctan2( -np.cos(delta)*np.sin(alpha-alphap), np.sin(delta)*np.cos(deltap)- \
                                 np.cos(delta)*np.sin(deltap)*np.cos(alpha-alphap) )

        theta = np.arcsin( limit((np.sin(delta)*np.sin(deltap)+np.cos(delta)*np.cos(deltap)*np.cos(alpha-alphap)),-1,1) )

        # Preparing the output
        nlon = phi*radeg
        nlat = theta*radeg

        # Ascending Node
        #  By default the origin of nlon is the point where the two equators
        #  cross the first time
        if anode is not None: nlon += anode

    # REVERSE
    else:
        phi = np.array(lon/radeg)
        theta = np.array(lat/radeg)

        # Ascending Node
        if anode is not None: phi = (lon-anode)/radeg

        # arg(x,y) but atan(y,x)
        alpha = alphap + np.arctan2( -np.cos(theta)*np.sin(phi-phip), np.sin(theta)*np.cos(deltap) - \
                                     np.cos(theta)*np.sin(deltap)*np.cos(phi-phip))
        delta = np.arcsin( np.sin(theta)*np.sin(deltap) + np.cos(theta)*np.cos(deltap)*np.cos(phi-phip) )

        # Preparing the output
        nlon = alpha*radeg
        nlat = delta*radeg

    # Want everything less than 360.0
    nlon = nlon % 360.0

    # Make negative points positive
    bd, = np.where(nlon < 0.0)
    if len(bd)>0:
        nlon[bd] = nlon[bd]+360.0

    return nlon, nlat


def rotsphcen(lon,lat,clon,clat,polar=False,gnomic=False,reverse=False):
    '''
    This is very similar to rotsph.pro except that the coordinates
    input are not for the north pole but for the new equator.
    Everything is in DEGREES.
    
    Parameters
    ----------
    lon       Array of longitudes to be rotated
    lat       Array of latitudes to be rotated
    clon      Longitude of the new EQUATOR in the old coordinate system
    clat      Latitude of the new EQUATOR in the old coordinate system
    /polar    Return polar coordinates (rad,phi) instead of LON/LAT.
    phi starts at North.
    /gnomic   Also do a gnomic (tangent plane) projection.
    /reverse  The reverse operation.  In that case (nlon,nlat) should be input
    as (lon,lat). E.g.

           rotsphcen,ra,dec,cra,cdec,nlon,nlat
           rotsphcen,nlon,nlat,cra,cdec,nra,ndec,/reverse
           
           (ra,dec) and (nra,ndec) should be identical to 1E-10.

    Returns
    -------
    nlon  Array of rotated longitudes.  If /polar then this is PHI
       the polar angle (measured from N toward E).
       
    nlat  Array of rotated latitudes.  If /polar then this is RAD
       the polar radial distance.
    '''

    radeg = 180.0/np.pi

    # NOT polar coordinates
    if (polar is False) and (gnomic is False):

        # Get coordinates for the north pole
        np_lon = np.array(clon)
        np_lat = np.array(clat+90.0)
        if (np_lat > 90.0):
            np_lon = np.array(clon+180.0)
            np_lon = np_lon % 360.0
            np_lat = 90.0-clat

        # Run rotsph.pro
        # NORMAL
        if reverse is False:
            nlon, nlat = rotsph(lon,lat,np_lon,np_lat,original=True)
        # REVERSE
        else:
            nlon, nlat = rotsph(lon,lat,np_lon,np_lat,reverse=True,original=True)

            # need to flip them around by 180 deg b/c the zero-point
            #  is set by the NP lon
            if ((clat+90.0) > 90.0):
                nlon = (nlon+180.0) % 360.0
                nlat = -nlat

        # Make the longitudes continuous
        nlon = (nlon+180.0) % 360.0
        nlon = nlon-180.0


    # POLAR or GNOMIC
    else:

        # Making polar coordinates

        #-----------------
        # NORMAL
        #------------------
        if reverse is False:
            # Run rotsph.pro and specify the center of the field (the origin) as the
            #  the Npole
            phi, theta = rotsph(lon,lat,clon[0],clat[0],original=True)
            # phi is now going clockwise and from South
            orig_phi = phi
            phi = -phi+180.0      # counterclockwise
            phi = phi % 360.0
            rad = 90.0-theta

            # Making gnomic projection
            if gnomic:
                # Scale the radius
                rad = radeg * np.cos(theta/radeg)/np.sin(theta/radeg)
                # Now convert from gnomic polar to X/Y
                # phi is from N toward E
                # x = R*sin(phi)
                # y = R*cos(phi)
                nlon = rad*np.sin(phi/radeg)
                nlat = rad*np.cos(phi/radeg)

            # Output polar coordinates
            if polar:
                nlon = phi
                nlat = rad

        #-----------------
        # REVERSE
        #-----------------
        else:

            # Polar
            if polar:
                phi = lon
                rad = lat
                theta = 90.0-rad

            # Removing gnomic projection
            if gnomic:
                # Now convert from X/Y to gnomic polar
                # phi is from N toward E
                # x = R*sin(phi)
                # y = R*cos(phi)
                #nlon = rad*sin(phi/radeg)
                #nlat = rad*cos(phi/radeg)
                rad = sqrt(lon**2.0+lat**2.0)
                phi = radeg*np.arctan2(lon,lat)      # in degrees
                # Scale the radius
                #rad = radeg * cos(theta/radeg)/sin(theta/radeg)
                theta = radeg*np.arctan(radeg/rad)   # in degrees

            # phi is now going clockwise and from South
            phi = -phi+180.0       # reverse phi

            #Run rotsph.pro and specify the center of the field (the origin) as the
            #  the Npole
            nlon, nlat = rotsph(phi,theta,clon[0],clat[0],reverse=True,original=True)

    return nlon, nlat



