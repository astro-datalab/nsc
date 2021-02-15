#!/usr/bin/env python
#
# PHOT.PY - SExtractor and DAOPHOT routines
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180823'  # yyyymmdd

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, Column
import time
import shutil
import subprocess
import logging
#from scipy.signal import convolve2d
from scipy.ndimage.filters import convolve
import astropy.stats
import struct
import tempfile
from dlnpyutils.utils import *

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


# Parse the DAOPHOT PSF profile errors
def parseprofs(lines):
    '''
    This parses the PSF profile errors output from the DAOPHOT PSF program.
    It returns a numpy structured array with ID, SIG, FLAG for the PSF stars.

    This is an example of the PDF profile error lines:
    1044  0.010      1039  0.010       304  0.013      1118  0.020       119  0.027   
     610  0.012       580  0.013       373  0.010       617  0.017      1087  0.027   
     696  0.010       229  0.012       170  0.016       628  0.018      1211  0.030 

    Parameters
    ----------
    lines : list
          The list of string lines from the DAOPHOT PSF program for the profile errors.

    Returns
    -------
    profs : numpy structured array
          The catalog containing ID, SIG, and FLAG columns for the PSF stars.

    Example
    -------

    .. code-block:: python

        profs = parseprofs(lines)

    '''

    # From psf.f
    # SIG is the root-mean-square scatter about the best-fitting analytic
    # function averaged over the central disk of radius FITRAD, expressed
    # as a fraction of the peak amplitude of the analytic model.

    dtype = np.dtype([('ID',int),('SIG',float),('FLAG',np.str_,10)])
    profs = np.zeros(len(lines)*5,dtype=dtype)
    profs['ID'] = -1
    cnt = 0
    for i in range(len(lines)):
        l = lines[i].rstrip()
        if l != "":
            # Loop through five columns
            for j in range(5):
                line1 = l[j*17:j*17+17]
                id1 = line1[0:7]
                sig1 = line1[7:14]
                flag1 = line1[14:17]
                if sig1 == " satura":
                    sig1 = 99.99
                    flag1 = "saturated"
                if id1.strip() != "":
                    profs[cnt]['ID'] = int(id1)
                    profs[cnt]['SIG'] = float(sig1)
                    profs[cnt]['FLAG'] = flag1.strip()
                    cnt = cnt + 1
    # Trimming any blank ones
    gd = (profs['ID'] > -1)
    profs = profs[gd]
    return profs


# Parse the DAOPHOT PSF parameter errors
def parsepars(lines):
    '''
    This parses the PSF parameter errors output from the DAOPHOT PSF program.
    It returns a list (one element per line) where each element constains
    a list of the 3-5 parameters.

    This is an example of lines of the PSF parameter errors:
           Chi    Parameters...
    >>   0.0319   1.79190   1.69498   
    >>   0.0382   1.21314   1.26585  -0.00693   
    >>   0.0215   1.62418   1.52379  -0.00385   
    >>   0.0196   1.66754   1.57059  -0.00304   
    >>   0.0543   1.41140   1.30613  -0.00290   
    >>   0.0197   1.68487   1.58727   0.68797  -0.00305   

    Parameters
    ----------
    lines : list
          The list of string lines from the DAOPHOT PSF program for the parameter errors.

    Returns
    -------
    out : list
        The list of lists containing the individual parameters.
    chi : list
        The list of chi values per line.

    Example
    -------

    .. code-block:: python

        out, chi = parseparse(lines)

    '''
    # Only want lines with ">>"
    #  somtimes theare is a  "Failed to converge." line
    lines2 = grep(lines,">>")
    if len(lines2)==0:
        print("No lines with >> found")
        return None, None
    out = []
    chi = []
    for i in range(len(lines2)):
        line1 = lines2[i].strip()
        if line1[0:2] == ">>": line1=line1[2:]  # strip leading >>
        line1.strip()
        arr = line1.split()                # split on whitespace
        if len(arr)>0:
            chi.append(float(arr[0]))
            chi.append(float(arr[0]))
            out.append(arr)
    return out, chi


# Write DAOPHOT apertures files
def aperswrite(filename=None,apertures=None):
    '''
    This program creates a DAOPHOT file with apertures with an array/list of apertures.
    The last two are assumed to be the inner and outer sky apertures.

    Parameters
    ----------
    filename : str
        The filename for the apertures.
    apertures : list or array
        The array of apertures.

    Returns
    -------
    Nothing is returned but the apertures file is created.

    Example
    -------

    .. code-block:: python

        aperswrite("photo.opt",apertures)

    '''
    # Not enough inputs
    if filename is None:
        print("No file name input")
        return
    if apertures is None:
        print("No apertures input")
        return

    # Make apertures file
    nap = len(apertures)
    if nap<3:
        print("Only "+str(nap)+" apertures input.  Need at least 3")
        return
    f = open(filename,'w')
    for i in range(nap-2):
        # use hexidecimal for aperture id, 2 digits, first starts with A
        id = hex(160+i+1)
        id = id[2:].capitalize()
        f.write("%2s = %7.4f\n" % (id,apertures[i]))
    f.write("IS = %7.4f\n" % apertures[nap-2])
    f.write("OS = %7.4f\n" % apertures[nap-1])
    f.close()


# Read DAOPHOT files
def daoread(fil):
    '''
    This program reads in DAOPHOT-style files and return an astropy table.
    The supported types are .coo, .lst, .ap (in development), and .als.

    Parameters
    ----------
    fil : str
        The filename of the DAOPHOT catalog file.

    Returns
    -------
    cat : astropy table
        The DAOPHOT catalog as an astropy table.

    Example
    -------

    Load an ALLSTAR catalog file:

    .. code-block:: python

        cat = daoread("image1.als")

    '''

    # Not enough inputs
    if fil is None:
        print("No file name input")
        return None
    # Make sure the file exists
    if os.path.exists(fil) is False:
        print(fil+" NOT found")
        return None
    lines = readlines(fil)
    nstars = len(lines)-3
    if nstars == 0:
        print("No stars in "+file)
        return None
    # Check header
    line2 = lines[1]
    nl = int(line2.strip().split(' ')[0])
    # NL  is a code indicating the file type:
    # NL = 3 a group file
    # NL = 2 an aperture photometry file
    # NL = 1 other (output from FIND, PEAK, or NSTAR) or ALLSTAR, DAOGROW
    # NL = 0 a file without a header
    
    # Check number of columns
    arr1 = lines[3].split()
    if len(arr1)==0: arr1 = lines[4].split()
    ncols = len(arr1)

    # NL = 1  coo file
    if (nl==1) & (ncols==7):
        dtype = np.dtype([('ID',long),('X',float),('Y',float),('MAG',float),('SHARP',float),('ROUND',float),('ROUND2',float)])
        cat = np.zeros(nstars,dtype=dtype)
        lengths = np.array([7,9,9,9,9,9,9])
        lo = np.concatenate((np.array([0]), np.cumsum(lengths[0:-1])))
        hi = lo+lengths
        names = cat.dtype.names
        for i in range(nstars):
            line1 = lines[i+3]
            for j in range(len(names)):
                cat[i][names[j]] = np.array(line1[lo[j]:hi[j]],dtype=dtype[names[j]])
    # NL = 1  tot file
    if (nl==1) & (ncols==9) & (arr1[-1].isdigit() is True):
        # NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
        #  1  2046  4094   117.7 38652.0   13.12    3.00    3.91    1.55    6.00
        #  
        #     11  454.570   37.310  13.9710   0.0084  164.683   14.040  -0.0690        6
        #     36  287.280   93.860  14.5110   0.0126  165.018   14.580  -0.0690        6
        dtype = np.dtype([('ID',long),('X',float),('Y',float),('MAG',float),('ERR',float),('SKY',float),('MAGFAP',float),('APCORR',float),('FINALAP',int)])
        cat = np.zeros(nstars,dtype=dtype)
        lengths = np.array([7,9,9,9,9,9,9,9,9])
        lo = np.concatenate((np.array([0]), np.cumsum(lengths[0:-1])))
        hi = lo+lengths
        names = cat.dtype.names
        for i in range(nstars):
            line1 = lines[i+3]
            for j in range(len(names)):
                cat[i][names[j]] = np.array(line1[lo[j]:hi[j]],dtype=dtype[names[j]])

    # NL = 1  als file
    elif (nl==1) & (ncols==9) & (arr1[-1].isdigit() is False):
        dtype = np.dtype([('ID',long),('X',float),('Y',float),('MAG',float),('ERR',float),('SKY',float),('ITER',float),('CHI',float),('SHARP',float)])
        cat = np.zeros(nstars,dtype=dtype)
        lengths = np.array([7,9,9,9,9,9,9,9,9])
        lo = np.concatenate((np.array([0]), np.cumsum(lengths[0:-1])))
        hi = lo+lengths
        names = cat.dtype.names
        for i in range(nstars):
            line1 = lines[i+3]
            for j in range(len(names)):
                cat[i][names[j]] = np.array(line1[lo[j]:hi[j]],dtype=dtype[names[j]])
    # NL = 2  aperture photometry
    elif nl==2:
        #
        #      1 1434.670   15.590   99.999   99.999   99.999   99.999   99.999
        #      1615.662 20.90  0.00  9.9999   9.9999   9.9999   9.9999   9.9999
        #
        #      2  233.850   18.420   99.999   99.999   99.999   99.999   99.999
        #      1613.601 20.96  0.02  9.9999   9.9999   9.9999   9.9999   9.9999
        #
        #  The columns are: ID, X, Y, Mag1, Mag2, etc..
        #                   Sky, St.Dev. of sky, skew of sky, Mag1err, Mag2err, etc.
        ncols = len(lines[4].split())
        naper = ncols-3   # apertures
        nstars = long((numlines(fil)-3.0)/3.0)  # stars
        dtype = np.dtype([('ID',long),('X',float),('Y',float),('SKY',float),('SKYSIG',float),('SKYSKEW',float),('MAG',float,naper),('ERR',float,naper)])
        cat = np.zeros(nstars,dtype=dtype)
        # for line 1
        lengths1 = np.concatenate([np.array([7,9,9]),np.zeros(naper,dtype=int)+9])
        dtype1 = np.concatenate([np.array([7,9,9]),np.zeros(naper,dtype=int)+9])
        lo1 = np.concatenate((np.array([0]), np.cumsum(lengths1[0:-1])))
        hi1 = lo1+lengths1
        names1 = ['ID','X','Y']+['MAG'+f for f in (np.arange(naper)+1).astype(str)]
        # for line 2
        lengths2 = np.concatenate([np.array([14,6,6]),np.zeros(naper,dtype=int)+9])
        lo2 = np.concatenate((np.array([0]), np.cumsum(lengths2[0:-1])))
        hi2 = lo2+lengths2
        names1 = ['SKY','SKYSIG','SKYSKEW']+['ERR'+f for f in (np.arange(naper)+1).astype(str)]
        for i in range(nstars):
            # line 1
            # ID, X, Y, Mag1, Mag2, etc.. 
            line1 = lines[i*3+4]
            cat[i]['ID'] = long(line1[lo1[0]:hi1[0]])
            cat[i]['X'] = float(line1[lo1[1]:hi1[1]])
            cat[i]['Y'] = float(line1[lo1[2]:hi1[2]])
            mag = np.zeros(naper,dtype=float)
            for j in range(naper):
                mag[j] = np.array(line1[lo1[j+3]:hi1[j+3]],dtype=float)
            cat[i]['MAG'] = mag
            # line 2
            # Sky, St.Dev. of sky, skew of sky, Mag1err, Mag2err, etc.  
            line2 = lines[i*3+5]
            cat[i]['SKY'] = float(line2[lo2[0]:hi2[0]])
            cat[i]['SKYSIG'] = float(line2[lo2[1]:hi2[1]])
            cat[i]['SKYSKEW'] = float(line2[lo2[2]:hi2[2]])
            err = np.zeros(naper,dtype=float)
            for j in range(naper):
                err[j] = np.array(line2[lo2[j+3]:hi2[j+3]],dtype=float)
            cat[i]['ERR'] = err
    # NL = 3  list
    elif nl==3:
        dtype = np.dtype([('ID',long),('X',float),('Y',float),('MAG',float),('ERR',float),('SKY',float)])
        cat = np.zeros(nstars,dtype=dtype)
        lengths = np.array([7,9,9,9,9,9,9])
        lo = np.concatenate((np.array([0]), np.cumsum(lengths[0:-1])))
        hi = lo+lengths
        names = cat.dtype.names
        for i in range(nstars):
            line1 = lines[i+3]
            for j in range(len(names)):
                cat[i][names[j]] = np.array(line1[lo[j]:hi[j]],dtype=dtype[names[j]])
    else:
        print("Cannot load this file")
        return None
    # Return as astropy Table
    return Table(cat)


# Make meta-data dictionary for an image:
def makemeta(fluxfile=None,header=None):
    '''
    This creates a meta-data dictionary for an exposure that is used by many
    of the photometry programs.  Either the filename or the header must be input.
    Note that sometimes in multi-extension FITS (MEF) files the information needed
    is both in the primary header and the extension header.  In that case it is best
    to combine them into one header and input that to makemeta().  This can easily
    be accomplished like this:
      
       head0 = fits.getheader("image1.fits",0)
       head = fits.getheader("image1.fits",1)
       head.extend(head0,unique=True)
       meta = makemeta(header=head)

    Parameters
    ----------
    fluxfile : str, optional
             The filename of the FITS image.
    header : str, optional
           The header of the image.

    Returns
    -------
    meta : astropy header
        The meta-data dictionary which is an astropy header with additional
        keyword/value pairs added.

    Example
    -------

    Create the meta-data dictionary for `image.fits`

    .. code-block:: python

        meta = makemeta("image.fits")

    Create the meta-data dictionary from `head`.

    .. code-block:: python

        meta = makemeta(header=head)

    '''

    # You generally need BOTH the PDU and extension header
    # To get all of this information

    if (fluxfile is None) & (header is None):
        print("No fluxfile or headerinput")
        return
    # Initialize meta using the header
    if fluxfile is not None:
        header = fits.getheader(fluxfile,0)
    meta = header

    #- INSTCODE -
    if "DTINSTRU" in meta.keys():
        if meta["DTINSTRU"] == 'mosaic3':
            meta["INSTCODE"] = 'k4m'
        elif meta["DTINSTRU"] == '90prime':
            meta["INSTCODE"] = 'ksb'
        elif meta["DTINSTRU"] == 'decam':
            meta["INSTCODE"] = 'c4d'
        else:
            print("Cannot determine INSTCODE type")
            return
    else:
        print("No DTINSTRU found in header.  Cannot determine instrument type")
        return

    #- RDNOISE -
    if "RDNOISE" not in meta.keys():
        # Check DECam style rdnoise
        if "RDNOISEA" in meta.keys():
            rdnoisea = meta["RDNOISEA"]
            rdnoiseb = meta["RDNOISEB"]
            rdnoise = (rdnoisea+rdnoiseb)*0.5
            meta["RDNOISE"] = rdnoise
        # Check other names
        if meta.get('RDNOISE') is None:
            for name in ['READNOIS','ENOISE']:
                if name in meta.keys(): meta['RDNOISE']=meta[name]
        # Bok
        if meta['INSTCODE'] == 'ksb':
            meta['RDNOISE']= [6.625, 7.4, 8.2, 7.1][meta['CCDNUM']-1]
        if meta.get('RDNOISE') is None:
            print('No RDNOISE found')
            return
    #- GAIN -
    if "GAIN" not in meta.keys():
        try:
            gainmap = { 'c4d': lambda x: 0.5*(x.get('GAINA')+x.get('GAINB')),
                        'k4m': lambda x: x.get('GAIN'),
                        'ksb': lambda x: [1.3,1.5,1.4,1.4][x.get['CCDNUM']-1] }  # bok gain in HDU0, use list here
            gain = gainmap[meta["INSTCODE"]](meta)
            meta["GAIN"] = gain
        except:
            gainmap_avg = { 'c4d': 3.9845419, 'k4m': 1.8575, 'ksb': 1.4}
            gain = gainmap_avg[meta["INSTCODE"]]
            meta["GAIN"] = gain
    #- CPFWHM -
    # FWHM values are ONLY in the extension headers
    cpfwhm_map = { 'c4d': 1.5 if meta.get('FWHM') is None else meta.get('FWHM')*0.27, 
                   'k4m': 1.5 if meta.get('SEEING1') is None else meta.get('SEEING1'),
                   'ksb': 1.5 if meta.get('SEEING1') is None else meta.get('SEEING1') }
    cpfwhm = cpfwhm_map[meta["INSTCODE"]]
    meta['CPFWHM'] = cpfwhm
    #- PIXSCALE -
    if "PIXSCALE" not in meta.keys():
        pixmap = { 'c4d': 0.27, 'k4m': 0.258, 'ksb': 0.45 }
        try:
            meta["PIXSCALE"] = pixmap[meta["INSTCODE"]]
        except:
            w = WCS(meta)
            meta["PIXSCALE"] = np.max(np.abs(w.pixel_scale_matrix))

    return meta


# Write SE catalog in DAO format
def sextodao(cat=None,meta=None,outfile=None,format="lst",naxis1=None,naxis2=None,saturate=None,rdnoise=None,gain=None,lowbad=None,thresh=None,logger=None):
    '''
    This writes out a Source Extractor catalog in a DAOPHOT format.

    Parameters
    ----------
    cat : numpy structured arrray or astropy Table format
        The Source Extractor catalog.
    meta : astropy header
         The image meta-data dictionary (naxis1, naxis2, saturate, rdnoise, gain, etc.).  The parameters
         can be input individually (see below).
    outfile : str
            The output filename.
    format : str, (lst, coo, ap, als)
           The output DAOPHOT format (lst, coo, ap, als).
    naxis1 : int, optional
           The X-dimensional size (in pixels) of the image.
    naxis2 : int, optional
           The Y-dimenaional size (in pixels) of the image.
    saturate : float, optional
           The saturation level of the image.
    rdnoise : float, optional
           The read noise of the image (in electrons).
    gain : float, optional
           The gain of the image (electrons/ADU).
    lowbad : float, optional
           The lower limit of the good range of values.
    thresh : float, optional
           The detection threshold.
    logger : logger object, optional
           The Logger to use for logging output.

    Returns
    -------
    Nothing is returned.  The catalog is written to `outfile`.

    Example
    -------

    .. code-block:: python

        sextodao(cat,meta,"cat.coo")

    '''

    if logger is None: logger = basiclogger('phot')   # set up basic logger if necessary
    # Not enough inputs
    if cat is None:
        logger.warning("No catalog input")
        return
    if meta is None:
        logger.warning("No image meta-data dictionary input")
        return
    if outfile is None:
        logger.warning("No outfile given")
        return
    # Delete outfile
    if os.path.exists(outfile): os.remove(outfile)

    # Get meta-data parameters, keyword inputs take priority over "meta"
    if naxis1 is None: naxis1=meta['NAXIS1']
    if naxis2 is None: naxis2=meta['NAXIS2']
    if saturate is None: saturate=meta['SATURATE']
    if rdnoise is None: rdnoise=meta['RDNOISE']
    if gain is None: gain=meta['GAIN']
    if lowbad is None:
        if meta.get('SKYMED') is not None:
            skymed = meta['SKYMED']
            skyrms = meta['SKYRMS']
            lowbad = skymed-7.*skyrms > 0.0
            thresh = skyrms*3.5
        else:
            logger.info("No sky value found in meta.  Using LOWBAD=1.0")
            lowbad = 1.0
    if thresh is None: thresh=20.0


    # Formats: coo, lst, ap, als

    # Header values:  this information comes from daophot2.pdf pg.69
    # NL: Originally meant "number of lines" but not anymore
    # NX: size of X-dimension of image in pixels
    # NY: size of Y-dimension of image in pixels
    # LOWBAD: lower good data limit, calculated by FIND
    # HIGHBAD: upper good data limit, specified in option file
    # THRESH: threshold calculated by FIND
    # AP1: radius (pixels) of the first aperture used by PHOTOMETRY
    # PH/ADU: gain in photons/ADU used when running FIND
    # RDNOISE: rdnoise (ADU) used when running FIND
    # FRAD: value of fitting radius

    # Go through the formats
    # "coo" file from FIND
    if format == "coo":

        #NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
        #  1  2046  4094  1472.8 38652.0   80.94    0.00    3.91    1.55    3.90
        # 
        #      1  1434.67    15.59   -0.045    0.313    0.873    1.218
        #      2   233.85    18.42   -0.018    0.218   -0.781    1.433
        #    ID      X         Y       MAG     SHARP    ROUND    ROUND2
        f = open(outfile,'w')
        # Header
        f.write(" NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n")
        f.write("  3 %5d %5d %7.1f %7.1f %7.2f %7.2f %7.2f %7.2f %7.2f\n" %
                (naxis1,naxis2,lowbad,saturate,thresh,3.0,gain,rdnoise/gain,3.9))
        f.write("\n")
        #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
        # Write the data
        for e in cat:
            f.write("%7d %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f\n" %
                    (e["NUMBER"],e["X_IMAGE"],e["Y_IMAGE"],e["MAG_AUTO"],0.6,0.0,0.0))
        f.close()

    # "lst" file from PICKPSF
    elif format == "lst":

        #NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
        # 3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90
        #
        #   318 1519.850  622.960   10.963    0.001    0.315
        #  1199 1036.580 2257.650   11.008    0.001    0.321
        #   ID     X        Y         MAG      ERR      SKY?
        f = open(outfile,'w')
        # Header
        f.write(" NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n")
        f.write("  3 %5d %5d %7.1f %7.1f %7.2f %7.2f %7.2f %7.2f %7.2f\n" %
                (naxis1,naxis2,lowbad,saturate,thresh,3.0,gain,rdnoise/gain,3.9))
        f.write("\n")
        #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
        # Write the data
        for e in cat:
            f.write("%7d %8.3f %8.3f %8.3f %8.3f %8.3f\n" %
                    (e["NUMBER"],e["X_IMAGE"]+1,e["Y_IMAGE"]+1,e["MAG_AUTO"],e["MAGERR_AUTO"],0.3))
        f.close()

    # "ap" file from PHOTOMETRY
    elif format == "ap":
        logger.warning(".ap files not supported yet")
        return

    # "als" file from ALLSTAR
    elif format == "als":

        # NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
        #  1  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90
        # 
        #      7  219.110   30.895   16.934   0.0935 1613.224       4.    0.872    0.040
        #     25 1396.437   62.936   12.588   0.0063 1615.938       4.    1.102   -0.042
        #    ID      X        Y       MAG      ERR     SKY        ITER     CHI     SHARP
        f = open(outfile,'w')
        # Header
        f.write(" NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n")
        f.write("  3 %5d %5d %7.1f %7.1f %7.2f %7.2f %7.2f %7.2f %7.2f\n" %
                (naxis1,naxis2,lowbad,saturate,thresh,3.0,gain,rdnoise/gain,3.9))
        f.write("\n")
        #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
        # Write the data
        for e in cat:
            f.write("%7d %8.3f %8.3f %8.3f %8.4f %8.3f %8.0f %8.3f %8.3f\n" %
                    (e["NUMBER"],e["X_IMAGE"]+1,e["Y_IMAGE"]+1,e["MAG_AUTO"],e["MAGERR_AUTO"],1500.0,1,1.0,0.0))
        f.close()

    # Not supported
    else:
        logger.warning(format+" NOT supported")
        return


# Run Source Extractor
#---------------------
def runsex(fluxfile=None,wtfile=None,maskfile=None,meta=None,outfile=None,configdir=None,logfile=None,logger=None):
    '''
    Run Source Extractor on an exposure.  The program is configured to work with files
    created by the NOAO Community Pipeline.

    Parameters
    ----------
    fluxfile : str
             The filename of the flux FITS image.
    wtfile : str
           The filename of the weight (1/variance) FITS image.
    maskfile : str
             The filename of the mask FITS image.
    meta : astropy header
         The meta-data dictionary for the exposure.
    outfile : str
            The output filename of the final catalog.
    configdir : str
              The directory that contains the Source Extractor configuration files.
              default.config, default.conv, default.nnw, default.param
    logfile : str, optional
            The name to use for the logfile.  If this is not input then the name will
            be the base name of `fluxfile` with the suffix "_sex.log".
    logger : logger object, optional
           The Logger to use for logging output.

    Returns
    -------
    cat : astropy Table
        The final Source Extractor catalog.
    maglim : float
        The magnitude limit of the exposure.

    The catalog is written to `outfile` and the output of Source Extractor to `logfile`.

    Example
    -------

    .. code-block:: python

        cat, maglim = runsex("flux.fits","wt.fits","mask.fits",meta,"cat.fits","/data/config/","sex.log")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running SExtractor --")

    # Not enough inputs
    if fluxfile is None:
        logger.warning("No fluxfile input")
        return
    if wtfile is None:
        logger.warning("No wtfile input")
        return
    if maskfile is None:
        logger.warning("No maskfile input")
        return
    if meta is None:
        logger.warning("No meta-data dictionary input")
        return
    if outfile is None:
        logger.warning("No outfile input")
        return
    if configdir is None:
        logger.warning("No configdir input")
        return

    # Check that necessary files exist
    for f in [fluxfile,wtfile,maskfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    base = os.path.basename(fluxfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if logfile is None: logfile=base+"_sex.log"

    # Working filenames
    sexbase = base+"_sex"
    sfluxfile = sexbase+".flux.fits"
    swtfile = sexbase+".wt.fits"
    smaskfile = sexbase+".mask.fits"

    if os.path.exists(outfile): os.remove(outfile)
    if os.path.exists(sfluxfile): os.remove(sfluxfile)
    if os.path.exists(swtfile): os.remove(swtfile)
    if os.path.exists(smaskfile): os.remove(smaskfile)
    if os.path.exists(logfile): os.remove(logfile)

    # Load the data
    flux,fhead = fits.getdata(fluxfile,header=True)
    wt,whead = fits.getdata(wtfile,header=True)
    mask,mhead = fits.getdata(maskfile,header=True)

    # 3a) Make subimages for flux, weight, mask

    # Turn the mask from integer to bitmask
    if ((meta["INSTCODE"]=='c4d') & (meta["plver"]>='V3.5.0')) | (meta["INSTCODE"]=='k4m') | (meta["INSTCODE"]=='ksb'):
         #  1 = bad (in static bad pixel mask) -> 1
         #  2 = no value (for stacks)          -> 2
         #  3 = saturated                      -> 4
         #  4 = bleed mask                     -> 8
         #  5 = cosmic ray                     -> 16
         #  6 = low weight                     -> 32
         #  7 = diff detect                    -> 64
         omask = mask.copy()
         mask *= 0
         nonzero = (omask>0)
         mask[nonzero] = 2**((omask-1)[nonzero])    # This takes about 1 sec
    # Fix the DECam Pre-V3.5.0 masks
    if (meta["INSTCODE"]=='c4d') & (meta["plver"]<'V3.5.0'):
      # --CP bit masks, Pre-V3.5.0 (PLVER)
      # Bit   DQ Type  PROCTYPE
      # 1  detector bad pixel          ->  1 
      # 2  saturated                   ->  4
      # 4  interpolated                ->  32
      # 16  single exposure cosmic ray ->  16
      # 64  bleed trail                ->  8
      # 128  multi-exposure transient  ->  0 TURN OFF
      # --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks
      #  1 = bad (in static bad pixel mask)
      #  2 = no value (for stacks)
      #  3 = saturated
      #  4 = bleed mask
      #  5 = cosmic ray
      #  6 = low weight
      #  7 = diff detect
      omask = mask.copy()
      mask *= 0     # re-initialize
      mask += (np.bitwise_and(omask,1)==1) * 1    # bad pixels
      mask += (np.bitwise_and(omask,2)==2) * 4    # saturated
      mask += (np.bitwise_and(omask,4)==4) * 32   # interpolated
      mask += (np.bitwise_and(omask,16)==16) * 16  # cosmic ray
      mask += (np.bitwise_and(omask,64)==64) * 8   # bleed trail

    # Mask out bad pixels in WEIGHT image
    #  set wt=0 for mask>0 pixels
    wt[ (mask>0) | (wt<0) ] = 0   # CP sets bad pixels to wt=0 or sometimes negative

    # Write out the files
    shutil.copy(fluxfile,sfluxfile)
    fits.writeto(swtfile,wt,header=whead,output_verify='warn')


    # 3b) Make SExtractor config files
    # Copy the default files
    shutil.copyfile(configdir+"default.conv","default.conv")
    shutil.copyfile(configdir+"default.nnw","default.nnw")
    shutil.copyfile(configdir+"default.param","default.param")

    # Read in configuration file and modify for this image
    lines = readlines(configdir+'default.config')

    # Gain, saturation, pixscale

    # Things to change
    # SATUR_LEVEL     59000.00         # level (in ADUs) at which arises saturation
    # GAIN            43.52             # detector gain in e-/ADU.
    # SEEING_FWHM     1.46920            # stellar FWHM in arcsec
    # WEIGHT_IMAGE  F4-00507860_01_comb.mask.fits

    filter_name = ''
    cnt = 0
    for l in lines:
        # CATALOG_NAME
        m = re.search('^CATALOG_NAME',l)
        if m != None:
            lines[cnt] = "CATALOG_NAME     "+outfile+"         # name of the output catalog\n"
        # FLAG_IMAGE
        m = re.search('^FLAG_IMAGE',l)
        if m != None:
            lines[cnt] = "FLAG_IMAGE     "+smaskfile+"         # filename for an input FLAG-image\n"
        # WEIGHT_IMAGE
        m = re.search('^WEIGHT_IMAGE',l)
        if m != None:
            lines[cnt] = "WEIGHT_IMAGE     "+swtfile+"  # Weight image name\n"
        # SATUR_LEVEL
        m = re.search('^SATUR_LEVEL',l)
        if m != None:
            lines[cnt] = "SATUR_LEVEL     "+str(meta["saturate"])+"         # level (in ADUs) at which arises saturation\n"
        # Gain
        m = re.search('^GAIN',l)
        if m != None:
            lines[cnt] = "GAIN            "+str(meta["gain"])+"            # detector gain in e-/ADU.\n"
        # SEEING_FWHM
        m = re.search('^SEEING_FWHM',l)
        if m != None:
            lines[cnt] = "SEEING_FWHM     "+str(meta["cpfwhm"])+"            # stellar FWHM in arcsec\n"
        # PHOT_APERTURES, aperture diameters in pixels
        m = re.search('^PHOT_APERTURES',l)
        if m != None:
            aper_world = np.array([ 0.5, 1.0, 2.0, 3.0, 4.0]) * 2  # radius->diameter, 1, 2, 4, 6, 8"
            aper_pix = aper_world / meta["pixscale"]
            lines[cnt] = "PHOT_APERTURES  "+', '.join(np.array(np.round(aper_pix,2),dtype='str'))+"            # MAG_APER aperture diameter(s) in pixels\n"            
        # Filter name
        m = re.search('^FILTER_NAME',l)
        if m != None:
            filter_name = (l.split())[1]
        cnt = cnt+1
    # Write out the new config file
    if os.path.exists("default.config"):
        os.remove("default.config")
    fo = open('default.config', 'w')
    fo.writelines(lines)
    fo.close()

    # Convolve the mask file with the convolution kernel to "grow" the regions
    # around bad pixels the SE already does to the weight map
    if (filter_name != ''):
        # Load the filter array
        f = open(filter_name,'r')
        linenum = 0
        for line in f:
            if (linenum == 1):
                shape = line.split(' ')[1]
                # Make it two pixels larger
                filter = np.ones(np.array(shape.split('x'),dtype='i')+2,dtype='i')
                #filter = np.zeros(np.array(shape.split('x'),dtype='i'),dtype='f')
            #if (linenum > 1):
            #    linedata = np.array(line.split(' '),dtype='f')
            #    filter[:,linenum-2] = linedata
            linenum += 1
        f.close()
        # Normalize the filter array
        #filter /= np.sum(filter)
        # Convolve with mask
        #filter = np.ones(np.array(shape.split('x'),dtype='i'),dtype='i')
        #mask2 = convolve2d(mask,filter,mode="same",boundary="symm")
        mask2 = convolve(mask,filter,mode="reflect")
        bad = ((mask == 0) & (mask2 > 0))
        newmask = np.copy(mask)
        newmask[bad] = 1     # mask out the neighboring pixels
        # Write new mask
        fits.writeto(smaskfile,newmask,header=mhead,output_verify='warn')

    # 3c) Run SExtractor
    try:
        # Save the SExtractor info to a logfile
        sf = open(logfile,'w')
        retcode = subprocess.call(["sex",sfluxfile,"-c","default.config"],stdout=sf,stderr=subprocess.STDOUT)
        sf.close()
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("SExtractor Execution failed:"+str(e))
        logger.error(e)
        raise

    # Check that the output file exists
    if os.path.exists(outfile) is True:
        # Load the catalog and keep it in memory for later use
        cat = Table.read(outfile,2)
        # How many sources were detected, final catalog file
        logger.info(str(len(cat))+" sources detected")
        logger.info("Final catalog is "+outfile)
        # Get the magnitude limit, use 90th percentile
        gdcat = (cat["MAG_AUTO"]<50)
        ngdcat = np.sum(gdcat)
        mag = cat["MAG_AUTO"][gdcat]
        mag_sorted = np.sort(mag)
        maglim = mag_sorted[int(np.round(0.90*ngdcat))]
        logger.info("Estimated magnitude limit = %6.2f mag" % maglim)
        # Get background value and RMS and add to meta
        plines = readlines(logfile)
        plines2 = grep(plines,'Background')
        arr = plines2[0].split()
        ind1 = grep(arr,'Background:',index=True)
        ind2 = grep(arr,'RMS',index=True)
        ind3 = grep(arr,'Threshold',index=True)
        background = np.float(arr[ind1[0]+1])
        rms = np.float(arr[ind2[0]+1])
        meta["SKYMED"] = (background,"Median sky background")
        meta["SKYRMS"] = (rms,"RMS of sky")

    # Delete temporary files
    os.remove(sfluxfile)
    os.remove(smaskfile)
    os.remove(swtfile)
    #os.remove("default.conv")

    return cat,maglim


# Determine seeing FWHM using SE catalog
#---------------------------------------
def sexfwhm(cat=None,logger=None):
    '''
    Determine the seeing FWHM using a Source Extractor catalog.

    Parameters
    ----------
    cat : astropy Table
        The Source Extractor catalog.

    Returns
    -------
    fwhm : float
         The seeing FWHM in arcsec.

    Example
    -------

    .. code-block:: python

        fwhm = sexfwhm(cat)

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    # Make sure we have the SE catalog
    if cat is None:
        logger.warning("No catalog input")
        return
    # Select good sources
    gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.05) & (cat['CLASS_STAR']>0.8) &
             (cat['FLAGS']==0) & (cat['IMAFLAGS_ISO']==0))
    ngdcat = np.sum(gdcat)
    # Not enough good source, remove FLAGS cut
    if (ngdcat<10):
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.05) & (cat['CLASS_STAR']>0.8) &
                 (cat['IMAFLAGS_ISO']==0))
        ngdcat = np.sum(gdcat)
    # Not enough good source, remove FLAGS/CLASS_STAR cuts
    if (ngdcat<10):
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.05) & (cat['IMAFLAGS_ISO']==0))
        ngdcat = np.sum(gdcat)
    # Not enough sources, lower thresholds
    if (ngdcat<10):
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.08))
        ngdcat = np.sum(gdcat)            
    medfwhm = np.median(cat[gdcat]['FWHM_WORLD']*3600.)
    logger.info('FWHM = %5.2f arcsec (%d sources)' % (medfwhm, ngdcat))

    return medfwhm


# Pick PSF candidates using SE catalog
#-------------------------------------
def sexpickpsf(cat=None,fwhm=None,meta=None,outfile=None,nstars=100,logger=None):
    '''
    Pick PSF stars using a Source Extractor catalog and output to a DAOPHOT-style file.

    Parameters
    ----------
    cat : astropy Table
        The Source Extractor catalog.
    fwhm : float
         The seeing FWHM of the exposure (in arcsec).
    meta : astropy dictionary
         The metal-data dictionary for the image.    
    outfile : str
           The filaname of the DAOPHOT-style lst file to write the PSF stars to.
    nstars : int, optional, default is 100
           The number of PSF stars to pick.
    logger : logging object
          The logger to use for logging information.

    Returns
    -------
    psfcat : astropy Table
         The table of PSF stars.

    The table of PSF stars is also written to `outfile`.

    Example
    -------

    .. code-block:: python

        psfcat = sexpickpsf(cat,fwhm,meta,"psfstars.lst",nstars=100)

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary

    # Make sure we have the SE catalog
    if cat is None:
        logger.warning("No catalog input")
        return
    # Make sure we have FWHM
    if fwhm is None:
        logger.warning("No FWHM input")
        return
    # Make sure we have meta
    if meta is None:
        logger.warning("No meta-data dictionary input")
        return
    # Make sure we have the output file
    if outfile is None:
        logger.warning("No outfile input")
        return

    # Select good sources
    gdcat1 = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.05) & (cat['CLASS_STAR']>0.8))
    ngdcat1 = np.sum(gdcat1)
    # Bright and faint limit, use 5th and 95th percentile
    minmag, maxmag = np.sort(cat[gdcat1]['MAG_AUTO'])[[int(np.round(0.05*ngdcat1)),int(np.round(0.95*ngdcat1))]]
    # Select stars with
    # -good FWHM values
    # -good clas_star values (unless FWHM too large)
    # -good mag range, bright but not too bright
    # -no flags set
    if fwhm<1.8:
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.1) & (cat['CLASS_STAR']>0.8) & 
                 (cat['FWHM_WORLD']*3600.>0.5*fwhm) & (cat['FWHM_WORLD']*3600.<1.5*fwhm) &
                 (cat['MAG_AUTO']>(minmag+1.0)) & (cat['MAG_AUTO']<(maxmag-0.5)) &
                 (cat['FLAGS']==0) & (cat['IMAFLAGS_ISO']==0))
        ngdcat = np.sum(gdcat)
    # Do not use CLASS_STAR if seeing bad, not as reliable
    else:
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.1) & 
                 (cat['FWHM_WORLD']*3600.>0.5*fwhm) & (cat['FWHM_WORLD']*3600.<1.5*fwhm) &
                 (cat['MAG_AUTO']>(minmag+1.0)) & (cat['MAG_AUTO']<(maxmag-0.5)) &
                 (cat['FLAGS']==0) & (cat['IMAFLAGS_ISO']==0))
        ngdcat = np.sum(gdcat)
    # No candidate, loosen cuts
    if ngdcat<10:
        logger.info("Too few PSF stars on first try. Loosening cuts")
        gdcat = ((cat['MAG_AUTO']< 50) & (cat['MAGERR_AUTO']<0.15) & 
                 (cat['FWHM_WORLD']*3600.>0.2*self.seeing) & (cat['FWHM_WORLD']*3600.<1.8*fwhm) &
                 (cat['MAG_AUTO']>(minmag+0.5)) & (cat['MAG_AUTO']<(maxmag-0.5)))
        ngdcat = np.sum(gdcat)
    # No candidates
    if ngdcat==0:
        logger.error("No good PSF stars found")
        raise

    # Candidate PSF stars, use only Nstars, and sort by magnitude
    si = np.argsort(cat[gdcat]['MAG_AUTO'])
    psfcat = cat[gdcat][si]
    if ngdcat>nstars: psfcat=psfcat[0:nstars]
    logger.info(str(len(psfcat))+" PSF stars found")

    # Output them in DAO format
    sextodao(psfcat,meta,outfile,format="lst")
    if os.path.exists(outfile) is False:
        logger.error("Output file "+outfile+" NOT found")
        raise

    return psfcat


    # Do we a need separate aperture photometry file?
    


# Make DAOPHOT option files
#--------------------------
def mkopt(base=None,meta=None,VA=1,LO=7.0,TH=3.5,LS=0.2,HS=1.0,LR=-1.0,HR=1.0,
          WA=-2,AN=-6,EX=5,PE=0.75,PR=5.0,CR=2.5,CE=6.0,MA=50.0,RED=1.0,WA2=0.0,
          fitradius_fwhm=1.0,HI=None,RD=None,GA=None,FW=None,logger=None):
    '''
    Create the DAOPHOT and ALLSTAR option files (.opt and .als.opt) for an exposure.

    Parameters
    ----------
    base : str
         The base name to use for the option files.  The DAOPHOT option file will
         be called `base`.opt and the ALLSTAR option file `base`.als.opt
    meta : astropy dictionary
         The metal-data dictionary for the image.    
    VA : int, default = 1
       The variable type of PSF to use.
       -1: Analytic PSF only
        0: Analytic PSF and look-up table of empirical corrections
        1: linear variations across the field
        2: quadratic variations across the field
    LO : float, default = 7.0
       Low good datum (7. works fine on most imags).
    TH : float, default = 3.5
       Threshold in sigma above the background (3.5 works fine).
    LS : float, default = 0.2
       Lower sharpness cutoff.
    HS : float, default = 1.0
       High sharpness cutoff.
    LR : float, default = -1.0
       Lower roundness cutoff.
    HR : float, default = 1.0
       High roundness cutoff.
    WA : int, default = -2
       Watch progress for DAOPHOT.  Determines what output is displayed.
    AN : int, default = -6
       Analytic model PSF.
        1: Gaussian (3 pararameters)
        2: Moffat function (3 parameters), beta=1.5
        3: Moffat function (3 parameters), beta=2.5
        4: Lorentz function (3 parameters)
        5: Penny function, Gauss+Lorentz (4 parameters), G+L are parallel
        6: Penny function (5 parameters), G and L can be in different directions
        A negative sign in front means to try all functions up to X and pick the best one.
    EX : int, default = 5
       Extra PSF cleaning passes.
    PE : float, default = 0.75
       Percent error due to the uncertainty in the fine-scale structure of the flat field.
    PR : float, default = 5.0
       Profile error due to the incompleteness of the PSF model.
    CR : float, default = 2.5
       Clipping range.  Used to remove outlier pixels. Parameter "a" in the formula given in
       Stetson 1987, PASP, 99, 191, section III.D.2.d "Resisting bad data".
    CE : float, default = 6.0
       Clipping exponent.  Parameter b in above clipping formula.
    MA : float, default = 50.0
       Maximum group size
    RED : float, default = 1.0
        Redetermine centroid (0 = no, 1 = yes).
    WA2 : float, default = 0.0
        Watch progress for ALLSTAR.      
    fitradius_fwhm : float, default = 1.0
        The fitting radius size in units of the seeing FWHM for the area to be fit.
    HI : float, optional
       High good datum.  Normally set by `saturate` from `meta`.
    RD : float, optional
       The read noise in electrons. Normally set by `rdnoise` from `meta`.
    GA : float, optional
       The gain in electrons/ADU. Normally set by `gain` from `meta`.
    FW : float, optional
       The seeing FWHM in pixels.  Normally set by `fwhm`/`pixscale` from `meta`.
    logger : logger object, optional
           The Logger to use for logging output.

    Returns
    -------
    Nothing is returned.  The DAOPHOT option file is written to `base`.opt and the ALLSTAR
    option file to `base`.als.opt.

    Example
    -------

    .. code-block:: python

        mkopt("image",meta)

    '''

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % MAKING THE OPT FILES
    #
    # (1) DAOPHOT parameters
    # 
    # LO    : Low good datum (7. works fine on most imags)
    # TH    : Threshold (3.5 works fine)
    # LS,HS : Low and high sharpness (default : 0.2 - 1.0)
    # LR,HR : Low roundness and high roundness (default : -1.0 - 1.0)
    # WA    : Watch progress
    # VA    : Variable PSF
    # AN    : Analytic model PSF
    # EX    : Extra PSF cleaning passes
    # PE    : Percent error
    # PR    : Profile error
    #
    # (2) ALLSTAR parameters
    # 
    # CR    : Clipping range (leave it)
    # CE    : Clipping exponent (leave it)
    # MA    : Maximum group size
    # RED   : Redetermine centroid (0 = no, 1 = yes)
    #
    # Frame-specific parameters.
    #
    # GA    : gain (e/ADU)
    # RD    : readout noise (e)
    # RE    : readout noise (ADU)
    # FW    : FWHM
    # HI    : hi good datum in ADU - saturation level
    # FI    : fitting radius
    # PS    : PSF radius
    # IS,OS : inner and outer sky annalus
    # VA  defined above
    #AN = -6     # It will try all PSF models (#1-6) and use the one with the lowest chi value
    #EX =  5     # extra PSF passes

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary

    optfile = base+".opt"
    alsoptfile = base+".als.opt"

    # Get frame specific parameters from meta if necessary
    if GA is None: GA = meta['gain']
    if RD is None: RD = meta['rdnoise']
    if FW is None: FW = meta['fwhm'] / meta['pixscale']
    if HI is None: HI = meta['saturate']


    # Calculating some things
    FW = np.min([ FW , 20 ])            # daophot won't accept anything higher than 20
    RE = RD/GA
    FI = np.min([ fitradius_fwhm*FW , 51 ])                  # daophot won't accept anything higher than 51
    PS = np.min([ (4.0*FW) , 51 ])       # daophot won't accept anything higher than 51
    IS = np.min([ (FI - 1.0) , 35 ])     # daophot won't accept anything higher than 35
    OS = np.min([ (PS + 1.0) , 100 ])    # daophot won't accept anything higher than 100

    # Writing the DAOPHOT parameter
    #------------------------------
    #
    # RE    : readout noise (ADU)
    # GA    : gain (e/ADU)
    # LO    : Low good datum (7. works fine on most imags)
    # HI    : hi good datum in ADU - saturation level
    # FW    : FWHM
    # TH    : Threshold (3.5 works fine)
    # LS,HS : Low and high sharpness (default : 0.2 - 1.0)
    # LR,HR : Low roundness and high roundness (default : -1.0 - 1.0)
    # WA    : Watch progress
    # FI    : fitting radius
    # PS    : PSF radius
    # VA    : Variable PSF
    # AN    : Analytic model PSF
    # EX    : Extra PSF cleaning passes
    # PE    : Percent error
    # PR    : Profile error

    outarr = [RE,GA,LO,HI,FW,TH,LS,HS,LR,HR,WA,FI,PS,VA,AN,EX,PE,PR]
    anotarr = ['RE','GA','LO','HI','FW','TH','LS','HS','LR','HR','WA','FI','PS','VA','AN','EX','PE','PR']
    nanot = len(anotarr)

    # Delete file if it exists
    if os.path.exists(optfile):
        os.remove(optfile)
    # Write opt file
    f = open(optfile,'w')
    for j in range(len(outarr)):
        if anotarr[j] == "HI":
            f.write("%2s = %8d\n" % (anotarr[j], outarr[j]))
        else:
            f.write("%2s = %8.2f\n" % (anotarr[j], outarr[j]))
    f.close()

    # Writing the ALLSTAR parameter file
    #-----------------------------------
    #
    # FI    : fitting radius
    # IS    :  ??
    # OS    :  ??
    # RED   : Redetermine centroid (0 = no, 1 = yes)
    # WA2   : Watch progress
    # PE    : Percent error
    # PR    : Profile error
    # CR    : Clipping range (leave it)
    # CE    : Clipping exponent (leave it)
    # MA    : Maximum group size

    outarr2 = [FI,IS,OS,RED,WA2,PE,PR,CR,CE,MA]
    anotarr2 = ['FI','IS','OS','RE','WA','PE','PR','CR','CE','MA']
    nanot2 = len(anotarr2)
    form = '(A5,F8.2)'

    # Delete file if it exists
    if os.path.exists(alsoptfile):
        os.remove(alsoptfile)
    # Write opt file
    f = open(alsoptfile,'w')
    for j in range(len(outarr2)):
        f.write("%2s = %8.2f\n" % (anotarr2[j], outarr2[j]))
    f.close()

    logger.info("Created "+optfile+" and "+alsoptfile)


# Make image ready for DAOPHOT
def mkdaoim(fluxfile=None,wtfile=None,maskfile=None,meta=None,outfile=None,logger=None):
    '''
    This constructs a FITS image that is prepared for DAOPHOT.
    This program was designed for exposures from the NOAO Community Pipeline.

    Parameters
    ----------
    fluxfile : str
             The filename of the flux FITS image.
    wtfile : str
           The filename of the weight (1/variance) FITS image.
    maskfile : str
             The filename of the mask FITS image.
    meta : astropy header
         The meta-data dictionary for the exposure.
    outfile : str
            The name of the output FITS file.

    Returns
    -------
    Nothing is returned.  The DAOPHOT-ready image is written to `outfile`.

    Example
    -------

    .. code-block:: python

        mkdaoim("flux.fits","wt.fits","mask.fits","image.fits")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary

    # Not enough inputs
    if fluxfile is None:
        logger.warning("No fluxfile input")
        return
    if wtfile is None:
        logger.warning("No wtfile input")
        return
    if maskfile is None:
        logger.warning("No maskfile input")
        return
    if meta is None:
        logger.warning("No meta-data dictionary input")
        return
    if outfile is None:
        logger.warning("No outfile input")
        return

    # Check that necessary files exist
    for f in [fluxfile,wtfile,maskfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    # Load the FITS files
    flux,fhead = fits.getdata(fluxfile,header=True)
    wt,whead = fits.getdata(wtfile,header=True)
    mask,mhead = fits.getdata(maskfile,header=True)

    # Set bad pixels to saturation value
    # --DESDM bit masks (from Gruendl):
    # BADPIX_BPM 1          /* set in bpm (hot/dead pixel/column)        */
    # BADPIX_SATURATE 2     /* saturated pixel                           */
    # BADPIX_INTERP 4
    #     /* interpolated pixel                        */
    # BADPIX_LOW     8      /* too little signal- i.e. poor read         */
    # BADPIX_CRAY   16      /* cosmic ray pixel                          */
    # BADPIX_STAR   32      /* bright star pixel                         */
    # BADPIX_TRAIL  64      /* bleed trail pixel                         */
    # BADPIX_EDGEBLEED 128  /* edge bleed pixel                          */
    # BADPIX_SSXTALK 256    /* pixel potentially effected by xtalk from super-saturated source */
    # BADPIX_EDGE   512     /* pixel flagged to exclude CCD glowing edges */
    # BADPIX_STREAK 1024    /* pixel associated with satellite (airplane/meteor) streak     */
    # BADPIX_FIX    2048    /* a bad pixel that was fixed                */
    # --CP bit masks, Pre-V3.5.0 (PLVER)
    # Bit   DQ Type  PROCTYPE
    # 1  detector bad pixel          InstCal
    # 1  detector bad pixel/no data  Resampled
    # 1  No data                     Stacked
    # 2  saturated                   InstCal/Resampled
    # 4  interpolated                InstCal/Resampled
    # 16  single exposure cosmic ray InstCal/Resampled
    # 64  bleed trail                InstCal/Resampled
    # 128  multi-exposure transient  InstCal/Resampled
    # --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks
    #  1 = bad (in static bad pixel mask)
    #  2 = no value (for stacks)
    #  3 = saturated
    #  4 = bleed mask
    #  5 = cosmic ray
    #  6 = low weight
    #  7 = diff detect
    # You can't have combinations but the precedence as in the order
    # of the list (which is also the order in which the processing
    # discovers them).  So a pixel marked as "bad" (1) won't ever be
    # flagged as "diff detect" (7) later on in the processing.
    #
    # "Turn off" the "difference image masking", clear the 8th bit
    # 128 for Pre-V3.5.0 images and set 7 values to zero for V3.5.0 or later.

    #logger.info("Turning off the CP difference image masking flags")
    if meta["plver"] > 0:      # CP data
        # V3.5.0 and on, Integer masks
        versnum = meta["plver"].split('.')
        if (versnum[0]>3) | ((versnum[0]==3) & (versnum[1]>=5)):
            bdpix = (mask == 7)
            nbdpix = np.sum(bdpix)
            if nbdpix > 0: mask[bdpix]=0

        # Pre-V3.5.0, Bitmasks
        else:
            bdpix = ( (mask & 2**7) == 2**7)
            nbdpix = np.sum(bdpix)                
            if nbdpix > 0: mask[bdpix]-=128   # clear 128

        logger.info("%d pixels cleared of difference image mask flag" % nbdpix)

    bdpix = (mask > 0.0)
    nbdpix = np.sum(bdpix)
    if nbdpix>0: flux[bdpix]=6e4
    logger.info("%d bad pixels masked" % nbdpix)

    fhead.append('GAIN',meta["GAIN"])
    fhead.append('RDNOISE',meta["RDNOISE"])

    # DAOPHOT can only handle BITPIX=16, 32, -32
    if fhead['BITPIX'] not in [16,32,-32]:
        logger.info("BITPIX="+str(fhead['BITPIX'])+" DAOPHOT can only handle 16,32,-32.  Changing to -32")
        flux = np.array(flux,dtype=np.float32)
        fhead['BITPIX'] = -32

    # Write new image
    logger.info("Wrote DAOPHOT-ready image to "+outfile)
    fits.writeto(outfile,flux,fhead,overwrite=True)


# DAOPHOT FIND detection
#-----------------------
def daofind(imfile=None,optfile=None,outfile=None,logfile=None,logger=None):
    '''
    This runs DAOPHOT FIND on an image.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    outfile : str, optional
            The output filename of the FIND catalog.  By default this is
            the base name of `imfile` with a ".coo" suffix.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".coo.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    cat : astropy table
        The DAOPHOT FIND catalog.

    The output catalog and logfile will also be created.

    Example
    -------

    .. code-block:: python

        cat = daofind("image.fits")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running DAOPHOT detection --")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return None

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if outfile is None: outfile = base+".coo"
    if logfile is None: logfile = base+".coo.log"
    scriptfile = base+".coo.sh"
    for f in [outfile,logfile,scriptfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,optfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tcoo",dir=".")
    os.close(tid)   # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    toptfile = tbase+".opt"
    toutfile = tbase+".coo"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daophot << END_DAOPHOT >> "+logfile+"\n" \
            "OPTIONS\n" \
            ""+toptfile+"\n" \
            "\n" \
            "ATTACH "+timfile+"\n" \
            "FIND\n" \
            "1,1\n" \
            ""+toutfile+"\n" \
            "y\n" \
            "EXIT\n" \
            "EXIT\n" \
            "END_DAOPHOT\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=False)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("DAOPHOT detection failed:"+str(e))
        logger.error(e)
        raise Exception("DAOPHOT failed")

    # Check that the output file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile]: os.remove(f)

        # Get info from the logfile
        if os.path.exists(logfile) is True:
            dlines = readlines(logfile)
            l1 = grep(dlines,"Sky mode and standard deviation")
            if len(l1)>0:
                logger.info(l1[0].strip())   # clip \n
                #l1 = l1[0]
                #lo = l1.find("=")
                #sky = np.array( l1[lo+1:].split('  '),dtype=float)
            l2 = grep(dlines,"Clipped mean and median")
            if len(l2)>0:
                logger.info(l2[0].strip())
                #l2 = l2[0]
                #lo = l2.find("=")
                #mnmed = np.array( l2[lo+2:].split(' '),dtype=float)
            # Number of sources
            l3 = grep(dlines," stars.")
            if len(l3)>0:
                logger.info(l3[0].rstrip().strip())
    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("Output not found")

    # Delete the script
    os.remove(scriptfile)

    # Load and return the catalog
    logger.info("Output file = "+outfile)
    return daoread(outfile)


# DAOPHOT aperture photometry
#----------------------------
def daoaperphot(imfile=None,coofile=None,apertures=None,outfile=None,optfile=None,apersfile=None,logfile=None,logger=None):
    '''
    This runs DAOPHOT aperture photometry on an image.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    coofile : str, optional
            The filename of the catalog of sources for which to obtain aperture photometry.
            By default it is assumed that this is the base name of `imfile` with a ".coo" suffix.
    apertures : list or array, optional
             The list of aperture to use.  The last two are used as the inner and outer sky radius.
             The default apertures are: apertures = [3.0, 6.0803, 9.7377, 15.5952, 19.7360, 40.0, 50.0]
    outfile : str, optional
            The output filename of the aperture photometry catalog.  By default this is
            the base name of `imfile` with a ".ap" suffix.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    apersfile : str, optional
              The file that will constrain the apertures used.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".ap.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    cat : astropy table
        The DAOPHOT aperture photometry catalog.
    maglim : float
        The magnitude limit of the exposure.

    The output catalog and logfile will also be created.

    Example
    -------

    .. code-block:: python

        cat, maglim = daoaperphot("image.fits","image.coo")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running DAOPHOT aperture photometry --")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return None

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if coofile is None: coofile = base+".coo"
    if outfile is None: outfile = base+".ap"
    if logfile is None: logfile = base+".ap.log"
    if apersfile is None: apersfile = base+".apers"
    scriptfile = base+".coo.sh"
    for f in [outfile,apersfile,logfile,scriptfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,optfile,coofile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tap",dir=".")
    os.close(tid)   # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    cooext = os.path.splitext(coofile)[1]
    tcoofile = tbase+cooext
    toptfile = tbase+".opt"
    tapersfile = tbase+".apers"
    toutfile = tbase+".ap"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)
    os.symlink(coofile,tcoofile)

    logger.info("coofile = "+coofile)

    # Make apertures file
    if apertures is None:
        # The last two are inner and outer sky apertures
        #apertures = [3.0, 3.7965, 4.8046, 6.0803, 7.6947, 9.7377, 12.3232, 15.5952, 19.7360, \
        #             24.9762, 31.6077, 40.0000, 50.0000]
        apertures = [3.000, 6.0803, 9.7377, 15.5952, 19.7360, 40.0000, 50.0000]
    aperswrite(tapersfile,apertures)
    #nap = len(apertures)
    #if nap<3:
    #    logger.warning("Only "+str(nap)+" apertures input.  Need at least 3")
    #    return None
    #f = open(tapersfile,'w')
    #for i in range(nap-2):
    #    # use hexidecimal for aperture id, 2 digits, first starts with A
    #    id = hex(160+i+1)
    #    id = id[2:].capitalize()
    #    f.write("%2s = %7.4f\n" % (id,apertures[i]))
    #f.write("IS = %7.4f\n" % apertures[nap-2])
    #f.write("OS = %7.4f\n" % apertures[nap-1])
    #f.close()

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daophot << END_DAOPHOT >> "+logfile+"\n" \
            "OPTIONS\n" \
            ""+toptfile+"\n" \
            "\n" \
            "ATTACH "+timfile+"\n" \
            "PHOTOMETRY\n" \
            ""+tapersfile+"\n" \
            " \n" \
            ""+tcoofile+"\n" \
            ""+toutfile+"\n" \
            "EXIT\n" \
            "EXIT\n" \
            "END_DAOPHOT\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

    # If PSF file exists temporarily move it out of the way
    if os.path.exists(base+".psf"):
        logger.info(base+".psf exists.  Temporarily moving it out of the way to perform aperture photometry.")
        psftemp = base+".psf.bak"
        if os.path.exists(psftemp): os.remove(psftemp)
        os.rename(base+".psf",psftemp)
        movedpsf = True
    else:
        movedpsf = False

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("DAOPHOT aperture photometry failed:"+str(e))
        logger.error(e)
        raise Exception("DAOPHOT failed")

    # Check that the output file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        if apersfile is not None: shutil.copyfile(tapersfile,apersfile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile,tcoofile,tapersfile]: os.remove(f)

        # Get info from the logfile
        maglim = None
        if os.path.exists(logfile):
            plines = readlines(logfile)
            l1 = grep(plines,"Estimated magnitude limit")
            if len(l1)>0:
                l1 = l1[0]
                l1 = l1[0:len(l1)-7]   # strip BELL at end \x07\n
                lo = l1.find(":")
                hi = l1.find("+-")
                maglim = np.float(l1[lo+1:hi])
                logger.info(l1.strip())   # clip leading/trailing whitespace
    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("Output not found")

    # Delete the script
    os.remove(scriptfile)

    # Move PSF file back
    if movedpsf is True: os.rename(psftemp,base+".psf")

    # Return the catalog
    logger.info("Output file = "+outfile)
    return daoread(outfile), maglim


# Pick PSF stars using DAOPHOT
#-----------------------------
def daopickpsf(imfile=None,catfile=None,maglim=None,outfile=None,nstars=100,optfile=None,logfile=None,logger=None):
    '''
    This runs DAOPHOT aperture photometry on an image.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    catfile : str
           The catalog file from which to pick PSF stars.  This is normally the .ap file.
    maglim : float
           The magnitude limit for this image.  
    nstars : int, optional, default = 100
           The number of PSF stars to pick.  
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    outfile : str, optional
            The output filename of the aperture photometry catalog.  By default this is
            the base name of `imfile` with a ".lst" suffix.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".lst.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    cat : astropy table
        The list of PSF stars.

    The output catalog and logfile will also be created.

    Example
    -------

    .. code-block:: python

        psfcat = daopickpsf("image.fits","image.coo",19.5,nstars=100)

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running DAOPHOT PICKPSF -- ")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return None
    # Make sure we have the catalog file name
    if catfile is None:
        logger.warning("No catalog filename input")
        return None

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if outfile is None: outfile = base+".lst"
    if logfile is None: logfile = base+".lst.log"
    scriptfile = base+".pickpsf.sh"
    for f in [outfile,logfile,scriptfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,catfile,optfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tlst",dir=".")
    os.close(tid)  # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    toptfile = tbase+".opt"
    catext = os.path.splitext(catfile)[1]
    tcatfile = tbase+catext
    toutfile = tbase+".lst"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)
    os.symlink(catfile,tcatfile)

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daophot << END_DAOPHOT >> "+logfile+"\n" \
            "OPTIONS\n" \
            ""+toptfile+"\n" \
            "\n" \
            "ATTACH "+timfile+"\n" \
            "PICKPSF\n" \
            ""+tcatfile+"\n" \
            ""+str(nstars)+","+str(maglim-1.0)+"\n" \
            ""+toutfile+"\n" \
            "EXIT\n" \
            "EXIT\n" \
            "END_DAOPHOT\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("DAOPHOT PICKPSF failed:"+str(e))
        logger.error(e)
        raise Exception("DAOPHOT failed")

    # Check that the output file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile,tcatfile]: os.remove(f)

        # Get info from the logfile
        if os.path.exists(logfile):
            plines = readlines(logfile)
            l1 = grep(plines,"suitable candidates were found.")
            if len(l1)>0:
                logger.info(l1[0].strip())   # clip \n
    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("DAOPHOT failed")

    # Delete the script
    os.remove(scriptfile)

    # Return the catalog
    logger.info("Output file = "+outfile)
    return daoread(outfile)


# Run DAOPHOT PSF
#-------------------
def daopsf(imfile=None,listfile=None,apfile=None,optfile=None,neifile=None,outfile=None,logfile=None,verbose=False,logger=None):
    '''
    This runs DAOPHOT PSF to create a .psf file.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    listfile : str
           The filename of the list of PSF stars.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    apfile : str, optional
           The filename of the aperture photometry file.  By default it is assumed
           that this is the base name of `imfile` with a ".ap" suffix.
    neifile : str, optional
            The output PSF stars and neighbors file.  By default this is the base name of `imfile`
            with a ".nei" suffix.
    outfile : str, optional
            The output filename of the aperture photometry catalog.  By default this is
            the base name of `imfile` with a ".psf" suffix.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".psf.log" suffix.
    verbose : bool, default is False
            Verbose output of the DAOPHOT PSF parameter errors and PSF star profile errors.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    pararr : list
           A list of lists giving the parameters for the various PSF parameter fits.
    parchi : list
           The array of chi values for the various parameter fits.
    profs : structured  numpy array
          The catalog of PSF star profiles giving ID, CHI and FLAG.

    The output catalog and logfile will be created.

    Example
    -------

    .. code-block:: python

        pararr, parchi, profs = daopsf("image.fits","image.lst")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running DAOPHOT PSF -- ")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return None
    # Make sure we have the list file name
    if listfile is None:
        logger.warning("No list filename input")
        return None

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if apfile is None: apfile = base+".ap"
    if outfile is None: outfile = base+".psf"
    if logfile is None: logfile = base+".psf.log"
    if neifile is None: neifile = base+".nei"
    scriptfile = base+".psf.sh"
    for f in [outfile,neifile,logfile,scriptfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,listfile,optfile,apfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return None

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tpsf",dir=".")
    os.close(tid)  # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    toptfile = tbase+".opt"
    tapfile = tbase+".ap"
    listext = os.path.splitext(listfile)[1]
    tlistfile = tbase+listext
    toutfile = tbase+".psf"
    tneifile = tbase+".nei"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)
    os.symlink(listfile,tlistfile)
    os.symlink(apfile,tapfile)

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daophot << END_DAOPHOT >> "+logfile+"\n" \
            "OPTIONS\n" \
            ""+toptfile+"\n" \
            "\n" \
            "ATTACH "+timfile+"\n" \
            "PSF\n" \
            ""+tapfile+"\n" \
            ""+tlistfile+"\n" \
            ""+toutfile+"\n" \
            "\n" \
            "EXIT\n" \
            "EXIT\n" \
            "END_DAOPHOT\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("DAOPHOT PSF failed:"+str(e))
        logger.error(e)
        raise Exception("DAOPHOT failed")

    # Check that the output file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        os.rename(tneifile,neifile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile,tlistfile,tapfile]: os.remove(f)        

        # Get info from the logfile
        if os.path.exists(logfile):
            plines = readlines(logfile)
            # Get parameter errors
            l1 = grep(plines,"Chi    Parameters",index=True)
            l2 = grep(plines,"Profile errors",index=True)
            l3 = grep(plines,"File with PSF stars and neighbors",index=True)
            if len(l1)>0:
                parlines = plines[l1[0]+1:l2[0]-1]
                pararr, parchi = parsepars(parlines)
                minchi = np.min(parchi)
                logger.info("Chi = "+str(minchi))
            # Get profile errors
            if len(l2)>0:
                proflines = plines[l2[0]+1:l3[0]-1]
                if verbose: logger.info(" ".join(proflines))
                profs = parseprofs(proflines)
                logger.info(str(len(profs))+" PSF stars used")
            else:
                logger.error("No DAOPHOT profile errors found in logfile")
                raise Exception("DAOPHOT problem")
    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("DAOPHOT output not found")

    # Delete the script
    os.remove(scriptfile)

    # Return the parameter and profile error information
    logger.info("Output file = "+outfile)
    return pararr, parchi, profs


# Subtract neighbors of PSF stars
#--------------------------------
def subpsfnei(imfile=None,listfile=None,photfile=None,outfile=None,optfile=None,psffile=None,
              nstfile=None,grpfile=None,logfile=None,logger=None):
    '''
    This subtracts neighbors of PSF stars so that an improved PSF can be made.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    listfile : str
           The filename of the list of PSF stars.
    photfile : str, optional
           The filename of the photometry file (normally the .nei aperture photometry file).
           By default it is assumed that this is the base name of `imfile` with a ".nei" suffix.
    outfile : str
            The FITS filename for the image with the neighbors subtracted.  By default this is
            the base name of `imfile` with a "a.fits" suffix.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    psffile : str, optional
           The name of the PSF file.  By default it is assumed that this is the base name of
           `imfile` with a ".psf" suffix.
    nstfile : str, optional
           The name of the output .nst file.
           By default it is assumed that this is the base name of `imfile` with a ".nst" suffix.
    grpfile : str, optional
           The name of the output .grp file that contains information on the groups of stars.
           By default it is assumed that this is the base name of `imfile` with a ".grp" suffix.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".subnei.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    Nothing is returned.  The subtracted image and logfile will be created.

    Example
    -------

    .. code-block:: python

        subpsfnei("image.fits","image.lst","image.nei","imagea.fits")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Subtracting PSF stars neighbors -- ")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return
    # Make sure we have the list file name
    if listfile is None:
        logger.warning("No list filename input")
        return
    # Make sure we have the subtracted image (output) file name
    if outfile is None:
        logger.warning("No subtracted image file name input")
        return

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if photfile is None: photfile = base+".nei"
    if outfile is None: outfile = base+"a.fits"
    if logfile is None: logfile = base+".subnei.log"
    if psffile is None: psffile = base+".psf"
    if nstfile is None: nstfile = base+".nst"
    if grpfile is None: grpfile = base+".grp"
    scriptfile = base+".subnei.sh"
    for f in [outfile,logfile,scriptfile,nstfile,grpfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,listfile,optfile,psffile,photfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tsubnei",dir=".")
    os.close(tid)   # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    toptfile = tbase+".opt"
    tphotfile = tbase+".ap"
    listext = os.path.splitext(listfile)[1]
    tlistfile = tbase+listext
    toutfile = tbase+"a.fits"
    tpsffile = tbase+".psf"
    tnstfile = tbase+".nst"
    tgrpfile = tbase+".grp"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)
    os.symlink(listfile,tlistfile)
    os.symlink(photfile,tphotfile)
    os.symlink(psffile,tpsffile)

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daophot << END_DAOPHOT >> "+logfile+"\n" \
            "OPTIONS\n" \
            ""+toptfile+"\n" \
            "\n" \
            "ATTACH "+timfile+"\n" \
            "GROUP\n" \
            ""+tphotfile+"\n" \
            ""+tpsffile+"\n" \
            "5.\n" \
            ""+tgrpfile+"\n" \
            "NSTAR\n" \
            ""+tpsffile+"\n" \
            ""+tgrpfile+"\n" \
            ""+tnstfile+"\n" \
            "SUBSTAR\n" \
            ""+tpsffile+"\n" \
            ""+tnstfile+"\n" \
            "y\n" \
            ""+tlistfile+"\n" \
            ""+toutfile+"\n" \
            "\n" \
            "EXIT\n" \
            "END_DAOPHOT\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("PSF star neighbor subtracting failed:"+str(e))
        logger.error(e)
        raise Exception("PSF subtraction failed")

    # Check that the output file exists
    if os.path.exists(toutfile):
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        os.rename(tnstfile,nstfile)
        os.rename(tgrpfile,grpfile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile,tlistfile,tphotfile]: os.remove(f)
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("PSF subtraction failed")

    # Delete the script
    os.remove(scriptfile)

    # Print final output filename
    logger.info("Output file = "+outfile)


# Create DAOPHOT PSF
#-------------------
def createpsf(imfile=None,apfile=None,listfile=None,psffile=None,doiter=True,maxiter=5,minstars=6,nsigrej=2,subneighbors=True,
              subfile=None,optfile=None,neifile=None,nstfile=None,grpfile=None,meta=None,logfile=None,verbose=False,logger=None):
    '''
    Iteratively create a DAOPHOT PSF for an image.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    apfile : str, optional
           The filename of the photometry file (normally the .ap aperture photometry file).
           By default it is assumed that this is the base name of `imfile` with a ".ap" suffix.
    listfile : str, optional
           The filename that will contain the final list of PSF stars.  By default this is the
           base name of `imfile` with a ".lst" suffix.
    psffile : str, optional
           The name of the PSF file.  By default it is assumed that this is the base name of
           `imfile` with a ".psf" suffix.
    doiter : bool, default is True
          Iteratively remove bad or suspect PSF stars and refit the PSF.
    maxiter : int, optional, default = 5
            The maximum number of iterations of removing suspect stars.
    minstars : int, optional, default = 6
            The minimum required stars for a PSF.
    nsigrej : float, optional, default = 2
            Reject stars with profile rms scatter higher than 2x the median.
    subneighbors : bool, optional, default = True
             Subtract stars neighboring the PSF stars and then refit the PSF.
    subfile : str, optional
            The FITS filename for the image with the neighbors subtracted.  By default this is
            the base name of `imfile` with a "a.fits" suffix.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".opt" suffix.
    neifile : str, optional
           The name of the output .nei file of PSF stars and neighbors.  By default is is assumed
           that this is the base name of `imfile` with a ".nei" suffix.
    nstfile : str, optional
           The name of the output .nst file created by NSTAR.
           By default it is assumed that this is the base name of `imfile` with a ".nst" suffix.
    grpfile : str, optional
           The name of the output .grp file that contains information on the groups of stars.
           By default it is assumed that this is the base name of `imfile` with a ".grp" suffix.
    meta : str, optional
           The meta-data dictionary for this image.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".subnei.log" suffix.
    verbose : bool, default is False
            Verbose output of the DAOPHOT PSF parameter errors and PSF star profile errors.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    Nothing is returned.  The PSF, subtracted image and logfile are created.

    Example
    -------

    .. code-block:: python

        createpsf("image.fits","image.ap","image.lst","image.psf")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Creating PSF Iteratively --")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if optfile is None: optfile = base+".opt"
    if listfile is None: listfile = base+".lst"
    if apfile is None: apfile = base+".ap"
    if psffile is None: psffile = base+".psf"
    if subfile is None: outfile = base+"a.fits"
    if logfile is None: logfile = base+".cpsf.log"
    if neifile is None: neifile = base+".nei"
    if nstfile is None: nstfile = base+".nsf"
    if grpfile is None: grpfile = base+".grp"
    for f in [outfile,logfile,psffile,nstfile,grpfile,neifile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,optfile,listfile,apfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return

    # Working list file
    wlistfile = listfile+"1"
    if os.path.exists(wlistfile): os.remove(wlistfile)
    shutil.copy(listfile,wlistfile)

    # Make copy of original PSF list
    if os.path.exists(listfile+".orig"): os.remove(listfile+".orig")
    shutil.copy(listfile,listfile+".orig")

    # Iterate
    #---------
    if doiter is False: maxiter=1
    iter = 1
    endflag = 0
    lastchi = 99.99
    dchi_thresh = 0.002
    while (endflag==0):
        logger.info("Iter = "+str(iter))

        # Run DAOPSF
        try:
            pararr, parchi, profs = daopsf(imfile,wlistfile,apfile,logger=logger)
            chi = np.min(parchi)
        except:
            logger.error("Failure in DAOPSF")
            raise

        # Check for bad stars
        nstars = len(profs)
        gdstars = (profs['FLAG'] != 'saturated')
        medsig = np.median(profs['SIG'][gdstars])
        bdstars = (profs['FLAG'] != '') | (profs['SIG']>nsigrej*medsig)
        nbdstars = np.sum(bdstars)
        # Make sure we have enough stars left
        if (nstars-nbdstars < minstars):
            nbdstars = nstars-minstars
            si = np.argsort(profs['SIG'])[::-1]
            bdstars = si[0:nbdstars]  # take the worse ones
        logger.info("  "+str(nbdstars)+" stars with flag or high sig")
        # Delete stars with flags or high SIG values
        if (nbdstars>0) & (nstars>minstars):
            listlines = readlines(wlistfile)
            # Read the list
            lstcat = daoread(wlistfile)
            # Match up with the stars we are deleting
            mid, ind1, ind2 = np.intersect1d(profs[bdstars]['ID'],lstcat['ID'],return_indices=True)
            # Remove the lines from listlines
            newlistlines = remove_indices(listlines,ind2+3)
            # Write new list
            writelines(wlistfile,newlistlines,overwrite=True)
            logger.info("  Removing IDs="+str(" ".join(profs[bdstars]['ID'].astype(str))))
            logger.info("  "+str(nbdstars)+" bad stars removed. "+str(nstars-nbdstars)+" PSF stars left")
        # Should we end
        if (iter==maxiter) | (nbdstars==0) | (nstars<=minstars) | (np.abs(lastchi-chi)<dchi_thresh): endflag=1
        iter = iter+1
        lastchi = chi

    # Subtract PSF star neighbors
    if subneighbors:
        subfile = base+"a.fits"
        try:
            subpsfnei(imfile,wlistfile,neifile,subfile,psffile=psffile,logger=logger)
        except:
            logger.error("Subtracting neighbors failed.  Keeping original PSF file")
        # Check that the subtracted image exist and rerun DAOPSF
        if os.path.exists(subfile):
            # Final run of DAOPSF
            logger.info("Final DAOPSF run")
            try:
                pararr, parchi, profs = daopsf(imfile,wlistfile,apfile,logger=logger)
                chi = np.min(parchi)
            except:
                logger.error("Failure in DAOPSF")
                raise

    # Put information in meta
    if meta is not None:
        meta['PSFCHI'] = (chi,"Final PSF Chi value")
        meta['PSFSTARS'] = (len(profs),"Number of PSF stars")

    # Copy working list to final list
    if os.path.exists(listfile): os.remove(listfile)
    shutil.move(wlistfile,listfile)
    logger.info("Final list of PSF stars in "+listfile+".  Original list in "+listfile+".orig")


# Run ALLSTAR
#-------------
def allstar(imfile=None,psffile=None,apfile=None,subfile=None,outfile=None,optfile=None,meta=None,logfile=None,logger=None):
    '''
    Run DAOPHOT ALLSTAR on an image.

    Parameters
    ----------
    imfile : str
           The filename of the DAOPHOT-ready FITS image.
    psffile : str, optional
           The name of the PSF file.  By default it is assumed that this is the base name of
           `imfile` with a ".psf" suffix.
    apfile : str, optional
           The filename of the photometry file (normally the .ap aperture photometry file).
           By default it is assumed that this is the base name of `imfile` with a ".ap" suffix.
    subfile : str, optional
            The FITS filename for the image with all stars subtracted.  By default this is
            the base name of `imfile` with a "s.fits" suffix.
    outfile : str, optional
            The file name of the final .als source catalog.
    optfile : str, optional
            The option file for `imfile`.  By default it is assumed that this is
            the base name of `imfile` with a ".als.opt" suffix.
    meta : str, optional
           The meta-data dictionary for this image.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".subnei.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    cat : astropy table
        The catalog of ALLSTAR sources.

    The PSF subtracted image and logfile will also be created.

    Example
    -------

    .. code-block:: python

        cat = allstar("image.fits","image.psf")

    '''

    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running ALLSTAR --")

    # Make sure we have the image file name
    if imfile is None:
        logger.warning("No image filename input")
        return

    # Set up filenames, make sure they don't exist
    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if psffile is None: psffile = base+".psf"
    if optfile is None: optfile = base+".als.opt"
    if apfile is None: apfile = base+".ap"
    if subfile is None: subfile = base+"s.fits"
    if outfile is None: outfile = base+".als"
    if logfile is None: logfile = base+".als.log"
    scriptfile = base+".als.sh"
    for f in [outfile,subfile,logfile,scriptfile]:
        if os.path.exists(f): os.remove(f)

    # Check that necessary files exist
    for f in [imfile,psffile,apfile,optfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tals",dir=".")
    os.close(tid)  # close open file
    tbase = os.path.basename(tfile)
    timfile = tbase+".fits"
    toptfile = tbase+".als.opt"
    tapfile = tbase+".ap"
    tpsffile = tbase+".psf"
    tsubfile = tbase+"s.fits"
    toutfile = tbase+".als"
    os.symlink(imfile,timfile)
    os.symlink(optfile,toptfile)
    os.symlink(apfile,tapfile)
    os.symlink(psffile,tpsffile)

    # Load the option file lines
    optlines = readlines(optfile)
    # Lines for the DAOPHOT ALLSTAR script
    lines = ["#!/bin/sh\n",
             "allstar << END_ALLSTAR >> "+logfile+"\n"]
    lines += optlines
    lines += ["\n",
              timfile+"\n",
              tpsffile+"\n",
              tapfile+"\n",
              toutfile+"\n",
              tsubfile+"\n",
              "EXIT\n",
              "EXIT\n",
              "END_ALLSTAR\n"]
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Copy option file to daophot.opt
    if os.path.exists("allstar.opt") is False: shutil.copyfile(optfile,"allstar.opt")

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=False)
        if retcode < 0:
            logger.warning("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.warning("ALLSTAR failed:"+str(e))
        logger.warning(e)
        raise Exception("ALLSTAR failed")

    # Check that the output file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        os.rename(tsubfile,subfile)
        # Remove the temporary links
        for f in [tfile,timfile,toptfile,tpsffile,tapfile]: os.remove(f)        
        # How many sources converged
        num = numlines(outfile)-3
        logger.info(str(num)+" stars converged")
        logger.info("Output file = "+outfile)
    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("ALLSTAR failed")

    # Delete the script
    os.remove(scriptfile)

    # Put information in the header
    if meta is not None:
        meta["NALLSTAR"] = (num,"Number of ALLSTAR converged sources")

    # Return the final catalog
    return daoread(outfile)


# Calculate aperture corrections
#-------------------------------
def daogrow(photfile,aperfile,meta,nfree=3,fixedvals=None,maxerr=0.2,logfile=None,logger=None):
    '''
    Run DAOGROW that calculates curve of growths using aperture photometry.

    Parameters
    ----------
    photfile : str
             The aperture photometry file.
    aperfile : str
             The file containing the apertures used for the aperture photometry.
    meta : astropy header
           The meta-data dictionary for the image.
    nfree : float, optional, default = 3
          The number of parameters to fit.  Max is 5.
    fixedvals : float, optional
          The values for the parameters that are fixed.  Shoul have 5-nfree elements.
          By default they are [1.03, 0.2, 0.1, 0.6, 0.0].
    maxerr : float, optional, default = 0.2
           The maximum error to allow in DAOGROW.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".gro.log" suffix.
    logger : logging object
           The logger to use for the logging information.

    Returns
    -------
    totcat : astropy table
           The aperture corrected photometry file (.tot).
    Also, the .tot and other DAOGROW files are created. 

    Example
    -------

    .. code-block:: python

        totcat = daogrow("im101a.ap","photo.opt",meta)

    '''
    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Running DAOGROW --")

    # Make sure we have apfile
    if photfile is None:
        logger.warning("No photfile input")
        return None
    # Make sure we have apfile
    if aperfile is None:
        logger.warning("No aperfile input")
        return None
    # Make sure we have the meta-data dictionary
    if meta is None:
        logger.warning("No meta input")
        return None

    # Checked number of elements for fixedvals
    if fixedvals is not None:
        if len(fixedvals) != 5-nfree:
            logger.warning("Fixedvals must have 5-nfree elements."+str(len(fixedvals))+" found.")
            return None

    # Check that necessary files exist
    for f in [photfile,aperfile]:
        if os.path.exists(f) is False:
            logger.warning(apfile+" NOT found")
            return None

    # Set up filenames, make sure they don't exist
    base = os.path.basename(photfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    if logfile is None: logfile = base+".gro.log"
    outfile = base+".tot"
    scriptfile = base+".gro.sh"
    for f in [logfile,scriptfile,outfile,base+".poi",base+".cur",base+".gro",base+".crl"]:
        if os.path.exists(f): os.remove(f)

    # Make temporary short filenames to DAOPHOT can handle them
    tid,tfile = tempfile.mkstemp(prefix="tcoo",dir=".")
    os.close(tid)   # close open file
    tbase = os.path.basename(tfile)
    tphotfile = tbase+".ap"
    taperfile = tbase+".opt"
    textfile = tbase+".ext"
    tinffile = tbase+".inf"
    toutfile = tbase+".tot"
    os.symlink(photfile,tphotfile)
    os.symlink(aperfile,taperfile)

    # Write the .inf and .ext files
    # F1-00507800_01a                11  04 51  1.900    30.000
    dateobs = meta['DATE-OBS']
    timearr = (dateobs.split('T')[1]).split(':')
    if meta.get('airmass') is not None:
        airmass = meta['airmass']
    else: 
        airmass = 1.0
    lines = " %-23s %9d %3d %2d %6.3f %9.3f\n" % (tbase,int(timearr[0]),int(timearr[1]),int(float(timearr[2])),airmass,meta['exptime'])
    writelines(tinffile,lines)
    # .ext just has the .ap filename
    writelines(textfile,tphotfile+"\n")

    # The fixed values for the other parameters that are fixed
    if fixedvals is None:
        allfixedvals = [1.03, 0.2, 0.1, 0.6, 0.0]
        fixedvals = allfixedvals[nfree:]

    # Lines for the DAOPHOT script
    lines = "#!/bin/sh\n" \
            "daogrow << DONE >> "+logfile+"\n" \
            ""+taperfile+"\n" \
            "\n" \
            ""+tinffile+"\n" \
            ""+textfile+"\n" \
            ""+str(nfree)+"\n" \
            ""+",".join(np.array(fixedvals).astype(str))+"\n" \
            ""+str(maxerr)+"\n" \
            "DONE\n"
    # Write the script
    f = open(scriptfile,'w')
    f.writelines(lines)
    f.close()
    os.chmod(scriptfile,509)

    # Run the script
    try:
        retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=False)
        if retcode < 0:
            logger.error("Child was terminated by signal"+str(-retcode))
        else:
            pass
    except OSError as e:
        logger.error("DAOGROW failed:"+str(e))
        logger.error(e)
        raise Exception("DAOGROW failed")

    # Check that the outfile file exists
    if os.path.exists(toutfile) is True:
        # Move output file to the final filename
        os.rename(toutfile,outfile)
        if os.path.exists(tbase+".poi"): os.rename(tbase+".poi",base+".poi")
        if os.path.exists(tbase+".cur"): os.rename(tbase+".cur",base+".cur")
        if os.path.exists(tbase+".gro"): os.rename(tbase+".gro",base+".gro")
        if os.path.exists(tbase+".crl"): os.rename(tbase+".crl",base+".crl")
        # Remove the temporary links
        for f in [tfile,tphotfile,taperfile,tinffile,textfile]: os.remove(f)

    # Failure
    else:
        logger.error("Output file "+outfile+" NOT Found")
        raise Exception("Output not found")

    # Delete the script
    os.remove(scriptfile)

    # Load and return the catalog
    logger.info("Output file = "+outfile)

    # Return the .tot catalog
    return daoread(outfile)


# Calculate aperture corrections
#-------------------------------
def apcor(imfile=None,listfile=None,psffile=None,meta=None,optfile=None,alsoptfile=None,logger=None):
    '''
    Calculate the aperture correction for an image.

    Parameters
    ----------
    imfile : str
           The filename of the PSF-neighbor-subtracted image.
    listfile : str
           The list of PSF stars.
    psffile : str, optional
           The name of the PSF file.
    meta : astropy headre
         The meta-data dictionary for the image.
    optfile : str, optional
            The DAOPHOT option file for `imfile`.
    alsoptfile : str, optional
            The ALLSTAR option file for `imfile`.
    logfile : str, optional
            The name of the logfile to constrain the output of the DAOPHOT FIND
            run.  By default this is the base name of `imfile` with a ".daogrow.log" suffix.
    logger : logging object
           The logger to use for the loggin information.

    Returns
    -------
    apcor : float
        The aperture correction in magnitudes.

    Example
    -------

    .. code-block:: python

        apcor = apcor("im101a.fits","im101.lst","im101.psf",meta,"im101.opt")

    '''
    if logger is None: logger=basiclogger('phot')   # set up basic logger if necessary
    logger.info("-- Calculating aperture correction --")

    # Make sure we have imfile
    if imfile is None:
        logger.warning("No image filename input")
        return
    # Make sure we have listfile
    if listfile is None:
        logger.warning("No listfile input")
        return
    # Make sure we have psffile
    if psffile is None:
        logger.warning("No psffile input")
        return
    # Make sure we have optfile
    if optfile is None:
        logger.warning("No optfile input")
        return
    # Make sure we have alsoptfile
    if alsoptfile is None:
        logger.warning("No alsoptfile input")
        return
    # Make sure we have meta
    if meta is None:
        logger.warning("No meta input")
        return


    # Check that necessary files exist
    for f in [imfile,listfile,psffile,optfile]:
        if os.path.exists(f) is False:
            logger.warning(f+" NOT found")
            return

    base = os.path.basename(imfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]

    # Step 1: Get aperture photometry for the PSF stars on the PSF-neighbor subtracted image
    logger.info("Getting aperture photometry for PSF stars")
    apertures = [3.0, 3.7965, 4.8046, 6.0803, 7.6947, 9.7377, 12.3232, 15.5952, 19.7360, \
                 24.9762, 31.6077, 40.0000, 50.0000]
    apersfile = base+".apers"
    apcat, maglim = daoaperphot(imfile,listfile,apertures,optfile=optfile,apersfile=apersfile,logger=logger)

    # Step 2: Get PSF photometry from the same image
    psfcat = allstar(imfile,psffile,base+".ap",optfile=alsoptfile,logger=logger)

    # Step 3: Run DAOGROW
    #  it creates a .tot, .cur, .poi files
    #  use .tot and .als files to calculate delta mag for each star (see mkdel.pro)
    #  and then a total aperture correction for all the stars.
    totcat = daogrow(base+".ap",apersfile,meta,logger=logger)
    # Check that the magnitudes arent' all NANs, this can sometimes happen
    if np.sum(np.isnan(totcat['MAG'])) > 0:
        logger.info("DAOGROW .tot file has NANs.  Trying 2 free parameters instead.")
        totcat = daogrow(base+".ap",apersfile,meta,nfree=2,logger=logger)
    if np.sum(np.isnan(totcat['MAG'])) > 0:
        logger.info("DAOGROW .tot file has NANs.  Trying 4 free parameters instead.")
        totcat = daogrow(base+".ap",apersfile,meta,nfree=4,logger=logger)

    # Step 4: Calculate median aperture correction
    totcat = daoread(base+".tot")
    # Match up with the stars we are deleting
    mid, ind1, ind2 = np.intersect1d(psfcat['ID'],totcat['ID'],return_indices=True)
    apcorr = np.median(psfcat[ind1]['MAG']-totcat[ind2]['MAG'])

    logger.info("aperture correction = %7.3f mag" % apcorr)

    return apcorr
