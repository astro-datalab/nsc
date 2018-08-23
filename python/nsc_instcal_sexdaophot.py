#!/usr/bin/env python
#
# NSC_INSTCAL_SEXDAOPHOT.PY -- Run SExtractor and DAOPHOT on an exposure
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180819'  # yyyymmdd

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table
import time
import shutil
import re
import subprocess
import glob
import logging
import socket
#from scipy.signal import convolve2d
from scipy.ndimage.filters import convolve
import astropy.stats
import struct

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# Get NSC directories
def getnscdirs(version=None):
    # Host
    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    # Version
    verdir = ""
    if version is not None:
       verdir = version if version.endswith('/') else version+"/"
    # on thing/hulk use
    if (host == "thing") | (host == "hulk"):
        basedir = "/dl1/users/dnidever/nsc/instcal/"+verdir
        tmproot = "/d0/dnidever/nsc/instcal/"+verdir+"tmp/"
    # on gp09 use
    if (host == "gp09") | (host == "gp08") | (host == "gp07") | (host == "gp06") | (host == "gp05"):
        basedir = "/net/dl1/users/dnidever/nsc/instcal/"+verdir
        tmproot = "/data0/dnidever/nsc/instcal/"+verdir+"tmp/"
    return basedir,tmproot

# Standard grep function that works on string list
def grep(lines,expr,index=False):
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

# Parse the DAOPHOT PSF profile errors
def parseprofs(lines):
    dtype = np.dtype([('id',int),('chi',float),('flag',np.str_,3)])
    profs = np.zeros(len(lines)*5,dtype=dtype)
    profs['id'] = -1
    cnt = 0L
    for i in range(len(lines)):
        l = lines[i].rstrip()
        if l != "":
            # Loop through five columns
            for j in range(5):
                line1 = l[j*17:j*17+17]
                id1 = line1[0:7]
                chi1 = line1[7:14]
                flag1 = line1[14:17]
                if id1.strip() != "":
                    profs[cnt]['id'] = int(id1)
                    profs[cnt]['chi'] = float(chi1)
                    profs[cnt]['flag'] = flag1.strip()
                    cnt = cnt + 1
    # Trimming any blank ones
    gd = (profs['id'] > -1)
    profs = profs[gd]
    return profs

# Parse the DAOPHOT PSF parameter errors
def parsepars(lines):
    out = []
    chi = []
    for i in range(len(lines)):
        line1 = lines[i].strip()
        if line1[0:2] == ">>": line1=line1[2:]  # strip leading >>
        line1.strip()
        arr = line1.split()                # split on whitespace
        if len(arr)>0:
            chi.append(float(arr[0]))
            out.append(arr)
    return out, chi

# Read DAOPHOT files
def daoread(fil):
    if os.path.exists(fil) is False:
        print(fil+" NOT found")
        return None
    f = open(fil,'r')
    lines = f.readlines()
    f.close()
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
    # NL = 1 other (output from FIND, PEAK, or NSTAR) or ALLSTAR
    # NL = 0 a file without a header
    
    # Check number of columns
    ncols = len(lines[3].split())

    # NL = 1  coo file
    if (nl==1) & (ncols==7):
        dtype = np.dtype([('id',long),('x',float),('y',float),('mag',float),('sharp',float),('round',float),('round2',float)])
        cat = np.zeros(nstars,dtype=dtype)
        lengths = np.array([7,9,9,9,9,9,9])
        lo = np.concatenate((np.array([0]), np.cumsum(lengths[0:-1])))
        hi = lo+lengths
        names = cat.dtype.names
        for i in range(nstars):
            line1 = lines[i+3]
            for j in range(len(names)):
                cat[i][names[j]] = np.array(line1[lo[j]:hi[j]],dtype=dtype[names[j]])
    # NL = 1  als file
    elif (nl==1) & (ncols==9):
        dtype = np.dtype([('id',long),('x',float),('y',float),('mag',float),('err',float),('sky',float),('iter',float),('chi',float),('sharp',float)])
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
        print("Reading aperture photometry files not supported yet.")
        return
    # NL = 3  list
    elif nl==3:
        dtype = np.dtype([('id',long),('x',float),('y',float),('mag',float),('err',float),('sky',float)])
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
        return
    return cat

# Remove indices from a list
def remove_indices(lst,index):
    newlst = []
    for i in range(len(lst)):
       if i not in index: newlst.append(lst[i])
    return newlst

# Little function used by numlines
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

# Read number of lines in a file
def numlines(fil):
    with open(fil, "r") as f:
        return (sum(bl.count("\n") for bl in blocks(f)))

    # Could also use this
    #count=0
    #for line in open(fil): count += 1

# Class to represent an exposure to process
class Exposure:

    def __init__(self,fluxfile,wtfile,maskfile,nscver="t3a"):
        self.fluxfile = fluxfile
        self.wtfile = wtfile
        self.maskfile = maskfile
        base = os.path.basename(fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        self.base = base
        self.nscver = nscver
        self.logfile = base+".log"
        self.logger = None
        self.chip = None

        if os.path.exists(fluxfile) is False:
            print(fluxfile+" NOT found")
            return
        if os.path.exists(wtfile) is False:
            print(wtfile+" NOT found")
            return
        if os.path.exists(maskfile) is False:
            print(maskfile+" NOT found")
            return

        # Set up logging to screen and logfile
        #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
        logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
        rootLogger = logging.getLogger()

        #logfile = tmpdir+"/"+base+".log"
        #fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, fileName))
        fileHandler = logging.FileHandler(self.logfile)
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        rootLogger.addHandler(consoleHandler)
        rootLogger.setLevel(logging.NOTSET)
        self.logger = rootLogger

        self.logger.info("Starting logfile at "+self.logfile)
        
        #rootLogger.info("  Temporary directory is: "+tmpdir)

        
    # Load chip
    def loadchip(self,extension,fluxfile="flux.fits",wtfile="wt.fits",maskfile="mask.fits"):
        # Load the data
        self.logger.info(" Loading chip "+str(extension))
        try:
            flux,fhead = fits.getdata(self.fluxfile,extension,header=True)
            fhead0 = fits.getheader(self.fluxfile,0)  # add PDU info
            fhead.extend(fhead0,unique=True)
            wt,whead = fits.getdata(self.wtfile,extension,header=True)
            mask,mhead = fits.getdata(self.maskfile,extension,header=True)
        except:
            self.logger.info("No extension "+str(extension))
            return
        # Write the data to the appropriate files
        if os.path.exists(fluxfile):
            os.remove(fluxfile)
        fits.writeto(fluxfile,flux,header=fhead,output_verify='warn')
        if os.path.exists(wtfile):
            os.remove(wtfile)
        fits.writeto(wtfile,wt,header=whead,output_verify='warn')
        if os.path.exists(maskfile):
            os.remove(maskfile)
        fits.writeto(maskfile,mask,header=mhead,output_verify='warn')        
        # Create the chip object
        self.chip = Chip(fluxfile,wtfile,maskfile)
        self.chip.bigextension = extension
        self.chip.nscver = self.nscver
        # Add logger information
        self.chip.logger = self.logger
        
    # setup
    def setup():
        pass
    
    # teardown
    def teardown():
        pass
        
        
# Class to represent a single chip of an exposure
class Chip:

    def __init__(self,fluxfile,wtfile,maskfile):
        self.fluxfile = fluxfile
        self.wtfile = wtfile
        self.maskfile = maskfile
        self.bigextension = None
        base = os.path.basename(fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        self.dir = os.path.abspath(os.path.dirname(fluxfile))
        self.base = base
        self.header = fits.getheader(fluxfile,0)
        self.sexfile = self.dir+"/"+self.base+"_sex.fits"
        self.daofile = self.dir+"/"+self.base+"_dao.fits"
        self.sexcatfile = None
        self.sexcat = None
        self.seeing = None
        self.daomaglim = None
        # Internal hidden variables
        self._rdnoise = None
        self._gain = None
        self._ccdnum = None
        self._pixscale = None
        self._saturate = None
        self._wcs = None
        self._exptime = None
        self._instrument = None
        self._plver = None
        self._cpfwhm = None
        # Logger
        self.logger = None

    
    def __repr__(self):
        return "Chip object"
        
        
    @property
    def rdnoise(self):
        # We have it already, just return it
        if self._rdnoise is not None:
            return self._rdnoise
        # Can't get rdnoise, no header yet
        if self.header is None:
            self.logger.warning("Cannot get RDNOISE, no header yet")
            return None
        # Check DECam style rdnoise
        if "RDNOISEA" in self.header.keys():
            rdnoisea = self.header["RDNOISEA"]
            rdnoiseb = self.header["RDNOISEB"]
            rdnoise = (rdnoisea+rdnoiseb)*0.5
            self._rdnoise = rdnoise
            return self._rdnoise
        # Get rdnoise from the header
        for name in ['RDNOISE','READNOIS','ENOISE']:
            # We have this key, set _rndoise and return
            if name in self.header.keys():
                self._rdnoise = self.header[name]
                return self._rdnoise
        self.logger.warning('No RDNOISE found')
        return None
            
    @property
    def gain(self):
        # We have it already, just return it
        if self._gain is not None:
            return self._gain
        try:
            gainmap = { 'c4d': lambda x: 0.5*(x.get('gaina')+x.get('gainb')),
                        'k4m': lambda x: x.get('gain'),
                        'ksb': lambda x: [1.3,1.5,1.4,1.4][ccdnum-1] }  # bok gain in HDU0, use list here
            gain = gainmap[self.instrument](self.header)
        except:
            gainmap_avg = { 'c4d': 3.9845419, 'k4m': 1.8575, 'ksb': 1.4}
            gain = gainmap_avg[self.instrument]
        self._gain = gain
        return self._gain
            
        ## Can't get gain, no header yet
        #if self.header is None:
        #    print("Cannot get GAIN, no header yet")
        ## Get rdnoise from the header
        #for name in ['GAIN','EGAIN']:
        #    # We have this key, set _gain and return
        #    if self.header.has_key(name):
        #        self._gain = self.header[name]
        #        return self._gain
        #print('No GAIN found')
        #return None
            
    @property
    def ccdnum(self):
        # We have it already, just return it
        if self._ccdnum is not None:
            return self._ccdnum
        # Can't get ccdnum, no header yet
        if self.header is None:
            self.logger.warning("Cannot get CCDNUM, no header yet")
            return None
        # Get ccdnum from the header
        # We have this key, set _rndoise and return
        if 'CCDNUM' in self.header.keys():
            self._ccdnum = self.header['CCDNUM']
            return self._ccdnum
        self.logger.warning('No CCDNUM found')
        return None
            
    @property
    def pixscale(self):
        # We have it already, just return it
        if self._pixscale is not None:
            return self._pixscale
        pixmap = { 'c4d': 0.27, 'k4m': 0.258, 'ksb': 0.45 }
        try:
            pixscale = pixmap[self.instrument]
            self._pixscale = pixscale
            return self._pixscale
        except:
            self._pixscale = np.max(np.abs(self.wcs.pixel_scale_matrix))
            return self._pixscale
            
    @property
    def saturate(self):
        # We have it already, just return it
        if self._saturate is not None:
            return self._saturate
        # Can't get saturate, no header yet
        if self.header is None:
            self.logger.warning("Cannot get SATURATE, no header yet")
            return None
        # Get saturate from the header
        # We have this key, set _saturate and return
        if 'SATURATE' in self.header.keys():
            self._saturate = self.header['SATURATE']
            return self._saturate
        self.logger.warning('No SATURATE found')
        return None
    
    @property
    def wcs(self):
        # We have it already, just return it
        if self._wcs is not None:
            return self._wcs
        # Can't get wcs, no header yet
        if self.header is None:
            self.logger.warning("Cannot get WCS, no header yet")
            return None
        try:
            self._wcs = WCS(self.header)
            return self._wcs
        except:
            self.logger.warning("Problem with WCS")
            return None
            
    @property
    def exptime(self):
        # We have it already, just return it
        if self._exptime is not None:
            return self._exptime
        # Can't get exptime, no header yet
        if self.header is None:
            self.logger.warning("Cannot get EXPTIME, no header yet")
            return None
        # Get rdnoise from the header
        # We have this key, set _rndoise and return
        if 'EXPTIME' in self.header.keys():
                self._exptime = self.header['EXPTIME']
                return self._exptime
        print('No EXPTIME found')
        return None

    @property
    def instrument(self):
        # We have it already, just return it
        if self._instrument is not None:
            return self._instrument
        # Can't get instrument, no header yet
        if self.header is None:
            self.logger.warning("Cannot get INSTRUMENT, no header yet")
            return None
        # instrument, c4d, k4m or ksb
        # DTINSTRU = 'mosaic3 '
        # DTTELESC = 'kp4m    '
        # Bok 90Prime data has
        if self.header.get("DTINSTRU") == 'mosaic3':
            self._instrument = 'k4m'
            return self._instrument
        elif self.header.get("DTINSTRU") == '90prime':
            self._instrument = 'ksb'
            return self._instrument
        else:
            self._instrument = 'c4d'
            return self._instrument

    @property
    def plver(self):
        # We have it already, just return it
        if self._plver is not None:
            return self._plver
        # Can't get plver, no header yet
        if self.header is None:
            self.logger.warning("Cannot get PLVER, no header yet")
            return None
        plver = self.header.get('PLVER')
        if plver is None:
            self._plver = 'V1.0'
        self._plver = plver
        return self._plver

    @property
    def cpfwhm(self):
        # We have it already, just return it
        if self._cpfwhm is not None:
            return self._cpfwhm
        # Can't get fwhm, no header yet
        if self.header is None:
            self.logger.warning("Cannot get CPFWHM, no header yet")
            return None
        # FWHM values are ONLY in the extension headers
        cpfwhm_map = { 'c4d': 1.5 if self.header.get('FWHM') is None else self.header.get('FWHM')*0.27, 
                       'k4m': 1.5 if self.header.get('SEEING1') is None else self.header.get('SEEING1'),
                       'ksb': 1.5 if self.header.get('SEEING1') is None else self.header.get('SEEING1') }
        cpfwhm = cpfwhm_map[self.instrument]
        self._cpfwhm = cpfwhm
        return self._cpfwhm


    # Write SE catalog in DAO format
    def sextodao(self,cat=None,outfile=None,format="lst"):

        # Not enough inputs
        if cat is None:
            self.logger.warning("No catalog input")
            return
        if outfile is None:
            self.logger.warning("No outfile given")
            return
        # Delete outfile
        if os.path.exists(outfile): os.remove(outfile)
    
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
                    (self.header['naxis1'],self.header['naxis2'],1000.0,self.saturate,100.0,3.0,self.gain,self.rdnoise/self.gain,3.9))
            #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
            # Write the data
            for e in cat:
                f.write("%7d %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f\n" %
                        (e["NUMBER"],e["X_IMAGE"]+1,e["Y_IMAGE"]+1,e["MAG_AUTO"],0.6,0.0,0.0))
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
                    (self.header['naxis1'],self.header['naxis2'],1000.0,self.saturate,100.0,3.0,self.gain,self.rdnoise/self.gain,3.9))
            #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
            # Write the data
            for e in cat:
                f.write("%7d %8.3f %8.3f %8.3f %8.3f %8.3f\n" %
                        (e["NUMBER"],e["X_IMAGE"]+1,e["Y_IMAGE"]+1,e["MAG_AUTO"],e["MAGERR_AUTO"],0.3))
            f.close()

        # "ap" file from PHOTOMETRY
        elif format == "ap":
            self.logger.warning(".ap files not supported yet")
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
                    (self.header['naxis1'],self.header['naxis2'],1000.0,self.saturate,100.0,3.0,self.gain,self.rdnoise/self.gain,3.9))
            #f.write("  3  2046  4094  1472.8 38652.0   80.94    3.00    3.91    1.55    3.90\n")
            # Write the data
            for e in cat:
                f.write("%7d %8.3f %8.3f %8.3f %8.4f %8.3f %8.0f %8.3f %8.3f\n" %
                        (e["NUMBER"],e["X_IMAGE"]+1,e["Y_IMAGE"]+1,e["MAG_AUTO"],e["MAGERR_AUTO"],1500.0,1,1.0,0.0))
            f.close()

        # Not supported
        else:
            self.logger.warning(format+" NOT supported")
            return


    # Run Source Extractor
    #---------------------
    def runsex(self,outfile=None):

        self.logger.info("-- Running SExtractor --")

        if outfile is None: outfile = "flux_sex.cat.fits"
        fluxfile = "flux_sex.fits"
        wtfile = "wt_sex.fits"
        maskfile = "mask_sex.fits"
        logfile = "flux_sex.log"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(fluxfile): os.remove(fluxfile)
        if os.path.exists(wtfile): os.remove(wtfile)
        if os.path.exists(maskfile): os.remove(maskfile)
        if os.path.exists(logfile): os.remove(logfile)

        # Load the data
        flux,fhead = fits.getdata(self.fluxfile,header=True)
        wt,whead = fits.getdata(self.wtfile,header=True)
        mask,mhead = fits.getdata(self.maskfile,header=True)

        # 3a) Make subimages for flux, weight, mask

        # Turn the mask from integer to bitmask
        if ((self.instrument=='c4d') & (self.plver>='V3.5.0')) | (self.instrument=='k4m') | (self.instrument=='ksb'):
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
        if (self.instrument=='c4d') & (self.plver<'V3.5.0'):
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
        shutil.copy(self.fluxfile,fluxfile)
        fits.writeto(wtfile,wt,header=whead,output_verify='warn')


        # 3b) Make SExtractor config files
        # Copy the default files
        basedir, tmpdir = getnscdirs(self.nscver)
        shutil.copyfile(basedir+"config/default.conv","default.conv")
        shutil.copyfile(basedir+"config/default.nnw","default.nnw")
        shutil.copyfile(basedir+"config/default.param","default.param")

        # Read in configuration file and modify for this image
        f = open(basedir+'config/default.config', 'r') # 'r' = read
        lines = f.readlines()
        f.close()

        # Gain, saturation, pixscale

        # Things to change
        # SATUR_LEVEL     59000.00         # level (in ADUs) at which arises saturation
        # GAIN            43.52             # detector gain in e-/ADU.
        # SEEING_FWHM     1.46920            # stellar FWHM in arcsec
        # WEIGHT_IMAGE  F4-00507860_01_comb.mask.fits

        filter_name = ''
        cnt = 0L
        for l in lines:
            # CATALOG_NAME
            m = re.search('^CATALOG_NAME',l)
            if m != None:
                lines[cnt] = "CATALOG_NAME     "+outfile+"         # name of the output catalog\n"
            # FLAG_IMAGE
            m = re.search('^FLAG_IMAGE',l)
            if m != None:
                lines[cnt] = "FLAG_IMAGE     "+maskfile+"         # filename for an input FLAG-image\n"
            # WEIGHT_IMAGE
            m = re.search('^WEIGHT_IMAGE',l)
            if m != None:
                lines[cnt] = "WEIGHT_IMAGE     "+wtfile+"  # Weight image name\n"
            # SATUR_LEVEL
            m = re.search('^SATUR_LEVEL',l)
            if m != None:
                lines[cnt] = "SATUR_LEVEL     "+str(self.saturate)+"         # level (in ADUs) at which arises saturation\n"
            # Gain
            m = re.search('^GAIN',l)
            if m != None:
                lines[cnt] = "GAIN            "+str(self.gain)+"            # detector gain in e-/ADU.\n"
            # SEEING_FWHM
            m = re.search('^SEEING_FWHM',l)
            if m != None:
                lines[cnt] = "SEEING_FWHM     "+str(self.cpfwhm)+"            # stellar FWHM in arcsec\n"
            # PHOT_APERTURES, aperture diameters in pixels
            m = re.search('^PHOT_APERTURES',l)
            if m != None:
                aper_world = np.array([ 0.5, 1.0, 2.0, 3.0, 4.0]) * 2  # radius->diameter, 1, 2, 4, 6, 8"
                aper_pix = aper_world / self.pixscale
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
            fits.writeto(maskfile,newmask,header=mhead,output_verify='warn')

        # 3c) Run SExtractor
        try:
            # Save the SExtractor info to a logfile
            sf = open(logfile,'w')
            retcode = subprocess.call(["sex",fluxfile,"-c","default.config"],stdout=sf,stderr=subprocess.STDOUT)
            sf.close()
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("SExtractor Execution failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # Load the catalog and keep it in memory for later use
            self.sexcatfile = outfile
            cat = Table.read(outfile,2)
            self.sexcat = cat
            # How many sources were detected, final catalog file
            head = fits.getheader(outfile,2)
            self.logger.info(str(head["naxis2"])+" sources detected")
            self.logger.info("Final catalog is "+outfile)

        # Delete temporary files
        os.remove(fluxfile)
        os.remove(maskfile)
        os.remove(wtfile)
        #os.remove("default.conv")


    # Determine FWHM using SE catalog
    #--------------------------------
    def sexfwhm(self):

        self.logger.info("-- Determining seeing FWHM using SExtractor catalog --")

        # Make sure we have the SE catalog
        if self.sexcatfile is None:
            self.logger.warning("No SE catalog found")
            return
        # Load the catalog if necessary
        if self.sexcat is None:
            self.sexcat = Table.read(self.sexcatfile,2)
        # Select good sources
        gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.05) & (self.sexcat['CLASS_STAR']>0.8))
        ngdcat = np.sum(gdcat)
        # CLASS_STAR is not as reliable if the seeing is bad
        if (ngdcat<10) & (self.cpfwhm>1.8):
            gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.05))
            ngdcat = np.sum(gdcat)            
        # Not enough sources, lower thresholds
        if (ngdcat<10):
            gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.08))
            ngdcat = np.sum(gdcat)            
        medfwhm = np.median(self.sexcat[gdcat]['FWHM_WORLD']*3600.)
        self.logger.info('  FWHM = %5.2f arcsec (%d sources)' % (medfwhm, ngdcat))
        self.seeing = medfwhm


    # Pick PSF candidates using SE catalog
    #-------------------------------------
    def sexpickpsf(self,nstars=100):

        self.logger.info("-- Picking PSF stars using SExtractor catalog --")

        # Make sure we have the SE catalog
        if self.sexcatfile is None:
            self.logger.warning("No SE catalog found")
            return
        # Load the catalog if necessary
        if self.sexcat is None:
            self.sexcat = Table.read(self.sexcatfile,2)
        # Make sure we have seeing calculated
        if self.seeing is None: self.sexfwhm()

        base = os.path.basename(self.sexfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        outfile = base+".lst"

        # Select good sources
        gdcat1 = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.05) & (self.sexcat['CLASS_STAR']>0.8))
        ngdcat1 = np.sum(gdcat1)
        # Bright and faint limit, use 5th and 95th percentile
        minmag, maxmag = np.sort(self.sexcat[gdcat1]['MAG_AUTO'])[[int(np.round(0.05*ngdcat1)),int(np.round(0.95*ngdcat1))]]
        # Select stars with
        # -good FWHM values
        # -good clas_star values (unless FWHM too large)
        # -good mag range, bright but not too bright
        if self.cpfwhm<1.8:
            gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.1) & (self.sexcat['CLASS_STAR']>0.8) & 
                     (self.sexcat['FWHM_WORLD']*3600.>0.5*self.seeing) & (self.sexcat['FWHM_WORLD']*3600.<1.5*self.seeing) &
                     (self.sexcat['MAG_AUTO']>(minmag+1.0)) & (self.sexcat['MAG_AUTO']<(maxmag-0.5)))
            ngdcat = np.sum(gdcat)
        # Do not use CLASS_STAR if seeing bad, not as reliable
        else:
            gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.1) & 
                     (self.sexcat['FWHM_WORLD']*3600.>0.5*self.seeing) & (self.sexcat['FWHM_WORLD']*3600.<1.5*self.seeing) &
                     (self.sexcat['MAG_AUTO']>(minmag+1.0)) & (self.sexcat['MAG_AUTO']<(maxmag-0.5)))
            ngdcat = np.sum(gdcat)
        # No candidate, loosen cuts
        if ngdcat<10:
            self.logger.info("Too few PSF stars on first try. Loosening cuts")
            gdcat = ((self.sexcat['MAG_AUTO']< 50) & (self.sexcat['MAGERR_AUTO']<0.15) & 
                     (self.sexcat['FWHM_WORLD']*3600.>0.2*self.seeing) & (self.sexcat['FWHM_WORLD']*3600.<1.8*self.seeing) &
                     (self.sexcat['MAG_AUTO']>(minmag+0.5)) & (self.sexcat['MAG_AUTO']<(maxmag-0.5)))
            ngdcat = np.sum(gdcat)
        # No candidates
        if ngdcat==0:
            self.logger.error("No good PSF stars found")
            raise

        # Candidate PSF stars, use only Nstars, and sort by magnitude
        psfcat = np.sort(self.sexcat[gdcat],order='MAG_AUTO')
        if ngdcat>nstars: psfcat=psfcat[0:nstars]

        # Output them in DAO format
        self.sextodao(psfcat,outfile,format="lst")
        if os.path.exists(outfile) is False:
            self.logger.error("Output file "+outfile+" NOT found")
            raise

        # Do we a need separate aperture photometry file?

    # Make DAOPHOT option files
    #--------------------------
    def mkopt(self,VA=1,LO=7.0,TH=3.5,LS=0.2,HS=1.0,LR=-1.0,HR=1.0,WA=-2,AN=-6,
                   EX=5,PE=0.75,PR=5.0,CR=2.5,CE=6.0,MA=50.0,RED=1.0,WA2=0.0,
                   fitradius_fwhm=1.0):

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

        self.logger.info("-- Creating DAOPHOT option file --")

        base = os.path.basename(self.daofile)
        dir = os.path.abspath(os.path.dirname(self.daofile))
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = dir+"/"+base+".opt"
        alsoptfile = dir+"/"+base+".als.opt"
        
        # Frame specific parameters
        GA = self.gain
        RD = self.rdnoise
        if self.seeing is not None:
            FW = self.seeing / self.pixscale
        else:
            self.logger.info("No FWHM using CPFWHM")
            FW = self.cpfwhm / self.pixscale
        HI = self.saturate


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

        self.logger.info(" Created "+optfile+" and "+alsoptfile)
        
        
    # Make image ready for DAOPHOT
    def mkdaoim(self):
        flux,fhead = fits.getdata(self.fluxfile,header=True)
        wt,whead = fits.getdata(self.wtfile,header=True)
        mask,mhead = fits.getdata(self.maskfile,header=True)

        self.logger.info("-- Creating DAOPHOT-ready image --")

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

        #self.logger.info("Turning off the CP difference image masking flags")
        if self.plver > 0:      # CP data
            # V3.5.0 and on, Integer masks
            versnum = self.plver.split('.')
            if (versnum[0]>3) | ((versnum[0]==3) & (versnum[1]>=5)):
                bdpix = (mask == 7)
                nbdpix = np.sum(bdpix)
                if nbdpix > 0: mask[bdpix]=0

            # Pre-V3.5.0, Bitmasks
            else:
                bdpix = ( (mask & 2**7) == 2**7)
                nbdpix = np.sum(bdpix)                
                if nbdpix > 0: mask[bdpix]-=128   # clear 128

            self.logger.info("%d pixels cleared of difference image mask flag" % nbdpix)

        bdpix = (mask > 0.0)
        nbdpix = np.sum(bdpix)
        if nbdpix>0: flux[bdpix]=6e4
        self.logger.info("%d bad pixels masked" % nbdpix)

        fhead.append('GAIN',self.gain)
        fhead.append('RDNOISE',self.rdnoise)

        # Write new image
        self.logger.info("Wrote DAOPHOT-ready image to "+self.daofile)
        fits.writeto(self.daofile,flux,fhead,overwrite=True)


    # DAOPHOT detection
    #----------------------
    def daodetect(self,imfile=None):

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        scriptfile = base+".coo.sh"
        outfile = base+".coo"
        logfile = base+".coo.log"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)
        # Image file
        if imfile is None:
            if os.path.exists(base+".fits") is False:
                self.logger.warning("No image file input and "+base+".fits NOT found")
                return
            imfile = base+".fits"

        # Lines for the DAOPHOT script
        lines = "#!/bin/sh\n" \
                "daophot << END_DAOPHOT >> "+logfile+"\n" \
                "OPTIONS\n" \
                ""+optfile+"\n" \
                "\n" \
                "ATTACH "+imfile+"\n" \
                "FIND\n" \
                "1,1\n" \
                ""+outfile+"\n" \
                "y\n" \
                "EXIT\n" \
                "EXIT\n" \
                "END_DAOPHOT\n"
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

        # Run the script
        self.logger.info("-- Running DAOPHOT detection --")
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=False)
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("DAOPHOT detection failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # Get info from the logfile
            if os.path.exists(logfile) is True:
                f = open(logfile,'r')
                dlines = f.readlines()
                f.close()
                l1 = grep(dlines,"Sky mode and standard deviation")
                if len(l1)>0:
                    self.logger.info(l1[0].strip())   # clip \n
                    #l1 = l1[0]
                    #lo = l1.find("=")
                    #sky = np.array( l1[lo+1:].split('  '),dtype=float)
                l2 = grep(dlines,"Clipped mean and median")
                if len(l2)>0:
                    self.logger.info(l2[0].strip())
                    #l2 = l2[0]
                    #lo = l2.find("=")
                    #mnmed = np.array( l2[lo+2:].split(' '),dtype=float)
                # Number of sources
                l3 = grep(dlines," stars.")
                if len(l3)>0:
                    self.logger.info(l3[0].rstrip().strip())
        # Failure
        else:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        os.remove(scriptfile)

    # DAOPHOT aperture photometry
    #----------------------------
    def daoaperphot(self,coofile=None,apertures=None,imfile=None,outfile=None):

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        # Image file
        if imfile is None: imfile = base+".fits"
        imbase = os.path.basename(imfile)
        imbase = os.path.splitext(os.path.splitext(imbase)[0])[0]
        apfile = imbase+".apers"
        scriptfile = imbase+".ap.sh"
        logfile = imbase+".ap.log"
        # Output file
        if outfile is None:
            outfile = imbase+".ap"
        if os.path.exists(apfile): os.remove(apfile)
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)
        # Detection/coordinate file
        if coofile is None:
            if os.path.exists(base+".coo") is False:
                self.logger.warning("No detection/coordinate file input and "+base+".coo NOT found")
                return
            coofile = base+".coo"

        # Make apertures file
        if apertures is None:
            # The last two are inner and outer sky apertures
            #apertures = [3.0, 3.7965, 4.8046, 6.0803, 7.6947, 9.7377, 12.3232, 15.5952, 19.7360, \
            #             24.9762, 31.6077, 40.0000, 50.0000]
            apertures = [3.000, 6.0803, 9.7377, 15.5952, 19.7360, 40.0000, 50.0000]
        nap = len(apertures)
        if nap<3:
            self.logger.warning("Only "+str(nap)+" apertures input.  Need at least 3")
            return
        f = open(apfile,'w')
        for i in range(nap-2):
            # use hexidecimal for aperture id, 2 digits, first starts with A
            id = hex(160+i+1)
            id = id[2:].capitalize()
            f.write("%2s = %7.4f\n" % (id,apertures[i]))
        f.write("IS = %7.4f\n" % apertures[nap-2])
        f.write("OS = %7.4f\n" % apertures[nap-1])
        f.close()

        # Lines for the DAOPHOT script
        lines = "#!/bin/sh\n" \
                "daophot << END_DAOPHOT >> "+logfile+"\n" \
                "OPTIONS\n" \
                ""+optfile+"\n" \
                "\n" \
                "ATTACH "+imfile+"\n" \
                "PHOTOMETRY\n" \
                ""+apfile+"\n" \
                " \n" \
                ""+coofile+"\n" \
                ""+outfile+"\n" \
                "EXIT\n" \
                "EXIT\n" \
                "END_DAOPHOT\n"
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

        # If PSF file exists temporarily move it out of the way
        if os.path.exists(base+".psf"):
            self.logger.info(base+".psf exists.  Temporarily moving it out of the way to perform aperture photometry.")
            psftemp = base+".psf.bak"
            if os.path.exists(psftemp): os.remove(psftemp)
            os.rename(base+".psf",psftemp)
            movedpsf = True
        else:
            movedpsf = False

        # Run the script
        self.logger.info("-- Running DAOPHOT aperture photometry --")
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("DAOPHOT aperture photometry failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # Get info from the logfile
            if os.path.exists(logfile):
                f = open(logfile,'r')
                plines = f.readlines()
                f.close()
                l1 = grep(plines,"Estimated magnitude limit")
                if len(l1)>0:
                    l1 = l1[0]
                    l1 = l1[0:len(l1)-7]   # strip BELL at end \x07\n
                    lo = l1.find(":")
                    hi = l1.find("+-")
                    maglim = np.float(l1[lo+1:hi])
                    self.daomaglim = maglim
                    self.logger.info(l1.strip())   # clip leading/trailing whitespace
        # Failure
        else:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        os.remove(scriptfile)

        # Move PSF file back
        if movedpsf is True: os.rename(psftemp,base+".psf")


    # Pick PSF stars using DAOPHOT
    #-----------------------------
    def daopickpsf(self,catfile=None,maglim=None,nstars=100):

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        scriptfile = base+".pickpsf.sh"
        outfile = base+".lst"
        logfile = base+".lst.log"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)

        # What detection/coordinate file
        if catfile is None:
            if os.path.exists(base+".ap") is False:
                self.logger.warning("No catalog file input and "+base+".ap NOT found")
                return
            catfile = base+".ap"

        # Magnitude limit
        if maglim is None:
            if self.daomaglim is None:
                self.logger.warning("No magnitude input and DAOMAGLIMIT not set yet")
                return
            maglim = self.daomaglim-1.0

        # Lines for the DAOPHOT script
        lines = "#!/bin/sh\n" \
                "daophot << END_DAOPHOT >> "+logfile+"\n" \
                "OPTIONS\n" \
                ""+optfile+"\n" \
                "\n" \
                "ATTACH "+base+".fits\n" \
                "PICKPSF\n" \
                ""+catfile+"\n" \
                ""+str(nstars)+","+str(maglim)+"\n" \
                ""+outfile+"\n" \
                "EXIT\n" \
                "EXIT\n" \
                "END_DAOPHOT\n"
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

        # Run the script
        self.logger.info("-- Running DAOPHOT PICKPSF -- ")
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("DAOPHOT PICKPSF failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # Get info from the logfile
            if os.path.exists(logfile):
                f = open(logfile,'r')
                plines = f.readlines()
                f.close()
                l1 = grep(plines,"suitable candidates were found.")
                if len(l1)>0:
                    self.logger.info(l1[0].strip())   # clip \n
        # Failure
        else:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        os.remove(scriptfile)


    # Run DAOPHOT PSF
    #-------------------
    def daopsf(self,listfile=None,apfile=None,imfile=None,verbose=False):

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        scriptfile = base+".psf.sh"
        outfile = base+".psf"
        logfile = base+".psf.log"
        neifile = base+".nei"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)
        if os.path.exists(neifile): os.remove(neifile)

        # Aperture photometry file
        if apfile is None:
            if os.path.exists(base+".ap") is False:
                self.logger.warning("No aperture photometry file input and "+base+".ap NOT found")
                return
            apfile = base+".ap"
        # List file
        if listfile is None:
            if os.path.exists(base+".lst") is False:
                self.logger.warning("No PSF candidates list input and "+base+".lst NOT found")
                return
            listfile = base+".lst"
        # Image file
        if imfile is None: imfile = base+".fits"

        self.logger.info("-- Running DAOPHOT PSF -- ")

        # Lines for the DAOPHOT script
        lines = "#!/bin/sh\n" \
                "daophot << END_DAOPHOT >> "+logfile+"\n" \
                "OPTIONS\n" \
                ""+optfile+"\n" \
                "\n" \
                "ATTACH "+imfile+"\n" \
                "PSF\n" \
                ""+apfile+"\n" \
                ""+listfile+"\n" \
                ""+outfile+"\n" \
                "\n" \
                "EXIT\n" \
                "EXIT\n" \
                "END_DAOPHOT\n"
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

        # Run the script
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("DAOPHOT PSF failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # Get info from the logfile
            if os.path.exists(logfile):
                f = open(logfile,'r')
                plines = f.readlines()
                f.close()
                # Get parameter errors
                l1 = grep(plines,"Chi    Parameters",index=True)
                l2 = grep(plines,"Profile errors",index=True)
                l3 = grep(plines,"File with PSF stars and neighbors",index=True)
                if len(l1)>0:
                    parlines = plines[l1[0]+1:l2[0]-1]
                    pararr, parchi = parsepars(parlines)
                    minchi = np.min(parchi)
                    self.logger.info("  Chi = "+str(minchi))
                # Get profile errors
                if len(l2)>0:
                    proflines = plines[l2[0]+1:l3[0]-1]
                    if verbose: self.logger.info(" ".join(proflines))
                    profs = parseprofs(proflines)
                else:
                    self.logger.error("No DAOPHOT profile errors found in logfile")
                    raise
        # Failure
        else:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        os.remove(scriptfile)

        return pararr, parchi, profs


    # Subtract neighbors of PSF stars
    #--------------------------------
    def subpsfnei(self,lstfile=None,photfile=None,outfile=None,psffile=None,imfile=None):

        # Set up filenames, make sure the output ones don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        scriptfile = base+".subnei.sh"
        logfile = base+".subnei.log"
        nstfile = base+".nst"
        grpfile = base+".grp"
        if imfile is None: imfile = base+".fits"
        if outfile is None: outfile = base+"a.fits"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)
        if os.path.exists(nstfile): os.remove(nstfile)
        if os.path.exists(grpfile): os.remove(grpfile)
        # Check that image file exists
        if os.path.exists(imfile) is False:
            self.logger.error(imfile+" NOT Found")
            return
        # Phot file (normally .nei file)
        if photfile is None: 
            if os.path.exists(base+".nei") is False:
                self.logger.warning("Photometry file not input and "+base+".nei NOT found")
                return
            photfile = base+".nei"
        # List file
        if lstfile is None: 
            if os.path.exists(base+".lst") is False:
                self.logger.warning("List file not input and "+base+".lst NOT found")
                return
            lstfile = base+".lst"
        # PSF file
        if psffile is None: 
            if os.path.exists(base+".psf") is False:
                self.logger.warning("PSF file not input and "+base+".psf NOT found")
                return
            psffile = base+".psf"

        self.logger.info("-- Subtracting PSF stars neighbors -- ")

        # Lines for the DAOPHOT script
        lines = "#!/bin/sh\n" \
                "daophot << END_DAOPHOT >> "+logfile+"\n" \
                "OPTIONS\n" \
                ""+optfile+"\n" \
                "\n" \
                "ATTACH "+imfile+"\n" \
                "GROUP\n" \
                ""+photfile+"\n" \
                ""+psffile+"\n" \
                "5.\n" \
                ""+grpfile+"\n" \
                "NSTAR\n" \
                ""+psffile+"\n" \
                ""+grpfile+"\n" \
                ""+nstfile+"\n" \
                "SUBSTAR\n" \
                ""+psffile+"\n" \
                ""+nstfile+"\n" \
                "y\n" \
                ""+lstfile+"\n" \
                ""+outfile+"\n" \
                "\n" \
                "EXIT\n" \
                "END_DAOPHOT\n"
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("daophot.opt") is False: shutil.copyfile(base+".opt","daophot.opt")

        # Run the script
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=True)
            if retcode < 0:
                self.logger.error("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.error("PSF star neighbor subtracting failed:"+str(e))
            self.logger.error(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is False:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        #os.remove(scriptfile)


    # Create DAOPHOT PSF
    #-------------------
    def createpsf(self,listfile=None,apfile=None,doiter=True,maxiter=5,minstars=6,subneighbors=True,verbose=False):

        self.logger.info("-- Creating PSF Iteratively --")

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".opt"
        scriptfile = base+".psf.sh"
        outfile = base+".psf"
        logfile = base+".psf.log"
        neifile = base+".nei"
        imfile = base+".fits"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)
        if os.path.exists(neifile): os.remove(neifile)

        # Aperture photometry file
        if apfile is None:
            if os.path.exists(base+".ap") is False:
                self.logger.warning("No aperture photometry file input and "+base+".ap NOT found")
                return
            apfile = base+".ap"
        # List file
        if listfile is None:
            if os.path.exists(base+".lst") is False:
                self.logger.warning("No PSF candidates list input and "+base+".lst NOT found")
                return
            listfile = base+".lst"
        # Working list file
        wlistfile = listfile+"1"
        if os.path.exists(wlistfile): os.remove(wlistfile)
        shutil.copy(listfile,wlistfile)

        # Make copy of original list
        if os.path.exists(listfile+".orig"): os.remove(listfile+".orig")
        shutil.copy(listfile,listfile+".orig")

        # Iterate
        #---------
        if doiter is False: maxiter=1
        iter = 1
        endflag = 0
        while (endflag==0):
            self.logger.info("Iter = "+str(iter))

            # Run DAOPSF
            try:
                pararr, parchi, profs = self.daopsf(wlistfile,apfile)
            except:
                self.logger.error("Failure in DAOPSF")
                raise

            # Check for bad stars
            nstars = len(profs)
            bdstars = (profs['flag'] != '')
            nbdstars = np.sum(bdstars)
            self.logger.info("  "+str(nbdstars)+" stars with flags")
            # Delete stars with flags
            if (nbdstars>0) & (nstars>minstars):
                f = open(wlistfile,'r')
                listlines = f.readlines()
                f.close()
                # Read the list
                lstcat = daoread(wlistfile)
                # Match up with the stars we are deleting
                mid, ind1, ind2 = np.intersect1d(profs[bdstars]['id'],lstcat['id'],return_indices=True)
                # Remove the lines from listlines
                newlistlines = remove_indices(listlines,ind2+3)
                # Write new list
                os.remove(wlistfile)
                f = open(wlistfile,'w')
                f.writelines(newlistlines)
                f.close()
                self.logger.info("  Removing IDs="+str(" ".join(profs[bdstars]['id'].astype(str))))
                self.logger.info("  "+str(nbdstars)+" bad stars removed. "+str(nstars-nbdstars)+" PSF stars left")

            # Should we end
            if (iter==maxiter) | (nbdstars==0) | (nstars<=minstars): endflag=1
            iter = iter+1

        # Subtract PSF star neighbors
        if subneighbors:
            subfile = base+"a.fits"
            #try:
            self.subpsfnei(wlistfile,neifile,subfile,imfile=imfile,psffile=outfile)
            #except:
            #    self.logger.error("Subtracting neighbors failed.  Keeping original PSF file")
            # Check that the subtracted image exist and rerun DAOPSF
            if os.path.exists(subfile):
                # Final run of DAOPSF
                self.logger.info("Final DAOPDF run")
                try:
                    pararr, parchi, profs = self.daopsf(wlistfile,apfile)
                except:
                    self.logger.error("Failure in DAOPSF")
                    raise
                # Get aperture photometry for PSF stars from subtracted image
                self.logger.info("Getting aperture photometry for PSF stars")
                apertures = [3.0, 3.7965, 4.8046, 6.0803, 7.6947, 9.7377, 12.3232, 15.5952, 19.7360, \
                             24.9762, 31.6077, 40.0000, 50.0000]
                self.daoaperphot(wlistfile,apertures,subfile)

        # Copy working list to final list
        if os.path.exists(listfile): os.remove(listfile)
        shutil.move(wlistfile,listfile)
        self.logger.info("Final list of PSF stars in "+listfile+".  Original list in "+listfile+".orig")

        
        
    # Run ALLSTAR
    #-------------
    def runallstar(self,psffile=None,apfile=None,subfile=None):

        # Set up filenames, make sure they don't exist
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        optfile = base+".als.opt"
        scriptfile = base+".als.sh"
        if apfile is None:
            if os.path.exists(base+".ap") is False:
                self.logger.warning("apfile not input and "+base+".ap NOT found")
                return
            apfile = base+".ap"
        if psffile is None:
            if os.path.exists(base+".psf") is False:
                self.logger.warning("psffile not input and "+base+".psf NOT found")
                return
            psffile = base+".psf"
        if subfile is None: subfile = base+"s.fits"
        outfile = base+".als"
        logfile = base+".als.log"
        if os.path.exists(outfile): os.remove(outfile)
        if os.path.exists(subfile): os.remove(subfile)
        if os.path.exists(logfile): os.remove(logfile)
        if os.path.exists(scriptfile): os.remove(scriptfile)

        # Load the option file lines
        f = open(optfile,'r')
        optlines = f.readlines()
        f.close()

        # Lines for the DAOPHOT ALLSTAR script
        lines = ["#!/bin/sh\n",
                 "allstar << END_ALLSTAR >> "+logfile+"\n"]
        lines += optlines
        lines += ["\n",
                  base+".fits\n",
                  psffile+"\n",
                  apfile+"\n",
                  outfile+"\n",
                  subfile+"\n",
                  "EXIT\n",
                  "EXIT\n",
                  "END_ALLSTAR\n"]
        # Write the script
        f = open(scriptfile,'w')
        f.writelines(lines)
        f.close()
        os.chmod(scriptfile,0775)

        # Copy option file to daophot.opt
        if os.path.exists("allstar.opt") is False: shutil.copyfile(optfile,"allstar.opt")

        # Run the script
        self.logger.info("-- Running ALLSTAR --")
        try:
            retcode = subprocess.call(["./"+scriptfile],stderr=subprocess.STDOUT,shell=False)
            if retcode < 0:
                self.logger.warning("Child was terminated by signal"+str(-retcode))
            else:
                pass
        except OSError as e:
            self.logger.warning("ALLSTAR failed:"+str(e))
            self.logger.warning(e)
            raise

        # Check that the output file exists
        if os.path.exists(outfile) is True:
            # How many sources converged
            num = numlines(outfile)-3
            self.logger.info(str(num)+" stars converged")
            self.logger.info("Final catalog is "+outfile)
        # Failure
        else:
            self.logger.error("Output file "+outfile+" NOT Found")
            raise

        # Delete the script
        os.remove(scriptfile)


    # Combine SE + DAOPHOT/ALLSTAR info into ONE file

        
    # Run DAOGROW to calculate aperture corrections
    #----------------------------------------------
    def daogrow(self):
        pass

    # Process a single chip
    #----------------------
    def process(self):
        self.runsex()
        self.mkopt()
        self.mkdaoim()
        self.daodetect()
        self.daoaperphot()
        self.daopickpsf()
        self.runpsf()
        self.runallstar()
