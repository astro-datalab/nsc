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
import time
import shutil
import re
import subprocess
import glob
import logging
import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve
import astropy.stats


# Represent an exposure to process
class Exposure:

    def __init__(self,fluxfile,wtfile,maskfile):
        self.fluxfile = fluxfile
        self.wtfile = wtfile
        self.maskfile = maskfile
        base = os.path.basename(fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        self.base = base
        self.logfile = base+".log"
        self.logger = None
        self.chip = None

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
        
        #rootLogger.info("Running SExtractor on "+base+" on host="+host)
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
        # Add logger information
        self.chip.logger = self.logger
        
    # setup
    def setup():
        pass
    
    # teardown
    def teardown():
        pass
        
        
# Represent a single chip of an exposure
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
        self.seeing = None
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
            print("Cannot get RDNOISE, no header yet")
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
        print('No RDNOISE found')
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
            print("Cannot get CCDNUM, no header yet")
        # Get ccdnum from the header
        # We have this key, set _rndoise and return
        if 'CCDNUM' in self.header.keys():
            self._ccdnum = self.header['CCDNUM']
            return self._ccdnum
        print('No CCDNUM found')
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
            print("Cannot get SATURATE, no header yet")
        # Get saturate from the header
        # We have this key, set _saturate and return
        if 'SATURATE' in self.header.keys():
            self._saturate = self.header['SATURATE']
            return self._saturate
        print('No SATURATE found')
        return None
    
    @property
    def wcs(self):
        # We have it already, just return it
        if self._wcs is not None:
            return self._wcs
        # Can't get wcs, no header yet
        if self.header is None:
            print("Cannot get WCS, no header yet")
        try:
            self._wcs = WCS(self.header)
            return self._wcs
        except:
            print("Problem with WCS")
            return None
            
    @property
    def exptime(self):
        # We have it already, just return it
        if self._exptime is not None:
            return self._exptime
        # Can't get exptime, no header yet
        if self.header is None:
            print("Cannot get EXPTIME, no header yet")
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
            print("Cannot get INSTRUMENT, no header yet")
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
            print("Cannot get PLVER, no header yet")
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
        # Can't get plver, no header yet
        if self.header is None:
            print("Cannot get CPFWHM, no header yet")
        # FWHM values are ONLY in the extension headers
        cpfwhm_map = { 'c4d': 1.5 if self.header.get('FWHM') is None else self.header.get('FWHM')*0.27, 
                       'k4m': 1.5 if self.header.get('SEEING1') is None else self.header.get('SEEING1'),
                       'ksb': 1.5 if self.header.get('SEEING1') is None else self.header.get('SEEING1') }
        cpfwhm = cpfwhm_map[self.instrument]
        self._cpfwhm = cpfwhm
        return self._cpfwhm

    
    # Make DAOPHOT option files:
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

        
    # Pick PSF stars
    def pickpsf():
        pass
        
    # Create DAOPHOT PSF
    def runpsf():
        pass
        
    # Subtract neighbors of PSF stars
    def subnei():
        pass
        
    # Run ALLSTAR
    def runallstar():
        pass
        
    # PSF star aperture photometry for calculating the aperture correction
    def brighstaraperphot():
        pass
        
    # Run DAOGROW to calculate aperture corrections
    def daogrow():
        pass

