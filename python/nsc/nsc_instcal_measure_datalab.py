#!/usr/bin/env python

#AUTHORS: David Nidever (original author)
#         david.nidever@montana.edi
#         Katie Fasbender (adapted for analysis on MSU Tempest Research Cluster)
#         katiefasbender@montana.edu
#
# NSC_INSTCAL_MEAS.PY -- Run SExtractor and DAOPHOT on an exposure from the
# NOIRLab Astro Data Archive (NOIRLab Source Catalog measurements procedure)



from __future__ import print_function

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, Column,vstack
import time
import shutil
import re
import subprocess
import logging
import socket
from scipy.ndimage.filters import convolve
import astropy.stats
import struct
from utils import *
from phot import *
import glob

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# Get NSC directories
def getnscdirs(version=None):
    # Host
    hostname = socket.gethostname()
    host = hostname.split('.')[0].strip()
    print("host = ",host)
    # Version
    verdir = ""
    if version is not None:
       verdir = version if version.endswith('/') else version+"/"
    # on gp07 use
    if (host == "gp09") | (host == "gp08") | (host == "gp07") | (host == "thing") | (host == "hulk"):
        basedir = "/net/dl2/kfas/nsc/instcal/"+verdir
        tmproot = "/net/dl2/kfas/nsc/instcal/"+verdir+"tmp/"
    # on tempest use
    else:
        basedir = "/home/x25h971/nsc/instcal/"+verdir
        tmproot = "/home/x25h971/nsc/instcal/"+verdir+"tmp/"
    return basedir,tmproot


# Class to represent an exposure to process
class Exposure:

    # Initialize Exposure object
    def __init__(self,fluxfile,wtfile,maskfile,nscversion="v4"): #"t3a"):
        # Check that the files exist
        if os.path.exists(fluxfile) is False:
            print(fluxfile+" NOT found")
            return
        if os.path.exists(wtfile) is False:
            print(wtfile+" NOT found")
            return
        if os.path.exists(maskfile) is False:
            print(maskfile+" NOT found")
            return
        # Setting up the object properties
        self.origfluxfile = fluxfile
        self.origwtfile = wtfile
        self.origmaskfile = maskfile
        self.fluxfile = None      # working files in temp dir
        self.wtfile = None        # working files in temp dir
        self.maskfile = None      # working files in temp dir
        base = os.path.basename(fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        self.base = base
        self.nscversion = nscversion
        self.logfile = base+".log"
        self.logger = None
        self.origdir = None
        self.wdir = None     # the temporary working directory
        self.outdir = None
        self.chip = None


        # Get instrument
        head0 = fits.getheader(fluxfile,0)
        if head0["DTINSTRU"] == 'mosaic3':
            self.instrument = 'k4m'
        elif head0["DTINSTRU"] == '90prime':
            self.instrument = 'ksb'
        elif head0["DTINSTRU"] == 'decam':
            self.instrument = 'c4d'
        else:
            print("Cannot determine instrument type")
            return
        # Get number of extensions
        hdulist = fits.open(fluxfile)
        nhdu = len(hdulist)
        hdulist.close()
        self.nexten = nhdu
        # Get night
        dateobs = head0.get("DATE-OBS")
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        self.night = night
        # Output directory
        basedir,tmpdir = getnscdirs(nscversion)
        self.outdir = basedir+self.instrument+"/"+self.night+"/"+self.base+"/"

        
    # Setup
    def setup(self):
        #print("nscversion = ",version)
        basedir,tmproot = getnscdirs(version)
        # Prepare temporary directory
        tmpcntr = 1#L 
        tmpdir = tmproot+self.base+"."+str(tmpcntr)
        print("temp dir = ",tmpdir)
        while (os.path.exists(tmpdir)):
            tmpcntr = tmpcntr+1
            tmpdir = tmproot+self.base+"."+str(tmpcntr)
            if tmpcntr > 20:
                print("Temporary Directory counter getting too high. Exiting")
                sys.exit()
        os.mkdir(tmpdir)
        origdir = os.getcwd()
        self.origdir = origdir
        os.chdir(tmpdir)
        self.wdir = tmpdir

        # Set up logging to screen and logfile
        logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
        rootLogger = logging.getLogger()
        # file handler
        fileHandler = logging.FileHandler(self.logfile)
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)
        # console/screen handler
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        rootLogger.addHandler(consoleHandler)
        rootLogger.setLevel(logging.NOTSET)
        self.logger = rootLogger
        self.logger.info("Setting up in temporary directory "+tmpdir)
        self.logger.info("Starting logfile at "+self.logfile)

        # Copy over images from zeus1:/mss
        fluxfile = "bigflux.fits.fz"
        wtfile = "bigwt.fits.fz"
        maskfile = "bigmask.fits.fz"
        #self.logger.info("Copying InstCal images from mass store archive")
        self.logger.info("Copying InstCal images downloaded from Astro Data Archive")
        shutil.copyfile(self.origfluxfile,tmpdir+"/"+os.path.basename(self.origfluxfile))
        self.logger.info("  "+self.origfluxfile)
        if (os.path.basename(self.origfluxfile) != fluxfile):
            os.symlink(os.path.basename(self.origfluxfile),fluxfile)
        shutil.copyfile(self.origwtfile,tmpdir+"/"+os.path.basename(self.origwtfile))
        self.logger.info("  "+self.origwtfile)
        if (os.path.basename(self.origwtfile) != wtfile):
            os.symlink(os.path.basename(self.origwtfile),wtfile)
        shutil.copyfile(self.origmaskfile,tmpdir+"/"+os.path.basename(self.origmaskfile))
        self.logger.info("  "+self.origmaskfile)
        if (os.path.basename(self.origmaskfile) != maskfile):
            os.symlink(os.path.basename(self.origmaskfile),maskfile)

        # Set local working filenames
        self.fluxfile = fluxfile
        self.wtfile = wtfile
        self.maskfile = maskfile
        
        # Make final output directory
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)   # will make multiple levels of directories if necessary
            self.logger.info("Making output directory: "+self.outdir)


    # Load chip
    def loadchip(self,extension,fluxfile="flux.fits",wtfile="wt.fits",maskfile="mask.fits"):
        # Load the data
        self.logger.info(" Loading chip "+str(extension))
        # Check that the working files set by "setup"
        if (self.fluxfile is None) | (self.wtfile is None) | (self.maskfile is None):
            self.logger.warning("Local working filenames not set.  Make sure to run setup() first")
            return
        try:
            flux,fhead = fits.getdata(self.fluxfile,extension,header=True)
            fhead0 = fits.getheader(self.fluxfile,0)  # add PDU info
            fhead.extend(fhead0,unique=True)
            wt,whead = fits.getdata(self.wtfile,extension,header=True)
            mask,mhead = fits.getdata(self.maskfile,extension,header=True)
        except:
            self.logger.error("No extension "+str(extension))
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
        self.chip = Chip(fluxfile,wtfile,maskfile,self.base)
        self.chip.bigextension = extension
        self.chip.nscversion = self.nscversion
        self.chip.outdir = self.outdir
        # Add logger information
        self.chip.logger = self.logger


    # Process all chips
    def process(self):
        self.logger.info("-------------------------------------------------")
        self.logger.info("Processing ALL extension images")
        self.logger.info("-------------------------------------------------")

        # LOOP through the HDUs/chips
        #----------------------------
        #for i in [int(sys.argv[6])]: #ktedit:createpsf_test,  only analyze 1 chip!
        for i in range(1,self.nexten):
            t0 = time.time()
            self.logger.info(" ")
            self.logger.info("=== Processing subimage "+str(i)+" ===")
            # Load the chip
            self.loadchip(i)
            self.logger.info("CCDNUM = "+str(self.chip.ccdnum))
            # Process it
            nsrc=self.chip.process()
            # Clean up 
            self.chip.cleanup()
            t1_check=time.time()
            self.logger.info("dt = "+str(time.time()-t0)+" seconds")
            if 2==1:
                chiptimes = Table.read(basedir+'lists/nsc_dr3_chiptimes.fits')
                chiptimes.add_row([(str(self.fluxfile).strip().split('/')[-1]).split('.')[0],i,int(self.chip.ccdnum),int(nsrc),int(t1_check-t0)])
                chiptimes = Table(np.unique(chiptimes))
                chiptimes.write(basedir+'lists/nsc_dr3_chiptimes.fits',overwrite=True)
    
    # Teardown
    def teardown(self):
        # Delete files and temporary directory
        self.logger.info("Deleting files and temporary directory.")
        # Move the final log file
        shutil.move(self.logfile,self.outdir+self.base+".log")
        # Delete temporary files and directory
        tmpfiles = glob.glob("*")
        for f in tmpfiles: os.remove(f)
        os.rmdir(self.wdir)
        # CD back to original directory
        os.chdir(self.origdir)
        

    # RUN all steps to process this exposure
    def run(self):
        self.setup()
        self.process()
        self.teardown()

        
# Class to represent a single chip of an exposure
class Chip:

    def __init__(self,fluxfile,wtfile,maskfile,bigbase):
        self.fluxfile = fluxfile
        self.wtfile = wtfile
        self.maskfile = maskfile
        self.bigbase = bigbase
        self.bigextension = None
        base = os.path.basename(fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        self.dir = os.path.abspath(os.path.dirname(fluxfile))
        self.base = base
        self.meta = makemeta(header=fits.getheader(fluxfile,0))
        self.sexfile = self.dir+"/"+self.base+"_sex.fits"
        self.daofile = self.dir+"/"+self.base+"_dao.fits"
        self.sexcatfile = None
        self.sexcat = None
        self.seeing = None
        self.apcorr = None
        #-----------------------------------------------------------#ktedit:sex2 T
        # For the second run of SExtractor on the ALLSTAR PSF-subtracted file
        self.allsubfile=self.dir+"/"+self.base+"_daos.fits"
        self.smeta = None #in runsex(), define self.meta=makemeta(header=fits.getheader(allsubfile,0)) 
        self.sexcatfile2 = None
        self.sexcat2 = None
        #-----------------------------------------------------------#ktedit:sex2 B      
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
        self._daomaglim = None    # set by daoaperphot()
        self._sexmaglim = None    # set by runsex()
        #self._sexmaglim2 = None   #ktedit:sex2; set by runsex() when run on ALLSTAR PSF-subtracted image
        #self._daomaglim2 = None   #ktedit:sex2; set by daoaperphot() when run on 2nd SExtractor cat
        self.sexiter = 1          #ktedit:sex2; to keep track of which SExtractor run we're on
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
        if self.meta is None:
            self.logger.warning("Cannot get RDNOISE, no header yet")
            return None
        # Check DECam style rdnoise
        if "RDNOISEA" in self.meta.keys():
            rdnoisea = self.meta["RDNOISEA"]
            rdnoiseb = self.meta["RDNOISEB"]
            rdnoise = (rdnoisea+rdnoiseb)*0.5
            self._rdnoise = rdnoise
            return self._rdnoise
        # Get rdnoise from the header
        for name in ['RDNOISE','READNOIS','ENOISE']:
            # We have this key, set _rndoise and return
            if name in self.meta.keys():
                self._rdnoise = self.meta[name]
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
            gain = gainmap[self.instrument](self.meta)
        except:
            gainmap_avg = { 'c4d': 3.9845419, 'k4m': 1.8575, 'ksb': 1.4}
            gain = gainmap_avg[self.instrument]
        self._gain = gain
        return self._gain
            
        ## Can't get gain, no header yet
        #if self.meta is None:
        #    print("Cannot get GAIN, no header yet")
        ## Get rdnoise from the header
        #for name in ['GAIN','EGAIN']:
        #    # We have this key, set _gain and return
        #    if self.meta.has_key(name):
        #        self._gain = self.meta[name]
        #        return self._gain
        #print('No GAIN found')
        #return None
            
    @property
    def ccdnum(self):
        # We have it already, just return it
        if self._ccdnum is not None:
            return self._ccdnum
        # Can't get ccdnum, no header yet
        if self.meta is None:
            self.logger.warning("Cannot get CCDNUM, no header yet")
            return None
        # Get ccdnum from the header
        # We have this key, set _rndoise and return
        if 'CCDNUM' in self.meta.keys():
            self._ccdnum = self.meta['CCDNUM']
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
        if self.meta is None:
            self.logger.warning("Cannot get SATURATE, no header yet")
            return None
        # Get saturate from the header
        # We have this key, set _saturate and return
        if 'SATURATE' in self.meta.keys():
            self._saturate = self.meta['SATURATE']
            return self._saturate
        self.logger.warning('No SATURATE found')
        return None
    
    @property
    def wcs(self):
        # We have it already, just return it
        if self._wcs is not None:
            return self._wcs
        # Can't get wcs, no header yet
        if self.meta is None:
            self.logger.warning("Cannot get WCS, no header yet")
            return None
        try:
            self._wcs = WCS(self.meta)
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
        if self.meta is None:
            self.logger.warning("Cannot get EXPTIME, no header yet")
            return None
        # Get rdnoise from the header
        # We have this key, set _rndoise and return
        if 'EXPTIME' in self.meta.keys():
                self._exptime = self.meta['EXPTIME']
                return self._exptime
        print('No EXPTIME found')
        return None

    @property
    def instrument(self):
        # We have it already, just return it
        if self._instrument is not None:
            return self._instrument
        # Can't get instrument, no header yet
        if self.meta is None:
            self.logger.warning("Cannot get INSTRUMENT, no header yet")
            return None
        # instrument, c4d, k4m or ksb
        # DTINSTRU = 'mosaic3 '
        # DTTELESC = 'kp4m    '
        # Bok 90Prime data has
        if self.meta.get("DTINSTRU") == 'mosaic3':
            self._instrument = 'k4m'
            return self._instrument
        elif self.meta.get("DTINSTRU") == '90prime':
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
        if self.meta is None:
            self.logger.warning("Cannot get PLVER, no header yet")
            return None
        plver = self.meta.get('PLVER')
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
        if self.meta is None:
            self.logger.warning("Cannot get CPFWHM, no header yet")
            return None
        # FWHM values are ONLY in the extension headers
        cpfwhm_map = { 'c4d': 1.5 if self.meta.get('FWHM') is None else self.meta.get('FWHM')*0.27, 
                       'k4m': 1.5 if self.meta.get('SEEING1') is None else self.meta.get('SEEING1'),
                       'ksb': 1.5 if self.meta.get('SEEING1') is None else self.meta.get('SEEING1') }
        cpfwhm = cpfwhm_map[self.instrument]
        self._cpfwhm = cpfwhm
        return self._cpfwhm

    @property
    def maglim(self):
        # We have it already, just return it
        if self._daomaglim is not None:
            return self._daomaglim
        if self._sexmaglim is not None:
            return self._sexmaglim
        self.logger.warning('Maglim not set yet')
        return None


    # Write SE catalog in DAO format
    #-------------------------------
    #def sextodao(self,cat=None,outfile=None,format="coo"):
    def sextodao(self,cat=None,outfile=None,format="coo",meta=None): #ktedit:sex2
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        if outfile is None: outfile=daobase+".coo"
        if cat is None: cat=self.sexcat
        if meta is None: meta=self.meta                      #ktedit:sex2
        else: offs=0                                         #ktedit:sex2
        #sextodao(self.sexcat,self.meta,outfile=outfile,format=format,logger=self.logger)
        sextodao(cat,meta,outfile=outfile,format=format,logger=self.logger) #ktedit:sex2

    # Run Source Extractor
    #---------------------
    #def runsex(self,outfile=None):
    def runsex(self,dthresh=1.1,outfile=None): #ktedit:sex2
        #--------------------------------------------------------------------------------------------ktedit:sex2 T
        # if allsub=False, run SExtractor on the fluxfile.  Otherwise, run on the ALLSTAR PSF-subtracted image.
        if self.sexiter==1: 
            infile=self.fluxfile
            meta=self.meta
            sexcatfile="flux_sex.cat.fits"
            offset=0
        else:
            daobase = os.path.basename(self.daofile)
            daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]             
            infile = daobase+str(self.sexiter-1)+"s.fits"
            self.smeta=makemeta(header=fits.getheader(infile,0))  #should this be self.smeta?  probably.
            meta=self.smeta
            sexcatfile = "flux_sex"+str(self.sexiter)+".cat.fits"
            if self.sexcat is not None: offset=int(self.sexcat['NUMBER'][-1]) #ktedit:sex2
        #--------------------------------------------------------------------------------------------ktedit:sex2 B
        basedir, tmpdir = getnscdirs(self.nscversion)
        configdir = basedir+"config/"
        #sexcat, maglim = runsex(self.fluxfile,self.wtfile,self.maskfile,self.meta,sexcatfile,configdir,logger=self.logger)
        sexcat, maglim = runsex(infile,self.wtfile,self.maskfile,meta,sexcatfile,configdir,offset=offset,sexiter=self.sexiter,dthresh=dthresh,logger=self.logger) #ktedit:sex2
        #--------------------------------------------------------------------------------------------ktedit:sex2 T
        sexcat.add_column(np.repeat(self.sexiter,len(sexcat)),name="NDET_ITER") # keep track of what SExtractor iteration each source is from
        sexcat.add_column(np.zeros(len(sexcat)),name="REPEAT")                  # keep track of sources that were detected in multiple iterations
        # Rules of "REPEAT" column: 
        #  0 = source only detected once
        #  1 = source detected in multiple iterations (all iterations but last), will be removed from sexcat
        #  2 = source detected in multiple iterations (last iteration source was detected in)
        temp_sexcat=Table.read(sexcatfile,2)
        temp_sexcat.add_column(np.repeat(self.sexiter,len(temp_sexcat)),name="NDET_ITER")
        temp_sexcat.add_column(np.zeros(len(temp_sexcat)),name="REPEAT")
        temp_sexcat.write(sexcatfile,overwrite=True)
        # --If first SExtractor iteration, define cat
        if self.sexiter==1:
            self.sexcatfile = sexcatfile
            self.sexcat = sexcat
            self._sexmaglim = maglim
            # Set the FWHM as well
            fwhm = sexfwhm(sexcat,logger=self.logger)
            self.meta['FWHM'] = fwhm
        # --If 2nd+ SExtractor iteration, compare sources with
        # those from previous iteration and combine catalogs 
        else: 
            sexcat = vstack([self.sexcat,sexcat])
            # lastsex -> newsexcat, restsex -> prevsexcat   
            newsexcat = sexcat[sexcat['NDET_ITER']==self.sexiter]
            prevsexcat = sexcat[sexcat['NDET_ITER']==(self.sexiter-1)]
            for newsource in newsexcat:
                dpix = 2
                prevsexcat_close = prevsexcat[(prevsexcat['X_IMAGE']<(newsource['X_IMAGE']+dpix)) & (prevsexcat['X_IMAGE']>(newsource['X_IMAGE']-dpix)) & (prevsexcat['Y_IMAGE']<(newsource['Y_IMAGE']+dpix)) & (prevsexcat['Y_IMAGE']>(newsource['Y_IMAGE']-dpix))]
                if len(prevsexcat_close)>0:
                    for oldsource in prevsexcat_close:
                        d_btwn_centers = np.sqrt((newsource['X_IMAGE']-oldsource['X_IMAGE'])**2+(newsource['Y_IMAGE']-oldsource['Y_IMAGE'])**2)
                        if d_btwn_centers <= dpix:
                            old_repeat_index = int(np.where(sexcat['NUMBER']==oldsource['NUMBER'])[0])
                            new_repeat_index = int(np.where(sexcat['NUMBER']==newsource['NUMBER'])[0])
                            sexcat[old_repeat_index]['REPEAT'] = 1
                            sexcat[new_repeat_index]['REPEAT'] = 2
            self.sexcat = sexcat[sexcat['REPEAT']!=1]

        #--------------------------------------------------------------------------------------------ktedit:sex2 B


    # Determine FWHM using SE catalog
    #--------------------------------
    def sexfwhm(self):
        self.seeing = sexfwhm(self.sexcat)
        return self.seeing

    # Pick PSF candidates using SE catalog
    #-------------------------------------
    def sexpickpsf(self,nstars=100):
        base = os.path.basename(self.sexfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        fwhm = self.sexfwhm() if self.seeing is None else self.seeing
        psfcat = sexpickpsf(self.sexcat,fwhm,self.meta,base+".lst",nstars=nstars,logger=self.logger)

    # Make DAOPHOT option files
    #--------------------------
    #def mkopt(self,**kwargs):
    def mkopt(self):
        base = os.path.basename(self.daofile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        #mkopt(base,self.meta,logger=self.logger,**kwargs)
        mkopt(base,self.meta,logger=self.logger)
        
    # Make image ready for DAOPHOT
    def mkdaoim(self):
        mkdaoim(self.fluxfile,self.wtfile,self.maskfile,self.meta,self.daofile,logger=self.logger)

    # DAOPHOT detection
    #----------------------
    def daofind(self):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        cat = daofind(self.daofile,outfile=daobase+".coo",logger=self.logger)

    # DAOPHOT aperture photometry
    #----------------------------
    def daoaperphot(self):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        #----------------------------------------------------------------#ktedit:sex2 T
        imfile=self.daofile
        if self.sexiter==1:
            coofile=daobase+".coo"
            outfile=daobase+".ap"
        else:
            coofile=daobase+str(self.sexiter)+".coo"
            outfile=daobase+str(self.sexiter)+".ap"
        #----------------------------------------------------------------#ktedit:sex2 B
        #apcat, maglim = daoaperphot(self.daofile,daobase+".coo",outfile=daobase+".ap",logger=self.logger)
        apcat, maglim = daoaperphot(imfile,coofile,outfile=outfile,optfile=daobase+".opt",logger=self.logger) #ktedit:sex2
        #self._daomaglim = maglim
        if self.sexiter==1: self._daomaglim = maglim #ktedit:sex2
        #else: self._daomaglim2 = maglim            #ktedit:sex2

    # Pick PSF stars using DAOPHOT
    #-----------------------------
    def daopickpsf(self,maglim=None,nstars=100):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        if maglim is None: maglim=self.maglim
        psfcat = daopickpsf(self.daofile,daobase+".ap",maglim,daobase+".lst",nstars,logger=self.logger)

    # Run DAOPHOT PSF
    #-------------------
    def daopsf(self,verbose=False):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        psfcat = daopsf(self.daofile,daobase+".lst",outfile=daobase+".psf",verbose=verbose,logger=self.logger)

    # Subtract neighbors of PSF stars
    #--------------------------------
    def subpsfnei(self):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        psfcat = subpsfnei(self.daofile,daobase+".lst",daobase+".nei",daobase+"a.fits",logger=self.logger)

    # Create DAOPHOT PSF
    #-------------------
    def createpsf(self,listfile=None,apfile=None,doiter=True,maxiter=5,minstars=6,subneighbors=True,verbose=False):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        subit=createpsf(daobase+".fits",daobase+".ap",daobase+".lst",meta=self.meta,logger=self.logger)
        self.subiter=subit
        
    # Run ALLSTAR
    #-------------
    def allstar(self,psffile=None,apfile=None,subfile=None):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        #-------------------------------------------------------------------#ktedit:sex2 T
        imfile = daobase+".fits"
        meta=self.meta
        subfile = daobase+str(self.sexiter)+"s.fits"
        if self.sexiter==1: 
            apfile = daobase+".ap"
            outfile = daobase+".als"
        else: 
#            subfile = daobase+str(self.sexiter)+"s.fits"
            apfile = daobase+str(self.sexiter)+".ap"
            outfile = daobase+str(self.sexiter)+".als"
        #-------------------------------------------------------------------#ktedit:sex2 B
        #alscat = allstar(daobase+".fits",daobase+".psf",daobase+".ap",outfile=daobase+".als",meta=self.meta,logger=self.logger)
        alscat = allstar(imfile,daobase+".psf",apfile=apfile,subfile=subfile,outfile=outfile,optfile=daobase+".als.opt",meta=meta,logger=self.logger) #ktedit:sex2


    # Combine total + new SExtractor & ALLSTAR catalog files #ktedit:sex2; this function is new
    #-----------------------------------------------------
    def combine_cats(self,type):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]

        if type=="sexcat":
            file1 = daobase+".coo"                                      # total SExcat stored in this file
            file2 = daobase+str(self.sexiter)+".coo"                    # new SExcat stored in this file

            # must also combie .fits format of SE cat
            file3 = "flux_sex.cat.fits"
            file4 = "flux_sex"+str(self.sexiter)+".cat.fits"
            sexcat1 = Table.read(file3)
            sexcat2 = Table.read(file4)
            sexcat_total = vstack([sexcat1,sexcat2])
            sexcat_total.write(file3,overwrite=True)

        elif type=="alscat":
            file1 = daobase+".als"                                      # total ALLSTAR cat stored in this file
            file2 = daobase+str(self.sexiter)+".als"                    # new ALLSTAR cat stored in this file

        cat1 = readlines(file1)
        cat2 = readlines(file2)
        combined_cat = cat1+cat2[3:]                                    # combine the catalogs 
        writelines(file1,combined_cat,overwrite=True)


    # Get aperture correction
    #------------------------
    def getapcor(self):
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        apcorr = apcor(daobase+"a.fits",daobase+".lst",daobase+".psf",self.meta,optfile=daobase+'.opt',alsoptfile=daobase+".als.opt",logger=self.logger)
        self.apcorr = apcorr
        self.meta['apcor'] = (apcorr,"Aperture correction in mags")

    # Combine SE and DAOPHOT catalogs
    #--------------------------------
    def finalcat(self,outfile=None,both=True,sexdetect=True):
        # both       Only keep sources that have BOTH SE and ALLSTAR information
        # sexdetect  SE catalog was used for DAOPHOT detection list

        self.logger.info("--  Creating final combined catalog --")

        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        if outfile is None: outfile=self.base+".cat.fits"

        # Check that we have the SE and ALS information
        if (self.sexcat is None) | (os.path.exists(daobase+".als") is None):
            self.logger.warning("SE catalog or ALS catalog NOT found")
            return

        # Load ALS catalog
        als = Table(daoread(daobase+".als")) 
        nals = len(als)
        # Apply aperture correction
        if self.apcorr is None:
            self.logger.error("No aperture correction available")
            return
        als['MAG'] -= self.apcorr

        # Just add columns to the SE catalog
        ncat = len(self.sexcat)
        newcat = self.sexcat.copy()
        alsnames = ['X','Y','MAG','ERR','SKY','ITER','CHI','SHARP']
        newnames = ['XPSF','YPSF','MAGPSF','ERRPSF','SKY','ITER','CHI','SHARP','RAPSF','DECPSF']
        newtypes = ['float64','float64','float','float','float','float','float','float','float64','float64']
        nan = float('nan')
        newvals = [nan, nan, nan, nan ,nan, nan, nan, nan, nan, nan]
        # DAOPHOT detection list used, need ALS ID
        if not sexdetect:
            alsnames = ['ID']+alsnames
            newnames = ['ALSID']+newnames
            newtypes = ['int32']+newtypes
            newvals = [-1]+newvals
        newcols = []
        for n,t,v in zip(newnames,newtypes,newvals):
            col = Column(name=n,length=ncat,dtype=t)
            col[:] = v
            newcols.append(col)
        newcat.add_columns(newcols)
        #self.logger.info("allstar catalog has "+str(nals)+" lines and "+str(len(np.unique(als['ID'])))+" unique ids") #ktedit:sex2
        #self.logger.info("newcat has "+str(len(newcat))+" lines and "+str(len(np.unique(newcat['NUMBER'])))+" unique starids") #ktedit:sex2
        # Match up with IDs if SE list used by DAOPHOT
        if sexdetect:
            mid, ind1, ind2 = np.intersect1d(newcat["NUMBER"],als["ID"],return_indices=True)
            for id1,id2 in zip(newnames,alsnames):
                newcat[id1][ind1] = als[id2][ind2]
            # Only keep sources that have SE+ALLSTAR information
            #  trim out ones that don't have ALS
            if (both is True) & (nals<ncat): newcat = newcat[ind1]
            #self.logger.info("newcat has "+str(len(newcat))+" lines") #ktedit:sex2

        # Match up with coordinates, DAOPHOT detection list used
        else:
            print("Need to match up with coordinates")

            # Only keep sources that have SE+ALLSTAR information
            #  trim out ones that don't have ALS
            if (both is True) & (nals<ncat): newcat = newcat[ind1]

        # Add RA, DEC
        r,d = self.wcs.all_pix2world(newcat["XPSF"],newcat["YPSF"],1)
        newcat['RAPSF'] = r
        newcat['DECPSF'] = d        
        #self.logger.info("length of final catalog = "+str(len(newcat))) #ktedit:sex2

        # Write to file
        self.logger.info("Final catalog = "+outfile)
        fits.writeto(outfile,None,self.meta,overwrite=True)  # meta in PDU header
        #  append the table in extension 1
        hdulist = fits.open(outfile)
        hdu = fits.table_to_hdu(newcat)
        hdulist.append(hdu)
        hdulist.writeto(outfile,overwrite=True)
        hdulist.close()
        #newcat.write(outfile,overwrite=True)
        #fits.append(outfile,0,self.meta)  # meta is header of 2nd extension


    # Process a single chip
    #----------------------
    def process(self):

        # Set up SE iteration-------------------------------------------------------------------------------
        sexiter_endflag=0                       
        while sexiter_endflag==0:

            # For every iteration, run Source Extractor-----------------------------------------------------
            self.logger.info("-- SExtractor run "+str(self.sexiter)+" --")

            # determine DETECT_THRESH value from SE iteration
            if self.sexiter==1: sex_dt=1.7
            else: sex_dt=1.1

            self.runsex(dthresh=sex_dt)

            # get the info for this iteration's catalog
            nowcat = self.sexcat[self.sexcat['NDET_ITER']==self.sexiter]  # cat for current SE iteration
            nowcat = nowcat[nowcat['REPEAT']==0]                          # select only newly detected sources
            self.logger.info(str(len(nowcat))+" new sources detected")
            nowcat_sn5 = nowcat[1/nowcat['MAGERR_AUTO']>=5]               # select only the entries with SN>=5
            if self.sexiter==1: 
                ogcat = nowcat                                            #to check the current cat against the first one
                ogcat_sn5 = nowcat_sn5
            #self.logger.info("sex"+str(self.sexiter)+" stats: mean"+str(np.mean(1/nowcat['MAGERR_AUTO']))+" median"+str(np.median(1/nowcat['MAGERR_AUTO'])))

            # Perform aperture photometry and PSF fitting with DAOPHOT--------------------------------------
            self.logger.info("-- Getting ready to run DAOPHOT --")

            # for first iteration only, make DAO-ready files
            if self.sexiter==1:
                self.mkopt()
                self.mkdaoim()
            
            # convert SE cat to DAO format
            #self.daodetect()
            # Create DAOPHOT-style coo file
            # Need to use SE positions
            if self.sexiter==1: sdao_ofile="flux_dao.coo"          #ktedit:sex2; select aperphot output filename
            else: sdao_ofile="flux_dao"+str(self.sexiter)+".coo"   #ktedit:sex2
            self.sextodao(outfile=sdao_ofile)

            self.daoaperphot()

            # for first iteration only, fit PSF 
            if self.sexiter==1:
                self.daopickpsf()   
                self.createpsf()

            # combine SE cats, run ALLSTAR, combine ALLSTAR cats--------------------------------------------
            if self.sexiter>1: self.combine_cats(type="sexcat")           
            self.allstar()
            if self.sexiter>1: self.combine_cats(type="alscat")                      

            # Check to see if we've run enough SExtractor iterations----------------------------------------
            # Requirements to end:
            # - at least 2 iterations, no more than 4
            # - median S/N of latest SE cat <=5
            # - #sources for which S/N>=5 latest cat is less than 25% #sources (also S/N>=5) in first cat
            if self.sexiter>2 and (self.sexiter==4 or (np.median(1/nowcat['MAGERR_AUTO'])<=5)  or len(nowcat_sn5)<(.25*len(ogcat_sn5))): sexiter_endflag=1
            self.sexiter+=1                                                          

        # Get aperture correction, create final cat from SE + ALLSTAR cats----------------------------------
        self.getapcor()
        self.finalcat()
        return(len(self.sexcat))


        # David's notes:------------------------------------------------------------------------------------

        # Do I need to rerun daoaperphot to get aperture
        # photometry at the FINAL allstar positions??
        
        # Is there a way to reduce the number of iterations needed to create the PSF?
        # what do the ?, * mean anyway?
        # maybe just remove the worse 10% of stars or something

        # Put all of the daophot-running into separate function (maybe separate module)
        # same for sextractor

        # Maybe make my own xmatch function that does one-to-one matching

    # Clean up the files
    #--------------------
    def cleanup(self):
        self.logger.info("Copying final files to output directory "+self.outdir)
        base = os.path.basename(self.fluxfile)
        base = os.path.splitext(os.path.splitext(base)[0])[0]
        daobase = os.path.basename(self.daofile)
        daobase = os.path.splitext(os.path.splitext(daobase)[0])[0]
        # Copy the files we want to keep
        # final combined catalog, logs
        outcatfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".fits"
        if os.path.exists(outcatfile): os.remove(outcatfile)
        shutil.copyfile("flux.cat.fits",outcatfile)
        # Copy DAOPHOT opt files
        outoptfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".opt"
        if os.path.exists(outoptfile): os.remove(outoptfile)
        shutil.copyfile(daobase+".opt",outoptfile)
        outalsoptfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".als.opt"
        if os.path.exists(outalsoptfile): os.remove(outalsoptfile)
        shutil.copyfile(daobase+".als.opt",outalsoptfile)
        # Copy DAOPHOT PSF star list
        outlstfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".psf.lst"
        if os.path.exists(outlstfile): os.remove(outlstfile)
        shutil.copyfile(daobase+".lst",outlstfile)
        # Copy DAOPHOT PSF file
        outpsffile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".psf"
        if os.path.exists(outpsffile): os.remove(outpsffile)
        shutil.copyfile(daobase+".psf",outpsffile)
        # Copy DAOPHOT .apers file??
        # Copy neighbor-subtracted images to output dir #ktedit
        #for i in range(0,int(self.subiter-1)):
        #    if int(i)==int(self.subiter-2): nsub=""
        #    else:nsub=str(i+1)
        #    outnsubfile=self.outdir+self.bigbase+"_"+str(self.ccdnum)+"_"+nsub+"a.fits"
        #    nsubfile=daobase+nsub+"a.fits"
        #    if os.path.exists(outnsubfile): os.remove(outnsubfile)
        #    shutil.copyfile(nsubfile,outnsubfile)
        # Copy daophot-ready image to output dir
        #outdimfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+"daoim.fits"
        #if os.path.exists(outdimfile): os.remove(outdimfile)
        #shutil.copyfile(daobase+".fits",outdimfile)
        # copy Allstar PSF subtracted files to output dir #ktedit
        #for i in range(1,int(self.sexiter)):
        #    outsubfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+"_"+str(i)+"s.fits"
        #    if os.path.exists(outsubfile): os.remove(outsubfile)
        #    shutil.copyfile(daobase+str(i)+"s.fits",outsubfile)
        # Copy SE config file
        outconfigfile = self.outdir+self.bigbase+"_"+str(self.ccdnum)+".sex.config"
        if os.path.exists(outconfigfile): os.remove(outconfigfile)
        shutil.copyfile("default.config",outconfigfile)
        # Copy SE segmentation files       #ktedit:sex2
        #for i in range(1,int(self.sexiter)):
        #    outsegfile=self.outdir+self.bigbase+"_"+str(self.ccdnum)+"_"+str(i)+"seg.fits"
        #    if os.path.exists(outsegfile): os.remove(outsegfile)
        #    shutil.copyfile("seg_"+str(i)+".fits",outsegfile)

        # Combine all the log files
        logfiles = glob.glob(base+"*.log")
        loglines = []
        for logfil in logfiles:
            loglines += ["==> "+logfil+" <==\n"]
            f = open(logfil,'r')
            lines = f.readlines()
            f.close()
            loglines += lines
            loglines += ["\n"]
        f = open(base+".logs","w")
        f.writelines("".join(loglines))
        f.close()
        outlogfile =  self.outdir+self.bigbase+"_"+str(self.ccdnum)+".logs"
        if os.path.exists(outlogfile): os.remove(outlogfile)
        shutil.copyfile(base+".logs",outlogfile)

        # Delete temporary directory/files
        self.logger.info("  Cleaning up")
        files1 = glob.glob("flux*")
        files2 = glob.glob("default*")        
        files = files1+files2+["flux.fits","wt.fits","mask.fits","daophot.opt","allstar.opt"]
        for f in files:
            if os.path.exists(f): os.remove(f)


# Main command-line program
if __name__ == "__main__":

    # Version
    verdir = ""
    if len(sys.argv) > 4:
    #if len(sys.argv) > 6: #ktedit:createpsf_test
       version = sys.argv[4]
       print("version = ",version)
       verdir = version if version.endswith('/') else version+"/"
    else: version = None
    # Get NSC directories
    basedir, tmpdir = getnscdirs(version)

    # Make sure the directories exist
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    t0 = time.time()
    print(sys.argv)

    # Not enough inputs
    n = len(sys.argv)
    if n < 4:
        print("Syntax - nsc_instcal_sexdaophot.py fluxfile wtfile maskfile version")
        sys.exit()

    # File names
    filedir ="/home/x25h971/nsc/instcal/v4/exposures/"
#    rawname = sys.argv[1]
    fluxfile = sys.argv[1]
    wtfile = sys.argv[2]
    maskfile = sys.argv[3]
    # Check that the files exist
    if os.path.exists(fluxfile) is False:
        print(fluxfile+" file NOT FOUND")
        sys.exit()
    if os.path.exists(wtfile) is False:
        print(wtfile+" file NOT FOUND")
        sys.exit()
    if os.path.exists(maskfile) is False:
        print(maskile+" file NOT FOUND")
        sys.exit()

    # Create the Exposure object
    exp = Exposure(fluxfile,wtfile,maskfile,nscversion=version)
    # Run
    exp.run()

    print("Total time = "+str(time.time()-t0)+" seconds")
