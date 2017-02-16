#!/usr/bin/env python

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

if __name__ == "__main__":


# Run SExtractor on one FULL DECam stacked image

    dir = "/datalab/users/dnidever/decamcatalog/instcal/"
    tmproot = "/data0/dnidever/decamcatalog/instcal/tmp/"

    t0 = time.time()

    print sys.argv

    # Not enough inputs
    n = len(sys.argv)
    if n < 4:
        print "Syntax - nsc_instcal.py fluxfile wtfile maskfile"
        sys.exit()

    #fluxfile = "/net/mss1/archive/pipeline/Q20160307/DECALS/201209/c4d_120914_040214_osi_r_a1.fits.fz"
    #wtfile = "/net/mss1/archive/pipeline/Q20160307/DECALS/201209/c4d_120914_040214_osw_r_a1.fits.fz"
    #maskfile = "/net/mss1/archive/pipeline/Q20160307/DECALS/201209/c4d_120914_040214_osd_r_a1.fits.fz"

    # SMASH 20160101
    # 178 c4d_160102_063033_ood_g_v1.fits.fz  00507860  Field33  2016-01-02T06:25:16.518568  g  267.000  4:57:28.15  -84:18:05.4
    # 179 c4d_160102_063033_ooi_g_v1.fits.fz  00507860  Field33  2016-01-02T06:25:16.518568  g  267.000  4:57:28.15  -84:18:05.4
    # 180 c4d_160102_063033_oow_g_v1.fits.fz  00507860  Field33  2016-01-02T06:25:16.518568  g  267.000  4:57:28.15  -84:18:05.4
    #fluxfile = "/net/mss1/archive/pipeline/Q20160107/DEC15B/20160101/c4d_160102_063033_ooi_g_v1.fits.fz"
    #wtfile = "/net/mss1/archive/pipeline/Q20160107/DEC15B/20160101/c4d_160102_063033_oow_g_v1.fits.fz"
    #maskfile = "/net/mss1/archive/pipeline/Q20160107/DEC15B/20160101/c4d_160102_063033_ood_g_v1.fits.fz"


    # File names
    fluxfile = sys.argv[1]
    wtfile = sys.argv[2]
    maskfile = sys.argv[3]
    # Check that the files exist
    if os.path.exists(fluxfile) == False:
        print fluxfile, "file NOT FOUND"
        sys.exit()
    if os.path.exists(wtfile) == False:
        print wtfile, "file NOT FOUND"
        sys.exit()
    if os.path.exists(maskfile) == False:
        print maskile, "file NOT FOUND"
        sys.exit()


    base = os.path.basename(fluxfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]

    # 1) Prepare temporary directory
    #---------------------------------
    #print "Step #1: Preparing temporary directory"
    tmpcntr = 1L
    tmpdir = tmproot+base+"."+str(tmpcntr)
    while (os.path.exists(tmpdir)):
        tmpcntr = tmpcntr+1
        tmpdir = tmproot+base+"."+str(tmpcntr)
        if tmpcntr > 10:
            print "Temporary Directory counter getting too high. Exiting"
            sys.exit()
    os.mkdir(tmpdir)
    origdir = os.getcwd()
    os.chdir(tmpdir)

    # Set up logging to screen and logfile
    #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    logfile = tmpdir+"/"+base+".log"
    #fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, fileName))
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    rootLogger.setLevel(logging.NOTSET)

    # Setting up logging
    ##logging.getLogger().addHandler(logging.StreamHandler())
    #logger = logging.getLogger('nsc_fullstack')
    #logfile = tmpdir+"/"+base+".log"
    #hdlr = logging.FileHandler(logfile)
    #formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    #hdlr.setFormatter(formatter)
    #logger.addHandler(hdlr)
    ##logger.setLevel(logging.WARNING)

    rootLogger.info("Running SExtractor on "+base)
    rootLogger.info("  Temporary directory is: "+tmpdir)

    # 2) Copy over images from zeus1:/mss
    #-------------------------------------
    rootLogger.info("Step #2: Copying InstCal images from mass store archive")
    shutil.copyfile(fluxfile,tmpdir+"/"+os.path.basename(fluxfile))
    rootLogger.info("  "+fluxfile)
    os.symlink(os.path.basename(fluxfile),"bigflux.fits.fz")
    shutil.copyfile(wtfile,tmpdir+"/"+os.path.basename(wtfile))
    rootLogger.info("  "+wtfile)
    os.symlink(os.path.basename(wtfile),"bigwt.fits.fz")
    shutil.copyfile(maskfile,tmpdir+"/"+os.path.basename(maskfile))
    rootLogger.info("  "+maskfile)
    os.symlink(os.path.basename(maskfile),"bigmask.fits.fz")

    #file = "/net/mss1/archive/pipeline/Q20160307/DECALS/201209/c4d_120914_040214_ose_r_a1.fits.fz"
    #filebase = os.path.basename(file)
    # this is trimming the fits.fz extension
    #shutil.copyfile(file,tmpdir+"/"+filebase)
    #print tmpdir+"/"+filebase
    ##cmd = 'rsync -av --password-file=/home/dnidever/.password dnidever@zeus1.sdm.noao.edu:'+file+' .'
    #cmd = 'scp -B dnidever@zeus1.sdm.noao.edu:'+file+' .'
    #spawn,cmd,out,errout
    # maybe use the DL SIA service

    # Get number of extensions
    hdulist = fits.open("bigflux.fits.fz")
    nhdu = len(hdulist)
    hdulist.close()

    # 3) Run Sextractor on all subimages
    rootLogger.info("Step #3: Running SExtractor on all subimages")
    head0 = fits.getheader("bigflux.fits.fz",0)
    #fwhmpix = head0.get("fwhm")
    #if fwhmpix is None:
    #    fwhm = 1.5
    #else:
    #    fwhm = fwhmpix*0.27
    dateobs = head0.get("DATE-OBS")
    night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]


    # Make final output directory
    if not os.path.exists("/datalab/users/dnidever/decamcatalog/instcal/"+night):
        os.mkdir("/datalab/users/dnidever/decamcatalog/instcal/"+night)
    if not os.path.exists("/datalab/users/dnidever/decamcatalog/instcal/"+night+"/"+base):
        os.mkdir("/datalab/users/dnidever/decamcatalog/instcal/"+night+"/"+base)
        rootLogger.info("  Making output directory: /datalab/users/dnidever/decamcatalog/instcal/"+night+"/"+base)

    # LOOP through the HDUs/chips
    for i in xrange(1,nhdu):
        rootLogger.info(" Processing subimage "+str(i))
        try:
            flux,fhead = fits.getdata("bigflux.fits.fz",i,header=True)
            wt,whead = fits.getdata("bigwt.fits.fz",i,header=True)
            mask,mhead = fits.getdata("bigmask.fits.fz",i,header=True)
        except:
            rootLogger.info("No extension "+str(i))

        # Use CCDNUM
        ccdnum = fhead['ccdnum']
        rootLogger.info("  CCDNUM = "+str(ccdnum))

        # FWHM values are ONLY in the extension headers
        fwhmpix = fhead.get("fwhm")
        if fwhmpix is None:
            fwhm = 1.5
        else:
            fwhm = fwhmpix*0.27

        # 3a) Make subimages for flux, weight, mask

        if os.path.exists("flux.fits"):
            os.remove("flux.fits")
        fits.writeto("flux.fits",flux,header=fhead,output_verify='warn')

        # I think we might need to set the weight map pixels 
        #  that are "bad" to 5e30
        # I think wt=0 for bad pixels
        if os.path.exists("wt.fits"):
            os.remove("wt.fits")
        fits.writeto("wt.fits",wt,header=whead,output_verify='warn')

        if os.path.exists("mask.fits"):
            os.remove("mask.fits")
        fits.writeto("mask.fits",mask,header=mhead,output_verify='warn')


        # 3b) Make SExtractor config files
        # Copy the default files
        shutil.copyfile(dir+"default.conv",tmpdir+"/default.conv")
        shutil.copyfile(dir+"default.nnw",tmpdir+"/default.nnw")
        shutil.copyfile(dir+"default.param",tmpdir+"/default.param")

        # fix gain, rdnoise in header
        # Read in the SExtractor config file and change gain
        f = open('/datalab/users/dnidever/decamcatalog/instcal/default.config', 'r') # 'r' = read
        lines = f.readlines()
        f.close()

        # Things to change
        # SATUR_LEVEL     59000.00         # level (in ADUs) at which arises saturation
        # GAIN            43.52             # detector gain in e-/ADU.
        # PIXEL_SCALE     0.26282645             # size of pixel in arcsec (0=use FITS WCS info).
        # SEEING_FWHM     1.46920            # stellar FWHM in arcsec
        # WEIGHT_IMAGE  F4-00507860_01_comb.mask.fits

        cnt = 0L
        for l in lines:
            # SATUR_LEVEL
            m = re.search('^SATUR_LEVEL',l)
            if m != None:
                lines[cnt] = "SATUR_LEVEL     59000.00         # level (in ADUs) at which arises saturation\n"
                #print "SATUR line ", cnt
            # Gain
            m = re.search('^GAIN',l)
            if m != None:
                lines[cnt] = "GAIN            1.0             # detector gain in e-/ADU.\n"
                #print "GAIN line ", cnt
            # PIXEL_SCALE
            m = re.search('^PIXEL_SCALE',l)
            if m != None:
                lines[cnt] = "PIXEL_SCALE     0.27             # size of pixel in arcsec (0=use FITS WCS info). \n"
                #print "PIXEL_SCALE line ", cnt
            # SEEING_FWHM
            m = re.search('^SEEING_FWHM',l)
            if m != None:
                lines[cnt] = "SEEING_FWHM     "+str(fwhm)+"            # stellar FWHM in arcsec\n"
                #print "FWHM line ", cnt
            # WEIGHT_IMAGE
            m = re.search('^WEIGHT_IMAGE',l)
            if m != None:
                lines[cnt] = "WEIGHT_IMAGE  wt.fits    # Weight image name.\n"
                #print "WEIGHT line ", cnt
            cnt = cnt+1
        # Write out the new config file
        if os.path.exists("default.config"):
            os.remove("default.config")
        fo = open('default.config', 'w')
        fo.writelines(lines)
        fo.close()

        # 3c) Run SExtractor
        #p = subprocess.Popen('sex', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rootLogger.info("  Running SExtractor")
        if os.path.exists("cat.fits"):
            os.remove("cat.fits")

        try:
            # Save the SExtractor info to a logfile
            slogfile = tmpdir+"/"+base+"_"+str(ccdnum)+".sex.log"
            sf = open(slogfile,'w')
            retcode = subprocess.call(["sex","flux.fits","-c","default.config"],stdout=sf,stderr=subprocess.STDOUT)
            sf.close()
            sf = open(slogfile,'r')
            slines = sf.readlines()
            sf.close()
            rootLogger.info(slines)
            #retcode = call("mycmd" + " myarg", shell=True)
            if retcode < 0:
                #rootLogger.info(>>sys.stderr, "Child was terminated by signal", -retcode)
                #rootLogger.info(sys.stderr+"Child was terminated by signal"+str(-retcode))
                rootLogger.info("Child was terminated by signal"+str(-retcode))
            else:
                #rootLogger.info(>>sys.stderr, "Child returned", retcode)
                #rootLogger.info(sys.stderr+"Child returned"+str(retcode))
                rootLogger.info("Child returned"+str(retcode))
        except OSError as e:
            #rootLogger.info(>>sys.stderr, "SExtractor Execution failed:", e)
            #rootLogger.info(sys.stderr+"SExtractor Execution failed:"+str(e))            
            rootLogger.info("SExtractor Execution failed:"+str(e))


        # Catch the output and put it in a logfile

        # 3d) Load the catalog (and logfile) and write final output file
        # Move the file to final locaation
        if os.path.exists("cat.fits"):
            outcatfile = "/datalab/users/dnidever/decamcatalog/instcal/"+night+"/"+base+"/"+base+"_"+str(ccdnum)+".fits"
            # Clobber if it already exists
            if os.path.exists(outcatfile):
                os.remove(outcatfile)
                rootLogger.info("  Copying final catalog to "+outcatfile)
            # Copy to final directory
            shutil.copyfile("cat.fits",outcatfile)
        else:
            rootLogger.info("  No output catalog")

        # 4) Delete temporary directory/files
        rootLogger.info("  Deleting subimages")
        if os.path.exists("check.fits"):
            os.remove("check.fits")
        if os.path.exists("flux.fits"):
            os.remove("flux.fits")
        if os.path.exists("wt.fits"):
            os.remove("wt.fits")
        if os.path.exists("mask.fits"):
            os.remove("mask.fits")

    # Move the log file
    #os.rename(logfile,"/datalab/users/dnidever/decamcatalog/"+night+"/"+base+"/"+base+".log")
    # The above rename gave an error on gp09, OSError: [Errno 18] Invalid cross-device link
    shutil.move(logfile,"/datalab/users/dnidever/decamcatalog/instcal/"+night+"/"+base+"/"+base+".log")

    # Delete temporary files and directory
    #rootLogger.info("Deleting all temporary files")
    tmpfiles = glob.glob("*")
    for f in tmpfiles:
        os.remove(f)
    os.rmdir(tmpdir)

    # CD back to original directory
    os.chdir(origdir)

    rootLogger.info(str(time.time()-t0)+" seconds")
