#!/usr/bin/env python

# Incomming command: 
#python /home/x25h971/nsc/instcal/v4/calibrate_exposure_tempest.py --expdir untarred_expdir? --refcat refcatfile --redo --gsynthphot True --psf True

# -------
# Imports
# -------

from argparse import ArgumentParser
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from glob import glob
import healpy as hp
import logging
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy import stats
import shutil
import socket
import subprocess
import time
#from utils_tempest import *

import warnings
warnings.resetwarnings()
warnings.filterwarnings('ignore',category=UserWarning,append=True)


from nsc import utils,query,modelmag
from calibrate_tempest import *

# ---------
# Main Code
# ---------

if __name__=="__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Calibrate the exposures from one HEALPix of the NSC')
    parser.add_argument('--version', type=str, nargs=1, help='Version number')
    parser.add_argument('--pix', type=str, nargs=1, help='HEALPix index')
    parser.add_argument('--nside', type=str, nargs=1, help='HEALPix NSIDE (ring ordering, please!)')
    parser.add_argument('--expdir',type=str,nargs=1,help='Uncompressed name of directory of the exposure you want to calibrate')
    parser.add_argument('--refcatfile', type=str, nargs=1, help="File with reference source catalog for this exposure's HEALPix")
    parser.add_argument('--gsynth', action='store_true',help='Photometric Calibration with Gaia Syntetic Photometry?')
    parser.add_argument('--psf',action='store_true',help='Astrometric Calibration with PSF coordinates?')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposures that were previously processed')
    args = parser.parse_args()

    # Start time, get hostname
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = (args.version[0])         # NSC version
    pix = int(args.pix[0])              # index of HEALPix to calibrate
    nside = int(args.nside[0])          # HEALPix NSIDE (ring ordering)
    #maxexp = int(args.maxexp[0])     # maximum number of exposures to maintain running at any time
    expdir = args.expdir[0]
    refcatfile = args.refcatfile[0]
    gsynth = args.gsynth                    # Photometric calibration with Gaia Synthetic Photometry?
    psf = args.psf                          # Astrometric calibration with PSF coordinates?
    redo = args.redo                         # if called, redo = True
    sleep_time = 10                          # seconds to sleep between checking job status
    #inputlist = args.list                    # list of exposures to analyze
    #if inputlist is not None:
    #    inputlist = inputlist[0]

    # Establish necessary directories
    # -tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"   # location of operations
    inddir = basedir+"healpix_indicators/"+str(nside)+"_"+str(int(pix))+"/"  # location of indicator
    #makedir(inddir)
    indfile = inddir+str(expdir.split("/")[-1].split(".tar.gz")[0])+".txt"

    # Log File
    #---------
    # Create Log file name;
    # format is nsc_combine_main.DATETIME.log
    ltime = time.localtime()
    # time.struct_time(tm_year=2019, tm_mon=7, tm_mday=22,
    #                  tm_hour=0, tm_min=30, tm_sec=20,
    #                  tm_wday=0, tm_yday=203, tm_isdst=1)
    smonth = str(ltime[1])
    if ltime[1]<10: smonth = '0'+smonth
    sday = str(ltime[2])
    if ltime[2]<10: sday = '0'+sday
    syear = str(ltime[0])[2:]
    shour = str(ltime[3])
    if ltime[3]<10: shour='0'+shour
    sminute = str(ltime[4])
    if ltime[4]<10: sminute='0'+sminute
    ssecond = str(int(ltime[5]))
    if ltime[5]<10: ssecond='0'+ssecond
    logtime = smonth+sday+syear+shour+sminute+ssecond
    logfile = basedir+'logs/calibrate_exposure_'+str(pix)+"_"+str(expdir.split("/")[-1])+'.'+logtime+'.log'

    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    rootLogger.setLevel(logging.NOTSET)

    rootLogger.info("Calibrating Exposure "+str(expdir))
    rootLogger.info("host = "+host)
    rootLogger.info("pix = "+str(pix))
    rootLogger.info("nside = "+str(nside))
    rootLogger.info("version = "+version)
    rootLogger.info("redo = "+str(redo))


    # Check for expdir status (tarred, zipped, extant?)
    # Unzip and untar the exposure repo if necessary
    expdir_base = expdir.split("/")[-1]
    expdir_nightfolder = "/".join(expdir.split("/")[:-1])
    runexp = False
    if os.path.exists(expdir+".tar.gz"): # if exp.tar.gz,
        os.chdir(expdir_nightfolder)
        rootLogger.info("tar -xzf "+expdir_base+".tar.gz") # -C /")
        os.system("tar -xzf "+expdir_base+".tar.gz") # -C /")
        os.chdir(basedir)
        runexp = True
    elif os.path.exists(expdir+".tar"): # elif exp.tar, (already unzipped)
        os.chdir(expdir_nightfolder)
        rootLogger.info("tar -xf "+expdir_base+".tar") # -C /")
        os.system("tar -xf "+expdir_base+".tar") # -C /")
        os.chdir(basedir)
        runexp = True
    elif os.path.exists(expdir): # elif exp (aleady unzipped and untarred)
        rootLogger.info("exposure "+expdir+" already untarred and unzipped")
        runexp = True
    else: # else if not there,
        rootLogger.info("exposure directory does not exist!")
        runexp = False

    # Calibrate the exposure
    if runexp:
        c = calibrate(expdir,refcatfile=refcatfile,logfilename=logfile,redo=redo,gsynthphot=gsynth,psf=psf,logger=rootLogger)
        # check to see if exposure has been calibrated
        if c:
            rootLogger.info("exposure "+str(expdir)+"complete!")
            # if so, write an indicator file
            with open(indfile,"w") as fl:
                fl.writelines("done! "+str(logtime))
                fl.close()
            exp_nightdir = "/".join(expdir.split("/")[:-1])
            exp_base = expdir.split("/")[-1]
            os.chdir(exp_nightdir)
            shutil.rmtree(exp_base)
            os.chdir(basedir)
        # re-tar and zip file
        #os.system()
    else: rootLogger.info("exposure can't be processed!")
