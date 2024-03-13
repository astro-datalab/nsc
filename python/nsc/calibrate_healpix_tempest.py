#!/usr/bin/env python


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
import socket
import subprocess
import time
#from utils_tempest import *

import warnings
warnings.resetwarnings()
warnings.filterwarnings('ignore',category=UserWarning,append=True)


from nsc import utils,query,modelmag
from slurm_funcs import *
from calibrate_tempest import *

# ---------
# Main Code
# ---------

if __name__=="__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Calibrate the exposures from one HEALPix of the NSC')
    parser.add_argument('--version', type=str, nargs=1,help='Version number')
    parser.add_argument('--pix',type=str,nargs=1,help='The HEALPix index, NSISDE of your choice!') 
    parser.add_argument('--nside', type=str, nargs=1,default=32,help='HEALPix NSIDE (ring ordering, please!)')
    parser.add_argument('--maxexp',type=str,nargs=1,default=10,help='Maximum number of exposures to calibrate at once (per job)')
    parser.add_argument('--partition',type=str,nargs=1,help='Partition to stick to')
    parser.add_argument('--gsynth', action='store_true',help='Photometric Calibration with Gaia Syntetic Photometry?')
    parser.add_argument('--psf',action='store_true',help='Astrometric Calibration with PSF coordinates?')
    parser.add_argument('-rexp','--redoexp', action='store_true',help='Redo exposures that were previously processed')
    args = parser.parse_args()

    # Start time, get hostname
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = (args.version[0])         # NSC version
    pix = int(args.pix[0])              # index of HEALPix to calibrate
    nside = int(args.nside[0])          # HEALPix NSIDE (ring ordering)
    maxexp = int(args.maxexp[0])     # maximum number of exposures to maintain running at any time
    gsynth = args.gsynth             # Photometric calibration with Gaia Synthetic Photometry?
    psf = args.psf                   # Astrometric calibration with PSF coordinates?
    redo = args.redoexp              # if called, redo exposures previously processed
    sleep_time = 10                  # seconds to sleep between checking job status

    # Establish necessary directories
    # -tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"   # location of operations
    inddir = basedir+"healpix_indicators/"+str(nside)+"_"+str(int(pix))
    makedir(inddir)

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
    logfile = basedir+'logs/calibrate_healpix_'+str(pix)+'.'+logtime+'.log'

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

    rootLogger.info("Calibrating HEALPix "+str(pix)+" of NSIDE="+str(nside))
    rootLogger.info(str(maxexp)+" exposures to transfer/process at once")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("redo = "+str(redo))

    calibrate_healpix(pix,version,nside=nside,maxexp=maxexp,logtime=logtime,gsynth=gsynth,psf=psf,redo=redo,logger=rootLogger)
