#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import time
import shutil
import re
import subprocess
import glob
import logging
import socket
from argparse import ArgumentParser
import nsc

if __name__ == "__main__":

    # Calibrate one full DECam/Mosaic3/Bok InstCal image

    parser = ArgumentParser(description='Calibrate one NSC InstCal image')
    parser.add_argument('expdir', type=str, nargs=1, help='Exposure directory')
    #parser.add_argument('inpref', type=str, nargs=1, help='Input reference')
    parser.add_argument('--eqnfile', type=str, nargs=1, default='', help='Modelmag equations filename')
    parser.add_argument('--redo', action='store_true', help='Redo this image')
    parser.add_argument('--selfcal', action='store_true', help='Self calibration')
    parser.add_argument('--saveref', action='store_true', help='Save the reference catalog')
    parser.add_argument('--ncpu', type=int, nargs=1, help='Number of cpus to use')
    args = parser.parse_args()

    t0 = time.time()

    # File names
    expdir = args.expdir
    eqnfile = args.eqnfile
     
    # Calibrate catalogs for one exposure
    dldir,mssdir,localdir = nsc.utils.rootdirs()
     
    # Make sure the directory exists 
    if os.path.exists(expdir) == False:
        raise ValueError(expdir+' NOT FOUND')
     
    t00 = time.time() 
     
    base = os.path.basename(expdir) 
    
    # Start the logfile 
    #------------------ 
    # format is delvered_forcebrick.brick.host.DATETIME.log
    host = socket.gethostname()
    hostname = host.split('.')[0]
    logtime = datetime.now().strftime("%Y%m%d%H%M%S") 
    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger() 
    while logger.hasHandlers(): # some existing loggers, remove them   
        logger.removeHandler(logger.handlers[0]) 
    logger = logging.getLogger()
    logtime = datetime.now().strftime("%Y%m%d%H%M%S")
    logfile = expdir+'/'+base+'_calib.log'
    if os.path.exists(logfile): os.remove(logfile)
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.NOTSET)


    nsc.calibrate.calibrate(expdir,eqnfile=args.eqnfile,redo=args.redo,
                            selfcal=args.selfcal,saveref=args.saveref,ncpu=args.ncpu)
