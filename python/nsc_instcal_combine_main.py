#!/usr/bin/env python

import os
import sys
import shutil
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table
from dlnpyutils import utils as dln, coords, job_daemon as jd
from astropy import units as u
from astropy.coordinates import SkyCoord
#import subprocess
import time
from argparse import ArgumentParser
import socket
import logging
import healpy as hp

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC Instcal Catalogs.')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    parser.add_argument('--hosts', type=str, nargs=1, help='Delimited list of hosts')
    parser.add_argument('-nm','--nmulti', type=int, nargs=1, default=15, help='Number of jobs')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)
    redo = args.redo
    nmulti = dln.first_el(args.nmulti)
    if args.hosts is not None:
        hosts = args.hosts[0].split(',')
    else:
        hosts = host
    nside = 128
    radeg = 180 / np.pi

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        #basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        basedir = "/net/dl2/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        #basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        basedir = "/net/dl2/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    if not os.path.exists(tmpdir): os.mkdir(tmpdir)
    if not os.path.exists(basedir+'combine'): os.mkdir(basedir+'combine/')
    if not os.path.exists(basedir+'combine/logs/'): os.mkdir(basedir+'combine/logs/')
    if not os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/'): os.mkdir(localdir+'dnidever/nsc/instcal/'+version+'/')
    if not os.path.exists(basedir+'logs'): os.mkdir(basedir+'logs/')
    # Hosts
    #if hosts is None: hosts = ['gp09','hulk','thing']
    if (host not in hosts):
        print('Current HOST='+host+' not in list of HOSTS = [ '+','.join(hosts)+' ] ')
        sys.exit()

    # Turn off mulit-threading
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    # Log file
    #------------------
    # format is nsc_combine_main.DATETIME.log
    ltime = time.localtime()
    # time.struct_time(tm_year=2019, tm_mon=7, tm_mday=22, tm_hour=0, tm_min=30, tm_sec=20, tm_wday=0, tm_yday=203, tm_isdst=1)
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
    logfile = basedir+'combine/logs/nsc_instcal_combine_main.'+logtime+'.log'

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

    rootLogger.info("Combining NOAO InstCal catalogs")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("nmulti = "+str(nmulti))
    rootLogger.info("hosts = "+','.join(np.array(hosts)))
    rootLogger.info("redo = "+str(redo))
    
    # Which healpix pixels have data
    listfile = basedir+'lists/nsc_instcal_combine_healpix_list.fits.gz'
    if not os.path.exists(listfile):
        rootLogger.info(listfile+' NOT FOUND.  Run nsc_instcal_combine_qacuts.pro/py')
        sys.exit()

    rootLogger.info("Reading list from "+listfile)
    healstr = fits.getdata(listfile,1)
    index = fits.getdata(listfile,2)
    upix = index['pix']
    npix = len(index)
    # Copy to local directory for faster reading speed
    if os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/'+os.path.basename(listfile)):
                 os.remove(localdir+'dnidever/nsc/instcal/'+version+'/'+os.path.basename(listfile))
    shutil.copyfile(listfile,localdir+'dnidever/nsc/instcal/'+version+'/'+os.path.basename(listfile))
    #file_copy,listfile,localdir+'dnidever/nsc/instcal/'+version+'/',/over

    # Create the commands
    allpix = upix.copy()
    allcmd = dln.strjoin("/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal_combine_cluster.py ",allpix.astype(np.str))
    allcmd = dln.strjoin(allcmd," "+version+" --nside "+str(nside))
    if redo: allcmd = dln.strjoin(allcmd,' -r')
    alldirs = np.zeros(npix,(np.str,200))
    alldirs[:] = tmpdir
    nallcmd = len(allcmd)        

    # Check what's been done already
    if not redo:
        rootLogger.info("Checking if any have already been done")
        exists = np.zeros(dln.size(allpix),bool)+False
        for ip,p in enumerate(allpix):
            outfile = basedir+'combine/'+str(p//1000)+'/'+str(p)+'.fits.gz'
            if os.path.exists(outfile): exists[ip]=True
        bd,nbd,gd,ngd = dln.where(exists,comp=True)
        if ngd==0:
            rootLogger.info('All pixels were previously completed. Nothing to run')
            sys.exit()
        if nbd>0:
            rootLogger.info(str(nbd)+' pixels previously completed.  Removing them. '+str(ngd)+' left.')
            allpix = allpix[gd]
            allcmd = allcmd[gd]
            alldirs = alldirs[gd]
            nallcmd = len(allcmd)

    rootLogger.info(str(nallcmd)+' healpix to process')

    ## Only keep MC region
    #ra,dec = hp.pix2ang(nside,allpix,lonlat=True)
    #coords = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    #lmc = SkyCoord(ra=81.9*u.degree,dec=-69.866*u.degree)
    #lrad = lmc.separation(coords).deg
    #smc = SkyCoord(ra=13.183*u.degree,dec=-72.8283*u.degree)
    #srad = smc.separation(coords).deg
    ##gd,ngd = dln.where( (coords.galactic.b.deg<-10) & (lrad>6) & (srad>6) ((lrad>6) & (lrad<25)) | ((srad>6) & (srad<15)) )
    ##gd,ngd = dln.where( (coords.galactic.b.deg<-10) & (lrad>6) & (srad>6) & ((lrad<25) | (srad<15)) )
    #gd,ngd = dln.where( (coords.galactic.b.deg<-10) & ((lrad<25) | (srad<15)) )
    #rootLogger.info('Only processing '+str(ngd)+' Magellanic Clouds HEALPix')
    #allpix = allpix[gd]
    #allcmd = allcmd[gd]
    #alldirs = alldirs[gd]
    #nallcmd = ngd

    # RANDOMIZE
    np.random.seed(0)
    rnd = np.argsort(np.random.rand(nallcmd))
    allpix = allpix[rnd]
    allcmd = allcmd[rnd]
    alldirs = alldirs[rnd]

    # Parcel out the jobs
    nhosts = dln.size(hosts)
    torun = np.arange(nallcmd)
    nperhost = int(np.ceil(nallcmd/nhosts))
    for i in range(nhosts):
        if host==hosts[i]: torun=torun[i*nperhost:(i+1)*nperhost]
    ntorun = len(torun)
    pix = allpix[torun]
    cmd = allcmd[torun]
    dirs = alldirs[torun]
    rootLogger.info('Running '+str(len(torun))+' on '+host)

    ## Check what's been done already
    #if not redo:
    #    exists = np.zeros(dln.size(pix),bool)+False
    #    for ip,p in enumerate(pix):
    #        outfile = basedir+'combine/'+str(p//1000)+'/'+str(p)+'.fits.gz'
    #        if os.path.exists(outfile): exists[ip]=True
    #    bd,nbd,gd,ngd = dln.where(exists,comp=True)
    #    if ngd==0:
    #        rootLogger.info('All pixels were previously completed. Nothing to run')
    #        sys.exit()
    #    if nbd>0:
    #        rootLogger.info(str(nbd)+' pixels previously completed.  Removing them. '+str(ngd)+' left.')
    #        pix = pix[gd]
    #        cmd = cmd[gd]
    #        dirs = dirs[gd]

    # Saving the structure of jobs to run
    runfile = basedir+'lists/nsc_instcal_combine_main.'+host+'.'+logtime+'_run.fits'
    rootLogger.info('Writing running information to '+runfile)
    runstr = np.zeros(len(cmd),dtype=np.dtype([('pix',int),('host',(np.str,20))]))
    runstr['pix'] = pix
    runstr['host'] = host
    Table(runstr).write(runfile)

    # Now run the combination program on each healpix pixel
    import pdb; pdb.set_trace()
    a = input("Press RETURN to start")
    jobs = jd.job_daemon(cmd,dirs,hyperthread=True,prefix='nsccmb',nmulti=nmulti)

    # Save the jobs
    Table(jobs).write(basedir+'lists/nsc_instcal_combine_main.'+host+'.'+logtime+'_jobs.fits')
