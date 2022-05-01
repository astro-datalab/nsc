#!/usr/bin/env python

import os
import sys
import shutil
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table
from dlnpyutils import utils as dln, coords, bindata, dbutils, job_daemon as jd
from astropy import units as u
from astropy.coordinates import SkyCoord
#import subprocess
import time
from argparse import ArgumentParser
import socket
import logging
import healpy as hp

# Driver for nsc_instcal_measure_update.py to update OBJECTIDs in exposure measurement catalogs
if __name__ == "__main__":
    parser = ArgumentParser(description='Update objectIDs in exposure measurement catalogs')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    parser.add_argument('--hosts', type=str, nargs=1, help='Delimited list of hosts')
    parser.add_argument('-l','--list', type=str, nargs=1, default='', help='List of exposures to run')
    parser.add_argument('-nm','--nmulti', type=int, nargs=1, default=12, help='Number of jobs')
    #parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)
    #redo = args.redo
    redo = False
    nmulti = dln.first_el(args.nmulti)
    if args.hosts is not None:
        hosts = args.hosts[0].split(',')
    else:
        hosts = [host]
    inplistfile = dln.first_el(args.list)
    if inplistfile == '': inplistfile = None
    nside = 128

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        basedir = "/net/dl2/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        basedir = "/net/dl2/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    if not os.path.exists(basedir+'combine'): os.mkdir(basedir+'combine/')
    if not os.path.exists(basedir+'combine/logs/'): os.mkdir(basedir+'combine/logs/')
    if not os.path.exists(basedir+'logs'): os.mkdir(basedir+'logs/')
    # Hosts
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
    logfile = basedir+'combine/logs/nsc_instcal_measure_update_main.'+logtime+'.log'

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

    rootLogger.info("Updating objectID in NSC exposure measurement catalogs")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("nmulti = "+str(nmulti))
    rootLogger.info("hosts = "+','.join(np.array(hosts)))
    #rootLogger.info("redo = "+str(redo))

    # Exposures to run
    #listfile = basedir+'lists/nsc_instcal_combine_healpix_list.db'
    #listfile = basedir+'lists/nsc_instcal_combine_healpix_list.fits.gz'
    listfile = basedir+'lists/nsc_'+version+'_exposures.fits.gz'
    if inplistfile is None:
        rootLogger.info("Reading list from "+listfile)

        #hlist = dbutils.query(listfile,'hlist',cols='FILE')
        #nlist = len(hlist)
        #healstr = fits.getdata(listfile,1)
        #expdir = [os.path.dirname(f) for f in healstr['FILE']]
        #expdir = np.unique(expdir)
        liststr = fits.getdata(listfile,1)
        expdir = liststr['EXPDIR']
        # change /dl1/users/dnidever/ to /dl2/dnidever/
        expdir = expdir.replace('/dl1/users/dnidever/','/dl2/dnidever/')
        # Trim trailing /
        #expdir = expdir.rstrip('/')  # this doesn't work for some reason
        expdir = np.char.array([a.rstrip('/') for a in expdir])
        nexpdir = len(expdir)
    else:
        rootLogger.info("Reading list from "+inplistfile)        
        expdir = dln.readlines(inplistfile)
        nexpdir = len(expdir)

    # Create the commands
    allexpdir = expdir.copy()
    allcmd = dln.strjoin("/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal_measure_update.py ",allexpdir)
    alldirs = np.zeros(nexpdir,(np.str,200))
    alldirs[:] = tmpdir
    nallcmd = len(allcmd)        

    # Check what's been done already
    check = False
    #if not redo:
    if check:
        rootLogger.info("Checking if any have already been done")
        exists = np.zeros(dln.size(allexpdir),bool)+False
        for ip,p in enumerate(allexpdir):
            base = os.path.basename(allexpdir[ip])
            outfile = allexpdir[ip]+'/'+base+'.updated'
            if os.path.exists(outfile): exists[ip]=True
            if ip % 1000 == 0: print(str(ip)+' '+str(p))
        bd,nbd,gd,ngd = dln.where(exists,comp=True)
        if ngd==0:
            rootLogger.info('All exposures were previously updated. Nothing to run')
            sys.exit()
        if nbd>0:
            rootLogger.info(str(nbd)+' exposures previously updated.  Removing them from the list. '+str(ngd)+' left.')
            allexpdir = allexpdir[gd]
            allcmd = allcmd[gd]
            alldirs = alldirs[gd]
            nallcmd = len(allcmd)

    rootLogger.info(str(nallcmd)+' exposures to process')

    # RANDOMIZE
    np.random.seed(0)
    rnd = np.argsort(np.random.rand(nallcmd))
    allexpdir = allexpdir[rnd]
    allcmd = allcmd[rnd]
    alldirs = alldirs[rnd]

    # Parcel out the jobs
    nhosts = dln.size(hosts)
    torun = np.arange(nallcmd)
    nperhost = int(np.ceil(nallcmd/nhosts))
    for i in range(nhosts):
        if host==hosts[i]: torun=torun[i*nperhost:(i+1)*nperhost]
    ntorun = len(torun)
    expdir = allexpdir[torun]
    cmd = allcmd[torun]
    dirs = alldirs[torun]
    rootLogger.info('Running '+str(len(torun))+' on '+host)

    # Saving the structure of jobs to run
    runfile = basedir+'lists/nsc_instcal_measure_update_main.'+host+'.'+logtime+'_run.fits'
    rootLogger.info('Writing running information to '+runfile)
    runstr = np.zeros(len(cmd),dtype=np.dtype([('expdir',(np.str,200)),('host',(np.str,20))]))
    runstr['expdir'] = expdir
    runstr['host'] = host
    Table(runstr).write(runfile)

    # Now run the measurement update program on the exposures
    import pdb; pdb.set_trace()
    a = input("Press RETURN to start")
    jobs = jd.job_daemon(cmd,dirs,hyperthread=True,prefix='nscmeasupd',nmulti=nmulti)

    # Save the jobs
    Table(jobs).write(basedir+'lists/nsc_instcal_measure_update_main.'+host+'.'+logtime+'_jobs.fits')

    rootLogger.info('dt='+str(time.time()-t0)+' sec')    
