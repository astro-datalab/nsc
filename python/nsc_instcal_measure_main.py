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


# Run Source Extractor on many NSC exposures
if __name__ == "__main__":
    parser = ArgumentParser(description='Run NSC Instcal Measurement.')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    parser.add_argument('--hosts', type=str, nargs=1, help='Delimited list of hosts')
    parser.add_argument('-nm','--nmulti', type=int, nargs=1, default=30, help='Number of jobs')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposure that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=70000, help='The maximum number of exposures to attempt to process per host')
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
    maxjobs = args.maxjobs
    nside = 128
    radeg = 180 / np.pi
    t0 = time.time()

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmpdir = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    if not os.path.exists(tmpdir): os.mkdir(tmpdir)
    subdirs = ['logs','c4d','k4m','ksb']
    for sub in subdirs:
        if not os.path.exists(basedir+sub): os.mkdir(basedir+sub)

    # Hosts
    if (host not in hosts):
        print('Current HOST='+host+' not in list of HOSTS = [ '+','.join(np.atleast_1d(hosts))+' ] ')
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
    logfile = basedir+'combine/logs/nsc_instcal_measure_main.'+logtime+'.log'

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

    rootLogger.info("Creating NOAO InstCal Measurement catalogs")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("nmulti = "+str(nmulti))
    rootLogger.info("hosts = "+','.join(np.atleast_1d(hosts)))
    rootLogger.info("redo = "+str(redo))
    
    # Loading the lists
    list1 = fits.getdata(basedir+'/lists/decam_instcal_list.fits.gz',1)
    list2 = fits.getdata(basedir+'/lists/mosaic3_instcal_list.fits.gz',1)
    list3 = fits.getdata(basedir+'/lists/bok90prime_instcal_list.fits.gz',1)
    lstr = dln.concatenate([list1,list2,list3])    # concatenate
    nlstr = dln.size(lstr)
    del(list1,list2,list3)
    rootLogger.info(str(nlstr)+' InstCal images')

    # Putting them in RANDOM but REPEATABLE order
    rootLogger.info('RANDOMIZING WITH SEED=1')
    rnd = np.arange(nlstr)
    np.random.seed(1)
    np.random.shuffle(rnd)
    lstr = lstr[rnd]

    gdexp = np.arange(nlstr)
    ngdexp = nlstr

    # Check the exposures
    rootLogger.info('Checking on the exposures')
    dtype_expstr = np.dtype([('instrument',np.str,100),('fluxfile',np.str,100),('wtfile',np.str,100),('maskfile',np.str,100),('allexist',bool),
                             ('outfile',np.str,100),('done',bool),('locked',bool),('torun',bool),('cmd',np.str,100),
                             ('cmddir',np.str,100),('submitted',bool)])
    expstr = np.zeros(ngdexp,dtype=dtype_expstr)
    for i in range(ngdexp):
        if i % 500 == 0: rootLogger.info(i)

        instrument = lstr['INSTRUMENT'][gdexp[i]].strip()
        if type(instrument) is bytes: instrument=instrument.decode()
        fluxfile = lstr['FLUXFILE'][gdexp[i]].strip()
        if type(fluxfile) is bytes: fluxfile=fluxfile.decode()
        wtfile = lstr['WTFILE'][gdexp[i]].strip()
        if type(wtfile) is bytes: wtfile=wtfile.decode()
        maskfile = lstr['MASKFILE'][gdexp[i]].strip()
        if type(maskfile) is bytes: maskfile=maskfile.decode()
        fdir,base = os.path.split(fluxfile)

        # Change the root directory name
        #  /net/mss1/blah/blah/
        fluxfile = mssdir+fluxfile[10:]
        wtfile = mssdir+wtfile[10:]
        maskfile = mssdir+maskfile[10:]

        expstr['instrument'][i] = instrument
        expstr['fluxfile'][i] = fluxfile
        expstr['wtfile'][i] = wtfile
        expstr['maskfile'][i] = maskfile

        #import pdb; pdb.set_trace()

        # Check if the output already exists.
        dateobs = lstr['DATE_OBS'][gdexp[i]]
        if type(dateobs) is np.bytes_: dateobs=dateobs.decode()
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        baseroot = base[0:base.find('.fits.fz')]
        #outfile = dldir+'users/dnidever/decamcatalog/instcal/'+night+'/'+baseroot+'/'+baseroot+'_'+strtrim(1,2)+'.fits'
        outfile = basedir+instrument+'/'+night+'/'+baseroot+'/'+baseroot+'_1.fits'
        expstr['outfile'][i] = outfile

        # Do all three files exist?
        #if file_test(fluxfile) eq 1 and file_test(wtfile) eq 1 and file_test(maskfile) eq 1 then expstr[i].allexist=1
        expstr['allexist'][i] = True    # THIS TAKES TOO LONG!!!
        # Does the output file exist
        #if file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1 then expstr[i].done = 1
        #expstr[i].done = 0
        expstr['done'][i] = False
        if os.path.exists(outfile): expstr['done'][i] = True

        # Not all three files exist
        if expstr['allexist'][i] is False:
            if silent is False: rootLogger.info('Not all three flux/wt/mask files found for '+fluxfile)
            continue

        # Already done
        if (expstr['done'][i] is True) and redo is False:
            if silent is False: rootLogger.info(outfile+' EXISTS and --redo NOT set')
            continue

        #lock = djs_lockfile(outfile)
        lockfile = outfile+'.lock'
        #testlock = file_test(lockfile)
        testlock = False

        # No lock file
        if (testlock is False or unlock is true):
            #dum = djs_lockfile(outfile)  ; this is slow
            #if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
            #if testlock eq 0 then touchzero,outfile+'.lock'  ; this is fast
            expstr['cmd'][i] = '/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal_measure.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version
            expstr['cmddir'][i] = tmpdir
            expstr['torun'][i] = True
        # Lock file exists
        else:
            expstr['locked'][i] = True
            expstr['torun'][i] = False
            if silent is False: rootLogger.info('Lock file exists '+outfile+'.lock')

    # Parcel out the jobs
    nhosts = dln.size(hosts)
    torun,nalltorun = dln.where((expstr['torun'] == True) & (expstr['done'] == False))
    nperhost = int(np.ceil(nalltorun/nhosts))
    for i in range(nhosts):
        if host==hosts[i]: torun=torun[i*nperhost:(i+1)*nperhost]
    ntorun = len(torun)

    if ntorun == 0:
        rootLogger.info('No exposures to process.')
        sys.exit()

    # Pick the jobs to run
    # MAXJOBS
    if ntorun > maxjobs:
        rootLogger.info('More jobs than MAXJOBS.  Cutting down to '+str(maxjobs)+' jobs')
        expstr['submitted'][torun[0:maxjobs]] = True
    else:
        expstr['submitted'][torun] = True
    tosubmit, = np.where(expstr['submitted'] == True)
    ntosubmit = len(tosubmit)
    rootLogger.info(str(ntosubmit)+' jobs to submit')
    cmd = expstr[tosubmit]['cmd']
    cmddir = expstr[tosubmit]['cmddir']


    # Lock the files that will be submitted
    dolock = False
    if dolock is True:
        print('Locking files to be submitted')
        for i in range(ntosubmit):
            outfile = expstr['outfile'][tosubmit[i]]
            if os.path.exists(os.path.dirname(outfile)) is False: os.mkdir(os.path.dirname(outfile))  # make directory 
            lockfile = outfile+'.lock'
            if os.path.exists(lockfile) is False: dln.touch(lockfile)
            expstr['locked'][tosubmit[i]] = True

    # Saving the structure of jobs to run
    runfile = basedir+'lists/nsc_instcal_measure_main.'+hostname+'.'+logtime+'_run.fits'
    rootLogger.info('Writing running information to '+runfile)
    if os.path.exists(runfile): os.delete(runfile)
    Table(expstr).write(runfile)


    # Now run measurement on each exposure
    import pdb; pdb.set_trace()
    a = input("Press RETURN to start")
    jobs = jd.job_daemon(cmd,cmddirs,hyperthread=True,prefix='nscmeas',wait=5,nmulti=nmulti)

    # Save the jobs
    Table(jobs).write(basedir+'lists/nsc_instcal_measure_main.'+host+'.'+logtime+'_jobs.fits')

    # Unlocking files
    if dolock is True:
        rootLogger.info('Unlocking processed files')
        lockfiles = dln.strjoin(expstr[tosubmit]['outfile'],'.lock')
        for f in lockfiles:
            if os.path.exists(f): os.delete(f)

    rootLogger.info('dt='+str(time.time()-t0)+' sec')
