#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table
from dlnpyutils import utils as dln, coords, job_daemon
#import subprocess
import time
from argparse import ArgumentParser
import socket

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC Instcal Catalogs.')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = args.version
    redo = args.redo
    nmulti = args.nmulti
    nside = 128
    radeg = 180 / np.pi

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        basedir = "/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    if ~os.path.exists(tmpdir): os.mkdir(tmpdir)
    if ~os.path.exists(basedir+'combine'): os.mkdir(basedir+'combine/')
    if ~os.path.exists(basedir+'combine/logs/'): os.mkdir(basedir+'combine/logs/')
    if ~os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/'): os.mkdir(localdir+'dnidever/nsc/instcal/'+version+'/')
    if ~os.path.exists(basedir+'logs'): os.mkdir(basedir+'logs/')
    # Hosts
    if hosts is None: hosts = ['gp09','hulk','thing']
    if ~(host in hosts):
      print('Current HOST='+host+' not in list of HOSTS = [ '+','.join(hosts)+' ] ')
      sys.exit()


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
    #JOURNAL,logfile

    print('Reading list from '+listfile)
    healstr = fits.getdata(listfile,1)
    index = fits.getdata(listfile,2)
    upix = index['pix']
    npix = len(index)
    # Copy to local directory for faster reading speed
    file_copy,listfile,localdir+'dnidever/nsc/instcal/'+version+'/',/over

    # Create the commands
    allpix = upix.copy()
    allcmd = "/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal_combine.py "+str(allpix)+" "+version+" --nside "+str(nside)
    if redo: allcmd += ' -r'
    alldirs = np.zeros(npix,(np.str,200))
    alldirs[:] = tmpdir
    nallcmd = len(allcmd)

    # RANDOMIZE
    np.random.seed(0)
    rnd = np.argsort(np.random.rand(npix))
    allpix = allpix[rnd]
    allcmd = allcmd[rnd]
    alldirs = alldirs[rnd]

    # Parcel out the jobs
    nhosts = len(hosts)
    torun = np.arange(nallcmd)
    nperhost = int(np.ceil(nallcmd/nhosts))
    for i in range(nhosts):
        if host==hosts[i]: torun=torun[i*nperhost:(i+1)*nperhost]
    ntorun = len(torun)
    pix = allpix[torun]
    cmd = allcmd[torun]
    dirs = alldirs[torun]
    print('Running '+str(len(torun))+' on '+hostname)

    # Saving the structure of jobs to run
    runfile = dir+'lists/nsc_instcal_combine_main.'+hostname+'.'+logtime+'_run.fits'
    print('Writing running information to '+runfile)
    runstr = np.zeros(len(cmd),dtype=np.dtype([('pix',int),('host',(np.str,20))]))
    runstr['pix'] = pix
    runstr['host'] = hostname
    Table(runstr).write(runfile)

    # Now run the combination program on each healpix pixel
    a = raw_input("Press RETURN to start")
    jobs = job_daemon(cmd,cmddir,hyperthread=True,prefix='nsccmb',nmulti=nmulti)

    # End logfile
    #------------
    #JOURNAL
