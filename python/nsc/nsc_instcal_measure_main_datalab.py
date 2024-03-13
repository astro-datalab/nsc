#!/usr/bin/env python

# AUTHORS: David Nidever 
#          david.nidever@montana.edu
#          Katie Fasbender
#          katiefasbender@montana.edu

# nsc_meas_wrapper.py will run the NSC Measurements process on all exposures in a given list on tempest.montana.edu


# This script writes a "job_name.sh" file for x exposures, y times every z seconds.
# What that means:
# - In each "job_name.sh" file, "x" exposures will be analyzed
# - "y" number of job files will be submitted in a batch before sleeping for "z" seconds.
# x = 1, y = numjobs, z = 60



#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from argparse import ArgumentParser
from astropy.table import Table,Column
from astropy.io import fits
from dlnpyutils import utils as dln, coords, job_daemon as jd
import logging
import numpy as np
import os
import socket
import subprocess
import sys
import time

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def makedir(dir):
    '''makes a directory with name "dir" if it does not exist
    Arguments:
    ----------
    dir (str)
            directory name
    Returns:
    --------
        None; directory "dir" is created if not already there
    '''
    if not os.path.exists(dir):
        os.mkdir(dir)

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Run NSC Instcal Measurement.')
    parser.add_argument('--version', type=str, nargs=1, help='Version number')
    parser.add_argument('--host',type=str,nargs=1,help='Delimited list of hosts to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposure that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=70000, help='Max number of jobs')
    parser.add_argument('--nmulti', type=int, nargs=1, default=10, help='Number of jobs')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list of exposures to use')
    #parser.add_argument('--chip',type=int,nargs=1,default=1,help='chip to analyze')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)
    redo = args.redo
    print(redo,type(redo))
    ##nmulti = dln.first_el(args.nmulti)
    hosts=args.host[0].split(',')
    nhost=len(hosts) #npar -> nhost
    maxjobs = args.maxjobs
    nmulti = args.nmulti[0]
    inputlist = args.list
    if inputlist is not None:
        inputlist = inputlist[0]
    nside = 128
    radeg = 180 / np.pi
    #chip = args.chip[0]
    t0 = time.time()

    # Establish necessary directories
    if host == "gp07" or host == "gp09":
        basedir = "/net/dl2/kfas/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/net/dl2/kfas/nsc/instcal/"+version+"/"
        tmpdir = localdir+"tmp/"
    elif host == "thing" or host == "hulk":
        basedir = "/net/dl2/kfas/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/net/dl2/kfas/nsc/instcal/"+version+"/"
        tmpdir = localdir+"tmp/"

    makedir(tmpdir)
    subdirs = ['logs','c4d','k4m','ksb']
    for sub in subdirs:
        makedir(basedir+sub)

    # Log File 
    #---------
    # Create Log file name;
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
    logfile = basedir+'logs/nsc_instcal_measure_main.'+logtime+'.log'

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
    ##rootLogger.info("hosts = "+','.join(np.atleast_1d(hosts)))
    rootLogger.info("redo = "+str(redo))

    # Loading the exposure list(s)
    #-----------------------------
    if inputlist is None:
        rootLogger.info('Using input lists:')
        rootLogger.info('  '+basedir+'/lists/decam_instcal_list.fits.gz')
        rootLogger.info('  '+basedir+'/lists/mosaic3_instcal_list.fits.gz')
        rootLogger.info('  '+basedir+'/lists/bok90prime_instcal_list.fits.gz')
        list1 = fits.getdata(basedir+'/lists/decam_instcal_list.fits.gz',1)
        list2 = fits.getdata(basedir+'/lists/mosaic3_instcal_list.fits.gz',1)
        list3 = fits.getdata(basedir+'/lists/bok90prime_instcal_list.fits.gz',1)
        lstr = dln.concatenate([list1,list2,list3])    # concatenate
        del(list1,list2,list3)
    else:
        rootLogger.info('Using input list: '+inputlist)
        lstr = fits.getdata(inputlist,1) # r for random (about to do)
        # Check that it has all the columns that we need
        needcols = ['INSTRUMENT','FLUXFILE','MASKFILE','WTFILE','DATE_OBS']
        for n in needcols:
            if n not in lstr.dtype.names:
                raise ValueError('Column '+n+' not in '+inputlist)
    nlstr = dln.size(lstr)
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
    #--------------------
    rootLogger.info('Checking on the exposures')
    dtype_expstr = np.dtype([('instrument',str,100),('rawname',str,100),('fluxfile',str,100),('wtfile',str,100),
                             ('maskfile',str,100),('outfile',str,150),('done',bool),('torun',bool),('cmd',str,1000),
                             ('cmddir',str,1000),('submitted',bool),('allexist',bool)])
    expstr = np.zeros(ngdexp,dtype=dtype_expstr)  # string array for exposure info
    for i in range(ngdexp):
        if i % 1000 == 0: rootLogger.info(str(i)+"/"+str(ngdexp))

        # Format exposure info for string array
        instrument = lstr['INSTRUMENT'][gdexp[i]].strip()
        if type(instrument) is bytes: instrument=instrument.decode()
        rawname = lstr['RAWNAME'][gdexp[i]].strip()
        if type(rawname) is bytes: rawname=rawname.decode()
        fluxfile = lstr['FLUXFILE'][gdexp[i]].strip()
        if type(fluxfile) is bytes: fluxfile=fluxfile.decode()
        wtfile = lstr['WTFILE'][gdexp[i]].strip()
        if type(wtfile) is bytes: wtfile=wtfile.decode()
        maskfile = lstr['MASKFILE'][gdexp[i]].strip()
        if type(maskfile) is bytes: maskfile=maskfile.decode()
        fdir,base = os.path.split(fluxfile)

        # Change the root directory name to reflect host repo structure
        # format on tempest will be basedir+/exposures/
        lo = fluxfile.find('/mss1/')
        fluxfile = mssdir+fluxfile[(lo+6):]
        #fluxfile = fluxfile.split('/')[-1]
        lo = wtfile.find('/mss1/')
        wtfile = mssdir+wtfile[(lo+6):]
        #wtfile = wtfile.split('/')[-1]
        lo = maskfile.find('/mss1/')
        maskfile = mssdir+maskfile[(lo+6):]
        #maskfile = maskfile.split('/')[-1]

        expstr['instrument'][i] = instrument
        expstr['rawname'][i] = rawname
        expstr['fluxfile'][i] = fluxfile
        expstr['wtfile'][i] = wtfile
        expstr['maskfile'][i] = maskfile

        # Check if the output already exists.
        dateobs = lstr['DATE_OBS'][gdexp[i]]
        if type(dateobs) is np.bytes_: dateobs=dateobs.decode()
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        baseroot = base[0:base.find('.fits.fz')]
        outfile = basedir+instrument+'/'+night+'/'+baseroot+'/'+baseroot+'.log'
        expstr['outfile'][i] = outfile
        #print("outfile = ",outfile)

        # Does the output file exist?
        expstr['done'][i] = False
        if os.path.exists(outfile): 
            expstr['done'][i] = True
            if i%100==0:print(outfile)
        if expstr['fluxfile'][i]==expstr['wtfile'][i]: expstr['done'][i] = True

        # If no outfile exists or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            #if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
            expstr['cmd'][i] = 'python '+basedir+'nsc_instcal_meas.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version
            expstr['cmddir'][i] = tmpdir
            expstr['torun'][i] = True
        # If outfile exists and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False
            #rootLogger.info(outfile+' EXISTS and --redo NOT set')

    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(expstr['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    rootLogger.info(str(ntorun)+" exposures")
    if ntorun == 0:
        rootLogger.info('No exposures to process.')
        sys.exit()
    nperhost = int(np.ceil(nalltorun/nhost))
    for i in range(nhost):
        if host==hosts[i]: torun=torun[i*nperhost:(i+1)*nperhost]
    rootLogger.info(str(ntorun)+' exposures to process on '+str(nhost)+' hosts, and '+str(nperhost)+' jobs per host.')

    # Pick the jobs to run
    if ntorun>maxjobs:
        rootLogger.info("More jobs than MAXJOBS. Cutting down to "+str(maxjobs)+" jobs")
        expstr['submitted'][torun[0:maxjobs]] = True
    else:
        expstr['submitted'][torun] = True
    tosubmit, = np.where(expstr['submitted']==True)
    ntosubmit = len(tosubmit)
    rootLogger.info(str(ntosubmit)+" jobs to submit")
    cmd = expstr[tosubmit]['cmd']
    cmddir = expstr[tosubmit]['cmddir']

    # Saving the structure of jobs to run
    runfile = basedir+'/nsc_instcal_measure_main.'+logtime+'_run.fits'
    rootLogger.info("Writing running information to "+runfile)
    if os.path.exists(runfile): os.remove(runfile)
    Table(expstr).write(runfile)

    # Now run measurement on each exposure
    #import pdb; pdb.set_trace()
    #a = input("Press RETURN to start")
    jobs = jd.job_daemon(cmd,cmddir,hyperthread=True,prefix='nscmeas',waittime=5,nmulti=nmulti)

    # Save jobs
    Table(jobs).write(basedir+'lists/nsc_instcal_measure_main.'+host+'.'+logtime+'_jobs.fits')

    rootLogger.info("Done")
