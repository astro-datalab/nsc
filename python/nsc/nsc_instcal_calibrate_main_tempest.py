#!/usr/bin/env python

# AUTHORS: David Nidever
#          david.nidever@montana.edu
#          Katie Fasbender
#          katiefasbender@montana.edu

# nsc_instcal_calibrate_main_tempest.py will run the NSC Calibration process on all exposures
# in a given list on tempest.montana.edu

# This script writes a "job_name.sh" file for an exposure, maintaining "maxjobs" jobs
# running across defined slurm partitions on tempest.
# What that means:
# - Each "job_name.sh" script written by THIS script includes the command to
#   run the mmt process on an exposure from provided NSC exposure list (--list)
# -User defins the slurm partitions (--partition) & how many exposures to process at once (--maxjobs)
# - This script will cycle through the exposure list, submitting jobs to the slurm queue until
#   all exposures are analyzed or something goes wrong.



#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from argparse import ArgumentParser
from astropy.table import Table,Column
from astropy.io import fits
from dlnpyutils import utils as dln, coords
import glob
import logging
import numpy as np
import os
import shutil
import socket
import subprocess
import sys
import time

from slurm_funcs import *

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Run NSC Instcal Calibration Process.')
    parser.add_argument('--version', type=str, nargs=1,help='Version number')
    parser.add_argument('--partitions',type=str,nargs=1,help='Comma-delimited list of partitions to divide HEALPix jobs between')
    parser.add_argument('--maxhp', type=int, nargs=1, default=1,help='Max number of HEALPix to process at any given time (1 job per HP)') 
    parser.add_argument('--maxexp',type=str,nargs=1,default=10,help='Maximum number of exposures to calibrate at once (maxexp exposures per job, each a task)')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list of HEALPix to use')
    parser.add_argument('--nside', type=str, nargs=1,help='HEALPix NSIDE (ring ordering, please!)')
    parser.add_argument('--gsynth', action='store_true',help='Photometric Calibration with Gaia Syntetic Photometry?')
    parser.add_argument('--psf',action='store_true',help='Astrometric Calibration with PSF coordinates?')
    parser.add_argument('-r','--redo', action='store_true',help='Redo HEALPix that were previously processed')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)     # NSC version
    partitions=args.partitions[0].split(',') # the slurm partitions to submit jobs to
    npar=len(partitions)                     # number of slurm partitions
    maxhp = int(args.maxhp[0])                  # maximum number of HEALPix to maintain running at any time
    maxexp = int(args.maxexp[0])                # maximum number of exposures to maintain running at any time
    nside = int(args.nside[0])                  # HEALPix NSIDE, ring ordering
    nchan = maxhp//npar                      # number of job channels per partition cpar -> nchan
    chans = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,nchan)],maxhp)
    nchantot = len(chans)
    #print("job channels = ",chans)
    gsynth = args.gsynth                     # Photometric calibration with Gaia Synthetic Photometry?
    if gsynth: gsth = " --gsynth "
    else: gsth = ""
    psf = args.psf                           # Astrometric calibration with PSF coordinates?
    if psf: ps = " --psf"
    else: ps = ""
    redo = args.redo                         # if called, redo = True
    if redo: rdo = " --redo"
    else: rdo = ""
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of HEALPix-labeled exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories
    # -tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"   # location of operations
    outfiledir = basedir+"outfiles/"                     # a place for the job files
    inddir = basedir+"healpix_indicators/"               # a place for files that indiate a complete healpix (all exposures run)
    refcatdir = basedir+"refcats/"                       # a place to save the refcat for each healpix
    makefile(inddir)
    makefile(refcatdir)
    makedir(outfiledir)

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
    logfile = basedir+'logs/nsc_instcal_calibrate_main.'+logtime+'.log'

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

    rootLogger.info("Calibrating NOIRLab Measurement Catalogs by HEALPix (NSIDE = "+str(nside)+")")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("partitions = "+str(partitions))
    rootLogger.info("redo = "+str(redo))

    # Loading the exposure list(s)
    #-----------------------------
    if inputlist is None:
        rootLogger.info('Using input list:')
        rootLogger.info('  '+basedir+'/lists/nsc_calibrate_healpix_list.fits.gz')
        lstr = fits.getdata(basedir+'/lists/nsc_calibrate_healpix_list.fits.gz')

    else:
        rootLogger.info('Using input list: '+inputlist)
        lstr = fits.getdata(inputlist,1) # r for random (about to do)
        # Check that it has all the columns that we need
        needcols = ['filter','instrument','expdir','blackmore_repo','pix']
        for n in needcols:
            if n not in lstr.dtype.names:
                raise ValueError('Column '+n+' not in '+inputlist)
    nlstr = dln.size(lstr)
    rootLogger.info(str(nlstr)+' measured exposures')
    unique_hpix = np.unique(lstr['pix'])
    nhpix = len(unique_hpix)
    rootLogger.info(str(nhpix)+' HEALPix, NSIDE='+str(nside))

    # Check the HEALPix
    #--------------------
    rootLogger.info('Checking on the HEALPix')
    dtype_hpstr = np.dtype([('pix',int),('partition',str,100),('submitted',bool),
                             ('hp_done',bool),('torun',bool),
                             ('jobname',str,100),('jobid',str,100),('jobstatus',str,100),
                             ('cmd',str,1000),('cputime',str,100),('maxrss',str,100)])
    hpstr = np.zeros(nhpix,dtype=dtype_hpstr)  # stuctured array for exposure info
    hpstr = Table(hpstr)
    hpstr['pix'] = unique_hpix
    hpstr['submitted'] = False
    hpstr['jobid']="-99.99"
    hpstr['hp_done'] = False
    for i in range(nhpix):
        if i % 50000 == 0: rootLogger.info("HEALPix count "+str(i)+"/"+str(nhpix))
        pix = int(unique_hpix[i]) # HEALPix index
        part = chans[i%nchantot]
        hpstr['pix'][i] = pix
        hpstr['partition'][i] = part
        # Check indicator directory to see if all exposures for the healpix have been calibrated
        # (if the healpix indicator file exists, the healpix has been completed!)
        hp_inddir = inddir+str(nside)+"_"+str(int(pix))+"/"
        makedir(hp_inddir)
        indfile = hp_inddir+str(int(pix))+".txt"
        if os.path.exists(indfile):  hpstr['done'][i] = True
        # If HEALPix is not completed or yes redo: run it!
        if (hpstr['hp_done'][i]==False) or (redo):
            hpstr['cmd'][i] = 'python '+basedir+'calibrate_healpix_tempest.py --version '+str(version)+' --pix '+str(pix)+' --nside '+str(nside)+' --maxexp '+str(maxexp)+' --partition '+part.split("_")[0]+gsth+ps+rdo
            hpstr['torun'][i] = True
        # else, no run
        else: hpstr['torun'][i] = False

    # Parcel out jobs
    #----------------
    # Define HEALPix to run & total #jobs/partition
    torun,nalltorun = dln.where(hpstr['torun'] == True)    # Total number of jobs to run (# HEALPix)
    hpstr['torun'][torun] = True
    ntorun = len(torun)
    rootLogger.info(str(ntorun)+" HEALPix")
    if ntorun == 0:
        rootLogger.info('No HEALPix to process.')
        sys.exit()
    njchan = 0 # number of jobs per channel (divide evenly)
    if ntorun!=0: njchan = (ntorun//nchan)+1
    rootLogger.info(str(ntorun)+' HEALPix (NSIDE='+str(nside)+') to process, '+str(maxhp)+' at a time on '+str(npar)+' slurm partitions.')


    # Start submitting jobs
    #----------------------
    runfile = basedir+'lists/nsc_instcal_calibrate_main.'+logtime+'_run.fits'
    Table(hpstr).write(runfile)
    jsub = 0          # one for each HEALPix to analyze, increment after submitting exposure processing job
    jcomp = 0         # one for each HEALPix to analyze, increment after a job has finished running
    pc = 0            # for looping through job channel list

    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#
    # For each channel, check status of last job & submit a new one accordingly
    while jcomp<ntorun:   # -> until we have nexp exposures to download...
        part = chans[pc]      # select channel
        rootLogger.info("Checking partition "+part+" queue, "+str(ntorun-jcomp)+" out of "+str(ntorun)+" HEALPix left")
        # - Check last HEALPix processing (proc) job submitted (l)
        #              (job status?)
        # - Check next HEALPix up in the queue (n)
        #               (last job not running/pending?) = r(eady)
        #              (all files exist?)  = e(xist)
        partition_ind = set(np.where(hpstr['partition'][torun]==part)[0])  # indices for this partition
        unsubmitted_ind = set(np.where(hpstr['submitted'][torun]==0)[0])   # indices for unsubmitted HPix
        submitted_ind = set(np.where(hpstr['submitted'][torun]==1)[0])     # indices for submitted HPix
        sub = list(partition_ind & submitted_ind)                          # indices of HPix jobs submitted to channel (last_sub)
        unsub = list(partition_ind & unsubmitted_ind)                      # indices of HPix jobs unsubmitted to channel (this_sub)

        #--------------------------------------------------------------#
        # -> get index  & info of LAST HPix job submitted to channel (l)
        if len(sub)>0:
            sub = np.sort(sub)
            lsub = [sub[-1]] # index of LAST EXPOSURE submitted
            lsub_jname = hpstr['jobname'][torun[lsub]][0]
            if maxhp<5:                                                  # let rest if not looping through many jobs
                rootLogger.info("wait a sec before checking last proc job...")
                time.sleep(sleep_time)
            lsub_jstat,lsub_jid = sacct_cmd(lsub_jname,["state","jobid"])  # check the status of last proc job
            hpstr['jobstatus'][torun[lsub]] = lsub_jstat
            hpstr['jobid'][torun[lsub]] = lsub_jid
            if lsub_jstat=="RUNNING" or lsub_jstat=="REQUEUED" or lsub_jstat=="PENDING": # If last job is still going, don't submit another.
                r==False
            else:                                                                        # Else, allow next one up to be submitted
                r==True
                jcomp+=1
                # check for healpix completeness indicator
                rootLogger.info("HEALPix "+str(hpstr[torun[lsub]]['pix'])+" processing over! Check out next pix up in "+str(part)+" partition")
        else: r = True                                                                    # If no job submitted to partition, go ahead and submit one.
        #-----------------------------------------------------------------#
        # -> get index & info of NEXT HPix to check & submit to channel (n)
        if len(unsub)>0:
            nsub = [np.sort(unsub)[0]] # index of NEXT HPix to be submitted
            rootLogger.info("Checking out next pix up, "+str(hpstr['pix'][torun[nsub]])+"...")
            if r:
                rootLogger.info("-- Submitting next job up to "+str(part)+" channel --")
                pix = hpstr['pix'][torun[nsub]][0]
                cmd = hpstr['cmd'][torun[nsub]][0]
                partition = hpstr['partition'][torun[nsub]][0].split("_")[0]
                rootLogger.info("--Submit Processing Job for HEALPix "+str(pix)+"--")
                # -> Write exposure processing job script to file
                job_name = 'nsc_calpix_'+str(pix)+'_'+str(logtime)+'_'+str(jsub)
                job_file = write_jscript(job_name,jsub,partition,[cmd],outfiledir,"katiefasbender@montana.edu",cpupt=2,mempc="4G",rtime="7-00:00:00",parallel=False)
                # -> Submit exposure processing job to slurm queue
                os.system("sbatch "+str(job_file))
                hpstr['submitted'][torun[nsub]] = True
                jsub+=1
                rootLogger.info("Job "+job_name+"  submitted to "+partition+" partition")
                hpstr['jobname'][torun[nsub]] = job_name
            else: rootLogger.info("Partition "+part+" has a job still running, or no HPix up, moving to next channel.")
        else: rootLogger.info("No more hpix to submit in this partition")
        #------------------------------#
        # - Set up for next partition -
        pc+=1
        pc = pc-(pc//nchantot)*nchantot # stay within range of partitions list
        # - Save job structure -
        Table(hpstr).write(runfile,overwrite=True)
        if ntorun<10 or maxhp<5:
           print("sleeping for ten sec before checking next partition...")
           time.sleep(sleep_time)
        else: time.sleep(sleep_time//10)

    Table(hpstr).write(runfile,overwrite=True)
    rootLogger.info("job structure saved to "+runfile+", HP processing complete!")
