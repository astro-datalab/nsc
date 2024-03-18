#!/usr/bin/env python

# AUTHORS: David Nidever <david.nidever@montana.edu>
#          Katie Fasbender <katiefasbender@montana.edu>

# nsc_instcal_measure_main.py will run the NSC Measurements process on all exposures
# in a given list on whatever server.  If slurm is available, it should be used, otherwise
# the script will use dlnpyutils jobdaemon.py to submit and maintain "jobs".

# If using slurm.......
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
from dlnpyutils import utils as dln, coords, job_daemon as jd
import glob
import healpy as hp
import logging
import numpy as np
import os
import re
import shutil
import socket
import subprocess
import sys
import time

from slurm_funcs import *

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def complete_job(struct,run_ind,lsub_ind,lsub_jname,base_dir,exp_dir,logger):
    ''' after an exposure has been measured and the job is completed, delete the
        exposure flux,wt,maskfiles if we're not on dl servers, and delete the
        untarred repository. 
    Arguments
    ---------
    struct :job structure
    run_ind
    lsub_ind
    lsub_jname
    base_dir
    exp_dir
    logger
    '''
    # -> get info about completed job
    lj_cputime,lj_maxrss = sacct_cmd(lsub_jname,["cputimeraw","maxrss"],c=True)
    struct['done'] = True
    struct['cputime'][run_ind[lsub_ind]] = lj_cputime
    struct['maxrss'][run_ind[lsub_ind]] = lj_maxrss
    # -> delete uncompressed folder
    outdir = struct['outdirectory'][run_ind[lsub_ind]][0].split(".tar.gz")[0]
    print("outdir = ",outdir,os.path.exists(outdir))
    if os.path.exists(outdir):
        logger.info("Deleting uncompressed output directory for exposure "+str(struct['fluxfile'][run_ind[lsub_ind]]))
        out_files = glob.glob(outdir+"/*")
        for f in out_files:
            try: os.remove(f)
            except: logger.info("File "+f+" not found in exposure repo.")
        shutil.rmtree(struct['outdirectory'][run_ind[lsub_ind]][0].split(".tar.gz")[0])
    # -> remove downloaded files!                         XX
    logger.info("Removing downloaded files for exposure "+str(struct['fluxfile'][run_ind[lsub_ind]]))
    if os.path.exists(exp_dir+struct['fluxfile'][run_ind[lsub_ind]][0].strip()):
        os.remove(exp_dir+struct['fluxfile'][run_ind[lsub_ind]][0].strip())
    if os.path.exists(exp_dir+struct['maskfile'][run_ind[lsub_ind]][0].strip()):
        os.remove(exp_dir+struct['maskfile'][run_ind[lsub_ind]][0].strip())
    if os.path.exists(exp_dir+struct['wtfile'][run_ind[lsub_ind]][0].strip()):
        os.remove(exp_dir+struct['wtfile'][run_ind[lsub_ind]][0].strip())
    # -> remove leftover temporary directories            XX
    odir = base_dir+"tmp/"+(struct['fluxfile'][run_ind[lsub_ind]][0].strip()).split(".")[0]+"."
    for itt in [1,2,3]:                             #     XX
        if os.path.exists(odir+str(itt)): shutil.rmtree(odir+str(itt),ignore_errors=True)
    # -> check for actual existence of compressed folder
    if os.path.exists(struct['outdirectory'][run_ind[lsub_ind]][0].strip()): return True
    else: return False

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Run NSC Instcal Measurement Process.')
    parser.add_argument('--version', type=str, nargs=1, help='Version number')
    parser.add_argument('--host',type=str,nargs=1,help='Hostname, or delimited list of hosts to divide jobs between if on datalab')    
    parser.add_argument('--partition',type=str,nargs=1,default="None",help='Delimited list of partitions to divide jobs between, if using slurm')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposures that were previously processed')
    parser.add_argument('--nmulti', type=str, nargs=1, default=0, help='Number of jobs, if on datalab servers')
    parser.add_argument('--maxjobs', type=str, nargs=1, default=1, help='Max number of exposures to process at any given time')
    parser.add_argument('--nside',type=str,nargs=1,default=32,help='HEALPix NSIDE (ring ordering, please!) to sort list for processing order, 0 for random order')    
    parser.add_argument('--nexp',type=str,nargs=1,default=0,help='Number of exposures to download at once, 0 if not downloading')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list of exposures to use')
    args = parser.parse_args()

    # Start time
    t0 = time.time()

    # Inputs
    version = dln.first_el(args.version)     # NSC version (like "v4")
    hosts = args.host[0].split(',')          # host(s) (if multiple, we're using datalab)
    nhost = len(hosts)
    if nhost==1:
        host = hosts[0]
    else:
        hostname = socket.gethostname()
        host = hostname.split('.')[0]
    partitions=args.partition[0].split(',')  # the slurm partitions to submit jobs to, if we're using slurm
    redo = args.redo                         # if called, redo = True
    npar = len(partitions)                     # number of slurm partitions
    nmulti = int(args.nmulti[0])                  # for running on datalab
    maxjobs = int(args.maxjobs[0])           # maximum number of jobs to maintain running at any time
    nchan = maxjobs//npar                    # number of job channels per partition, if using slurm
    nside = int(args.nside[0])                  # HEALPix nside (ns of cca version)
    nexp = int(args.nexp[0])                 # number of exposures to download at once, if not on datalab servers
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories
    if host=="tempest_katie" or host=="tempest_group":
        if host=="tempest_katie":
            basedir = "/home/x25h971/nsc/instcal/"+version+"/"       # location of operations
        elif host=="tempest_group":
            basedir = "/home/group/davidnidever/nsc/instcal/"+version+"/"
        scriptdir = basedir                                          # where scripts will be
        expdir = basedir+"exposures/"                                # where exposures will be downloaded
        unavail_dir = basedir+"unavails/"                            # where we keep a log of unavailable exposures
        x_dir = basedir+"xvers/"                                     # where we keep a log of exposures with version "vX" naming
        makedir(expdir)
        makedir(unavail_dir)
        makedir(x_dir)
        unavails = subprocess.getoutput("ls -ltrh "+basedir+"/unavails/").split("\n")
        unavails = np.array([i.split(" ")[-1].split("_unavails.txt") for i in unavails])
        xs = subprocess.getoutput("ls -ltrh "+basedir+"/xvers/").split("\n")
        xs = np.array([i.split(" ")[-1].split("_xvers.txt")[0] for i in xs])
    elif host=="cca":
        basedir = os.getcwd()+"/"+version+"/"
        scriptdir = basedir
        expdir = basedir+"exposures/"
    elif host=="gp09" or host=="gp07":
        basedir = "/net/dl2/kfas/nsc/instcal/"+version+"/"
        scriptdir = basedir # was localdir
        expdir = "/net/mss1/"
    else:
        print("we don't support this host yet, sorry!!!!!! <3")
        sys.exit(0)
    tmpdir = basedir+"tmp/"
    makedir(tmpdir)
    outfiledir = basedir+"outfiles/"
    makedir(outfiledir)
    subdirs = ['logs','c4d','k4m','ksb']
    for sub in subdirs:
        makedir(basedir+sub)


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

    rootLogger.info("Creating NOIRLab InstCal Measurement catalogs")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    rootLogger.info("partitions = "+str(partitions))
    rootLogger.info("nmulti = "+str(nmulti))
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
        needcols = ['INSTRUMENT','FLUXFILE','MASKFILE','WTFILE','DATE_OBS','RA','DEC']
        for n in needcols:
            if n not in lstr.dtype.names:
                raise ValueError('Column '+n+' not in '+inputlist)
    nlstr = dln.size(lstr)
    rootLogger.info(str(nlstr)+' InstCal images')

    # Putting them in RANDOM but REPEATABLE order
    if nside==0:
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
    if nexp>0:
        dtype_expstr = np.dtype([('instrument',str,100),('ring'+str(nside),int),('rawname',str,100),('fluxfile',str,100),('wtfile',str,100),
                                 ('maskfile',str,100),('outfile',str,150),('logfile',str,150),('partition',str,100),
                                 ('done',bool),('torun',bool),('submitted',bool),('jobname',str,100),('jobid',str,100),('jobstatus',str,100),
                                 ('cmd',str,1000),('cputime',str,100),('maxrss',str,100),('outdirectory',str,500),('exp',int),
                                 ('ur',bool),("x",bool),('esubmitted',bool),('exp_jobname',str,100),('exp_jobid',str,100),('exp_jobstatus',str,100),
                                 ('exp_cmd',str,1000),('exp_cputime',str,100),('exp_maxrss',str,100)])
    else:
        dtype_expstr = np.dtype([('instrument',str,100),('ring'+str(ns),int),('rawname',str,100),('fluxfile',str,100),('wtfile',str,100),
                                 ('maskfile',str,100),('outfile',str,150),('logfile',str,150),('partition',str,100),
                                 ('done',bool),('torun',bool),('submitted',bool),('jobname',str,100),('jobid',str,100),('jobstatus',str,100),
                                 ('cmd',str,1000),('cputime',str,100),('maxrss',str,100),('outdirectory',str,500),('exp',int)])
    expstr = np.zeros(ngdexp,dtype=dtype_expstr)  # string array for exposure info
    for i in range(ngdexp):
        if i % 50000 == 0: rootLogger.info("exposure "+str(i)+"/"+str(ngdexp))

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
        if nexp>0:
            fluxfile = fluxfile.split('/')[-1]
            wtfile = wtfile.split('/')[-1]
            maskfile = maskfile.split('/')[-1]
        elif host=="gp07" or host=="gp09":
            lo = fluxfile.find('/mss1/')
            fluxfile = expdir+fluxfile[(lo+6):]
            lo = wtfile.find('/mss1/')
            wtfile = expdir+wtfile[(lo+6):]
            lo = maskfile.find('/mss1/')
            maskfile = expdir+maskfile[(lo+6):]
        expstr['instrument'][i] = instrument
        expstr['rawname'][i] = rawname
        expstr['fluxfile'][i] = fluxfile
        expstr['wtfile'][i] = wtfile
        expstr['maskfile'][i] = maskfile
        if nside>0: expstr['ring'+str(nside)][i] = hp.ang2pix(nside,lstr['RA'][gdexp[i]],lstr['DEC'][gdexp[i]],lonlat=True)

        # Check if the output already exists.
        dateobs = lstr['DATE_OBS'][gdexp[i]]
        if type(dateobs) is np.bytes_: dateobs=dateobs.decode()
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        baseroot = base[0:base.find('.fits.fz')]

        # Check to see if file is unavailable or has x version:
        if host=="tempest_katie" or host=="tempest_group":
            xvers = False
            if rawname in unavails: expstr['ur'][i] = True
            if rawname in xs:
                expstr['x'][i] = True
                xvers = True
                vers = baseroot.split(".")[0].split("_")[-1]
                versx = re.split('(\d+)',vers)[0]+"x"
                baseroot = baseroot.split(vers)[0]+versx+baseroot.split(vers)[-1]

        # Get output file information
        outroot =      basedir+instrument+'/'+night+'/'+baseroot
        outfile =      outroot+"/"+baseroot+'_1.fits' # outfile for first chip
        outlogfile =   outroot+"/"+baseroot+'.log' # logfile for complete exposure
        outdirectory = outroot+'.tar.gz'         # tarred and zipped exposure directory
        expstr['outfile'][i] = outfile
        expstr['logfile'][i] = outlogfile
        expstr['outdirectory'][i] = outdirectory

        # Does the tarred & zipped output exposure directory exist? (new)
        if os.path.exists(outdirectory): expstr['done'][i] = True
        ## Do the exposure output & logfiles files exist? (old)
        #expstr['done'][i] = False
        #if (os.path.exists(outfile) and os.path.exists(outlogfile)) or (expstr['fluxfile'][i]==expstr['wtfile'][i]):
        #    expstr['done'][i] = True

        # If exposure is completed or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            expstr['cmd'][i] = 'python '+scriptdir+'nsc_instcal_measure.py --fluxfile '+fluxfile+' --wtfile '+wtfile+' --maskfile '+maskfile+' --version '+version+' --host '+host
            if host=="tempest_group" or host=="tempest_katie": expstr['exp_cmd'][i] = 'python '+scriptdir+'get_exposures.py '+rawname+' '+fluxfile+' '+wtfile+' '+maskfile+' '+basedir+' '+expdir
            expstr['torun'][i] = True
        # If exposure is completed and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False

    # Sort by healpix if not random
    if nside!=0: expstr = np.sort(expstr,order="pix"+str(ns))

    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(expstr['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    rootLogger.info(str(ntorun)+" exposures")
    if ntorun == 0:
        rootLogger.info('No exposures to process.')
        sys.exit()

    # ------------------------------------------------------------------------------------------------------------------------
    # If no slurm, use jobdaemon ---------------------------------------------------------------------------------------------
    if partitions[0]=="None":

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
        cmddir = Column(np.repeat(tmpdir,len(expstr[tosubmit])))
    
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
        
    # ------------------------------------------------------------------------------------------------------------------------
    # If slurm is available, use it! -----------------------------------------------------------------------------------------
    else:

        njchan = 0 #number of jobs per channel (divide evenly) njpar -> njchan
        if ntorun!=0: njchan = (ntorun//nchan)+1
        rootLogger.info(str(ntorun)+' exposures to process, '+str(maxjobs)+' at a time on '+str(npar)+' slurm partitions.')
    
        # Split exposures evenly among defined "partitions" as run with nsc_meas_wrapper
        chans = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,nchan)],maxjobs) # partions -> chans
        nchantot = len(chans)
        #print("job channels = ",chans)
        partis = np.array([chans[(i-maxjobs*(i//maxjobs))] for i in range(ntorun)])
        expstr['partition'][torun] = partis

        # Start submitting jobs
        #----------------------
        runfile = basedir+'lists/nsc_instcal_measure_main.'+logtime+'_run.fits'
        expstr['jobid']="-99.99"
        Table(expstr).write(runfile)
        
        # --------------------------------------------------------------------------------------------------------------------
        # If not downloading exposures ---------------------------------------------------------------------------------------
        if nexp==0:

            jsub = 0          # one for each exposure to analyze, increment after submitting exposure processing job
            jcomp = 0         # one for each exposure to analyze, increment after a job has finished running
            pc = 0  # for looping through job channel list
        
            #-------------------------------------------------------------------------#
            #-------------------------------------------------------------------------#
            # For each channel, check status of last job & submit a new one accordingly
            while jcomp<ntorun:   # -> until we have nexp exposures to download...
                part = chans[pc]      # select channel
                rootLogger.info("Checking partition "+part+" queue, "+str(ntorun-jcomp)+" out of "+str(ntorun)+" jobs left")
                # - Check last exposure processing (proc) job submitted (l)
                #              (job status?)
                # - Check next exposure up in the queue (n)
                #               (last job not running/pending?) = r(eady)
                #              (all files exist?)  = e(xist)
                partition_ind = set(np.where(expstr['partition'][torun]==part)[0])  # indices for this partition
                unsubmitted_ind = set(np.where(expstr['submitted'][torun]==0)[0])   # indices for unsubmitted exposures (to get tsub for this submisison round)
                submitted_ind = set(np.where(expstr['submitted'][torun]==1)[0])     # indices for submitted exposures (to get lsub for last submission round)
                sub = list(partition_ind & submitted_ind)          # indices of exp_proc jobs submitted to channel (last_sub)
                unsub = list(partition_ind & unsubmitted_ind)                       # indices of exp_proc jobs unsubmitted to channel (this_sub)
        
                #---------------------------------------------------------------------------------#
                # -> get index  & info of LAST EXPOSURE who's proc job was submitted to channel (l)
                if len(sub)>0:
                    sub = np.sort(sub)
                    lsub = [sub[-1]] # index of LAST EXPOSURE submitted
                    lsub_jname = expstr['jobname'][torun[lsub]][0]
                    if maxjobs<5:                                                  # let rest if not looping through many jobs
                        rootLogger.info("wait a sec before checking last proc job...")
                        time.sleep(sleep_time)
                    lsub_jstat,lsub_jid = sacct_cmd(lsub_jname,["state","jobid"])  # check the status of last proc job
                    expstr['jobstatus'][torun[lsub]] = lsub_jstat
                    expstr['jobid'][torun[lsub]] = lsub_jid
                    if lsub_jstat=="RUNNING" or lsub_jstat=="REQUEUED" or lsub_jstat=="PENDING": # If last exposure job is still going, don't submit another.
                        r==False
                    else:                                                                        # Else, allow next one up to be submitted
                        r==True
                        if lsub_jstat=="COMPLETED":                                              # if last exposure proc job is marked as done, clean up
                            jcomp+=1
                            compstat = complete_job(expstr,torun,lsub,lsub_jname,basedir,expdir,rootLogger) # remove all that leftover exposure stuff
                            rootLogger.info("exposure "+str(expstr[torun[lsub]]['fluxfile'])+" processing completed! Submit next exposure up in "+str(part)+" partition")
                        else: rootLogger.info("exposure "+str(expstr[torun[lsub]]['fluxfile'])+" processing stopped! Submit next exposure up in "+str(part)+" partition")
                else: r = True                                                                    # If no job submitted to partition, go ahead and submit one.
                #---------------------------------------------------------------------#
                # -> get index & info of NEXT EXPOSURE to check & submit to channel (n)
                if len(unsub)>0:
                    nsub = [np.sort(unsub)[0]] # index of NEXT EXPOSURE submitted
                    rootLogger.info("Checking out next exposure up, "+str(expstr['fluxfile'][torun[nsub]])+"...")
                    # do all (n)'s exposure files exist?          = e
                    exp_files = [expstr['fluxfile'][torun[nsub]],expstr['wtfile'][torun[nsub]],expstr['maskfile'][torun[nsub]]]
                    exp_bools = [os.path.exists(expdir+exp_files[i][0].strip()) for i in range(len(exp_files))]
                    if False in exp_bools: e = False
                    else: e = True
                    if r and e: # if channel is ready for a new job to be submitted AND all files exist, submit the job!
                        rootLogger.info("-- Submitting next job up in processing queue to "+str(part)+" channel --")
                        cmd = expstr['cmd'][torun[nsub]]
                        partition = expstr['partition'][torun[nsub]].split("_")[0]
                        rootLogger.info("--Submit Processing Job for Exposure "+str(expstr['fluxfile'][torun[nsub]].split("/")[-1])+"--")
                        # -> Write exposure processing job script to file
                        job_name = 'nsc_meas_'+str(logtime)+'_'+str(jsub) #jsub, not jbsb!!!!!
                        job_file = write_jscript(job_name,jsub,partition,[cmd],outfiledir,"david.nidever@montana.edu")
                        # -> Submit exposure processing job to slurm queue
                        os.system("sbatch "+str(job_file))
                        expstr['submitted'][torun[nsub]] = True
                        jsub+=1
                        rootLogger.info("Job "+job_name+"  submitted to "+partition+" partition")
                        expstr['jobname'][torun[nsub]] = job_name
                    else: rootLogger.info("Partition "+part+" has a job still running, or no exposures up, moving to next channel.")
                #------------------------------#
                # - Set up for next partition -
                pc+=1
                pc = pc-(pc//nchantot)*nchantot # stay within range of partitions list
                # - Save job structure -
                Table(expstr).write(runfile,overwrite=True)
                if ntorun<10 or maxjobs<5:
                   print("sleeping for ten sec before checking next partition...")
                   time.sleep(sleep_time)
                else: time.sleep(sleep_time//10)

        # --------------------------------------------------------------------------------------------------------------------
        # If downloading exposures -------------------------------------------------------------------------------------------
        else:

            jb = 0            # one for each exposure to analyze, increment after adding exposure to queuee
            jsub = 0          # one for each exposure to analyze, increment after submitting exposure processing job
            jexp = 0          # one for each exposure download job submitted
            jcomp = 0         # one for each exposure to analyze, increment after a job has finished running
            #ljinds = []       # hold exposures for download
            part_counter = 0  # for looping through job channel list
            jb_flag = 0       # jb_flag = 1 when all exposures processing jobs that have been submitted are done running
            while jb_flag==0: # jb_flag = 1 when (jcomp=jsub)
        
                #-------------------------------------------------------#
                # -- Start New Batch of Jobs for Exposure Downloading --#
                #-------------------------------------------------------#
                if jb<ntorun: rootLogger.info("--Starting New Batch of Jobs ("+str(nexp)+" exposures)--")
                dload_inds = [] # indices of exposures to be submitted in the next downloading (dload) job
                proc_inds = []  # indices of exposures to be submitted in the next round of processing (proc) jobs
                #print("New start, dload_inds = ",dload_inds,", proc_inds = ",proc_inds)
                #----------------------------------------------------------------------------------
                #---------------------------------------------------------#
                # -- Check Exposure Processing Stream, looping through -- #
                # -- channels until nexp exposures are retrieved       -- #
                #---------------------------------------------------------#
                # For each channel, until nexp exposures are in dload queue:
                pc = part_counter     # channel counter
                nji = len(dload_inds) # exposure download counter
                while nji<nexp:       # -> until we have nexp exposures to download...
                    part = chans[pc]
                    rootLogger.info("Checking partition "+part+" queue")
                    # - Check last exposure processing job submitted (l)
                    #              (job status?)
                    # - Check next exposure up in the queue (n)
                    #              (exposures unreleased?)                   = nu
                    #              (been submitted as an exposure download?) = d
                    #              (exposure download job done?)             = x
                    #              (all files exist?)                        = e
                    partition_ind = set(np.where(expstr['partition'][torun]==part)[0])  # indices for this partition
                    unsubmitted_ind = set(np.where(expstr['submitted'][torun]==0)[0])   # indices for unsubmitted exposures (to get tsub for this submisison round)
                    submitted_ind = set(np.where(expstr['submitted'][torun]==1)[0])     # indices for submitted exposures (to get lsub for last submission round)
                    esubmitted_ind = set(np.where(expstr['esubmitted'][torun]==1)[0])   # indices for exposure-download-submitted exposures
                    released_ind = set(np.where(expstr['ur'][torun]!=1)[0])             # indices for released exposures
                    sub = list((partition_ind & submitted_ind) & released_ind)          # indices of exp_proc jobs submitted to channel (last_sub)
                    unsub = list(partition_ind & unsubmitted_ind)                       # indices of exp_proc jobs unsubmitted to channel (this_sub)
                    dloading = list(set(unsub) & esubmitted_ind)                        # indices of downloading exposures (unsubmitted for proc, but submitted for dload)
                    #print("sub, unsub, dloading = ",sub,unsub,dloading)
                    # get index of last exposure who's proc job was submitted to channel (l)
                    if len(sub)>0:
                        sub = np.sort(sub)
                        lsub = [sub[-1]]
                        #lsub = [sub[np.where(expstr['jobid'][torun[sub]] == max(expstr['jobid'][torun[sub]]))[0][0]]] # get by finding highest job id of exposure proc jobs submitted to channel (the latest one submitted)
                        #print("lsub = ",lsub)
                        lsub_jname = expstr['jobname'][torun[lsub]][0]
                        lsub_jstat_logged = expstr['jobstatus'][torun[lsub]][0].strip()
                        #print("lsub_jname = ",lsub_jname,lsub_jstat_logged)
                        # status of (l)'s proc job?                   = lsub_jstat
                        if lsub_jstat_logged!="COMPLETED" and lsub_jstat_logged!="FAILED": # if last job hasnt already been confirmed as complete
                            if maxjobs<5:
                                print("wait a sec before checking last proc job...")       #     check the status of last proc job
                                time.sleep(10)
                            lsub_jstat,lsub_jid = sacct_cmd(lsub_jname,["state","jobid"])
                            expstr['jobstatus'][torun[lsub]] = lsub_jstat
                            expstr['jobid'][torun[lsub]] = lsub_jid
                            if lsub_jstat=="COMPLETED":                                    #     and if  (l) is completed:   c
                                jcomp+=1
                                compstat = complete_job(expstr,torun,lsub,lsub_jname,
                                                        basedir,expdir,rootLogger)         #         remove all that leftover exposure shit
                                rootLogger.info("exposure "+str(expstr[torun[lsub]]['fluxfile'])+" processing completed!")
                        else:                                                          # else, if the last job is complete,
                            lsub_jstat = lsub_jstat_logged                             # just grab the necessary info
                            lsub_jid = expstr['jobid'][torun[lsub]]
                        #print("last exp proc job in this partition is ",lsub_jstat)
                    else: lsub = False
                    # get index of next exposure to check & submit to channel (n)
                    if len(unsub)>0:
                        nsub = [np.sort(unsub)[0]]
                        # exposures unreleased?                       = nu
                        if expstr[torun[nsub]]['ur']==1 or os.path.exists(basedir+"unavails/"+str(expstr['rawname'][torun[nsub]])+"_unavail.txt"):
                            nu = True
                            expstr[torun[nsub]]['ur'] = 1
                            expstr[torun[nsub]]['submitted'] = 1
                            expstr[torun[nsub]]['esubmitted'] = 1
                        else: nu = False
                        # has (n) been submitted in an exp dload job? = d
                        if nsub in dloading: d=True
                        else: d=False
                        # is that exposure dload job complete?        = x
                        if d==True:
                            nexp_jname = expstr['exp_jobname'][torun[nsub]][0]
                            #print("next exposure up's dload job is being checked...",nexp_jname)
                            nexp_jstat,nexp_jid = sacct_cmd(nexp_jname,["state","jobid"])
                            if nexp_jstat != "RUNNING": x=True
                            else: x=False
                        else: x=False
                        # do all (n)'s exposure files exist?          = e
                        exp_files = [expstr['fluxfile'][torun[nsub]],expstr['wtfile'][torun[nsub]],expstr['maskfile'][torun[nsub]]]
                        expx_files = []
                        for expp in exp_files:
                            expp = expp[0]
                            #print("expp = ",expp)
                            vers = expp.split(".")[0].split("_")[-1]
                            if vers in ['MT1', 'a1', 'lg9', 'ls10', 'ls9', 'v1', 'v2', 'v3', 'v4']:
                                versx = re.split('(\d+)',vers)[0]+"x"
                                exppx = expp.split(vers)[0]+versx+expp.split(vers)[-1]
                                #print("old base = ",expp," and new base = ",exppx)
                                expx_files.append(exppx)
                        #print("exp_files = ",exp_files)
                        exp_bools = [os.path.exists(expdir+exp_files[i][0].strip()) for i in range(len(exp_files))]
                        expx_bools = [os.path.exists(expdir+expx_files[i][0].strip()) for i in range(len(expx_files))]
                        if (False in exp_bools) and (False in expx_bools): e=False
                        elif (False in exp_bools) and (False not in expx_bools):
                            e = True
                            expstr['x'][torun[nsub]] = True
                        else: e=True
                    else:
                        nsub=False
                # - Do the appropriate shit with that information:
                   #if (n) and no (l): #------------------------------------case at start of run
                    if nsub and (not lsub):
                        #print("nsub but no lsub")
                        rootLogger.info("Checking out next exposure up, "+str(expstr['fluxfile'][torun[nsub]])+"...")
                        if e:                                               # if all (n)'s files exist:              e(existance of exposure files check)
                            rootLogger.info("All files exist, no need to download!")
                            if nsub not in proc_inds:                       #     if (n) not in proc                 -  (p*)
                                proc_inds.append(nsub[0])                      #         add (n) to proc queue          -
                                jb+=1                                       #         increment jb                   -
                                rootLogger.info("Adding new exposure "+str(expstr['fluxfile'][torun[nsub]])+" to processing queue")
                            else:                                           #     else, if (n) in proc:
                                rootLogger.info("full channel loop, moving on")
                                nji=nexp                                    #                                        x(exit channel loop)
                        elif (not e) and (not nu):                          # elif (n)'s missing files & not unreleased:
                            if (not d) and (nsub not in dload_inds):        #     if (n) not sub for or in dload:    -- (d*)
                                dload_inds.append(nsub[0])                     #        add (n) to dload queue          --
                                nji+=1                                      #        increment nji                   --
                                rootLogger.info("Adding new exposure "+str(expstr['fluxfile'][torun[nsub]])+" to download queue")
                            elif (not d) and (nsub in dload_inds):          #     elif (n) not submitted for but in dload
                                rootLogger.info("full channel loop, moving on")
                                nji=nexp                                    #                                        x
                            else:                                           #     else, if (n) has been sub for dload:
                                if x and (not nu):                          #         if that job is done (& files not there):
                                    expstr['esubmitted'][torun[nsub]]=0     #             re-assign exp for dload
                                    rootLogger.info("Exposure "+str(expstr['fluxfile'][torun[nsub]])+" download failed, re-queue")
                                else:                                       #         else, if that job is not done: it's either still rolling or not submitted yet, so I think we can skip on over to the next partition?
                                    rootLogger.info("exposure "+str(expstr['fluxfile'][torun[nsub]])+" still downloading or unreleased")
                        elif (not e) and nu:                                               # else, if (n) unreleased, skip to next channel
                            rootLogger.info("exposure "+str(expstr['fluxfile'][torun[nsub]])+" unreleased")
                   #elif (l) and (n): #-------------------------------------case while running
                    elif nsub and lsub:
                        #print("lsub and nsub")
                        rootLogger.info("Checking out next exposure up, "+str(expstr['fluxfile'][torun[nsub]])+"...")
                        if e:                                               # if all (n)'s files exist:              e
                            rootLogger.info("All files exist, no need to download!")
                            if (lsub_jstat=="RUNNING") or (lsub_jstat=="PENDING") or (lsub_jstat=="REQUEUED"):
                                                                            #     if (l) is run/pend/rq:             r(un check)
                                rootLogger.info("Exposure "+str(expstr[torun[lsub]]['fluxfile'])+" (job id="+lsub_jid
                                    +", job name="+lsub_jname+") is still "+lsub_jstat+", check next partition")
                                if maxjobs<5: time.sleep(10)        #         if more exposures, next channel
                            else:                                           #     else, if l is not running/pending/requeued:
                                #print("last proc job status = ",lsub_jstat)
                                if nsub not in proc_inds:                   #         if (n) not in proc:            -  (p*)
                                    proc_inds.append(nsub[0])                  #             add (n) to proc            -
                                    jb+=1                                   #             increment jb               -
                                    rootLogger.info("Adding new exposure "+str(expstr['fluxfile'][torun[nsub]])+" to processing queue")
                                else:                                       #         else, if (n) in proc:
                                    rootLogger.info("full channel loop, moving on")
                                    nji=nexp                                #                                        x
                        elif (not e) and (not nu):                          # elif (n)'s files don't exist and n isn't unreleased:
                            if (not d) and (nsub not in dload_inds):        #     if (n) not sub for or in dload:    -- (d*)
                                dload_inds.append(nsub[0])                     #         add (n) to dload               --
                                nji+=1                                      #         increment nji                  --
                                rootLogger.info("Adding new exposure "+str(expstr['fluxfile'][torun[nsub]])+" to download queue")
                            elif (not d) and (nsub in dload_inds):          #     elif (n) not sub for but in dload:
                                rootLogger.info("full channel loop, moving on")
                                nji = nexp                                  #                                        x
                            else:                                           #     else, if (n) has been submitted for dload:
                                if x and (not nu):                          #         if that job is done (& files not there):
                                    expstr['esubmitted'][torun[nsub]]=0      #             re-assign exp for dload
                                    rootLogger.info("Exposure "+str(expstr['fluxfile'][torun[nsub]])+" download failed, re-queue")
                                else:                                       #         else, if that job is not done: it's either still rolling or not submitted yet, so I think we can skip on over to the next partition?
                                    rootLogger.info("exposure "+str(expstr['fluxfile'][torun[nsub]])+" still downloading or unreleased")
                        elif (not e) and nu:                                # else, add nothing and skip to next partition
                            rootLogger.info("exposure "+str(expstr['fluxfile'][torun[nsub]])+" unreleased")
                   #elif (l) and no (n): #----------------------------------case at end of run
                    elif (not nsub) and lsub:
                        #print("lsub but not nsub")
                        if (lsub_jstat=="RUNNING") or (lsub_jstat=="PENDING") or (lsub_jstat=="REQUEUED"):
                                                                           # if (l) is run/pend/rq:                 r(un check)
                            rootLogger.info("Exposure "+str(expstr[torun[lsub]]['fluxfile'])+" (job id="+lsub_jid
                                +", job name="+lsub_jname+") is still "+lsub_jstat+", check next partition")
                            if maxjobs<5: time.sleep(10)           #     if more exposures, next channel
                        else:                                              # else,
                            if jcomp==jsub:                                #     if all exposure proc jobs are done:
                                rootLogger.info("All exposures processed, ending now!")
                                jb_flag = 1                                #                                         f(inished processing exposures)
                                nji = nexp                                 #                                         x
                                sys.exit(0)
                            # else:                                        #     else, skip to the next partition
                # - Set up for next partition
                    print("checked partition, proc_inds = ",proc_inds,", dload_inds = ",dload_inds)
                    pc+=1
                    pc = pc-(pc//nchantot)*nchantot # stay within range of partitions list
                    if len(proc_inds)==maxjobs or len(dload_inds)==nexp: nji = nexp
                    #print("nji = ",nji,", jb = ",jb,", jcomp = ",jcomp)
                    if maxjobs<5:
                        print("sleeping for ten...")
                        time.sleep(10)
                #----------------------------------------------------------------------------------
                #------------------------------------------------
                # -- Update Exposure Processing Job Streams  -- #
                # -- (Submit jobs in proc queue)             -- #
                #------------------------------------------------
                # - Loop through exposures in processing queue, Write/Submit jobs
                if len(proc_inds)>0:
                    rootLogger.info("-- Submitting jobs in processing queue...")
                    for jbsb in proc_inds:
                        if os.path.exists(basedir+"unavails/"+str(expstr['rawname'][torun[jbsb]])+"_unavail.txt"):
                            expstr['submitted'][torun[jbsb]] = True
                            expstr['ur'][torun[jbsb]] = True
                            rootLogger.info("Exposure unreleased, skip")
                        else:
                            cmd = expstr['cmd'][torun[jbsb]]
                            if expstr['x'][torun[jbsb]]: cmd = cmd+" --x"
                            partition = expstr['partition'][torun[jbsb]].split("_")[0]
                            rootLogger.info("--Submit Processing Job for Exposure "+str(expstr['fluxfile'][torun[jbsb]].split("/")[-1])+"--")
                            # -> Write exposure processing job script to file
                            job_name = 'nsc_meas_'+str(logtime)+'_'+str(jsub) #jsub, not jbsb!!!!!
                            job_file = write_jscript(job_name,jsub,partition,[cmd],outfiledir,email="katiefasbender@montana.edu",cpupt=2,mempc="4G",rtime="2-00:00:00")
                            #job_file=write_jscript(job_name,partition,cmd,outfiledir,single=True)
                            # -> Submit exposure processing job to slurm queue
                            os.system("sbatch "+str(job_file))
                            expstr['submitted'][torun[jbsb]] = True
                            jsub+=1
                            rootLogger.info("Job "+job_name+"  submitted to "+partition+" partition")
                            expstr['jobname'][torun[jbsb]] = job_name
                #----------------------------------------------------------------------------------
                #-------------------------------------------------#
                # -- Check and Update Exposure Download Stream -- #
                # -- (Submit jobs in the dload queue)          -- #
                #-------------------------------------------------#
                # - Get last exposure job to be submitted (le)
                esub_ind = np.where(expstr['esubmitted'][torun]==1)[0]   # indices for exposure-dload-submitted exposures
                if len(list(esub_ind))>0:                                # if there were exposure dload jobs submitted, get indices of latest one
                    lesub = list(esub_ind[np.where(expstr['exp'][torun[esub_ind]] == max(expstr['exp'][torun[esub_ind]]))[0]]) # get by finding highest job id of exp dload jobs submitted (the latest one submitted)
                    #print("lesub = ",lesub)
                    # get status of that job
                    lesub_jname = expstr['exp_jobname'][torun[lesub]][0].strip()
                    if maxjobs<5: time.sleep(5)
                    lesub_jstat,lesub_jid = sacct_cmd(lesub_jname,["state","jobid"]) # check last exposure dload job
                    expstr['exp_jobid'][torun[lesub]] = lesub_jid
                    expstr['exp_jobstatus'][torun[lesub]] = lesub_jstat
                    #print("last exp download job = ",lesub_jstat)
                else: lesub = False
                # - Get next exposures up in dload queue (ne)
                if len(dload_inds)>0:
                    nesub = list(dload_inds)
                    # - Write exp_job script to file
                    exp_job_name = 'exp_dload_'+str(logtime)+'_'+str(jexp) # set up an exp dload job
                    expstr['exp_jobname'][torun[dload_inds]] = exp_job_name
                    #print(" -> Exposure names",expstr['fluxfile'][torun[tjob_inds]])
                    exp_cmds = list(expstr['exp_cmd'][torun[dload_inds]])
                    exp_job_file = write_jscript(exp_job_name,jexp,"priority",exp_cmds,outfiledir,email="katiefasbender@montana.edu",cpupt=2,mempc="2G",rtime="00:10:00",parallel=True)
                    #exp_job_file = write_jscript(exp_job_name,"priority",exp_cmds,outfiledir,single=False)
                else: nesub = False
                #print("nesub,lesub = ",nesub,lesub)
                # - Do the relevant shit with that information:
                # if (le): -------------------------------------------------------# if there was a dload job submitted:
                if lesub:                                                  # if (le):
                    # while (le) is still running:......
                    while (lesub_jstat=="RUNNING" or lesub_jstat=="PENDING" or lesub_jstat=="REQUEUED"):
                        rootLogger.info("waiting for last exposure download job to finish...")
                        expstr['exp_jobstatus'][torun[lesub]] = lesub_jstat
                        rootLogger.info("Job id="+lesub_jid+" is still "+lesub_jstat+", sleepin for "+str(10)+" sec...")
                        time.sleep(10)
                        lesub_jstat = sacct_cmd(lesub_jname,["state"]).strip()
                    # if last exposure dload job completed, get some info about it
                    expstr['exp_jobstatus'][torun[lesub]] = lesub_jstat        #     wait until (le) is done, then get info about it
                    if lesub_jstat=="COMPLETED":
                        le_info = sacct_cmd(lesub_jname,["cputimeraw","maxrss"],c=True,m=True)
                        expstr['exp_cputime'][torun[lesub]] = le_info[0]
                        expstr['exp_maxrss'][torun[lesub]] = le_info[1]
                        rootLogger.info("Exposure download job completed, go ahead!")
                #----------------------------------------------------------------------------------
                # -- Check status of last exp_job, wait until not runnning/pending/requeued --
                    # -- Check contents of exposure repo, sleep until empty enough --
                    #exp_rest = 0
                    #while exp_rest==0:
                        #rootLogger.info("----checking contents of exposure repository----")
                        #exp_mem = int((subprocess.getoutput("du "+expdir+" --si -hb")).split('\t')[0])
                        #if exp_mem>92274688000: #if /exposures repo has more than 86GB of files, sleep and check again.
                            #rootLogger.info("----/exposures repo has "+str(exp_mem)+" bytes of files, sleep for a few seconds----")
                            #time.sleep(20)
                        #else: exp_rest = 1
                #----------------------------------------------------------------------------------
                # if (ne): -------------------------------------------------------# if there is a new dload job:
                if nesub:                                                  # if (ne):
                    # -- Submit new exp_job script to slurm queue --
                    rootLogger.info("--Submitting new Exposure Download job--")
                    os.system("sbatch "+str(exp_job_file))                 #     submit (ne)
                    expstr['esubmitted'][torun[nesub]] = True
                    expstr['exp'][torun[nesub]] = jexp
                    jexp+=1
                    rootLogger.info("Job "+exp_job_name+" submitted")
                    if 1==1: #maxjobs<10: 
                        print("sleepin for ten...[sec]...")
                        time.sleep(10)
                #----------------------------------------------------------------------------------
                part_counter = pc  # save where we left off looping through the partition channels
                # -- Save job structure --
                Table(expstr).write(runfile,overwrite=True)
                if ntorun<10 or maxjobs<5:
                   print("sleeping for ten sec before checking next partition...")
                   time.sleep(sleep_time)
                else: time.sleep(sleep_time//10)
