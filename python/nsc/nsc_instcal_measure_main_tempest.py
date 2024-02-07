#!/usr/bin/env python

# AUTHORS: David Nidever
#          david.nidever@montana.edu
#          Katie Fasbender
#          katiefasbender@montana.edu

# nsc_instcal_meas_main.py will run the NSC Measurements process on all exposures
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


def write_jscript(job_name,partition,cmds,dir,single=True):
    '''writes a SLURM job script to "job_name.sh"
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    partition (str)
            node/partition the job will run on
    cmds (str or str list)
            python command(s) to run 
    dir (str)
            directory for output files
    single (bool, default True)
            if you have a list of commands, set single=False
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
    job_file = dir+job_name+".sh"
    job_num = int(job_name.split("_")[-1])
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        if partition=="priority": fh.writelines("#SBATCH --account=priority-davidnidever\n") #specify account
        fh.writelines("#SBATCH --job-name="+job_name+"\n")       # job name
        fh.writelines("#SBATCH --output="+dir+job_name+".out\n") # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+dir+job_name+".err\n")  # error file
        fh.writelines("#SBATCH --partition="+partition+"\n")     # queue partition to run the job in
        #fh.writelines("#SBATCH --nodes=1\n")
        if single==True: fh.writelines("#SBATCH --ntasks=1\n")   # for running in parallel
        else: fh.writelines("#SBATCH --ntasks="+str(len(cmds))+"\n")
        fh.writelines("#SBATCH --cpus-per-task 2\n")             # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem-per-cpu=9G\n")           # memory, set --mem with care!!!!! refer to hyalite quickstart guide
        if single==True: fh.writelines("#SBATCH --time=36:00:00\n")
        else: fh.writelines("#SBATCH --time=00:05:00\n")         # Maximum job run time
        if job_num % 1000 == 0:
            fh.writelines("#SBATCH --mail-user katiefasbender@montana.edu\n")
            fh.writelines("#SBATCH --mail-type BEGIN\n")
        fh.writelines("module load Anaconda3\n")                 # load anaconda, needed for running python on hyalite/tempest
        fh.writelines("module load GCC\n")                       # load GCC module
        fh.writelines("source activate $HOME/condaenv/\n")       # activate conda environment
        if single: fh.writelines(cmds+"\n")                       # write python command to analyze exposure
        else:                                                    # OR, write python cmds to download exposures
            for cmd in cmds:
                fh.writelines("srun -n 1 "+cmd+" & \n")
            fh.writelines("wait\n")
        fh.writelines("conda deactivate")
    return job_file

def sacct_cmd(job_name,outputs,c=False,m=False):
    '''parses the output of a sacct command, returning specified information
    Arguments:
    ----------
    job_name (str)
            you know what this is
    outputs (str list)
            a list of information to get with the sacct command
            see sacct manual page for options
    c (bool)
            if job is known to be completed, default=False
    m (bool)
            if multiple tasks in the job, default=False
    Returns:
    --------
    outputs (list)
            a list of the sacct outputs specified
    '''
    if len(outputs)>1: spref = "sacct -n -P --delimiter=',' --format "
    else: spref = "sacct -n -X --format "
    scommand = (''.join([spref]+[i+"," for i in outputs]))[:-1]+" --name "+job_name
    if c and not m: jindex=1
    elif c and m: jindex=2
    else: jindex=0
    job_info = (subprocess.getoutput(scommand).split("\n")[jindex]).split(",")
    jinfo = [i.strip() for i in job_info]
    print("sacct cmd = ",scommand)
    print("sacct output = ",jinfo)
    if len(outputs)>1 and len(jinfo)>1: return(jinfo)
    elif len(outputs)>1 and len(jinfo)<=1: return(["" for ou in outputs])
    else: return(jinfo[0])

def complete_job(struct,run_ind,lsub_ind,lsub_jname,base_dir,exp_dir,logger):
    '''
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
    parser.add_argument('--partition',type=str,nargs=1,help='Delimited list of partitions to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposures that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=1, help='Max number of exposures to process at any given time')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list of exposures to use')
    parser.add_argument('--nexp',type=str,nargs=1,default=3,help='Number of exposures to download at once')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)     # NSC version
    redo = args.redo                         # if called, redo = True
    partitions=args.partition[0].split(',')  # the slurm partitions to submit jobs to
    npar=len(partitions)                     # number of slurm partitions
    maxjobs = int(args.maxjobs[0])           # maximum number of jobs to maintain running at any time
    nchan = maxjobs//npar                    # number of job channels per partition cpar -> nchan
    nexp = int(args.nexp[0])                 # number of exposures to download at once
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories
    # -tempest
    ltype = inputlist.split("_")[1]
    if ltype=="a": basedir = "/home/x25h971/nsc/instcal/"+version+"/"    # location of operations
    else: basedir = "/home/group/davidnidever/nsc/instcal/"+version+"/"  # location of extra output
    scriptdir = "/home/x25h971/nsc/instcal/"+version+"/"      # location of scripts
    expdir = basedir+"exposures/"                        # where exposures will be
    outfiledir = basedir+"outfiles/"
    makedir(expdir)
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
                             ('maskfile',str,100),('outfile',str,150),('logfile',str,150),('partition',str,100),
                             ('done',bool),('torun',bool),('ur',bool),('submitted',bool),
                             ('jobname',str,100),('jobid',str,100),('jobstatus',str,100),
                             ('cmd',str,1000),('cputime',str,100),('maxrss',str,100),('esubmitted',bool),
                             ('exp_jobname',str,100),('exp_jobid',str,100),('exp_jobstatus',str,100),
                             ('exp_cmd',str,1000),('exp_cputime',str,100),('exp_maxrss',str,100),('outdirectory',str,500),('exp',int)])
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
        # format on tempest will be basedir+/exposures/
        fluxfile = fluxfile.split('/')[-1]
        wtfile = wtfile.split('/')[-1]
        maskfile = maskfile.split('/')[-1]

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
        outroot =      basedir+instrument+'/'+night+'/'+baseroot
        outfile =      outroot+"/"+baseroot+'_1.fits' # outfile for first chip
        outlogfile =   outroot+"/"+baseroot+'.log' # logfile for complete exposure
        outdirectory = outroot+'.tar.gz'         # tarred and zipped exposure directory
        expstr['outfile'][i] = outfile
        expstr['logfile'][i] = outlogfile
        expstr['outdirectory'][i] = outdirectory

        ## Do the exposure output & logfiles files exist? (old)
        #expstr['done'][i] = False
        #if (os.path.exists(outfile) and os.path.exists(outlogfile)) or (expstr['fluxfile'][i]==expstr['wtfile'][i]):
        #    expstr['done'][i] = True
        #
        # Does the tarred & zipped output exposure directory exist? (new)
        if os.path.exists(outdirectory): expstr['done'][i] = True

        # If exposure is completed or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            #if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
            expstr['cmd'][i] = 'python '+scriptdir+'nsc_instcal_measure_tempest.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version+' '+ltype
            expstr['exp_cmd'][i] = 'python '+scriptdir+'get_exposures.py '+rawname+' '+fluxfile+' '+wtfile+' '+maskfile+' '+expdir+' '+ltype
            expstr['torun'][i] = True
        # If exposure is completed and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False

    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(expstr['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    rootLogger.info(str(ntorun)+" exposures")
    if ntorun == 0:
        rootLogger.info('No exposures to process.')
        sys.exit()
    njchan = 0 #number of jobs per channel (divide evenly) njpar -> njchan
    if ntorun!=0: njchan = (ntorun//nchan)+1
    rootLogger.info(str(ntorun)+' exposures to process, '+str(maxjobs)+' at a time on '+str(npar)+' slurm partitions.')

    # Split exposures evenly among defined "partitions" as run with nsc_meas_wrapper
    chans = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,nchan)],maxjobs) # partions -> chans
    nchantot = len(chans)
    #print("job channels = ",chans)
    expstr['partition'][torun] = [chans[(i-maxjobs*(i//maxjobs))] for i in range(ntorun)]


    # Start submitting jobs
    #----------------------
    runfile = basedir+'lists/nsc_instcal_measure_main.'+logtime+'_run.fits'
    expstr['jobid']="-99.99"
    Table(expstr).write(runfile)
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
                #print("exp_files = ",exp_files)
                exp_bools = [os.path.exists(expdir+exp_files[i][0].strip()) for i in range(len(exp_files))]
                if False in exp_bools: e=False
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
                    partition = expstr['partition'][torun[jbsb]].split("_")[0]
                    rootLogger.info("--Submit Processing Job for Exposure "+str(expstr['fluxfile'][torun[jbsb]].split("/")[-1])+"--")
                    # -> Write exposure processing job script to file
                    job_name = 'nsc_meas_'+str(logtime)+'_'+str(jsub) #jsub, not jbsb!!!!!
                    job_file=write_jscript(job_name,partition,cmd,outfiledir,single=True)
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
            exp_job_file=write_jscript(exp_job_name,"priority",exp_cmds,outfiledir,single=False)
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
           print("wait a sec...")
           time.sleep(15)
        else: time.sleep(1)
