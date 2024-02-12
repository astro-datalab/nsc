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
import healpy as hp
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

    
def write_jscript(job_name,jobnum,partition,cmds,dir,email,cpupt=2,mempc="9G",rtime="36:00:00",parallel=False):
    ''' Creates SLURM job script to <dir>/<job_name>.sh to run on <partition>,
        & writes appropriate lines to the file.  Formatted for use on Tempest.
        - Lines starting with #SBATCH are read by Slurm. 
        - Lines starting with ## are comments.
        - All other lines are read by the shell.
    ----------
    Arguments:
    ----------
    job_name (str)
            name of job, to be used for naming job script, 
            output, & error files
    job_num (int)
            an index to count the number of jobs you're submitting
    partition (str)
            node/partition the job will run on
    cmds (str list)
            python command(s) to run 
    dir (str)
            directory for job files (.sh, .out, .err)
    email (str)
            If you're sending yourself email updates (advised)
    cpupt (int)
            CPUs per task (logical CPUs!)
    mempc (str)
            memory per CPU, with Byte unit (B,KB,MB,GB...)
    rtime (str)
            Maximum runtime for job, in format "hh:mm:ss"
    parallel (bool)
            True if running several commands in parallel
            as 1 job (by defining --ntasks>1)
    --------
    Returns:
    --------
    job_file (str)
            job filename the job script is written to <dir>/<job_name>.sh
    '''
    # - Setup -
    job_file = dir+job_name+".sh"
    job_num = int(job_name.split("_")[-1]) # I include an indexing number in each job name for email updates
    ncmd = len(cmds)
    # - Create job file (<dir>/<job_name>.sh) & write appropriate lines -
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")                            # Tell computer what kind of file this is
        if partition=="priority":#----------------------------------If priority partition, use our group priority
            fh.writelines("#SBATCH --account=priority-davidnidever\n")     # account name (priority-davidnidever)
        fh.writelines("#SBATCH --partition="+partition+"\n")      # Slurm partition to run the job in
        fh.writelines("#SBATCH --job-name="+job_name+"\n")        # Job name
        fh.writelines("#SBATCH --output="+dir+job_name+".out\n")  # Job output file (<dir>/<job_name>.out)
        fh.writelines("#SBATCH --error="+dir+job_name+".err\n")   # Job error file (<dir>/<job_name>.err)
        fh.writelines("#SBATCH --cpus-per-task "+str(cpupt)+"\n") # Number of cores to allocate per task
        fh.writelines("#SBATCH --mem-per-cpu="+mempc+"\n")        # Memory per CPU, set with care!!!!!!
        fh.writelines("#SBATCH --time="+rtime+"\n")               # Maximum job run time
        # -- Set up for parallel tasks, if applicable --
        if ncmd==1 or (ncmd>1 and not parallel):#-------------------If one cmd, or multiple cmds running in sequence
            fh.writelines("#SBATCH --ntasks=1\n")
            cmdterm0 = ""
            cmdterm2=""
        elif (ncmd>1 and parallel):#--------------------------------If multiple cmds running in parallel
            fh.writelines("#SBATCH --ntasks="+str(len(cmds))+"\n")
            cmdterm0 = "srun -n 1 "
            cmdterm2 = " & "
        # -- Send me an email when every 1000th job begins --
        if job_num % 1000 == 0:
            fh.writelines("#SBATCH --mail-user "+email+"\n")
            fh.writelines("#SBATCH --mail-type BEGIN\n")
        # -- Set up environment to run python on tempest --
        fh.writelines("module load Anaconda3\n")
        fh.writelines("source activate $HOME/condaenv/\n")
        # -- Write python command(s) --
        for cmd in cmds:
            fh.writelines(cmdterm0+cmd+cmdterm2+"\n")
        if parallel: # don't finish job until all commands have run!
            fh.writelines("wait\n")
        fh.writelines("conda deactivate")
    return job_file # returns name of job script (<dir>/<job_name>.sh)


def sacct_cmd(job_name,fields,c=False,m=False):
    '''parses the output of a sacct command, returning specified information
    Arguments:
    ----------
    job_name (str)
            you know what this is
    fields (str list)
            a list of information to get with the sacct command
            see sacct manual page for options
    c (bool)
            True if job is known to be completed, default=False
    m (bool)
            True if multiple tasks in the job, default=False
    Returns:
    --------
    outputs (list)
            a list of the sacct outputs specified
    '''
    nfields = len(fields)      # Number of outputs to expect
    # -- Setup sacct command format --
    if nfields>1: spref = "sacct -n -P --delimiter=',' --format "
    else: spref = "sacct -n -X --format "
    scommand = spref+(','.join(fields))+" --name "+job_name
    # -- Which line of the sacct cmd output will contain our information --
    if c and m: line_index=2
    elif c and not m: line_index=1
    else: line_index=0
    # -- Run sacct command --
    job_info = (subprocess.getoutput(scommand).split("\n")[line_index]).split(",")
    outputs = [i.strip() for i in job_info] #split into individual field outputs
    print("sacct cmd = ",scommand)
    print("sacct output = ",outputs)
    if nfields>1 and len(outputs)>1: return(outputs)
    elif nfields>1 and len(outputs)<=1: return(["" for ou in fields])
    else: return(outputs[0])

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
    parser.add_argument('--nside',type=str,nargs=1,default=32,help='HEALPix NSIDE to sort list for processing order')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)     # NSC version
    redo = args.redo                         # if called, redo = True
    partitions=args.partition[0].split(',')  # the slurm partitions to submit jobs to
    ns = int(args.nside)                          # HEALPix nside
    npar=len(partitions)                     # number of slurm partitions
    maxjobs = int(args.maxjobs[0])           # maximum number of jobs to maintain running at any time
    nchan = maxjobs//npar                    # number of job channels per partition cpar -> nchan
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories
    # -tempest
    basedir = os.getcwd()+"/"+version+"/"    # location of operations
    scriptdir = basedir                      # location of scripts
    expdir = basedir+"exposures/"            # where exposures will be
    outfiledir = basedir+"outfiles/"         # wher job files will go
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
        needcols = ['INSTRUMENT','FLUXFILE','MASKFILE','WTFILE','DATE_OBS','RA','DEC']
        for n in needcols:
            if n not in lstr.dtype.names:
                raise ValueError('Column '+n+' not in '+inputlist)
    nlstr = dln.size(lstr)
    rootLogger.info(str(nlstr)+' InstCal images')
    gdexp = np.arange(nlstr)
    ngdexp = nlstr

    # Check the exposures
    #--------------------
    rootLogger.info('Checking on the exposures')
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
        # format on tempest will be basedir+/exposures/
        fluxfile = fluxfile.split('/')[-1]
        wtfile = wtfile.split('/')[-1]
        maskfile = maskfile.split('/')[-1]

        expstr['instrument'][i] = instrument
        expstr['ring'+str(ns)][i] = hp.ang2pix(ns,lstr['RA'][gdexp[i]],lstr['DEC'][gdexp[i]],lonlat=True)
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
            expstr['cmd'][i] = 'python '+scriptdir+'nsc_instcal_measure_cca.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version
            expstr['torun'][i] = True
        # If exposure is completed and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False
    # 
    expstr = np.sort(expstr,order="pix"+str(ns))
    
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
    partis = np.array([chans[(i-maxjobs*(i//maxjobs))] for i in range(ntorun)])
    expstr['partition'][torun] = partis


    # Start submitting jobs
    #----------------------
    runfile = basedir+'lists/nsc_instcal_measure_main.'+logtime+'_run.fits'
    expstr['jobid']="-99.99"
    Table(expstr).write(runfile)
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
                job_name = 'nsc_meas_'+str(logtime)+'_'+str(nsub) #jsub, not jbsb!!!!!
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
        
        
