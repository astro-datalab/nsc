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
        if single==True: fh.writelines("#SBATCH --ntasks=1\n")   # for running in parallel
        else: fh.writelines("#SBATCH --ntasks="+str(len(cmds))+"\n")
        fh.writelines("#SBATCH --cpus-per-task 2\n")             # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem-per-cpu=3000\n")           # memory, set --mem with care!!!!! refer to hyalite quickstart guide
        if single==True: fh.writelines("#SBATCH --time=48:00:00\n")
        else: fh.writelines("#SBATCH --time=00:05:00\n")         # Maximum job run time
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

def sacct_cmd(job_name,outputs):
    '''parses the output of a sacct command, returning specified information
    Arguments:
    ----------
    job_name (str)
            you know what this is
    outputs (str list)
            a list of information to get with the sacct command
            see sacct manual page for options
    Returns:
    --------
    outputs (list)
            a list of the sacct outputs specified
    '''
    if len(outputs)>1: spref = "sacct -n -P --delimiter=',' --format "
    else: spref = "sacct -n -X --format "
    scommand = (''.join([spref]+[i+"," for i in outputs]))[:-1]+" --name "+job_name
    job_info = (subprocess.getoutput(scommand).split("\n")[0]).split(",")
    jinfo = [i.strip() for i in job_info]
    if len(outputs)>1 and len(jinfo)>1: return(jinfo)
    elif len(outputs)>1 and len(jinfo)<=1: return(["" for ou in outputs])
    else: return(jinfo[0])

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
    nexp = int(args.nexp)                   # number of exposures to download at once
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories 
    # -tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"   # location of operations
    expdir = basedir+"exposures/"                        # where exposures will be
    outfiledir = basedir+"outfiles/"                     # a place for the job files
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
                             ('cmd',str,1000),('cputime',str,100),('maxrss',str,100),
                             ('exp_jobname',str,100),('exp_jobid',str,100),('exp_jobstatus',str,100),
                             ('exp_cmd',str,1000),('exp_cputime',str,100),('exp_maxrss',str,100)])
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
        outfile = basedir+instrument+'/'+night+'/'+baseroot+'/'+baseroot+'_1.fits' # outfile for first chip
        outlogfile = basedir+instrument+'/'+night+'/'+baseroot+'/'+baseroot+'.log' # logfile for complete exposure
        expstr['outfile'][i] = outfile
        expstr['logfile'][i] = outlogfile

        # Does the output file exist?
        expstr['done'][i] = False
        if (os.path.exists(outfile) and os.path.exists(outlogfile)) or (expstr['fluxfile'][i]==expstr['wtfile'][i]):
            expstr['done'][i] = True

        # If no outfile exists or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            #if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
            expstr['cmd'][i] = 'python '+basedir+'nsc_instcal_meas.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version
            expstr['exp_cmd'][i] = 'python '+basedir+'get_exposures.py '+rawname+' '+fluxfile+' '+wtfile+' '+maskfile+' '+expdir
            expstr['torun'][i] = True
        # If outfile exists and no redo:
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
    if ntorun!=0: njpar = (ntorun//nchan)+1
    rootLogger.info(str(ntorun)+' exposures to process, '+str(maxjobs)+' at a time on '+str(npar)+' slurm partitions.')

    # Split exposures evenly among defined "partitions" as run with nsc_meas_wrapper
    chans = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,nchan)],maxjobs) # partions -> chans
    nchans = len(chans)
    print("job channels = ",chans)
    expstr['partition'][torun] = [chans[(i-maxjobs*(i//maxjobs))] for i in range(ntorun)]


    # Start submitting jobs
    #----------------------
    runfile = basedir+'lists/nsc_instcal_measure_main.'+logtime+'_run.fits'
    expstr['jobid']="-99.99"
    Table(expstr).write(runfile)
    jb_flag = 0
    jb = 0            # one for each exposure to analyze
    ljinds = []       # hold exposures for download
    part_counter = 0  # for looping through channel list
    while jb_flag==0: # jb_flag = 1 when (jb = ntorun)

        # -- Start New Batch of Jobs for Downloading --
        rootLogger.info("--Starting New Batch of Jobs (3 exposures)--")
        tjob_inds = []      # hold expstr indices for exposures to be included in This batch (get using tsub)
        ljob_inds = ljinds  # hold expstr indices for exposures included in Last batch (get using lsub)
        #print("lastjob_inds = ",ljob_inds)

        # -- Loop through channels util nexp exposures are retrieved --
        pc = part_counter
        #print("pc = ",pc)
        nji = len(tjob_inds)
        while nji<nexp:
            part = chans[pc]
            rootLogger.info("Checking partition "+part+" queue")
            partition_ind = set(np.where(expstr['partition'][torun]==part)[0]) # indices for this partition
            unsubmitted_ind = set(np.where(expstr['submitted'][torun]==0)[0])  # indices for unsubmitted exposures (to get tsub for this submisison round)
            submitted_ind = set(np.where(expstr['submitted'][torun]==1)[0])    # indices for submitted exposures (to get lsub for last submission round)

            #------------------------------------------------------------------------------
            # Check Last Exposure Submitted to Channel (update tjob_inds & ljob_inds)
            #------------------------------------------------------------------------------
            # - Get index, id, & status of last exposure submitted to this channel
            last_sub = list(partition_ind & submitted_ind)
            if len(last_sub)==0: lsub=np.sort(list(partition_ind))[0]
            else: lsub=np.sort(last_sub)[-1]
            lj_name = expstr[torun[lsub]]['jobname']
            if (len(last_sub)>1 and expstr[torun[lsub]]['ur']!=1): lj_status,lj_id = sacct_cmd(lj_name,["state","jobid"])
            elif (expstr[torun[lsub]]['ur']!=1): lj_status,lj_id = ["UR","-99.99"] # last exposure is unreleased
            else: lj_status,lj_id = ["NONE","-99.99"]                              # no jobs have been submitted
            expstr[torun[lsub]]['jobstatus'] = lj_status
            expstr[torun[lsub]]['jobid'] = lj_id
            # - If last exposure is still running: skip to next partition!
            if (lj_status=="RUNNING" or lj_status=="PENDING" or lj_status=="REQUEUED"):
                rootLogger.info("Exposure "+str(expstr[torun[lsub]]['fluxfile'])
                    +" (job id="+lj_id+", job name="+lj_name
                    +") is still "+lj_status+", check next partition")
                time.sleep(1)
            # - Else: add next exposure to batch
            else:
                # -> if last exposure was completed, get some info about it
                if lj_status=="COMPLETED":
                    lj_cputime,lj_maxrss = sacct_cmd(lj_name,["cputimeraw","maxrss"])
                    expstr['done'] = True
                    expstr['cputime'][torun[lsub]] = lj_cputime
                    expstr['maxrss'][torun[lsub]] = lj_maxrss
                # -> remove downloaded files!
                if lj_status!="NONE":
                    if os.path.exists(expdir+expstr['fluxfile'][torun[lsub]].strip()):
                        os.remove(expdir+expstr['fluxfile'][torun[lsub]].strip())
                    if os.path.exists(expdir+expstr['maskfile'][torun[lsub]].strip()):
                        os.remove(expdir+expstr['maskfile'][torun[lsub]].strip())
                    if os.path.exists(expdir+expstr['wtfile'][torun[lsub]].strip()):
                        os.remove(expdir+expstr['wtfile'][torun[lsub]].strip())
                    os.system("rm -rf "+basedir+"/tmp/"+(expstr['fluxfile'][torun[lsub]].strip()).split(".")[0]+".*")
                # -> get indices & info of next exposure to submit in channel queue
                this_sub = list(partition_ind & unsubmitted_ind)
                if len(this_sub)!=0:
                    tsub = np.sort(this_sub)[0] # jbsub -> nesub
                    # check to see whether all exposure files exist
                    exp_files = [expstr['fluxfile'][torun[tsub]],expstr['wtfile'][torun[tsub]],expstr['maskfile'][torun[tsub]]]
                    exp_bools = [os.path.exists(expdir+exp_files[i].strip()) for i in range(len(exp_files))]
                    #print("exp_files = ",exp_files)
                    #print("exp_bools = ",exp_bools)
                    # --> If exposure not already in lineup, add to appropriate batch
                    if (tsub not in tjob_inds) and (tsub not in ljob_inds):
                        if (False in exp_bools): tjob_inds.append(tsub) # add exposure to this batch (submitted for processing next round)
                        else: ljob_inds.append(tsub)                # add exposure to last batch (about to be submitted for processing)
                        rootLogger.info("exposure "+expstr['fluxfile'][torun[tsub]]+" added to queue")
                        jb+=1
                    else: time.sleep(1)
            # -Set up for next partition
            pc+=1
            pc = pc-(pc//nchans)*nchans # stay within range of partitions list
            nji = len(tjob_inds) # check to see if we've reached nexp exposures
            #print("nji = ",nji,", jb = ",jb)
            if jb==ntorun-1:
                nji=nexp  # to keep within range of exposure list
            #------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------
        # Update Exposure Downloading Job Stream
        #----------------------------------------------------------------------------------
        # -- Write exp_job script to file --
        #print("Indices of next 3 exposures to download: = ",tjob_inds)
        exp_job_name = 'exp_dload_'+str(logtime)+'_'+str(jb//nexp)
        expstr['exp_jobname'][torun[tjob_inds]] = exp_job_name
        #print(" -> Exposure names",expstr['fluxfile'][torun[tjob_inds]])
        exp_cmds = list(expstr['exp_cmd'][torun[tjob_inds]])
        exp_job_file=write_jscript(exp_job_name,"priority",exp_cmds,outfiledir,single=False)

        # -- Check status of last exp_job, wait until not runnning/pending/requeued --
        #print("again, last job indices = ",ljob_inds)
        if len(ljob_inds)!=0: # if there was a previously submitted download job,
            lexps = ljob_inds[:3] # check the exposures being downloaded for completeness
            rootLogger.info("checking on status of last exposure-download job")
            le_name = expstr['exp_jobname'][torun[lexps]][0]
            #print("le_name = ",le_name)
            time.sleep(5)
            # check last exposure dload job
            le_stat,le_id = sacct_cmd(le_name,["state","jobid"])
            #print("last exp_download jobID = ",le_id,", status = ",le_stat)
            expstr['exp_jobid'][torun[lexps]] = le_id
            while (le_stat=="RUNNING" or le_stat=="PENDING" or le_stat=="REQUEUED"):
                expstr['exp_jobstatus'][torun[lexps]] = le_stat
                rootLogger.info("Job id="+le_id+" is still "+le_stat+", sleepin for "+str(sleep_time))
                time.sleep(sleep_time)
                le_stat = sacct_cmd(le_name,["state"]).strip()
                #print("le_stat = ",le_stat)
            # if last exposure dload job completed, get some info about it
            expstr['exp_jobstatus'][torun[lexps]] = le_stat
            if le_stat=="COMPLETED":
                le_info = sacct_cmd(le_name,["cputimeraw","maxrss"])
                expstr['exp_cputime'][torun[lexps]] = le_info[0]
                expstr['exp_maxrss'][torun[lexps]] = le_info[1]
                rootLogger.info("Exposure download job completed, go ahead!")

        # -- Check contents of exposure repo, sleep until empty enough --

        # -- Submit new exp_job script to slurm queue --
        rootLogger.info("--Submitting new Exposure Download job--")
        os.system("sbatch "+str(exp_job_file))
        rootLogger.info("Job "+exp_job_name+" submitted")
        #----------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------
        # Update Exposure Processing Job Streams
        #----------------------------------------------------------------------------------
        # -- Loop through last batch's exposures, Write/Submit jobs --
        for jbsb in ljob_inds:
            # - Check to see if the exposures were correctly downloaded
            # - if they have, create and submit the job!
            jb_files = [expstr['fluxfile'][torun[jbsb]],expstr['wtfile'][torun[jbsb]],expstr['maskfile'][torun[jbsb]]]
            jb_bools = [os.path.exists(expdir+jb_files[i].strip()) for i in range(len(jb_files))]
            if False not in jb_bools: 
                cmd = expstr['cmd'][torun[jbsb]]
                partition = expstr['partition'][torun[jbsb]].split("_")[0]
                rootLogger.info("--Submit Processing Job for Exposure "+jb_files[0].split("/")[-1])
                # -> Write exposure processing job script to file
                job_name = 'nsc_meas_'+str(logtime)+'_'+str(jb)
                job_file=write_jscript(job_name,partition,cmd,outfiledir,single=True)
                # -> Submit exposure processing job to slurm queue
                os.system("sbatch "+str(job_file))
                expstr['submitted'][torun[jbsb]] = True
                rootLogger.info("Job "+job_name+"  submitted to "+part+" partition")
                expstr['jobname'][torun[jbsb]] = job_name
            # - if not, ?????????????? just fuckin skip i guess 
            elif os.path.exists(basedir+"unavails/"+str(expstr['rawname'][torun[jbsb]])+"_unavail.txt"):
                expstr['submitted'][torun[jbsb]] = True
                rootLogger.info("Exposure unreleased, skip")
                expstr['ur'][torun[jbsb]] = True
            else: rootLogger.info("Exposure files not all downloaded, skip for now")
        #----------------------------------------------------------------------------------

        # -- Save job structure --
        part_counter = pc  # save where we left off looping through the partition channels
        ljinds = tjob_inds # put exposures in the next batch to be submitted for processing        
        Table(expstr).write(runfile,overwrite=True)
        if jb==ntorun: jb_flag=1
