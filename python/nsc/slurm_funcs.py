#!/usr/bin/env python
# Author: Katie Fasbender
#         katiefasbender@montana.edu

#--------
# Imports
#--------
import numpy as np
import os
import subprocess


#----------
# Functions
#----------
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

def write_jscript(job_name,jobnum,partition,cmds,dir,email,cpupt=2,mempc="1G",rtime="02:00:00",parallel=False,dependent=[]):
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
    dependent (str list)
            Job ID of each job that must finish before this one may begin
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
        # -- If job is dependent on completion of another job, make it so! -- 
        if len(dependent)!=0:
            jids = ":".join(dependent)
            fh.writelines("#SBATCH --dependency=afterany:"+jids+"\n") 
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

#----------
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
