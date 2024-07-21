import os
import numpy as np
from pwd import getpwuid
from grp import getgrgid
import string
import random
from datetime import datetime,timezone
import traceback
from astropy.table import Table
from astropy.io import fits
from dlnpyutils import utils as dln
import subprocess

# Some slurm helper functions

def genkey(n=20):
    characters = string.ascii_lowercase + string.digits
    key =  ''.join(random.choice(characters) for i in range(n))
    return key

# submitting tempest measure jobs
def submit(tasks,label,nodes=1,cpus=64,version='v4',account='priority-davidnidever',
           partition='priority',staggertime=60,host='tempest_group',shared=True,
           walltime='12-00:00:00',notification=False,memory=7500,numpy_num_threads=2,
           precommands=None,nosubmit=False,slurmroot='/tmp',slurmdir=None,
           verbose=True,logger=None):
    """
    Submit a bunch of jobs

    tasks : table
      Table with the information on the tasks.  Must have columns of:
        cmd, outfile, errfile, dir (optional)

    """

    if logger is None:
        logger = dln.basiclogger()

    username = getpwuid(os.getuid())[0]
    if slurmdir is None:
        slurmdir = os.path.join(slurmroot,username,'slurm')
        if os.path.exists(slurmdir)==False:
            os.makedirs(slurmdir)

    # Generate unique key
    key = genkey()
    if verbose:
        logger.info('key = '+key)
    # make sure it doesn't exist yet

    # job directory
    jobdir = os.path.join(slurmdir,label,key)
    if os.path.exists(jobdir)==False:
        os.makedirs(jobdir)

    # Figure out number of tasks per cpu
    ntasks = len(tasks)

    # Add column to tasks table
    if isinstance(tasks,Table)==False:
        tasks = Table(tasks)
    tasks['task'] = -1
    tasks['node'] = -1
    
    # Parcel out the tasks to the nodes
    #  -loop over all the nodes until we've
    #    exhausted all of the tasks
    count = 0
    while (count<ntasks):
        for i in range(nodes):
            node = i+1
            if count>=ntasks: break
            tasks['task'][count] = count+1
            tasks['node'][count] = node
            count += 1

    # Make the main slurm script and the scripts for each node
            
    # Node loop
    node_index = dln.create_index(tasks['node'])    
    tasknum = 0
    inventory = []
    for i in range(len(node_index['value'])):
        ind = node_index['index'][node_index['lo'][i]:node_index['hi'][i]+1]
        ind = np.sort(ind)
        nind = len(ind)
        node = node_index['value'][i]
        nodefile = 'node{:02d}.slurm'.format(node)
        nodename = 'node{:02d}'.format(node)

        # Separate the tasks for each node
        nodetasks = tasks[ind]
        nodetasksfile = os.path.join(jobdir,'node{:02d}_tasks.fits'.format(node))
        if verbose:
            logger.info('Writing node tasks to '+nodetasksfile)
        nodetasks.write(nodetasksfile)
        
        lines = []
        # Create the lines
        lines = []
        lines += ['#!/bin/bash']
        lines += ['# Auto-generated '+datetime.now().ctime()+' -- '+label+' ['+nodefile+']']
        if account is not None:
            lines += ['#SBATCH --account='+account]
        if partition is not None:
            lines += ['#SBATCH --partition='+partition]
        lines += ['#SBATCH --nodes=1']
        lines += ['#SBATCH --ntasks='+str(cpus)]
        lines += ['#SBATCH --mem-per-cpu='+str(memory)]
        lines += ['#SBATCH --cpus-per-task=1']
        lines += ['#SBATCH --time='+walltime]
        lines += ['#SBATCH --job-name='+label]
        lines += ['#SBATCH --output='+label+'_%j.out']
        lines += ['#SBATCH --err='+label+'_%j.err']
        lines += ['# ------------------------------------------------------------------------------']
        lines += ['export OMP_NUM_THREADS=2']
        lines += ['export OPENBLAS_NUM_THREADS=2']
        lines += ['export MKL_NUM_THREADS=2']
        lines += ['export VECLIB_MAXIMUM_THREADS=2']
        lines += ['export NUMEXPR_NUM_THREADS=2']
        lines += ['# ------------------------------------------------------------------------------']
        lines += ['export CLUSTER=1']
        lines += [' ']
        # nscjob_manager tasksfile version --host host --staggertime 60 --njobs 64
        cmd = 'nscjob_manager '+nodetasksfile+' '+version+' --host '+host
        cmd += ' --staggertime '+str(staggertime)+' --njobs '+str(cpus)
        nodemanagerscript = os.path.join(jobdir,'node{:02d}_manager.sh'.format(node))
        dln.writelines(nodemanagerscript,[cmd])
        lines += ['source '+nodemanagerscript+' &']
        lines += ['wait']
        lines += ['echo "Done"']
        if verbose:
            logger.info('Writing '+os.path.join(jobdir,nodefile))
        dln.writelines(os.path.join(jobdir,nodefile),lines)


    # Create the "master" slurm file
    masterfile = label+'.slurm'
    lines = []
    lines += ['#!/bin/bash']
    lines += ['# Auto-generated '+datetime.now().ctime()+' ['+masterfile+']']
    if account is not None:
        lines += ['#SBATCH --account='+account]
    if partition is not None:
        lines += ['#SBATCH --partition='+partition]
    lines += ['#SBATCH --nodes=1']
    #lines += ['#SBATCH --ntasks='+str(nproc)]
    lines += ['#SBATCH --mem-per-cpu='+str(memory)]
    lines += ['#SBATCH --cpus-per-task=1']
    lines += ['#SBATCH --time='+walltime]
    lines += ['#SBATCH --array=1-'+str(nodes)]
    lines += ['#SBATCH --job-name='+label]
    lines += ['#SBATCH --output='+label+'_%A[%a].out']
    lines += ['#SBATCH --err='+label+'_%A[%a].err']
    lines += ['# ------------------------------------------------------------------------------']
    lines += ['export OMP_NUM_THREADS=2']
    lines += ['export OPENBLAS_NUM_THREADS=2']
    lines += ['export MKL_NUM_THREADS=2']
    lines += ['export VECLIB_MAXIMUM_THREADS=2']
    lines += ['export NUMEXPR_NUM_THREADS=2']
    lines += ['# ------------------------------------------------------------------------------']
    lines += ['SBATCH_NODE=$( printf "%02d']
    lines += ['" "$SLURM_ARRAY_TASK_ID" )']
    lines += ['source '+jobdir+'/node${SBATCH_NODE}.slurm']
    if verbose:
        logger.info('Writing '+os.path.join(jobdir,masterfile))
    dln.writelines(os.path.join(jobdir,masterfile),lines)

    # Write the number of tasks
    dln.writelines(os.path.join(jobdir,label+'.ntasks'),ntasks)

    # Write the inventory file
    dln.writelines(os.path.join(jobdir,label+'_inventory.txt'),inventory)

    # Write the tasks list
    tasks.write(os.path.join(jobdir,label+'_tasks.fits'),overwrite=True)
    # Write the list of logfiles
    dln.writelines(os.path.join(jobdir,label+'_logs.txt'),list(tasks['outfile']))

    if nosubmit==False:
        # Now submit the job
        logger.info('Submitting '+os.path.join(jobdir,masterfile))
        # Change to the job directory, because that's where the outputs will go
        curdir = os.path.abspath(os.curdir)
        os.chdir(jobdir)
        try:
            res = subprocess.check_output(['sbatch',os.path.join(jobdir,masterfile)])
            success = True

            if type(res)==bytes: res = res.decode()
            res = res.strip()  # remove \n
            if verbose:
                logger.info(res)
            # Get jobid
            #  Submitted batch job 5937773 on cluster notchpeak
            res = res.split('\n')[-1]
            jobid = res.split()[3]
            if verbose:
                logger.info('jobid = '+jobid)            
        except:
            logger.info('Submitting job to SLURM failed with sbatch.')
            success = False
            tb = traceback.format_exc()
            logger.info(tb)
            jobid = -1

        return slurmdir,key,jobid

    return slurmdir,key
