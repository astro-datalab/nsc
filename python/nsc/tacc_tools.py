import os
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.io import fits
#from dlnpyutils import utils as dln
import shutil
import time
from datetime import datetime
import subprocess
from . import utils

def make_transfer_list(n=10000,checkprev=True):
    """
    Make a list of exposures to transfer from NOIRLab to TACC.
    """

    listdir = '/net/dl2/dnidever/nsc/instcal/v4/lists/'
    tab = Table.read(listdir+'decam_instcal_list_exptime10sec_20240727_left2transfer.fits.gz')
    #tab = Table.read(listdir+'decam_instcal_list_exptime10sec_20240727_left.fits.gz')
    #tab['fluxfile'] = np.char.array(tab['fluxfile']).astype(str).replace('/net/archive','/net/mss1/archive')
    #tab['wtfile'] = np.char.array(tab['wtfile']).astype(str).replace('/net/archive','/net/mss1/archive')
    #tab['maskfile'] = np.char.array(tab['maskfile']).astype(str).replace('/net/archive','/net/mss1/archive')
    #tab = Table.read(listdir+'decam_instcal_list_exptime10sec_20240714.fits.gz')
    #tab = Table.read(listdir+'r16avails_decam_instcal_list.fits.gz')

    ## Remove exposures that are done
    #done = utils.readlines(listdir+'/exposures_done_corral_20240714.txt')
    #done_exposure = [os.path.basename(d) for d in done]
    #_,ind1,ind2 = np.intersect1d(tab['base'],done_exposure,return_indices=True)
    #tab.remove_rows(ind1)

    print('Making TACC image transfer list')

    # Checking previous lists
    if checkprev:
        # Check any existing lists
        files = glob(listdir+'transfer*list_*.lst')
        files.sort()
        print('Found',len(files),'previous lists')

        # Load the previous lists
        prevlines = []
        for i in range(len(files)):
            print(files[i])
            lines = utils.readlines(files[i])
            prevlines += lines

        # Match them to FLUXFILE
        _,ind1,ind2 = np.intersect1d(prevlines,tab['fluxfile'],return_indices=True)
        if len(ind1)>0:
            print(len(ind1),' exposures in previous lists')
            # Delete them from the list
            del tab[ind2]
        print(len(tab),' exposures left')

    if len(tab)<n:
        print('Only',len(tab),' remain')
        n = len(tab)

    # Start the list of files
    print('Making list of',n,'exposures')
    lines = []
    for i in range(n):
        fluxfile = tab['fluxfile'][i]
        wtfile = tab['wtfile'][i]
        maskfile = tab['maskfile'][i]
        print(i,os.path.basename(fluxfile))
        if os.path.exists(fluxfile) and os.path.exists(wtfile) and os.path.exists(maskfile):
            lines += [fluxfile,wtfile,maskfile]

    # Write the list to a file
    tstamp = datetime.now().strftime('%Y%m%d%H%M%S')
    outfile = 'transfer'+str(n)+'list_'+tstamp+'.lst'
    utils.writelines(listdir+outfile,lines)
    print('List written to '+listdir+outfile)

def make_transfer_list_tempest(n=10000,checkprev=True):
    """
    Make a list of exposures to transfer from TACC to Tempest.
    """

    listdir = '/scratch1/09970/dnidever/nsc/instcal/v4/lists/'
    #listdir = '/home/x51j468/group/nsc/instcal/v4/lists/'
    #tab = Table.read(listdir+'decam_instcal_list_exptime10sec_20240727_left.fits.gz')
    tab = Table.read(listdir+'decam_instcal_list_exptime10sec_corrupted_fwhm2_20250119_left.fits.gz')

    ## Remove exposures that are done
    #done = utils.readlines(listdir+'/exposures_done_corral_20240714.txt')
    #done_exposure = [os.path.basename(d) for d in done]
    #_,ind1,ind2 = np.intersect1d(tab['base'],done_exposure,return_indices=True)
    #tab.remove_rows(ind1)

    print('Making TACC image transfer list')

    # Checking previous lists
    if checkprev:
        # Check any existing lists
        files = glob(listdir+'transfer*list_*.lst')
        files.sort()
        print('Found',len(files),'previous lists')

        # Load the previous lists
        prevlines = []
        for i in range(len(files)):
            print(files[i])
            lines = utils.readlines(files[i])
            prevlines += lines

        # Match them to FLUXFILE
        _,ind1,ind2 = np.intersect1d(prevlines,tab['fluxfile'],return_indices=True)
        if len(ind1)>0:
            print(len(ind1),' exposures in previous lists')
            # Delete them from the list
            del tab[ind2]
        print(len(tab),' exposures left')

    if len(tab)<n:
        print('Only',len(tab),' remain')
        n = len(tab)

    # Fix image filenames on TACC
    basedir = '/scratch1/09970/dnidever/nsc/instcal/v4/'
    imlines = utils.readlines(basedir+'images/allimages.lst')
    imlines = [basedir+'images/'+f[1:] for f in imlines]  # make absolute
    imlines = np.array(imlines)
    imbase = np.array([os.path.basename(f) for f in imlines])
    tab['fluxfile'] = tab['fluxfile'].astype((str,300))
    fluxfiles = tab['fluxfile']
    fluxbase = np.array([os.path.basename(f) for f in fluxfiles])
    _,ind1,ind2 = np.intersect1d(fluxbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' fluxfile names')
        tab['fluxfile'][ind1] = imlines[ind2]
    tab['wtfile'] = tab['wtfile'].astype((str,300))
    wtfiles = tab['wtfile']
    wtbase = np.array([os.path.basename(w) for w in wtfiles])
    _,ind1,ind2 = np.intersect1d(wtbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' wtfile names')
        tab['wtfile'][ind1] = imlines[ind2]
    tab['maskfile'] = tab['maskfile'].astype((str,300))
    maskfiles = tab['maskfile']
    maskbase = np.array([os.path.basename(m) for m in maskfiles])
    _,ind1,ind2 = np.intersect1d(maskbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' maskfile names')
        tab['maskfile'][ind1] = imlines[ind2]

    # Start the list of files
    print('Making list of',n,'exposures')
    lines = []
    for i in range(n):
        fluxfile = tab['fluxfile'][i]
        wtfile = tab['wtfile'][i]
        maskfile = tab['maskfile'][i]
        # Fix paths for TACC
        print(i,os.path.basename(fluxfile))
        if os.path.exists(fluxfile) and os.path.exists(wtfile) and os.path.exists(maskfile):
            lines += [fluxfile,wtfile,maskfile]

    # Write the list to a file
    tstamp = datetime.now().strftime('%Y%m%d%H%M%S')
    outfile = 'transfer'+str(n)+'list_'+tstamp+'.lst'
    utils.writelines(listdir+outfile,lines)
    print('List written to '+listdir+outfile)

def reorganize_files(stagedate):
    """
    Reorganize images transferred to TACC.
    """

    staging_dir = '/scratch1/09970/dnidever/nsc/instcal/v4/staging/'
    image_dir = '/scratch1/09970/dnidever/nsc/instcal/v4/images/'

    print('Checking staging directory '+os.path.join(staging_dir,stagedate))

    files = glob(os.path.join(staging_dir,stagedate,'*.fits*'))
    files.sort()
    print(len(files),'files found')

    # Move files
    for i in range(len(files)):
        # move file
        base = os.path.basename(files[i])
        print(i,base)
        src = files[i]
        if base[:3] == 'c4d':
            instrument = base.split('_')[0]
            night = '20'+base.split('_')[1]
        else:
            try:
                head = fits.getheader(files[i],0)
            except:
                print('Problem reading',files[i],'skipping')
                continue
            dateobs = head['date-obs']
            instrument = 'c4d'  # assume it's decam
            year = dateobs[:4]
            month = dateobs[5:7]
            day = dateobs[8:10]
            night = year+month+day
        year = night[:4]
        outdir = os.path.join(image_dir,instrument,year,night)
        if os.path.exists(outdir)==False:
            os.makedirs(outdir)
        dst = os.path.join(outdir,base)
        shutil.move(src,dst)

    #import pdb; pdb.set_trace()


def reorganize_files_tempest():
    """
    Reorganize images transferred to TACC.
    """

    staging_dir = '/home/group/davidnidever/nsc/instcal/v4/staging/'
    image_dir = '/home/group/davidnidever/nsc/instcal/v4/images/'

    print('Checking staging directory '+staging_dir)

    files = glob(os.path.join(staging_dir,'*.fits*'))
    files.sort()
    print(len(files),'files found')

    # Move files
    for i in range(len(files)):
        # move file
        base = os.path.basename(files[i])
        print(i,base)
        src = files[i]
        if base[:3] == 'c4d':
            instrument = base.split('_')[0]
            night = '20'+base.split('_')[1]
        else:
            try:
                head = fits.getheader(files[i],0)
            except:
                print('Problem reading',files[i],'skipping')
                continue
            dateobs = head['date-obs']
            instrument = 'c4d'  # assume it's decam
            year = dateobs[:4]
            month = dateobs[5:7]
            day = dateobs[8:10]
            night = year+month+day
        year = night[:4]
        outdir = os.path.join(image_dir,instrument,year,night)
        if os.path.exists(outdir)==False:
            os.makedirs(outdir)
        dst = os.path.join(outdir,base)
        shutil.move(src,dst)

    #import pdb; pdb.set_trace()

def measure_status():
    """
    Check how many exposures have been successfully processed with measurement.
    """

    basedir = '/scratch1/09970/dnidever/nsc/instcal/v4/'
    corraldir = '/corral/projects/NOIRLab/nsc/instcal/v4/'
    imagedir = '/scratch1/09970/dnidever/nsc/instcal/v4/images/'

    # Might be faster to just search for tgz files
    # lfind tgz > tgzfiles

    tab = Table.read(basedir+'lists/decam_instcal_list_exptime10sec_20240714.fits.gz')
    instrument = 'c4d'

    tab['corraldone'] = False
    tab['scratchdone'] = False
    tab['haveimages'] = False

    for i in range(len(tab)):
        base = tab['base'][i]
        dateobs = tab['date_obs'][i]
        night = dateobs[:4]+dateobs[5:7]+dateobs[8:10]
        year = dateobs[:4]
        # Check scratch output files
        outdir = os.path.join(basedir,instrument,year,night,base)
        measfile = os.path.join(outdir,base+'_meas.fits')
        tgzfile = os.path.join(outdir,base+'.tgz')
        if os.path.exists(outdir) and os.path.exists(measfile) and os.path.exists(tgzfile):
            tab['scratchdone'][i] = True
        # Check corral output files
        coutdir = os.path.join(corraldir,instrument,year,night,base)
        cmeasfile = os.path.join(coutdir,base+'_meas.fits')
        ctgzfile = os.path.join(coutdir,base+'.tgz')
        if os.path.exists(coutdir) and os.path.exists(cmeasfile) and os.path.exists(ctgzfile):
            tab['corraldone'][i] = True
        # Check images
        fluxfile = os.path.join(imagedir,instrument,year,night,os.path.basename(tab['fluxfile'][i]))
        wtfile = os.path.join(imagedir,instrument,year,night,os.path.basename(tab['wtfile'][i]))
        maskfile = os.path.join(imagedir,instrument,year,night,os.path.basename(tab['maskfile'][i]))
        if os.path.exists(fluxfile) and os.path.exists(wtfile) and os.path.exists(maskfile):
            tab['haveimages'][i] = True

        print(i,base,tab['corraldone'][i],tab['haveimages'][i],tab['scratchdone'][i])

    corraldone = np.sum(tab['corraldone'])
    haveimages = np.sum(tab['haveimages'])
    scratchdone = np.sum(tab['scratchdone'])
    done = np.sum(tab['corraldone'] | tab['scratchdone'])
    print('corral done = ',corraldone)
    print('scratch done = ',scratchdone)
    print('done = ',done)
    print('have images = ',haveimages)

    #import pdb; pdb.set_trace()

    return tab


    #yeardir = glob(basedir+'20??')
    #expdir = []
    #count = 0
    #for y in range(len(yeardir)):
    #    nightdir = glob(yeardir[y]+'/20??????')
    #    for i in range(len(nightdir)):
    #        edir = glob(nightdir[i]+'/*')
    #        edir = [e for e in edir if os.path.isdir(e)]
    #        print(i,nightdir[i],len(edir))
    #        for j in range(len(edir)):
    #            base = os.path.basename(edir[j])
    #            tarfile = edir[j]+'/'+base+'.tgz'
    #            measfile = edir[j]+'/'+base+'_meas.fits'
    #            headfile = edir[j]+'/'+base+'_header.fits'
    #            if os.path.exists(tarfile) and os.path.exists(measfile) and os.path.exists(headfile):
    #                print(count,j,edir[j],'good')
    #                expdir.append(edir[j])
    #            else:
    #                print(count,j,edir[j],'bad')
    #            count += 1
    #print(len(expdir),' exposures successfully completed measurement')
    #return expdir

def slurmsummary(skey,clobber=False):
    """ Get summary information for a slurm measure job """
    slurmdir = '/scratch1/09970/dnidever/dnidever/slurm/measure'
    sdir = slurmdir+'/'+skey
    measdir = '/scratch1/09970/dnidever/nsc/instcal/v4/c4d/'
    print(sdir)
    if os.path.exists(sdir)==False:
        raise Exception(sdir+' not found')
    tasksfile = sdir+'/measure_tasks.fits'
    if os.path.exists(tasksfile)==False:
        raise FileNotFoundError(tasksfile)
    tasks = Table.read(tasksfile)
    ntasks = len(tasks)
    print(ntasks,'tasks')
    #logsfile = sdir+'/measure_logs.txt'
    #if os.path.exists(logsfile)==False:
    #    raise FileNotFoundError(logsfile)
    #logfiles = utils.readlines(logsfile)
    #print(len(logfiles),'tasks')
    errfile = glob(sdir+'/measure-*.err')
    if len(errfile)>0:
        errfile = errfile[0]
    errmtime = os.path.getmtime(errfile)
    jobid = errfile.split('-')[-1][:-4].strip()
    print('JobID =',jobid)
    sumfile = sdir+'/'+skey+'_'+jobid+'_summary.fits'
    if os.path.exists(sumfile) and clobber==False:
        print(sumfile,'already exists and clobber not set')
        return
    outfile = glob(sdir+'/measure-*.out')
    if len(outfile)>0:
        outfile = outfile[0]
    outmtime = os.path.getmtime(outfile)
    outlines = utils.readlines(outfile)
    # Get "running" and "completed" lines
    rlines = utils.grep(outlines,'running')
    clines = utils.grep(outlines,'completed')
    # Get information for each task
    dt = [('logfile',str,200),('base',str,50),('exists',bool),
          ('ctime',float),('mtime',float),('size',float),('jobjobid',int),
          ('jobstarted',bool),('jobtaskid',int),('jobcompleted',bool),('jobtruncated',bool),
          ('jobelapsed',float),('slurmstart',float),('slurmend',float),('state',str,20),
          ('measfile',str,200),('done',bool)]
    info = np.zeros(ntasks,dtype=np.dtype(dt))
    info['jobtaskid'] = -1
    for i in range(ntasks):
        info['logfile'][i] = tasks['outfile'][i]
        info['base'][i] = tasks['name'][i]
        info['exists'][i] = os.path.exists(tasks['outfile'][i])
        if info['exists'][i]:
            info['size'][i] = os.path.getsize(tasks['outfile'][i])
        info['jobjobid'][i] = i+1
        rline = utils.grep(rlines,' job '+str(i+1)+' ')
        if len(rline)>0:
            taskid = rline[0].split()[2]
            info['jobstarted'][i] = True
            info['jobtaskid'][i] = taskid
            cline = utils.grep(clines,' Job '+str(i+1)+' ')
            if len(cline)>0:
                info['jobcompleted'][i] = len(cline)>0
                info['jobelapsed'][i] = cline[0].split()[-2]
        if info['jobstarted'][i] and info['jobcompleted'][i]==False:
            info['jobtruncated'][i] = True
        if info['exists'][i] and (info['size'][i]>0):
            info['ctime'][i] = os.path.getctime(tasks['outfile'][i])
            info['mtime'][i] = os.path.getmtime(tasks['outfile'][i])
            measfile = tasks['dir'][i]+'/'+tasks['name'][i]+'_meas.fits'
            info['measfile'][i] = measfile
            info['done'][i] = os.path.exists(measfile)
    # Get slurm job related information
    res = subprocess.run(['sacct','-j',jobid,'--format','JobID,JobName,Start,End,State'],capture_output=True)
    out = res.stdout.decode()
    lines = out.split('\n')
    line = lines[2]
    starttimestamp = line.split()[2]
    endtimestamp = line.split()[3]
    starttime = datetime.fromisoformat(starttimestamp).timestamp()
    endtime =  datetime.fromisoformat(endtimestamp).timestamp()
    state = line.split()[4]
    info['slurmstart'] = starttime
    info['slurmend'] = endtime
    info['state'] = state
    print(starttimestamp,endtimestamp,state)
    # Check logfile mtime against the slurm job endtime
    # to see if the tasks were 
    nrun = np.sum(info['exists'])
    print(nrun,'tasks were run')
    ndone = np.sum(info['done'])
    print(ndone,'tasks finished successfully with meas.fits files')
    ntruncated = np.sum(info['jobtruncated'])
    print(ntruncated,'tasks truncated')
    print('Saving summary to',sumfile)
    Table(info).write(sumfile,overwrite=True)
