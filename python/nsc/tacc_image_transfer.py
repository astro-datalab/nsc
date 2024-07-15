import os
import numpy as np
from glob import glob
from astropy.table import Table
from dlnpyutils import utils as dln
import shutil
import time
from datetime import datetime

def make_transfer_list(n=20000):
    """
    Make a list of exposures to transfer from NOIRLab to TACC.
    """

    listdir = '/net/dl2/dnidever/nsc/instcal/v4/lists/'
    tab = Table.read(listdir+'r16avails_decam_instcal_list.fits.gz')

    print('Making TACC image transfer list')

    # Check any existing lists
    files = glob(listdir+'transfer*list_*.lst')
    files.sort()
    print('Found',len(files),'previous lists')

    # Load the previous lists
    prevlines = []
    for i in range(len(files)):
        print(files[i])
        lines = dln.readlines(files[i])
        prevlines += lines

    # Match them to FLUXFILE
    _,ind1,ind2 = np.intersect1d(prevlines,tab['FLUXFILE'],return_indices=True)
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
        fluxfile = tab['FLUXFILE'][i]
        wtfile = tab['WTFILE'][i]
        maskfile = tab['MASKFILE'][i]
        print(i,os.path.basename(fluxfile))
        if os.path.exists(fluxfile) and os.path.exists(wtfile) and os.path.exists(maskfile):
            lines += [fluxfile,wtfile,maskfile]

    # Write the list to a file
    tstamp = datetime.now().strftime('%Y%m%d%H%M%S')
    outfile = 'transfer'+str(n)+'list_'+tstamp+'.lst'
    dln.writelines(listdir+outfile,lines)
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
        instrument = base.split('_')[0]
        night = '20'+base.split('_')[1]
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

    basedir = '/scratch1/09970/dnidever/nsc/instcal/v4/c4d/'
    #basedir = '/corral/projects/NOIRLab/nsc/instcal/v4/c4d/'

    # Might be faster to just search for tgz files
    # lfind tgz > tgzfiles

    yeardir = glob(basedir+'20??')
    expdir = []
    count = 0
    for y in range(len(yeardir)):
        nightdir = glob(yeardir[y]+'/20??????')
        for i in range(len(nightdir)):
            edir = glob(nightdir[i]+'/*')
            edir = [e for e in edir if os.path.isdir(e)]
            print(i,nightdir[i],len(edir))
            for j in range(len(edir)):
                base = os.path.basename(edir[j])
                tarfile = edir[j]+'/'+base+'.tgz'
                measfile = edir[j]+'/'+base+'_meas.fits'
                headfile = edir[j]+'/'+base+'_header.fits'
                if os.path.exists(tarfile) and os.path.exists(measfile) and os.path.exists(headfile):
                    print(count,j,edir[j],'good')
                    expdir.append(edir[j])
                else:
                    print(count,j,edir[j],'bad')
                count += 1
    print(len(expdir),' exposures successfully completed measurement')
    return expdir
