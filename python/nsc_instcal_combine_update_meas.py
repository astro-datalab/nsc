#!/usr/bin/env python

# Update measurement catalogs using the broken up measid/objectid lists
# from nsc_instcal_combine_breakup_idstr.py

import os
import sys
import numpy as np
import time
from dlnpyutils import utils as dln, db
from astropy.table import Table
from astropy.io import fits
import sqlite3
import socket
from argparse import ArgumentParser
import logging
from glob import glob
import subprocess

def exposure_update(exposure,redo=False):
    """ Update the measurement table using the broken up measid/objectid lists."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    iddir = '/data0/dnidever/nsc/instcal/v3/idstr/'
    version = 'v3'

    # Load the exposures table
    print('Loading exposure table')
    expcat = fits.getdata('/net/dl2/dnidever/nsc/instcal/'+version+'/lists/nsc_v3_exposure_table.fits.gz',1)

    # Make sure it's a list
    if type(exposure) is str: exposure=[exposure]

    # Match exposures to exposure catalog
    eind1,eind2 = dln.match(expcat['EXPOSURE'],exposure)

    print('Updating measid for '+str(len(exposure))+' exposures')

    # Loop over files
    for i,exp in enumerate(exposure):
        print(str(i+1)+' '+exp)

        t0 = time.time()
        instcode = expcat['INSTRUMENT'][eind1[i]]
        dateobs = expcat['DATEOBS'][eind1[i]]
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        expdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/'+instcode+'/'+night+'/'+exp
        edir = iddir+instcode+'/'+night+'/'+exp+'/'   # local directory for ID files
        #outdir = edir
        outdir = expdir

        # Check that the directory exists
        if os.path.exists(expdir) is False:
            print(expdir+' NOT FOUND')
            continue

        # Check output file
        measfile = outdir+'/'+exp+'_meas.fits'
        if (os.path.exists(measfile+'.gz')) & (redo is False):
            print(measfile+'.gz already exists.  Skipping')
            continue

        # Log file
        #------------------
        # format is EXPOSURE_measure_update.DATETIME.log
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
        logfile = outdir+'/'+exp+'_measure_update.'+logtime+'.log'
        if os.path.exists(logfile): os.remove(logfile)

        # Set up logging to screen and logfile
        logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
        if logging.getLogger().hasHandlers() is True: rootLogger.handlers=[]   # remove all handlers
        rootLogger = logging.getLogger()
        fileHandler = logging.FileHandler(logfile)
        fileHandler.setFormatter(logFormatter)
        rootLogger.addHandler(fileHandler)
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        rootLogger.addHandler(consoleHandler)
        rootLogger.setLevel(logging.NOTSET)

        rootLogger.info('Adding objectID for measurement catalogs for exposure = '+exp)
        rootLogger.info("expdir = "+expdir)
        rootLogger.info("host = "+host)
        rootLogger.info(" ")

        #  Load the exposure and metadata files
        metafile = expdir+'/'+exp+'_meta.fits'
        meta = Table.read(metafile,1)
        nmeta = len(meta)
        chstr = Table.read(metafile,2)
        rootLogger.info('KLUDGE!!!  Changing /dl1 filenames to /dl2 filenames')
        cols = ['EXPDIR','FILENAME','MEASFILE']
        for c in cols:
            f = np.char.array(chstr[c]).decode()
            f = np.char.array(f).replace('/dl1/users/dnidever/','/dl2/dnidever/')
            chstr[c] = f
        nchips = len(chstr)

        # Get "good" chips, astrometrically calibrated
        gdch,ngdch,bdch,nbdch = dln.where(chstr['NGAIAMATCH']>0,comp=True)
        if nbdch>0:
            rootLogger.info(str(nbdch)+' chips were not astrometrically calibrated')

        measdtype = np.dtype([('MEASID', 'S50'), ('OBJECTID', 'S50'), ('EXPOSURE', 'S50'), ('CCDNUM', '>i2'), ('FILTER', 'S2'), ('MJD', '>f8'), ('X', '>f4'),
                              ('Y', '>f4'), ('RA', '>f8'), ('RAERR', '>f4'), ('DEC', '>f8'), ('DECERR', '>f4'), ('MAG_AUTO', '>f4'), ('MAGERR_AUTO', '>f4'),
                              ('MAG_APER1', '>f4'), ('MAGERR_APER1', '>f4'), ('MAG_APER2', '>f4'), ('MAGERR_APER2', '>f4'), ('MAG_APER4', '>f4'),
                              ('MAGERR_APER4', '>f4'), ('MAG_APER8', '>f4'), ('MAGERR_APER8', '>f4'), ('KRON_RADIUS', '>f4'), ('ASEMI', '>f4'), ('ASEMIERR', '>f4'),
                              ('BSEMI', '>f4'), ('BSEMIERR', '>f4'), ('THETA', '>f4'), ('THETAERR', '>f4'), ('FWHM', '>f4'), ('FLAGS', '>i2'), ('CLASS_STAR', '>f4')])

        # Load and concatenate the meas catalogs
        chstr['MEAS_INDEX'] = -1   # keep track of where each chip catalog starts
        count = 0
        meas = Table(data=np.zeros(int(np.sum(chstr['NMEAS'][gdch])),dtype=measdtype))
        rootLogger.info('Loading and concatenating the chip measurement catalogs')
        for j in range(ngdch):
            jch = gdch[j]
            meas1 = Table.read(chstr['MEASFILE'][jch].strip(),1)   # load chip meas catalog
            nmeas1 = len(meas1)
            meas[count:count+nmeas1] = meas1
            chstr['MEAS_INDEX'][jch] = count
            count += nmeas1
        measid = np.char.array(meas['MEASID']).strip().decode()
        nmeas = len(meas)
        rootLogger.info(str(nmeas)+' measurements')


        # Look for the id files
        files = glob(edir+exp+'__*.npy')
        nfiles = len(files)
        rootLogger.info(str(nfiles)+' ID files to load')

        # Loop over ID files and load them up
        df = np.dtype([('measid',np.str,50),('objectid',np.str,50)])
        idcat = np.zeros(10000,dtype=df)
        count = 0
        for k in range(nfiles):
            idcat1 = np.load(files[k])
            nidcat1 = len(idcat1)
            # Add more elements
            if count+nidcat1 > len(idcat):
                idcat = dln.add_elements(idcat,np.maximum(100000,nidcat1))
            # Stuff in the data
            idcat[count:count+nidcat1] = idcat1
            count += nidcat1
        # Trim extra elements
        if len(idcat)>count: idcat=idcat[0:count]
        rootLogger.info('IDs for '+str(len(idcat))+' measurements')

        # Match up with measid
        ind1,ind2 = dln.match(idcat['measid'],measid)
        nmatch = len(ind1)
        rootLogger.info('Matches for '+str(nmatch)+' measurements')
        if nmatch>0:
            meas['OBJECTID'][ind2] = idcat['objectid'][ind1] 

        # Checking for missing objectid
        ind,nind = dln.where(np.char.array(meas['OBJECTID']).strip().decode() == '')
        # There can be missing/orphaned measurements at healpix boundaries in crowded
        # regions when the DBSCAN eps is different.  But there should be very few of these.
        # At this point, let's allow this to pass
        if nind>0:
            rootLogger.info('WARNING: '+str(nind)+' measurements are missing OBJECTIDs')
        if ((nmeas>=20000) & (nind>20)) | ((nmeas<20000) & (nind>3)):
            rootLogger.info('More missing OBJECTIDs than currently allowed.')
            dln.writelines(outdir+'/'+exp+'_meas.ERROR','')
            continue

        # Output the updated measurement catalog
        #  Writing a single FITS file is much faster than many small ones
        # could put it in /data0 but db01 won't be able to access that
        rootLogger.info('Writing final measurement catalog to '+measfile)
        meas.write(measfile,overwrite=True)
        if os.path.exists(measfile+'.gz'): os.remove(measfile+'.gz')
        ret = subprocess.call(['gzip',measfile])    # compress final catalog

        # Update the meta file as well, need to update the /dl2 filenames
        metafile = outdir+'/'+exp+'_meta.fits'        
        rootLogger.info('Updating meta file '+metafile)
        meta.write(metafile,overwrite=True)
        hdulist = fits.open(metafile)
        hdu = fits.table_to_hdu(chstr)
        hdulist.append(hdu)
        hdulist.writeto(metafile,overwrite=True)
        hdulist.close()

        # Create a file saying that the files were updated okay.
        #dln.writelines(expdir+'/'+exp+'_meas.updated','')
        dln.writelines(outdir+'/'+exp+'_meas.updated','')
        # Remove meas.ERROR, if it exists
        if os.path.exists(outdir+'/'+exp+'_meas.ERROR'):
            os.remove(outdir+'/'+exp+'_meas.ERROR')
        
        rootLogger.info('dt = '+str(time.time()-t0)+' sec.')

    print('dt = %6.1f sec.' % (time.time()-t00))

if __name__ == "__main__":
    parser = ArgumentParser(description='Update measid in exposure.')
    parser.add_argument('exposure', type=str, nargs=1, help='Exposure name')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this exposure')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    exposure = args.exposure[0]
    redo = args.redo

    # Input is a list
    if exposure[0]=='@':
        listfile = exposure[1:]
        if os.path.exists(listfile): 
            exposure = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    # Update the measurement files
    exposure_update(exposure,redo=redo)
