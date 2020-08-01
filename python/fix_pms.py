#!/usr/bin/env python

# Update measurement catalogs using the broken up measid/objectid lists
# from nsc_instcal_combine_breakup_idstr.py

import os
import sys
import numpy as np
import shutil
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
import healpy as hp

def fix_pms(pix):
    """ Correct the proper motions in the healpix object catalogs."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    version = 'v3'

    # Make sure it's a list
    if type(pix) is str: pix=[pix]

    print('Correcting PMs for '+str(len(pix))+' HEALPix')

    # Loop over files
    for i in range(len(pix)):
        t0 = time.time()
        pix1 = pix[i]
        print(str(i+1)+' '+str(pix1))

        hdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/combin/'+str(int(pix1)//1000))+'/'
        objfile = hdir+str(pix1)+'.fits.gz'

        # Check that the object file exists
        if os.path.exists(objfile) is False:
            print(objfile+' NOT FOUND')
            return

        # Check fixed file
        fixedfile = hdir+str(pix1)+'.fixed'
        if os.path.exists(fixedfile) is True:
            print(str(pix1)+' already fixed')
            continue

        # Copy the object file to a temporary file while we're working on it
        tempfile = objfile = hdir+str(pix1)+'_temp.fits.gz'
        if os.path.exists(tempfile): os.remove(tempfile)
        shutil.copyfile(objfile,tempfile)

        # Load the object file
        meta = fits.getdata(objfile,1)
        obj = fits.getdata(objfile,2)
        nobj = len(obj)

        # Save the current PM values
        dtype = np.dtype([('objid',np.str,100),('pmra',float),('pmraerr',float)),('pmdec',float),('pmdecerr',float)])
        old = np.zeros(nobj,dtype)
        old['objid'] = obj['objid']
        old['pmra'] = obj['pmra']
        old['pmraerr'] = obj['pmraerr']
        old['pmdec'] = obj['pmdec']
        old['pmdecerr'] = obj['pmdecerr']

        # Correct the proper motions
        #  negative values are 1.05612245
        #  positive values are 0.94387755

        # Save the old pms
        oldfile = hdir+str(pix1)+'_oldpms.fits'
        print('Saving old PMs to '+oldfile)
        Table(old).write(oldfile)

        # Save the new version of obj
        # Write the output file
        outfile = hdir+str(pix1)+'.fits'
        print('Writing combined catalog to '+outfile)
        if os.path.exists(outfile): os.remove(outfile)
        Table(meta).write(outfile)               # first, summary table
        #  append other fits binary tables
        hdulist = fits.open(outfile)
        hdu = fits.table_to_hdu(Table(obj))        # second, catalog
        hdulist.append(hdu)
        hdulist.writeto(outfile,overwrite=True)
        hdulist.close()
        if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
        ret = subprocess.call(['gzip',outfile])    # compress final catalog

        # Check that the output file looks okay, has the right extensiosn and Nrows
        hd1 = fits.getheader(outfile+'.gz',1)
        hd2 = fits.getheader(outfile+'.gz',2)
        # Problems
        if hd2['NAXIS2'] != nobj:
            print('problem with output file.  Restoring original file')
            os.remove(objfile)
            shutil.move(tempfile,objfile)
        # Everything okay
        else:
            # Delete temporary file
            os.remove(tempfile)
            # Write fixed file
            dln.touch(fixedfile)


    print('dt = %6.1f sec.' % (time.time()-t00))

if __name__ == "__main__":
    parser = ArgumentParser(description='Fix pms in healpix object catalogs.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    pix = args.pix[0]

    # Input is a list
    if pix[0]=='@':
        listfile = pix[1:]
        if os.path.exists(listfile): 
            pix = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    # Fix the pms in healpix object catalogs
    fix_pms(pix)
