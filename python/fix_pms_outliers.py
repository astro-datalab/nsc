#!/usr/bin/env python

# Update measurement catalogs using the broken up measid/objectid lists
# from nsc_instcal_combine_breakup_idstr.py

import os
import sys
import numpy as np
import shutil
import time
from dlnpyutils import utils as dln, coords, db
from astropy.table import Table
from astropy.io import fits
import sqlite3
import socket
from argparse import ArgumentParser
import logging
from glob import glob
#import subprocess
#import healpy as hp
#import tempfile
#import psycopg2 as pq
#import psutil
from dl import queryClient as qc


def fix_pms(objectid):
    """ Correct the proper motions in the healpix object catalog."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    version = 'v3'
    radeg = np.float64(180.00) / np.pi

    meas = qc.query(sql="select * from nsc_dr2.meas where objectid='"+objectid+"'",fmt='table')
    nmeas = len(meas)
    mnra = np.median(meas['ra'].data)
    mndec = np.median(meas['dec'].data)

    lim = 20.0  # 50.0
    gd, = np.where( (np.abs(meas['ra'].data-mnra)/np.cos(np.deg2rad(mndec))*3600 < lim) &
                     (np.abs(meas['dec'].data-mndec)*3600 < lim))
    ngd = len(gd)
    nbd = nmeas-ngd
    print('bad measurements '+str(nbd))
    if nbd==0:
        return None
    meas = meas[gd]

    raerr = np.array(meas['raerr']*1e3,np.float64)    # milli arcsec
    ra = np.array(meas['ra'],np.float64)
    ra -= np.mean(ra)
    ra *= 3600*1e3 * np.cos(mndec/radeg)     # convert to true angle, milli arcsec
    t = meas['mjd'].copy()
    t -= np.mean(t)
    t /= 365.2425                          # convert to year
    # Calculate robust slope
    try:
        pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
    except:
        print('problem')
        import pdb; pdb.set_trace()

    decerr = np.array(meas['decerr']*1e3,np.float64)   # milli arcsec
    dec = np.array(meas['dec'],np.float64)
    dec -= np.mean(dec)
    dec *= 3600*1e3                         # convert to milli arcsec
    # Calculate robust slope
    try:
        pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
    except:
        print('problem')
        import pdb; pdb.set_trace()

    out = [pmra,pmraerr,pmdec,pmdecerr]

    return out

if __name__ == "__main__":
    parser = ArgumentParser(description='Fix pms in healpix object catalogs.')
    parser.add_argument('catfile', type=str, nargs=1, help='Catalog file')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    if os.path.exists(args.catfile[0]) is False:
        raise ValueError(args.catfile[0]+' NOT FOUND')

    # Save the corrected file
    catfile = args.catfile[0]
    outdir = os.path.dirname(catfile)
    outbase = os.path.basename(catfile)
    outfile = outdir+'/'+outbase.replace('.fits.gz','_corrected.fits')

    cat = Table.read(catfile)
    # make sure column names are lowercase
    colnames = cat.colnames
    for n in colnames:
        cat[n].name = n.lower()
    print(cat.colnames)
    objectid = cat['id'].data
    nobj = len(objectid)


    cat['old_pmra'] = cat['pmra'].copy()
    cat['old_pmraerr'] = cat['pmraerr'].copy()
    cat['old_pmdec'] = cat['pmdec'].copy()
    cat['old_pmdecerr'] = cat['pmdecerr'].copy()


    # Fix the pms in healpix object catalogs
    for i in range(nobj):
        objid = objectid[i].astype(str).strip()
        print(str(i+1)+' '+objid)
        out = fix_pms(objid)
        if out is not None:
            print('  OLD: %10.2f %10.2f' % (cat['pmra'][i],cat['pmdec'][i]))
            print('  NEW: %10.2f %10.2f' % (out[0],out[2]))
            cat['pmra'][i] = out[0]
            cat['pmraerr'][i] = out[1]
            cat['pmdec'][i] = out[2]
            cat['pmdecerr'][i] = out[3]


    # Save the corrected file
    print('Saving corrected file to '+outfile)
    cat.write(outfile)
