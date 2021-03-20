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
from sklearn import linear_model, datasets

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
    #if nbd==0:
    #    return None
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
        #pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)

        # Run RANSAC
        ransac = linear_model.RANSACRegressor()
        ransac.fit(t.reshape(-1,1), ra)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        gdmask = inlier_mask
        pmra_ransac = ransac.estimator_.coef_[0]
        print('  ransac '+str(np.sum(inlier_mask))+' inliers   '+str(np.sum(outlier_mask))+' outliers')

        # Robust, weighted linear with with INLIERS
        #pmra_coef, pmra_coeferr = dln.poly_fit(t[gdmask],ra[gdmask],1,sigma=raerr[gdmask],robust=True,error=True)
        pmra_coef, pmra_coeferr = dln.poly_fit(t,ra,1,sigma=raerr,robust=True,error=True)
        pmra = pmra_coef[0]
        pmraerr = pmra_coeferr[0]
        radiff = ra-dln.poly(t,pmra_coef)
        radiff -= np.median(radiff)
        rasig = dln.mad(radiff)
        # Reject outliers
        gdsig = (np.abs(radiff) < 2.5*rasig) | (np.abs(radiff) < 2.5*raerr)
        print('  '+str(nmeas-np.sum(gdsig))+' 2.5*sigma clip outliers rejected')
        if np.sum(gdsig) < nmeas:
            pmra_coef, pmra_coeferr = dln.poly_fit(t[gdsig],ra[gdsig],1,sigma=raerr[gdsig],robust=True,error=True)
            pmra = pmra_coef[0]
            pmraerr = pmra_coeferr[0]
            rasig = dln.mad(ra-dln.poly(t,pmra_coef))
    except:
        print('problem')
        import pdb; pdb.set_trace()

    decerr = np.array(meas['decerr']*1e3,np.float64)   # milli arcsec
    dec = np.array(meas['dec'],np.float64)
    dec -= np.mean(dec)
    dec *= 3600*1e3                         # convert to milli arcsec
    # Calculate robust slope
    try:
        #pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)

        # Run RANSAC
        ransac = linear_model.RANSACRegressor()
        ransac.fit(t.reshape(-1,1), dec)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        gdmask = inlier_mask
        pmdec_ransac = ransac.estimator_.coef_[0]
        print('  ransac '+str(np.sum(inlier_mask))+' inliers   '+str(np.sum(outlier_mask))+' outliers')

        # Robust, weighted linear with with INLIERS
        #pmdec_coef, pmdec_coeferr = dln.poly_fit(t[gdmask],dec[gdmask],1,sigma=decerr[gdmask],robust=True,error=True)
        pmdec_coef, pmdec_coeferr = dln.poly_fit(t,dec,1,sigma=decerr,robust=True,error=True)
        pmdec = pmdec_coef[0]
        pmdecerr = pmdec_coeferr[0]
        decdiff = dec-dln.poly(t,pmdec_coef)
        decdiff -= np.median(decdiff)
        decsig = dln.mad(decdiff)
        # Reject outliers
        gdsig = (np.abs(decdiff) < 2.5*decsig) | (np.abs(decdiff) < 2.5*decerr)
        print('  '+str(nmeas-np.sum(gdsig))+' 2.5*sigma clip outliers rejected')
        if np.sum(gdsig) < nmeas:
            pmdec_coef, pmdec_coeferr = dln.poly_fit(t[gdsig],dec[gdsig],1,sigma=decerr[gdsig],robust=True,error=True)
            pmdec = pmdec_coef[0]
            pmdecerr = pmdec_coeferr[0]
            decsig = dln.mad(dec-dln.poly(t,pmdec_coef))            

    except:
        print('problem')
        import pdb; pdb.set_trace()

    out = [pmra,pmraerr,pmra_ransac,rasig,pmdec,pmdecerr,pmdec_ransac,decsig]

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
    objectid = cat['id'].data
    nobj = len(objectid)


    cat['old_pmra'] = cat['pmra'].copy()
    cat['old_pmraerr'] = cat['pmraerr'].copy()
    cat['old_pmdec'] = cat['pmdec'].copy()
    cat['old_pmdecerr'] = cat['pmdecerr'].copy()
    cat['pmra_ransac'] = 999999.
    cat['pmdec_ransac'] = 999999.
    cat['rasig'] = 999999.
    cat['decsig'] = 999999.


    # Fix the pms in healpix object catalogs
    for i in range(nobj):
        objid = objectid[i].astype(str).strip()
        print(str(i+1)+' '+objid)
        out = fix_pms(objid)
        if out is not None:
            print('  OLD:    %10.2f %10.2f' % (cat['pmra'][i],cat['pmdec'][i]))
            print('  NEW:    %10.2f %10.2f' % (out[0],out[4]))
            print('  RANSAC: %10.2f %10.2f' % (out[2],out[6]))
            cat['pmra'][i] = out[0]
            cat['pmraerr'][i] = out[1]
            cat['pmra_ransac'][i] = out[2]
            cat['rasig'][i] = out[3]
            cat['pmdec'][i] = out[4]
            cat['pmdecerr'][i] = out[5]
            cat['pmdec_ransac'][i] = out[6]
            cat['decsig'][i] = out[7]

    # Save the corrected file
    print('Saving corrected file to '+outfile)
    cat.write(outfile,overwrite=True)
