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
import tempfile
import psycopg2 as pq

def get_meas(pix):
    """ Get the measurements for a particular healpix."""
    # objid, ra, raerr, dec, decerr, mjd

    connection = pq.connect(user="dlquery",host="db01.datalab.noao.edu",
                            password="",port = "5432",database = "tapdb")
    cur = connection.cursor()

    cmd = """SELECT m.objectid,m.ra,m.raerr,m.dec,m.decerr,m.mjd from nsc_dr2.meas as m join
             nsc_dr2.object as obj on m.objectid=obj.objectid where obj.pix={0};""".format(pix)

    cur.execute(cmd)
    data = cur.fetchall()
    # Convert to numpy structured array
    dtype = np.dtype([('objectid',np.str,50),('ra',np.float64),('raerr',float),
                      ('dec',np.float64),('decerr',float),('mjd',np.float64)])
    meas = np.zeros(len(data),dtype=dtype)
    meas[...] = data
    del(data)

    cur.close()
    connection.close()

    return meas


def fix_pms(pix):
    """ Correct the proper motions in the healpix object catalog."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    version = 'v3'
    radeg = np.float64(180.00) / np.pi

    hdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/combine/'+str(int(pix)//1000)+'/'
    objfile = hdir+str(pix)+'.fits.gz'
    outfile = hdir+str(pix)+'_pmcorr.fits'
    
    print('Correcting proper motions for '+str(pix))

    # Check that the object file exists
    if os.path.exists(objfile) is False:
        print(objfile+' NOT FOUND')
        return

    # Check fixed file  
    if os.path.exists(outfile+'.gz') == True:
        print(str(pix)+' already fixed')
        return

    # Load the object file
    meta = fits.getdata(objfile,1)
    obj = fits.getdata(objfile,2)
    nobj = len(obj)
    print(str(nobj)+' objects')

    # Get the measurements
    meas = get_meas(pix)
    nmeas = len(meas)
    print('Retrieved '+str(nmeas)+' measurements')

    idindex = dln.create_index(meas['objectid'])
    # Not all matched
    if len(idindex['value']) != nobj:
        print('Number of unique OBJECTIDs in object and meas catalogs do not match')
        return
    ind1,ind2 = dln.match(obj['objectid'],idindex['value'])
    # Not all matched
    if len(ind1) != nobj:
        print('Some objects are missing measurements')
        return
    # sort by object index
    si = np.argsort(ind1)
    ind1 = ind1[si]
    ind2 = ind2[si]

    # Loop over
    for i in range(nobj):
        if (i % 1000)==0: print(i)
        # Calculate the proper motions
        mind = idindex['index'][idindex['lo'][ind2[i]]:idindex['hi'][ind2[i]]+1]
        cat1 = meas[mind]
        ncat1 = len(cat1)
        if ncat1>1:
            raerr = np.array(cat1['raerr']*1e3,np.float64)    # milli arcsec
            ra = np.array(cat1['ra'],np.float64)
            ra -= np.mean(ra)
            ra *= 3600*1e3 * np.cos(obj['dec'][i]/radeg)     # convert to true angle, milli arcsec
            t = cat1['mjd'].copy()
            t -= np.mean(t)
            t /= 365.2425                          # convert to year
            # Calculate robust slope
            try:
                pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
            except:
                print('problem')
                import pdb; pdb.set_trace()
            obj['pmra'][i] = pmra                 # mas/yr
            obj['pmraerr'][i] = pmraerr           # mas/yr

            decerr = np.array(cat1['decerr']*1e3,np.float64)   # milli arcsec
            dec = np.array(cat1['dec'],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            # Calculate robust slope
            try:
                pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
            except:
                print('problem')
                import pdb; pdb.set_trace()
            obj['pmdec'][i] = pmdec               # mas/yr
            obj['pmdecerr'][i] = pmdecerr         # mas/yr

    # Save the new version of obj
    # Write the output file
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
    if type(pix) is not list: pix=[pix]
    npix = len(pix)
    print('Correcting PMs for '+str(npix)+' HEALPix')
    for i in range(len(pix)):
        fix_pms(pix[i])


