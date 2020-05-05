#!/usr/bin/env python

# Update the measurement catalog with the objectID

import os
import sys
import numpy as np
import time
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from dlnpyutils import utils as dln
import shutil
import sqlite3
from glob import glob
import socket
from argparse import ArgumentParser

def querydb(dbfile,table='meas',cols='rowid,*',where=None):
    """ Query database table """
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()
    cmd = 'SELECT '+cols+' FROM '+table
    if where is not None: cmd += ' WHERE '+where
    cur.execute(cmd)
    data = cur.fetchall()
    db.close()

    # Return results
    return data


def readidstrdb(dbfile,where=None):
    """ Get data from IDSTR database"""
    data = querydb(dbfile,table='idstr',cols='*',where=where)
    # Put in catalog
    dtype_idstr = np.dtype([('measid',np.str,200),('exposure',np.str,200),('objectid',np.str,200),('objectindex',int)])
    cat = np.zeros(len(data),dtype=dtype_idstr)
    cat[...] = data
    del data
    return cat

def measurement_update(expdir):

    t0 = time.time()

    # Get version number from exposure directory
    lo = expdir.find('nsc/instcal/')
    dum = expdir[lo+12:]
    version = dum[0:dum.find('/')]
    cmbdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/'
    edir = '/net/dl1/users/dnidever/nsc/instcal/'+version+'/'
    nside = 128

    # Check if output file already exists
    base = os.path.basename(expdir)

    print('Adding objectID for measurement catalogs for exposure = '+base)

    #  Load the exposure and metadata files
    metafile = expdir+'/'+base+'_meta.fits'
    meta = Table.read(metafile,1)
    nmeta = len(meta)
    chstr = Table.read(metafile,2)
    #chstr['measfile'] = str(chstr['measfile'])
    nchips = len(chstr)

    measdtype = np.dtype([('MEASID', 'S50'), ('OBJECTID', 'S50'), ('EXPOSURE', 'S50'), ('CCDNUM', '>i2'), ('FILTER', 'S2'), ('MJD', '>f8'), ('X', '>f4'),
                          ('Y', '>f4'), ('RA', '>f8'), ('RAERR', '>f4'), ('DEC', '>f8'), ('DECERR', '>f4'), ('MAG_AUTO', '>f4'), ('MAGERR_AUTO', '>f4'),
                          ('MAG_APER1', '>f4'), ('MAGERR_APER1', '>f4'), ('MAG_APER2', '>f4'), ('MAGERR_APER2', '>f4'), ('MAG_APER4', '>f4'),
                          ('MAGERR_APER4', '>f4'), ('MAG_APER8', '>f4'), ('MAGERR_APER8', '>f4'), ('KRON_RADIUS', '>f4'), ('ASEMI', '>f4'), ('ASEMIERR', '>f4'),
                          ('BSEMI', '>f4'), ('BSEMIERR', '>f4'), ('THETA', '>f4'), ('THETAERR', '>f4'), ('FWHM', '>f4'), ('FLAGS', '>i2'), ('CLASS_STAR', '>f4')])

    # Load and concatenate the meas catalogs
    chstr['MEAS_INDEX'] = 0   # keep track of where each chip catalog starts
    count = 0
    meas = Table(data=np.zeros(int(np.sum(chstr['NMEAS'])),dtype=measdtype))
    print('Loading and concatenating the chip measurement catalogs')
    for i in range(nchips):
        meas1 = Table.read(chstr['MEASFILE'][i].strip(),1)   # load chip meas catalog
        nmeas1 = len(meas1)
        meas[count:count+nmeas1] = meas1
        chstr['MEAS_INDEX'][i] = count
        count += nmeas1
    measid = np.char.array(meas['MEASID']).strip().decode()
    print(str(len(meas))+' measurements')

    # Get the OBJECTID from the combined healpix file IDSTR structure
    #  remove any sources that weren't used

    # Figure out which healpix this figure overlaps
    pix = hp.ang2pix(nside,meas['RA'],meas['DEC'],lonlat=True)
    upix = np.unique(pix)
    npix = len(upix)
    print(str(npix)+' HEALPix to query')

    # Loop over the HEALPix pixels
    ntotmatch = 0
    for i in range(npix):
        fitsfile = cmbdir+'combine/'+str(int(upix[i])//1000)+'/'+str(upix[i])+'.fits.gz'
        dbfile = cmbdir+'combine/'+str(int(upix[i])//1000)+'/'+str(upix[i])+'_idstr.db'
        if os.path.exists(dbfile):
            # Read meas id information from idstr database for this expoure
            idstr = readidstrdb(dbfile,where="exposure=='"+base+"'")
            nidstr = len(idstr)
            idstr_measid = np.char.array(idstr['measid']).strip()
            idstr_objectid = np.char.array(idstr['objectid']).strip()
            ind1,ind2 = dln.match(idstr_measid,measid)
            nmatch = len(ind1)
            if nmatch>0:
                meas['OBJECTID'][ind2] = idstr_objectid[ind1]
                ntotmatch += nmatch
            print(str(i+1)+' '+str(upix[i])+' '+str(nmatch))
        else:
            print(str(i+1)+' '+dbfile+' NOT FOUND.  Checking for high-resolution database files.')
            # Check if there are high-resolution healpix idstr databases
            hidbfiles = glob(cmbdir+'combine/'+str(int(upix[i])//1000)+'/'+str(upix[i])+'_n*_*_idstr.db')
            nhidbfiles = len(hidbfiles)
            if os.path.exists(fitsfile) & (nhidbfiles>0):
                print('Found high-resolution HEALPix IDSTR files')
                for j in range(nhidbfiles):
                    dbfile1 = hidbfiles[j]
                    dbbase1 = os.path.basename(dbfile1)
                    idstr = readidstrdb(dbfile1,where="exposure=='"+base+"'")
                    nidstr = len(idstr)
                    idstr_measid = np.char.array(idstr['measid']).strip()
                    idstr_objectid = np.char.array(idstr['objectid']).strip()
                    ind1,ind2 = dln.match(idstr_measid,measid)
                    nmatch = len(ind1)
                    if nmatch>0:
                        meas['OBJECTID'][ind2] = idstr_objectid[ind1]
                        ntotmatch += nmatch
                    print('  '+str(j+1)+' '+dbbase1+' '+str(upix[i])+' '+str(nmatch))

    # Only keep sources with an objectid
    ind,nind = dln.where(meas['OBJECTID'] == '')
    if nind > 0:
        raise ValueError(str(nind)+' measurements are missing OBJECTIDs')


    # Output the updated catalogs
    print('Updating measurement catalogs')
    for i in range(nchips):
        measfile1 = chstr['MEASFILE'][i].strip()
        lo = chstr['MEAS_INDEX'][i]
        hi = lo+chstr['NMEAS'][i]
        meas1 = meas[lo:hi]
        meta1 = Table.read(measfile1,2)        # load the meta extensions
        # Copy as a backup
        if os.path.exists(measfile1+'.bak'): os.remove(measfile1+'.bak')
        shutil.copyfile(measfile1,measfile1+'.bak')
        # Write new catalog
        meas1.write(measfile1,overwrite=True)  # first, measurement table
        # append other fits binary table
        hdulist = fits.open(measfile1)
        hdu = fits.table_to_hdu(meta1)         # second, catalog
        hdulist.append(hdu)
        hdulist.writeto(measfile1,overwrite=True)
        hdulist.close()
        # Create a file saying that the file was successfully updated.
        dln.writelines(measfile1+'.updated','')
        # Delete backups
        if os.path.exists(measfile1+'.bak'): os.remove(measfile1+'.bak')

    # Create a file saying that the files were updated okay.
    dln.writelines(expdir+'/'+base+'_meas.updated','')

    print('dt = ',str(time.time()-t0)+' sec.')


if __name__ == "__main__":
    parser = ArgumentParser(description='Update NSC exposure measurement catalogs with OBJECTID.')
    parser.add_argument('expdir', type=str, nargs=1, help='Exposure directory')
    #parser.add_argument('-r','--redo', action='store_true', help='Redo this exposure catalog')
    #parser.add_argument('-v','--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    expdir = args.expdir[0]

    # Check if the exposure has already been updated
    base = os.path.basename(expdir)
    if os.path.exists(expdir+'/'+base+'_meas.updated'):
        print(expdir+' has already been updated')
        sys.exit()

    measurement_update(expdir)
