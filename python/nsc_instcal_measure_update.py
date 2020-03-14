#!/usr/bin/env python

# Update the measurement catalog with the objectID

import os
import numpy as np
import time
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from dlnpyutils import utils as dln
#from nsc_instcal_combine_cluster import readidstrdb
import shutil
import sqlite3

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
    #NSC_ROOTDIRS,dldir,mssdir,localdir
    cmbdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/'
    edir = '/net/dl1/users/dnidever/nsc/instcal/'+version+'/'
    nside = 128
    radeg = 180 / np.pi

    # Check if output file already exists
    base = os.path.basename(expdir)
    #outfile = expdir+'/'+base+'_1_meas.fits'
    #if file_test(outfile) eq 0 then begin
    #  print,outfile,' NOT FOUND'
    #  return
    #endif

    print('Adding objectID for measurement catalogs for exposure = '+base)

    #  Load the exposure and metadata files
    metafile = expdir+'/'+base+'_meta.fits'
    meta = fits.getdata(metafile,1)
    nmeta = len(meta)
    #chstr = fits.getdata(metafile,2)
    chstr = Table.read(metafile,2)
    #chstr['measfile'] = str(chstr['measfile'])
    nchips = len(chstr)

    measdtype = np.dtype([('MEASID', 'S50'), ('OBJECTID', 'S50'), ('EXPOSURE', 'S50'), ('CCDNUM', '>i2'), ('FILTER', 'S2'), ('MJD', '>f8'), ('X', '>f4'),
                          ('Y', '>f4'), ('RA', '>f8'), ('RAERR', '>f4'), ('DEC', '>f8'), ('DECERR', '>f4'), ('MAG_AUTO', '>f4'), ('MAGERR_AUTO', '>f4'),
                          ('MAG_APER1', '>f4'), ('MAGERR_APER1', '>f4'), ('MAG_APER2', '>f4'), ('MAGERR_APER2', '>f4'), ('MAG_APER4', '>f4'),
                          ('MAGERR_APER4', '>f4'), ('MAG_APER8', '>f4'), ('MAGERR_APER8', '>f4'), ('KRON_RADIUS', '>f4'), ('ASEMI', '>f4'), ('ASEMIERR', '>f4'),
                          ('BSEMI', '>f4'), ('BSEMIERR', '>f4'), ('THETA', '>f4'), ('THETAERR', '>f4'), ('FWHM', '>f4'), ('FLAGS', '>i2'), ('CLASS_STAR', '>f4')])

    # Load the catalogs
    #chstr = np.lib.recfunctions.append_fields(chstr,'meas_index',np.zeros(len(chstr),int),usemask=False,asrecarray=True)
    #chstr = dln.addcatcols(chstr,np.dtype([('MEAS_INDEX',int)]))
    chstr['MEAS_INDEX'] = 0
    count = 0
    meas = None
    for i in range(nchips):
        #import pdb; pdb.set_trace()
        cat1 = Table.read(chstr['MEASFILE'][i].strip(),1)
        ncat1 = len(cat1)
        if meas is None:
            #dtype = cat1.dtype
            meas = Table(data=np.zeros(int(np.sum(chstr['NMEAS'])),dtype=measdtype))
        meas[count:count+ncat1] = cat1
        chstr['MEAS_INDEX'][i] = count
        count += ncat1
    #meas['MEASID'] = meas['MEASID'
    #import pdb; pdb.set_trace()  

    measid = np.char.array(meas['MEASID']).strip().decode()
    print(str(len(meas))+' measurements')

    # Get the OBJECTID from the combined healpix file IDSTR structure
    #  remove any sources that weren't used

    # Figure out which healpix this figure overlaps

    pix = hp.ang2pix(nside,meas['RA'],meas['DEC'],lonlat=True)
    #theta = (90-meas.dec)/radeg
    #phi = meas.ra/radeg
    #ANG2PIX_RING,nside,theta,phi,pix
    upix = np.unique(pix)
    npix = len(upix)
    print(str(npix)+' HEALPix to query')

    # Load the healpix list
    #listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_instcal_combine_healpix_list.fits.gz'
    #if file_test(listfile) eq 0 then begin
    #  print,listfile,' NOT FOUND'
    #  return
    #endif
    #healstr = MRDFITS(listfile,1,/silent)  # takes ~20s to load
    #healstr.file = strtrim(healstr.file,2)
    #healstr.base = strtrim(healstr.base,2)
    #ind = where(healstr.base eq base,nind)  # takes ~2s
    #upix = healstr[ind].pix

    # Loop over the pixels
    ntotmatch = 0
    for i in range(npix):
        dbfile = cmbdir+'combine/'+str(int(upix[i])//1000)+'/'+str(upix[i])+'_idstr.db'
        if os.path.exists(dbfile):
            #import pdb; pdb.set_trace()
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
            print(str(i+1)+' '+dbfile+' NOT FOUND')


    # Only keep sources with an objectid
    ind,nind = dln.where(meas['OBJECTID'] == '')
    if nind > 0:
        raise ValueError('some objects are missing IDs')
    #if nind gt 0 then begin
    #  print,'Keeping ',strtrim(nind,2),' of ',strtrim(ncat,2),' sources'
    #  meas = meas[ind]
    #endif else begin
    #  print,'No sources to keep'
    #  return
    #endelse

    import pdb; pdb.set_trace()

    # Output
    print('Updating measurement catalogs')
    for i in range(nchips):
        measfile1 = chstr['MEASFILE'][i].strip()
        lo = chstr['MEAS_INDEX'][i]
        hi = lo+chstr['NMEAS'][i]
        meas1 = meas[lo:hi]
        meta1 = Table.read(measfile1,2)        # load the meta extensions
        # Copy as a backup
        #if os.path.exists(measfile1+'.bak'): os.remove(measfile1+'.bak')
        #shutil.copyfile(measfile1,measfile1+'.bak')
        # meas1.write(measfile1,overwrite=True)  # first, measurement table
        #  append other fits binary table
        #hdulist = fits.open(measfile1)
        #hdu = fits.table_to_hdu(meta1)        # second, catalog
        #hdulist.append(hdu)
        #hdulist.writeto(measfile1,overwrite=True)
        #hdulist.close()

        #MWRFITS,meas1,chstr[i].measfile,/create
        #MWRFITS,meta1,chstr[i].measfile,/silent


    print('dt = ',str(time.time()-t0)+' sec.')
