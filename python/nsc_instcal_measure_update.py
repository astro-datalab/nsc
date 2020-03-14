#!/usr/bin/env python

# Update the measurement catalog with the objectID

import os
import numpy as np
import time
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from dlnpyutils import utils as dln

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

    # Load the catalogs
    #chstr = np.lib.recfunctions.append_fields(chstr,'meas_index',np.zeros(len(chstr),int),usemask=False,asrecarray=True)
    #chstr = dln.addcatcols(chstr,np.dtype([('MEAS_INDEX',int)]))
    chstr['MEAS_INDEX'] = 0
    count = 0
    meas = None
    for i in range(nchips):
        #import pdb; pdb.set_trace()
        cat1 = fits.getdata(chstr['MEASFILE'][i].strip(),1)
        ncat1 = len(cat1)
        if meas is None:
            dtype = cat1.dtype
            meas = np.zeros(int(np.sum(chstr['NMEAS'])),dtype=dtype)
        meas[count:count+ncat1] = cat1
        chstr['MEAS_INDEX'][i] = count
        count += ncat1
    meas['MEASID'] = str(meas['MEASID'])

    # Get the OBJECTID from the combined healpix file IDSTR structure
    #  remove any sources that weren't used

    # Figure out which healpix this figure overlaps

    pix = hp.ang2pix(nside,meas['RA'],meas['DEC'],lonlat=True)
    #theta = (90-meas.dec)/radeg
    #phi = meas.ra/radeg
    #ANG2PIX_RING,nside,theta,phi,pix
    upix = np.uniq(pix)
    npix = len(upix)

    import pdb; pdb.set_trace()

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
    for i in range(npix):
        objfile = cmbdir+'combine/'+str(int(upix[i])/1000)+'/'+str(upix[i])+'.fits.gz'
        if os.path.exists(objfile):
            idstr = fits.getdata(objfile,3)
            #idstr['measid'] = str(idstr['measid'])
            nidstr = len(idstr)
            ind1,ind2 = dln.match(idstr['measid'],meas['measid'])
            nmatch = len(ind1)
            if nmatch>0:
                meas['objectid'][ind2] = str(idstr['objectid'][ind1])
            print(str(i+1)+' '+str(upix[i])+' '+str(nmatch))
        else:
            print(objfile+' NOT FOUND')


    # Only keep sources with an objectid
    ind,nind = dln.where(meas['objectid'] == '')
    if nind > 0:
        raise ValueError('some objects are missing IDs')
    #if nind gt 0 then begin
    #  print,'Keeping ',strtrim(nind,2),' of ',strtrim(ncat,2),' sources'
    #  meas = meas[ind]
    #endif else begin
    #  print,'No sources to keep'
    #  return
    #endelse

    # Output
    print('Updating measurement catalogs')
    for i in range(nchips):
        lo = chstr['meas_index'][i]
        hi = lo+chstr['nmeas'][i].nmeas-1
        meas1 = meas[lo:hi+1]
        meta1 = fits.getdata(chstr['measfile'][i],2)        # load the meta extensions
        #MWRFITS,meas1,chstr[i].measfile,/create
        #MWRFITS,meta1,chstr[i].measfile,/silent


    print('dt = ',str(time.time()-t0)+' sec.')
