#!/usr/bin/env python

import os
import time
import numpy as np
from glob import glob
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from scipy.optimize import curve_fit
from scipy import stats
import subprocess
import traceback
from . import utils,query,modelmag

def concatenate(expdir):
    """
    Combine multiple chip-level measurement files into a single multi-extension FITS file.
    """
    if os.path.exists(expdir)==False:
        print(expdir,'not found')
        return

    print('Concatenate FITS files for ',expdir)
    
    curdir = os.getcwd()
    
    # Is this a tar file
    if os.path.isdir(expdir)==False and expdir.endswith('.tar.gz') or expdir.endswith('.tgz'):
        print('This is a tar file.  Uncompressing.')
        tarfile = os.path.basename(expdir)
        base = tarfile.replace('.tgz','').replace('.tar.gz','')
        nightdir = os.path.abspath(os.path.dirname(expdir))
        expdir = nightdir+'/'+base
        os.chdir(nightdir)        
        res = subprocess.run(['tar','xvf',tarfile],capture_output=True)
        if res.returncode==0:
            print('success: removing',tarfile)
            os.remove(tarfile)
        else:
            print('problem untarring',tarfile)
            import pdb; pdb.set_trace()

    os.chdir(expdir)
    base = os.path.basename(expdir)
    outfile = base+'_meas.fits'
    if os.path.exists(outfile):
        print(outfile,'already exists')
        os.chdir(curdir)
        return
    fitsfiles1 = glob('*_?.fits')
    fitsfiles1.sort()
    fitsfiles2 = glob('*_??.fits')
    fitsfiles2.sort()
    fitsfiles = fitsfiles1+fitsfiles2
    # this assumes that it is a c4d exposure
    if len(fitsfiles)<59:
        print(len(fitsfiles),'fits files found. not enough.  skipping')
        os.chdir(curdir)
        return
    chdu = fits.HDUList()
    hhdu = fits.HDUList()
    hhdu.append(fits.PrimaryHDU())
    for i in range(len(fitsfiles)):
        hdu1 = fits.open(fitsfiles[i])
        tab1 = Table.read(fitsfiles[i])
        base1 = os.path.basename(fitsfiles[i])
        print('{:3d} {:s} {:8d}'.format(i+1,base1,len(tab1)))
        ccdnum = base1[:-5].split('_')[-1]
        newhdu = fits.table_to_hdu(Table(hdu1[1].data))
        newhdu.header['extname'] = ccdnum
        newhdu.header['ccdnum'] = ccdnum
        newhead = hdu1[0].header.copy()
        newhead['extname'] = ccdnum
        chdu.append(newhdu)
        hhdu.append(fits.ImageHDU(header=newhead))
        hdu1.close()
    print('Writing',outfile)
    chdu.writeto(outfile,overwrite=True)
    chdu.close()
    houtfile = base+'_header.fits'
    print('Writing',houtfile)
    hhdu.writeto(houtfile,overwrite=True)
    hhdu.close()
    # Confirm that it is there
    if os.path.exists(outfile) and os.path.exists(houtfile):
        # Delete fits files
        print('Deleting',len(fitsfiles),'individual fits files')
        for f in fitsfiles:
            if os.path.exists(f): os.remove(f)
        # Tar up the rest of the files
        # leave out meas.fits, log
        logfile = base+'.log'
        tarfile = base+'.tgz'
        allfiles = glob('*')
        allfiles = [f for f in allfiles if os.path.isfile(f)]
        allfiles.sort()
        if outfile in allfiles:
            allfiles.remove(outfile)
        if houtfile in allfiles:
            allfiles.remove(houtfile)
        if logfile in allfiles:
            allfiles.remove(logfile)
        if tarfile in allfiles:
            allfiles.remove(tarfile)
        cmd = ['tar','cvzf',tarfile]+allfiles
        print('tarring',len(allfiles),'files in',tarfile)
        res = subprocess.run(cmd,capture_output=True)
        if res.returncode==0:
            print('tar success')
            print('Deleting tarred files')
            for f in allfiles:
                if os.path.exists(f): os.remove(f)
    os.chdir(curdir)

def uncalcoord(ra,dec,chinfo):
    # Uncalibrate coordinates, only good measurements
    #  we need the original lon/lat values to get the right
    #  correction terms, iterate until convergence
    lon,lat = coords.rotsphcen(ra,dec,chinfo['cenra'],chinfo['cendec'],gnomic=True)
    flag = 0
    fra = ra.copy()
    fdec = dec.copy()
    last_fra = fra.copy()*0
    last_fdec = fra.copy()*0
    count = 0
    while (flag==0):
        flon,flat = coords.rotsphcen(fra,fdec,chinfo['cenra'],chinfo['cendec'],gnomic=True)
        lon2 = lon - utils.func_poly2d(flon,flat,*chinfo['racoef'])
        lat2 = lat - utils.func_poly2d(flon,flat,*chinfo['deccoef'])
        fra,fdec = coords.rotsphcen(lon2,lat2,chinfo['cenra'],chinfo['cendec'],reverse=True,gnomic=True)
        dra = np.abs(fra-last_fra)*3600
        ddec = np.abs(fdec-last_fdec)*3600        
        maxdiff = np.max(np.concatenate((dra,ddec)))
        last_fra = fra.copy()
        last_fdec = fdec.copy()
        if maxdiff < 1e-6 or count>10:
            flag = 1
        count += 1            
    return fra,fdec
                
def recreatemeas(calfile,metafile,outfile):
    """
    Recreate original measurement table using the calibrated table and
    meta-data file.

    Parameters
    ----------
    calfile : str
       Filaname of the calibrated measurement file.
    metafile : str
       Filename of the metadata file.
    outfile : str
       Filename for the output recreated measurement catalog.

    Returns
    -------
    The recreated original measurement file is written to disk to "outfile".

    Examples
    --------

    recreatemeas(calfile,metafile,outfile)

    """

    # Load meta-data
    mhdu = fits.open(metafile)
    expinfo = Table(mhdu[1].data)
    chinfo = None
    for i in range(len(mhdu)-2):
        if chinfo is None:
            chinfo = mhdu[i+2].data
        else:
            chinfo = np.concatenate((chinfo,mhdu[i+2].data))
    chinfo = Table(chinfo)
    mhdu.close()

    # Original measurement table format
    dt = [('NUMBER', '>i4'), ('X_IMAGE', '>f4'), ('Y_IMAGE', '>f4'), ('MAG_APER', '>f4', (5,)),
          ('MAGERR_APER', '>f4', (5,)), ('MAG_ISO', '>f4'), ('MAGERR_ISO', '>f4'),
          ('MAG_AUTO', '>f4'), ('MAGERR_AUTO', '>f4'), ('KRON_RADIUS', '>f4'),
          ('BACKGROUND', '>f4'), ('THRESHOLD', '>f4'), ('ISOAREA_IMAGE', '>i4'),
          ('ISOAREA_WORLD', '>f4'), ('ALPHA_J2000', '>f8'), ('DELTA_J2000', '>f8'),
          ('X2_WORLD', '>f8'), ('Y2_WORLD', '>f8'), ('XY_WORLD', '>f8'),
          ('A_WORLD', '>f4'), ('B_WORLD', '>f4'), ('THETA_WORLD', '>f4'),
          ('ELLIPTICITY', '>f4'), ('ERRX2_WORLD', '>f8'), ('ERRY2_WORLD', '>f8'),
          ('ERRXY_WORLD', '>f8'), ('ERRA_WORLD', '>f4'), ('ERRB_WORLD', '>f4'),
          ('ERRTHETA_WORLD', '>f4'), ('FWHM_WORLD', '>f4'), ('FLAGS', '>i2'),
          ('IMAFLAGS_ISO', '>i4'), ('NIMAFLAGS_ISO', '>i4'), ('CLASS_STAR', '>f4'),
          ('NDET_ITER', '>i8'), ('REPEAT', '>f8'), ('XPSF', '>f8'), ('YPSF', '>f8'),
          ('MAGPSF', '>f8'), ('ERRPSF', '>f8'), ('SKY', '>f8'), ('ITER', '>f8'),
          ('CHI', '>f8'), ('SHARP', '>f8'), ('RAPSF', '>f8'), ('DECPSF', '>f8')]
    
    # Loop over the chips
    hdu = fits.open(calfile)
    ohdu = fits.HDUList()
    for i in range(len(hdu)-1):
        tab1 = Table(hdu[i+1].data)
        zpoff = 2.5*np.log10(expinfo['exptime'][0]) + chinfo['zpterm'][i]
        meas1 = Table(np.zeros(len(tab1),dtype=np.dtype(dt)))
        # measid = instrument+'.'+str(expnum)+'.'+meas['ccdnum']+'.'+meas['number']
        number = [n.split('.')[-1] for n in tab1['measid']]
        meas1['NUMBER'][:] = number
        meas1['X_IMAGE'][:] = tab1['x']
        meas1['Y_IMAGE'][:] = tab1['y']
        # Uncalibrate, only good mags
        for c in ['mag_auto','mag_aper1','mag_aper2','mag_aper4','mag_aper6','mag_aper8','magpsf']:
            if c[:8]=='mag_aper':
                indx = {'1':0,'2':1,'4':2,'6':3,'8':4}[c[-1]]
                mag = tab1[c]
                gdmag = ((mag < 50) & np.isfinite(mag))
                meas1['MAG_APER'][:,indx] = np.nan   # bad by default
                meas1['MAG_APER'][gdmag,indx] = mag[gdmag]-zpoff
            else:
                mag = tab1[c]
                gdmag = ((mag < 50) & np.isfinite(mag))
                meas1[c.upper()][:] = np.nan   # bad by default
                meas1[c.upper()][gdmag] = mag[gdmag]-zpoff
        meas1['MAGERR_AUTO'][:] = tab1['magerr_auto']
        meas1['MAGERR_APER'][:,0] = tab1['magerr_aper1']        
        meas1['MAGERR_APER'][:,1] = tab1['magerr_aper2']
        meas1['MAGERR_APER'][:,2] = tab1['magerr_aper4']
        meas1['MAGERR_APER'][:,3] = tab1['magerr_aper6']
        meas1['MAGERR_APER'][:,4] = tab1['magerr_aper8']
        # uncalibrate coordinates
        #  SE coordinates
        fra,fdec = uncalcoord(tab1['ra'],tab1['dec'],chinfo[i])
        meas1['ALPHA_J2000'][:] = fra
        meas1['DELTA_J2000'][:] = fdec
        # PSF coordinates
        meas1['RAPSF'][:] = tab1['rapsf']
        meas1['DECPSF'][:] = tab1['decpsf']
        gdcoo, = np.where(np.isfinite(tab1['rapsf']) & np.isfinite(tab1['decpsf']))
        fra,fdec = uncalcoord(tab1['rapsf'][gdcoo],tab1['decpsf'][gdcoo],chinfo[i])
        meas1['RAPSF'][gdcoo] = fra
        meas1['DECPSF'][gdcoo] = fdec
        # Rest of the columns
        meas1['MAG_ISO'][:] = tab1['mag_iso']
        meas1['MAGERR_ISO'][:] = tab1['magerr_iso']
        meas1['KRON_RADIUS'][:] = tab1['kron_radius']
        meas1['BACKGROUND'][:] = tab1['background']
        meas1['THRESHOLD'][:] = tab1['threshold']
        meas1['ISOAREA_IMAGE'][:] = tab1['isoarea_image']
        meas1['ISOAREA_WORLD'][:] = tab1['isoarea_world']
        meas1['X2_WORLD'][:] = tab1['x2_world']
        meas1['Y2_WORLD'][:] = tab1['y2_world']
        meas1['XY_WORLD'][:] = tab1['xy_world']
        meas1['A_WORLD'][:] = tab1['asemi'] / 3600.0          # convert arcsec -> deg
        meas1['B_WORLD'][:] = tab1['bsemi'] / 3600.0          # convert arcsec -> deg
        meas1['THETA_WORLD'][:] = 90-tab1['theta']
        meas1['ELLIPTICITY'][:] = tab1['ellipticity']        
        meas1['ERRX2_WORLD'][:] = tab1['errx2_world']
        meas1['ERRY2_WORLD'][:] = tab1['erry2_world']
        meas1['ERRXY_WORLD'][:] = tab1['errxy_world']
        meas1['ERRA_WORLD'][:] = tab1['asemierr'] / 3600.0    # convert arcsec -> deg
        meas1['ERRB_WORLD'][:] = tab1['bsemierr'] / 3600.0    # convert arcsec -> deg
        meas1['ERRTHETA_WORLD'][:] = tab1['thetaerr']
        meas1['FWHM_WORLD'][:] = tab1['fwhm'] / 3600.0        # convert arcsec -> deg
        meas1['FLAGS'][:] = tab1['flags']
        meas1['IMAFLAGS_ISO'][:] = tab1['imaflags_iso']
        meas1['NIMAFLAGS_ISO'][:] = tab1['nimaflags_iso']
        meas1['CLASS_STAR'][:] = tab1['class_star']
        meas1['NDET_ITER'][:] = tab1['ndet_iter']
        meas1['REPEAT'][:] = tab1['repeat']
        meas1['XPSF'][:] = tab1['xpsf']
        meas1['YPSF'][:] = tab1['ypsf']
        meas1['ERRPSF'][:] = tab1['errpsf']
        meas1['SKY'][:] = tab1['skypsf']
        meas1['ITER'][:] = tab1['iter']
        meas1['CHI'][:] = tab1['chi']
        meas1['SHARP'][:] = tab1['sharp']

        # Append measurements
        ohdu.append(fits.table_to_hdu(meas1))
        
    # Write to new file
    ohdu.writeto(outfile,overwrite=True)
    ohdu.close()

def getzpterm(meas1,ref1,mmags,expinfo,chinfo,kind='modelmag'):
    """
    Obtain the photometric zero-point using matched observed and
    reference catalogs and a specific method.

    Parameters
    ----------
    meas1 : astropy table
       Observed measurement table.
    ref1  : astropy table
       Reference table matched to "meas1".
    expinfo : astropy table
       Meta-data table for the entire exposure.
    chinfo : astropy table
       The table of meta-data for each chip.
    kind : str, optional
       The method to use to determine the zeropoint: "modelmag",
         "gaiaxpsynth", or "selfcal".  Default is "modelmag".

    Returns
    -------
    expinfo : astropy table
       Meta-data table for the entire exposure updated with
         zero-point information.
    chinfo : astropy table
       The table of meta-data for each chip updated with
         zero-point information.

    Examples
    --------

    expinfo,chinfo = getzpterm(meas1,ref1,mmags,expinfo,chinfo,kind='modelmag')

    """

    medfwhm = expinfo['fwhm'][0]
    exptime = expinfo['exptime'][0]
    filt = expinfo['filter'][0]
    
    # Model Magnitudes
    #-----------------
    if kind=='modelmag':
        
        # Get the good sources 
        gdmeas, = np.where((meas1['imaflags_iso'] == 0) & (~((meas1['flags'] & 8) > 0)) & (~((meas1['flags'] & 16) > 0)) &
                          (meas1['magpsf'] < 50) &  (meas1['magerr_auto'] < 0.05) & (meas1['class_star'] > 0.8) &
                          (meas1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
        #  if the seeing is bad then class_star sometimes doesn't work well 
        if medfwhm > 1.8 and len(gdmeas) < 100:
            gdmeas, = np.where((meas1['imaflags_iso'] == 0) & (~((meas1['flags'] & 8) > 0)) & (~((meas1['flags'] & 16) > 0)) &
                              (meas1['magpsf'] < 50) &  (meas1['magerr_auto'] < 0.05) & 
                              (meas1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
        if len(gdmeas) > 0: 
            ref2 = ref1[gdmeas] 
            mmags2 = mmags[gdmeas,:] 
            meas2 = meas1[gdmeas]
            # Matched structure
            mag2 = meas2['magpsf'] + 2.5*np.log10(exptime)   # correct for the exposure time
            mstr = {'col':mmags2[:,2],'mag':mag2,'model':mmags2[:,0],
                    'err':mmags2[:,1],'ccdnum':meas2['ccdnum']}
            # Measure the zero-point
            expinfo,chinfo = fitzpterm(mstr,expinfo,chinfo)
            expinfo['zptype'] = 1
        else:
            logger.info('No good reference sources')

    # Gaia XP synthetic photometry
    #-----------------------------
    elif kind=='gaiaxpsynth':
        
        refmagcol = 'gsynth_'+filt.lower()+'mag'
        referrcol = 'e_gsynth_'+filt.lower()+'mag'
        # Get the good sources 
        gdmeas, = np.where((meas1['imaflags_iso'] == 0) & (~((meas1['flags'] & 8) > 0)) &
                          (~((meas1['flags'] & 16) > 0)) & (meas1['magpsf'] < 50) &
                          (meas1['errpsf'] < 0.05) & (meas1['class_star'] > 0.8) &
                          (meas1['fwhm_world']*3600 < 2*medfwhm) &
                          (np.abs(meas1['sharp']) < 1) & (meas1['chi']<5) &
                          (ref1['bp'] > 0) & (ref1['bp'] < 50) &
                          (ref1['rp'] > 0) & (ref1['rp'] < 50) &                          
                          (ref1[refmagcol] > 0) & (ref1[refmagcol] < 50))
        #  if the seeing is bad then class_star sometimes doesn't work well 
        if medfwhm > 1.8 and len(gdmeas) < 100:
            gdmeas, = np.where((meas1['imaflags_iso'] == 0) & (~((meas1['flags'] & 8) > 0)) &
                              (~((meas1['flags'] & 16) > 0)) &
                              (meas1['magpsf'] < 50) & (meas1['errpsf'] < 0.05) & 
                              (meas1['fwhm_world']*3600 < 2*medfwhm) &
                          (ref1['bp'] > 0) & (ref1['bp'] < 50) &
                          (ref1['rp'] > 0) & (ref1['rp'] < 50) &                          
                          (ref1[refmagcol] > 0) & (ref1[refmagcol] < 50))                              
        if len(gdmeas) > 0:
            ref2 = ref1[gdmeas] 
            mmags2 = mmags[gdmeas,:] 
            meas2 = meas1[gdmeas]
            # Matched structure
            mag2 = meas2['magpsf'] + 2.5*np.log10(exptime)   # correct for the exposure time
            mstr = {'col':ref2['bp']-ref2['rp'],'mag':mag2,'model':ref2[refmagcol],
                    'err':ref2[referrcol],'ccdnum':meas2['ccdnum']}
            # Measure the zero-point
            expinfo,chinfo = fitzpterm(mstr,expinfo,chinfo)
            expinfo['zptype'] = 2
        else:
            logger.info('No good reference sources')

    # Self-calibration
    #-----------------
    elif kind=='selfcal':
        expinfo,chinfo = selfcalzpterm(expdir,meas,expinfo,chinfo)
        expinfo['zptype'] = 3 

    else:
        raise ValueError(kind+' not supported')

    return expinfo,chinfo


def fitzpterm(mstr,expinfo,chinfo):
    # Fit the global and ccd zero-points 
     
    n = len(mstr['col']) 
    diff = mstr['model'] - mstr['mag']
    err = mstr['err']
    # Make a sigma cut 
    med = np.median(diff) 
    sig = dln.mad(diff) 
    gd, = np.where(np.abs(diff-med) < 3*sig) 
    x = np.zeros(n,float)
    zpterm,zptermerr1 = dln.wtmean(diff[gd],err[gd],error=True)
    zptermerr = dln.bootstrap(diff[gd],dln.wtmean,args=err[gd],indexargs=True)
    # Save in exposure table
    expinfo['zpterm'] = zpterm 
    expinfo['zptermerr'] = zptermerr 
    expinfo['zptermsig'] = sig 
    expinfo['nrefmatch'] = n 
    expinfo['ngoodrefmatch'] = len(gd)
     
    # Measure chip-level zpterm and variations 
    nchips = len(chinfo) 
    for i in range(nchips): 
        ind1,ind2 = dln.match(chinfo['ccdnum'][i],mstr['ccdnum'])
        nmatch = len(ind1)
        if nmatch >= 5: 
            gd1, = np.where(np.abs(diff[ind2]-med) < 3*sig) 
            if len(gd1) >= 5: 
                zpterm1, _ = dln.wtmean(diff[ind2[gd1]],err[ind2[gd1]],error=True)
                zptermerr1 = dln.bootstrap(diff[ind2[gd1]],dln.wtmean,args=err[ind2[gd1]],indexargs=True)
                chinfo['zpterm'][i] = zpterm1 
                chinfo['zptermerr'][i] = zptermerr1 
                chinfo['nrefmatch'][i] = nmatch 
     
    # Measure spatial variations 
    gdchip, = np.where(chinfo['zpterm'] < 1000) 
    if len(gdchip) > 1: 
        expinfo['zpspatialvar_rms'] = np.std(chinfo['zpterm'][gdchip]) 
        expinfo['zpspatialvar_range'] = np.max(chinfo['zpterm'][gdchip])-np.min(chinfo['zpterm'][gdchip]) 
        expinfo['zpspatialvar_nccd'] = len(gdchip)

    return expinfo,chinfo


def selfcalzpterm(expdir,cat,expinfo,chinfo,logger=None,silent=False):
     
    # Measure the zero-point for this exposure using other 
    # successfully-calibrated exposures as secondary references 
     
    dldir,mssdir,localdir = utils.rootdirs()

    if logger is None:
        logger = dln.basiclogger()
     
    logger.info('Self-calibrating')
    # Set zero-point column to bad by default 
    expinfo['zpterm'] = 999999. 
    expinfo['zptermerr'] = 999999. 
    expinfo['zptermsig'] = 999999. 
    expinfo['zpspatialvar_rms'] = 999999. 
    expinfo['zpspatialvar_range'] = 999999. 
    expinfo['zpspatialvar_nccd'] = 0 
    expinfo['nrefmatch'] = 0 
    expinfo['ngoodrefmatch'] = 0 
    chinfo['zpterm'] = 999999. 
    chinfo['zptermerr'] = 999999. 
    chinfo['nrefmatch'] = 0 
     
    # What instrument is this? 
    instrument = 'c4d'  # by default 
    if expdir.find('/k4m/') > -1:
        instrument = 'k4m' 
    if expdir.find('/ksb/') > -1:        
        instrument = 'ksb' 
    instfilt = instrument+'-'+str(expinfo['filter'][0]) 
    # version
    # get version number 
    lo = expdir.find('nsc/instcal/') 
    dum = expdir[lo+12:] 
    version = dum[0:dum.find('/')]
    if version == 'c4d' or version == 'k4m' or version == 'ksb': # version1 
        version = '' 
     
    # Central coordinates 
    cenra = expinfo['ra']
    cendec = expinfo['dec']
     
    # Read in the summary structure 
    #sumfile = dldir+'dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_calibrate.fits' 
    sumfile = dldir+'dnidever/nsc/instcal/'+version+'/lists/nsc_calibrate_summary.fits' 
    if os.path.exists(sumfile) == False: 
        if silent==False:
            logger.info(sumfile+' NOT FOUND')
        return 
    logger.info('Loading '+sumfile)
    suminfo = Table.read(sumfile,1)
    for n in suminfo.colnames: suminfo[n].name = n.lower()  # lowercase names
    allinstfilt = np.char.array(suminfo['instrument'].astype(str))+'-'+np.char.array(suminfo['filter'].astype(str))
    gdfilt, = np.where((allinstfilt == instfilt) & (suminfo['nrefmatch'] > 10) &
                       (suminfo['success'] == 1) & (suminfo['fwhm'] < 2.0)) 
    if len(gdfilt) == 0: 
        if silent==False:
            logger.info('No good '+instfilt+' exposures')
        return 
    sumfilt = suminfo[gdfilt]
     
    # Minimum matching radius 
    minradius = 0.43 
    if instrument == 'ksb': 
        minradius = np.maximum(minradius, 0.75)
    if instrument == 'c4d': 
        minradius = np.maximum(minradius, 1.1)
     
    # Calculate healpix number 
    nside = 32
    allpix = hp.ang2pix(nside,sumfilt['ra'],sumfilt['dec'],lonlat=True)
    cenpix = hp.ang2pix(nside,cenra,cendec,lonlat=True)
    ind1,ind2 = dln.match(cenpix,allpix)
    nmatch = len(ind1)
    if nmatch == 0: 
        if silent==False:
            logger.info('No good '+instfilt+' exposure at this position')
        return 
    sumfiltpos = sumfilt[ind2] 
    dist = coords.sphdist(sumfiltpos['ra'],sumfiltpos['dec'],cenra,cendec)
    gddist, = np.where(dist < minradius) 
    if len(gddist) == 0: 
        if silent==False:
            logger.info('No good '+instfilt+' exposures within %0.2f deg' % minradius)
        return 
    refexp = sumfiltpos[gddist] 
    nrefexp = len(refexp) 
    if silent==False:
        logger.info(str(len(gddist))+' good exposures within %0.2f deg' % minradius)
     
    # Pick the best ones 
    si = np.flip(np.argsort(refexp['zptermerr'])) 
    refexp = refexp[si] 
     
    # Loop over the reference exposures 
    #---------------------------------- 
    refcnt = 0
    nrefexpused = 0
    for i in range(nrefexp): 
        # Load the catalog 
        #catfile = str(refexp['expdir'][i]).strip()+'/'+str(refexp['base'][i]).strip()+'_cat.fits'
        catfile = str(refexp['expdir'][i]).strip()+'/'+str(refexp['base'][i]).strip()+'_meas.fits'
        #if catfile[0:4]=='/net' and dldir[0:4]=='/net':
        #    catfile = catfile[4:]
        if catfile[0:4] != '/net':
            catfile = '/net'+catfile
        if catfile.find('/dl1/users/') > -1:
            catfile = catfile.replace('/dl1/users/','/dl2/')
        cat1 = Table.read(catfile)
        # Only want good ones 
        gdcat1, = np.where((cat1['imaflags_iso'] == 0) & (((cat1['flags'] & 8)==0)) & (((cat1['flags'] & 16)==0)) & 
                           (cat1['cmag'] < 50) & (1.087/cat1['cerr'] > 10) & (cat1['fwhm_world']*3600 < 2.5) &  # Remove bad measurements 
                           (~((cat1['x_image'] > 1000) & (cat1['ccdnum'] == 31))))
        # X: 1-1024 okay 
        # X: 1025-2049 bad 
        # use 1000 as the boundary since sometimes there's a sharp drop 
        # at the boundary that causes problem sources with SExtractor 

         
        # Only want sources inside the radius of the exposure 
        if silent==False:
            logger.info(str(i+1)+' '+str(ngdcat1)+'   '+refexp['base'][i])
         
        # Some good sources 
        if ngdcat1 > 0: 
            nrefexpused += 1 
             
            # Convert to small catalog version 
            nnew = ngdcat1 
            dt = [('sourceid',(np.str,100)),('ra',float),('dec',float),('ndet',int),('cmag',float),('cerr',float)]
            newcat = np.zeros(nnew,dtype=np.dtype(dt))
            newcat = Table(newcat)
            for n in newcat.colnames:
                newcat[n] = cat[n][gdcat1]
             
            # First one, start REF structure 
            if len(ref) == 0: 
                # Start ref structure 
                ref = np.zeros(np.maximum(100000,nnew),dtype=np.dtype(dt))
                ref = Table(ref)
                nref = len(ref) 
                # Convert CMAG/CMAGERR to cumulative quantities 
                cmag = 2.5118864**newcat['cmag'] * (1.0/newcat['cerr']**2) 
                cerr = 1.0/newcat['cerr']**2 
                newcat['cmag'] = cmag 
                newcat['cerr'] = cerr 
                # Add to REF 
                ref[0:nnew] = newcat 
                refcnt += nnew 
                 
            # Second or later, add to them 
            else:                  
                # Crossmatch 
                ind1,ind2,dist = coords.xmatch(ref['ra'][0:refcnt],ref['dec'][0:refcnt],newcat['ra'],newcat['dec'],0.5)
                nmatch = len(ind1)
                if silent==False:
                    logger.info('  '+str(nmatch)+' crossmatches')
                 
                # Combine crossmatches 
                if nmatch > 0: 
                    # Cumulative sum of weighted average 
                    ref['cmag'][ind1] += 2.5118864**newcat['cmag'][ind2] * (1.0/newcat['cerr'][ind2]**2) 
                    ref['cerr'][ind1] += 1.0/newcat['cerr'][ind2]**2 
                    ref['ndet'][ind1] += 1 
                    # Remove 
                    if nmatch < len(newcat): 
                        newcat = np.delete(newcat,ind2)
                    else: 
                        newcat = None
                 
                # Add leftover ones 
                if len(newcat) > 0: 
                    if silent==False:
                        logger.info('  Adding '+str(len(newcat))+' leftovers')
                    # Convert CMAG/CMAGERR to cumulative quantities 
                    cmag = 2.5118864**newcat['cmag'] * (1.0/newcat['cerr']**2) 
                    cerr = 1.0/newcat['cerr']**2 
                    newcat['cmag'] = cmag 
                    newcat['cerr'] = cerr 
                    # Add more elements 
                    if refcnt+len(newcat) > nref: 
                        oldref = ref.copy()
                        nnewref = np.maximum(len(newcat),100000)
                        ref = np.zeros(nref+nnewref,dtype=np.dtype(dt))
                        ref = Table(ref)
                        ref[0:nref] = oldref 
                        nref = len(ref) 
                        del oldref
                    # Add new elements to REF 
                    ref[refcnt:refcnt+len(newcat)] = newcat 
                    refcnt += len(newcat) 
         
        # Do we have enough exposures 
        if nrefexpused >= (nrefexp < 10): 
            break

    # No good sources 
    if refcnt == 0: 
        if silent==False:
            logger.info('No good sources')
        return expinfo,chinfo

    # Remove extra elements 
    ref = ref[0:refcnt] 
    # Take average for CMAG and CERR 
    newflux = ref['cmag'] / ref['cerr']
    newmag = 2.5*np.log10(newflux) 
    newerr = np.sqrt(1.0/ref['cerr']) 
    ref['cmag'] = newmag 
    ref['cerr'] = newerr 
     
    # Cross-match with the observed data 
    ind1,ind2,dist = coords.xmatch(ref['ra'],ref['dec'],cat['ra'],cat['dec'],0.5)
    nmatch = len(ind1)
    if nmatch == 0: 
        if silent==False:
            logger.info('No reference sources match to the data')
        return expinfo,chinfo
    cat2 = cat[ind2] 
    ref2 = ref[ind1] 
    # Take a robust mean relative to CMAG 
    mag = cat2['mag_auto'] + 2.5*np.log10(expinfo['exptime'])  # correct for the exposure time 
    err = cat2['magerr_auto']
    model_mag = ref2['cmag']
    col2 = np.zeros(nmatch,float)
    mstr = {'col':col2,'mag':mag,'model':model_mag,'err':err,'ccdnum':cat2['ccdnum']} 
     
     
    # MEASURE THE ZERO-POINT 
    #----------------------- 
    # Same as from fitzpterm()
     
    # Fit the global and ccd zero-points 
    n = len(mstr['col']) 
    diff = mstr['model'] - mstr['mag']
    err = mstr['err']
    # Make a sigma cut 
    med = np.median(diff) 
    sig = dln.mad(diff) 
    gd, = np.where(np.abs(diff-med) < 3*sig) 
    ngd = len(gd)
    x = np.zeros(n,float)
    zpterm, _ = dln.wtmean(diff[gd],err[gd],error=True)
    zptermerr = dln.bootstrap(diff[gd],dln.wtmean,args=err[gd],indexargs=True)
    # Save in exposure structure 
    expinfo['zpterm'] = zpterm 
    expinfo['zptermerr'] = zptermerr 
    expinfo['zptermsig'] = sig 
    expinfo['nrefmatch'] = n 
    expinfo['ngoodrefmatch'] = ngd 
     
    # Measure chip-level zpterm and variations 
    nchips = len(chinfo) 
    for i in range(nchips): 
        ind1,ind2 = dln.match(chinfo['ccdnum'][i],mstr['ccdnum'])
        nmatch = len(ind1)
        if nmatch >= 5: 
            gd1, = np.where(np.abs(diff[ind2]-med) < 3*sig) 
            if len(gd1) >= 5: 
                zpterm1, _ = dln.wtmean(diff[ind2[gd1]],err[ind2[gd1]],error=True)
                zptermerr1 = dln.bootstrap(diff[ind2[gd1]],dln.wtmean,args=err[ind2[gd1]],indexargs=True)
                chinfo['zpterm'][i] = zpterm1 
                chinfo['zptermerr'][i] = zptermerr1 
                chinfo['nmatch'][i] = nmatch 
     
    # Measure spatial variations 
    gdchip, = np.where(chinfo['zpterm'] < 1000) 
    if len(gdchip) > 1: 
        expinfo['zpspatialvar_rms'] = np.std(chinfo['zpterm'][gdchip]) 
        expinfo['zpspatialvar_range'] = np.max(chinfo['zpterm'][gdchip]-np.min(chinfo['zpterm'][gdchip]))
        expinfo['zpspatialvar_nccd'] = len(gdchip)
     
    # Print out the results 
    #if not keyword_set(silent) then begin 
    #  logger.info('NMATCH = ',strtrim(expinfo['ngoodrefmatch'],2) 
    #  logger.info('ZPTERM=',stringize(expinfo['zpterm,ndec=4),'+/-',stringize(expinfo['zptermerr,ndec=4),'  SIG=',stringize(expinfo['zptermsig,ndec=4),'mag' 
    #  logger.info('ZPSPATIALVAR:  RMS=',stringize(expinfo['zpspatialvar_rms,ndec=3),' ',;           'RANGE=',stringize(expinfo['zpspatialvar_range,ndec=3),' NCCD=',strtrim(expinfo['zpspatialvar_nccd,2) 
    #endif 
     
    #dt = systime(1)-t0 
    #if not keyword_set(silent) then logger.info('dt = ',stringize(dt,ndec=2),' sec.' 
     
    return expinfo,chinfo

def loadheader(headfile):
    """
    Load header file
    """

    # FITS file
    if utils.file_isfits(headfile):
        hdu = fits.open(headfile)
        headdict = {}
        # Loop over the extendions
        for i in range(len(hdu)):
            head = hdu[i].header
            if i==0:
                # the main header is blank
                # use the next header which should have filter and other info
                headdict['main'] = hdu[1].header
            else:
                ccdnum = head['CCDNUM']
                headdict[ccdnum] = head
        hdu.close()
    # ASCII file
    else:
        headlines = dln.readlines(headfile)
        lo = dln.grep(headlines,'^SIMPLE  =',index=True)
        begind = dln.grep(headlines,'^XTENSION',index=True)
        begind = lo+begind
        endind = dln.grep(headlines,'^END',index=True)
        headdict = {}
        # Loop over the extendions
        for i in range(len(begind)):
            lines = headlines[begind[i]:endind[i]+1]
            head = fits.Header.fromstring('\n'.join(lines),sep='\n')
            if i==0:
                headdict['main'] = head
            else:
                ccdum = head['CCDNUM']
                headdict[ccdnum] = head
    return headdict
    
def calibrate(expdir,inpref=None,eqnfile=None,redo=False,selfcal=False,
              saveref=False,ncpu=1,photmethod=None,logger=None):
    """
    Perform photometry and astrometric calibration of an NSC exposure using
    external catalogs.

    Parameters
    ----------
    expdir : str
       The absolute path to the NSC exposure directory.
    inpref : astropy table, optional
       Input reference catalog.
    eqnfile : str, optional
       File that contains the model mag equations.
    redo : bool, optional
       Perform the calibration again on this exposure.  Default is False.
    selfcal : bool, optional
       Perform self-calibration instead of using external catalog.
         Default is False.
    saveref : bool, optional
       Save the reference catalog.  Default is False.
    ncpu : int, optional
       Number of cpus to use.  Default is 1.
    logger : logging object
       A logging object used for logging information.

    Returns
    -------
    A calibrated catalog is saved to the exposure directory

    Example
    -------

    calibrate(expdir,inpref,eqnfile)


    Written by D. Nidever in Nov 2017
    Translated to Python by D. Nidever, April 2022
    """
    
    # Calibrate catalogs for one exposure
    dldir,mssdir,localdir = utils.rootdirs()
    
    # Make sure the directory exists 
    if os.path.exists(expdir) == False:
        raise ValueError(expdir+' NOT FOUND')
     
    t00 = time.time() 

    if expdir[-1]=='/':
        expdir = expdir[0:-1]
    base = os.path.basename(expdir) 
    if logger is None:
        logger = dln.basiclogger()
    outfile = expdir+'/'+base+'_meta.fits' 
    # get version number 
    lo = expdir.find('nsc/instcal/') 
    dum = expdir[lo+12:] 
    version = dum[0:dum.find('/')]
    night = expdir.split('/')[-2]

    if photmethod is None:
        photmethod = 'gaiaxpsynth'
    
    logger.info('Calibrate catalogs for exposure '+base+' in '+expdir)
    logger.info('Version = '+version)
     
    # Check for output file 
    if os.path.exists(outfile) == True and redo==False:
        logger.info(outfile+' already exists and redo not set.')
        return

    # What instrument is this? 
    instrument = 'c4d'# by default
    if expdir.find('/k4m/') > -1:
        instrument = 'k4m'
    if expdir.find('/ksb/') > -1:
        instrument = 'ksb' 
    logger.info('This is a '+instrument+' exposure')
    
    # Model magnitude equation file
    if eqnfile is None:
        eqnfile = dldir+'dnidever/nsc/instcal/'+version+'/config/modelmag_equations.txt' 
    logger.info('Using model magnitude equation file '+eqnfile)
    if os.path.exists(eqnfile) == False: 
        raise ValueError(eqnfile+' NOT FOUND')
    eqnlines = dln.readlines(eqnfile)
    for l in eqnlines: logger.info(l)

    # Step 1. Read in the catalogs 
    #----------------------------- 
    logger.info('')
    logger.info('Step 1. Read in the catalogs')
    logger.info('-----------------------------')
    # tar file
    if expdir.endswith('.tar.gz'):
        measfile = None
        logger.info('tar file')
        exposure = os.path.basename(expdir).replace('.tar.gz','')
        result = subprocess.run(['tar','tvf',expdir],capture_output=True)
        res = result.stdout
        if type(res) is bytes: res = res.decode()
        res = res.split('\n')
        tarfiles = [d.split()[-1] for d in res if d!='']
        fitsfiles = [f for f in tarfiles if f.find('.fits')>-1]
        nchips = len(fitsfiles)
        logfile = exposure+'.log'
        logfiletest = logfile in ' '.join(tarfiles)
    # multi-extension meas.fits file
    elif os.path.exists(expdir+'/'+base+'_meas.fits'):
        measfile = expdir+'/'+base+'_meas.fits'
        mhdu = fits.open(measfile)
        nchips = len(mhdu)-1
        # leave the hdu open, we'll use it below
        catfiles = []
        logfile = os.path.join(*[expdir,base+'.log'])
        logfiletest = os.path.exists(logfile)
    # Individual measurement chip files
    else:
        measfile = None
        catfiles = []
        catfiles1 = glob(expdir+'/'+base+'_[1-9].fits')
        if len(catfiles1) > 0:
            catfiles += catfiles1
        catfiles2 = glob(expdir+'/'+base+'_[0-9][0-9].fits') 
        if len(catfiles2) > 0:
            catfiles += catfiles2
        catfiles.sort()
        ncatfiles = len(catfiles)
        if ncatfiles == 0:
            raise ValueError('No catalog files found')
        nchips = ncatfiles
        # check log file
        logfile = os.path.join(*[expdir,base+'.log'])
        logfiletest = os.path.exists(logfile)
    logger.info(str(nchips)+' catalogs found')
    
    if instrument=='c4d' and (nchips<59 or logfiletest==False):
        print('problems with the files')
        return
    
    # Check that this isn't a problematic Mosaic3 exposure 
    if expdir.find('/k4m/') > -1:
        dum = fits.getdata(catfiles[0],1) 
        head0 = dum['field_header_card']
        pixcnt = head0.get('PIXCNT*')
        if pixcnt is not None: 
            logger.info('This is a Mosaic3 exposure with pixel shift problems')
            return 
        wcscal = head0.get('WCSCAL') 
        if wcscal is not None and str(wcscal) == 'Failed': 
            logger.info('This is a Mosaic3 exposure with failed WCS calibration')
            return 

    # v4+ use separate header file
    if version >= 'v4':
        if measfile is not None:
            headfile = os.path.join(expdir,base+'_header.fits')            
        else:
            headfile = os.path.join(expdir,base+'.hdr')
            if os.path.exists(headfile)==False:
                headfile = os.path.join(dldir,'dnidever','nsc','instcal',version,
                                    'header',instrument,night,base+'.hdr')
        if os.path.exists(headfile)==False:
            raise ValueError(headfile+' not found')
        headdict = loadheader(headfile)
    else:
        headdict = None

    if version < 'v4':
        racol = 'alpha_j2000'
        deccol = 'delta_j2000'
        xcol = 'x_image'
        ycol = 'y_image'
    else:
        racol = 'rapsf'
        deccol = 'decpsf'
        xcol = 'xpsf'
        ycol = 'ypsf'

    # Figure out the number of sources 
    ncat = 0 
    ncatarr = np.zeros(nchips,int)
    for i in range(nchips):
        if measfile is None:
            head = fits.getheader(catfiles[i],1)
        else:
            head = mhdu[i+1].header
        ncatarr[i] = head['NAXIS2']
            
    ncat = int(np.sum(ncatarr)) 
    logger.info(str(ncat)+' total sources')
    if ncat == 0: 
        logger.info('No sources')
        return
    
    # Create structure, exten=1 has header now 
    ind, = np.where(ncatarr > 0)
    if measfile is not None:
        meas1 = mhdu[ind[0]+1].data
    else:
        meas1 = fits.getdata(catfiles[ind[0]],1)
    cdt = meas1.dtype.descr
    cdt += [('ccdnum',int),('ebv',float),('ra',float),('dec',float),('raerr',float),
            ('decerr',float),('measid',(str,50)),('filter',(str,50)),('mjd',float)]
    meas = np.zeros(ncat,dtype=np.dtype(cdt))
    meas = Table(meas)
    for c in meas.colnames: meas[c].name = c.lower()
    # Start the chips summary structure
    dt = [('expdir',(str,300)),('instrument',(str,10)),('filename',(str,300)),('measfile',(str,300)),
          ('ccdnum',int),('nsources',int),('nmeas',int),('cenra',float),('cendec',float),('ngaiamatch',int),
          ('ngoodgaiamatch',int),('rarms',float),('rastderr',float),('racoef',(float,4)),('decrms',float),
          ('decstderr',float),('deccoef',(float,4)),('vra',(float,4)),('vdec',(float,4)),
          ('zptype',int),('zpterm',float),('zptermerr',float),('nrefmatch',int),('depth95',float),('depth10sig',float)]
    chinfo = np.zeros(nchips,dtype=np.dtype(dt))
    chinfo = Table(chinfo)
    for c in chinfo.colnames:
        chinfo[c].name = c.lower()
    for c in ['cenra','cendec','rarms','rastderr','decrms','decstderr','zpterm','zptermerr']:
        chinfo[c] = 999999.
    chinfo['depth95'] = 99.99
    chinfo['depth10sig'] = 99.99
    chinfo['expdir'] = expdir 
    chinfo['instrument'] = instrument 
    # Load the files 
    cnt = 0
    for i in range(nchips):
        if measfile is not None:
            meas1 = Table(mhdu[i+1].data)
            hd = mhdu[i+1].header
            ccdnum = hd['ccdnum']
        else:
            dum = os.path.splitext(os.path.basename(catfiles[i]))[0].split('_')
            ccdnum = int(dum[-1])
            hd = fits.getheader(catfiles[i],1)
            meas1 = Table.read(catfiles[i],1)
        for c in meas1.colnames: meas1[c].name = c.lower()
        nmeas1 = hd['naxis2']   # handles empty catalogs 
         
        # Fix negative FWHM values 
        #  use A_WORLD and B_WORLD which are never negative 
        bdfwhm, = np.where(meas1['fwhm_world'] < 0.1)
        if len(bdfwhm) > 0: 
            meas['fwhm_world'][bdfwhm] = np.sqrt(meas1['a_world'][bdfwhm]**2+meas1['b_world'][bdfwhm]**2)*2.35 

        if measfile is None:
            chinfo['filename'][i] = catfiles[i] 
        chinfo['ccdnum'][i] = ccdnum 
        chinfo['nsources'][i] = nmeas1
        # Get the chip corners from full header
        if version < 'v4':
            dum = fits.getdata(catfiles[i],1)
            hd1 = fits.Header.fromstring('\n'.join(dum[0][0]),sep='\n')
        else:
            hd1 = headdict[int(ccdnum)]
        if instrument=='c4d':
            nx = 2046
            ny = 4094
        else:
            nx = hd1['NAXIS1']
            ny = hd1['NAXIS2']
        try:
            wcs1 = WCS(hd1)
        except:
            logger.info('Problem with WCS in header '+catfiles[i])    
            continue
        vcoo = wcs1.pixel_to_world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        chinfo['vra'][i] = vra 
        chinfo['vdec'][i] = vdec 
        if nmeas1 > 0: 
            temp = meas[cnt:cnt+nmeas1]
            for n in meas1.colnames:
                temp[n] = meas1[n]
            temp['ccdnum'] = ccdnum
            temp['ra'] = meas1['alpha_j2000']  # add these here in case the astrometric correction fails later on 
            temp['dec'] = meas1['delta_j2000']
            # Add coordinate uncertainties 
            #   sigma = 0.644*FWHM/SNR 
            #   SNR = 1.087/magerr 
            snr = 1.087/temp['magerr_auto']
            bderr, = np.where((temp['magerr_auto'] > 10) & (temp['magerr_iso'] < 10))
            if len(bderr) > 0: 
                snr[bderr] = 1.087/temp[bderr]['magerr_iso']
            bderr, = np.where((temp['magerr_auto'] > 10) & (temp['magerr_iso'] > 10))
            if len(bderr) > 0: 
                snr[bderr] = 1 
            coorderr = 0.664*(temp['fwhm_world']*3600)/snr 
            temp['raerr'] = coorderr 
            temp['decerr'] = coorderr 
            # Stuff into main structure 
            meas[cnt:cnt+nmeas1] = temp
            cnt += nmeas1
            cenra = np.mean(dln.minmax(meas1['alpha_j2000'])) 
            # Wrapping around RA=0 
            if dln.valrange(meas1['alpha_j2000']) > 100: 
                ra = meas1['alpha_j2000']
                bdra, = np.where(ra > 180) 
                if len(bdra) > 0: 
                    ra[bdra] -= 360 
                bdra2, = np.where(ra < -180) 
                if len(bdra2) > 0: 
                    ra[bdra2] += 360 
                cenra = np.mean(dln.minmax(ra)) 
                if cenra < 0: 
                    cenra += 360 
            chinfo['cenra'][i] = cenra 
            chinfo['cendec'][i] = np.mean(dln.minmax(meas1['delta_j2000']))
    
    if measfile is not None:
        mhdu.close()
            
    # Sort by ccdnum
    chinfo.sort('ccdnum')
    
    # Exposure level values 
    gdchip, = np.where((chinfo['nsources'] > 0) & (chinfo['cenra'] < 400))
    if len(gdchip) == 0: 
        logger.info('No good chip catalogs with good WCS.')
        return 
    # Central coordinates of the entire field 
    cendec = np.mean(dln.minmax(chinfo['cendec'][gdchip]))
    decrange = dln.valrange(chinfo['vdec'][gdchip])
    cenra = np.mean(dln.minmax(chinfo['cenra'][gdchip]))
    rarange = dln.valrange(chinfo['vra'][gdchip])*np.cos(np.deg2rad(cendec))

    # Wrapping around RA=0 
    if dln.valrange(dln.minmax(chinfo[gdchip]['cenra'])) > 100: 
        ra = chinfo['cenra'][gdchip]
        bdra, = np.where(ra > 180) 
        if len(bdra) > 0: 
            ra[bdra]-=360 
        bdra2, = np.where(ra < -180) 
        if len(bdra2) > 0: 
            ra[bdra2] += 360 
        cenra = np.mean(dln.minmax(ra)) 
        if cenra < 0: 
            cenra += 360 
        # use chip VRA to get RA range 
        vra = chinfo['vra'][gdchip]
        bdra, = np.where(vra > 180) 
        if len(bdra) > 0: 
            vra[bdra] -= 360 
        bdra2, = np.where(vra < -180) 
        if len(bdra2) > 0: 
            vra[bdra2] += 360 
        rarange = dln.valrange(vra)*np.cos(np.deg2rad(cendec))
        rawrap = True
    else:
        rawrap = False 

    logger.info('CENRA  = %.5f' % cenra)
    logger.info('CENDEC = %.5f' % cendec)
    coo = SkyCoord(ra=cenra,dec=cendec,unit='deg')
    glon = coo.galactic.l.deg
    glat = coo.galactic.b.deg
    logger.info('GLON = %.5f' % glon)
    logger.info('GLAT = %.5f' % glat)

    # Number of good sources 
    goodsources, = np.where((meas['imaflags_iso']==0) & (~((meas['flags']& 8)>0)) &
                            (~((meas['flags']&16)>0)) & (meas['mag_auto'] < 50))
    ngoodsources = len(goodsources)
    logger.info('GOOD SRCS = '+str(ngoodsources))
         
    # Measure median seeing FWHM 
    gdmeas, = np.where((meas['mag_auto'] < 50) & (meas['magerr_auto'] < 0.05) & (meas['class_star'] > 0.8))
    medfwhm = np.median(meas['fwhm_world'][gdmeas]*3600.) 
    logger.info('FWHM = %.2f arcsec' % medfwhm)
         
    # Load the logfile and get absolute flux filename
    loglines = dln.readlines(expdir+'/'+base+'.log')
    #ind = dln.grep(loglines,'Step #2: Copying InstCal images from mass store archive',index=True)
    ind = dln.grep(loglines,'Copying InstCal images',index=True)
    fline = loglines[ind[0]+1] 
    lo = fline.find('/archive')
    # make sure the mss1 directory is correct for this server 
    fluxfile = mssdir+str(fline[lo+1:]) 
    wline = loglines[ind[0]+2] 
    lo = wline.find('/archive') 
    wtfile = mssdir+str(wline[lo+1:]) 
    mline = loglines[ind[0]+3] 
    lo = mline.find('/archive') 
    maskfile = mssdir+str(mline[lo+1:])

    # Load the meta-data from the original header
    if headdict is None:
        head = fits.getheader(fluxfile,0)
    else:
        head = headdict['main']
    filterlong = head.get('filter')
    
    if filterlong is None:
        raise ValueError('no filter in header')
        dum = fits.getdata(catfiles[0],1)
        hd1 = dum['field_header_card']
        filterlong = str(hd1.get('filter'))
    if filterlong[0:2] == 'VR': 
        filt = 'VR' 
    else: 
        filt = filterlong[0:1] 
    if filterlong == 'bokr': 
        filt = 'r' 
    expnum = head.get('expnum')
    if instrument == 'ksb':  # Bok doesn't have expnum 
        #DTACQNAM= '/data1/batch/bok/20160102/d7390.0049.fits.fz' 
        dtacqnam = head.get('DTACQNAM')
        if dtacqnam is None:
            logger.info('Cannot create an EXPNUM for this Bok exposure')
            return 
        bokbase = os.path.basename(dtacqnam) 
        dum = bokbase.split('.')
        boknight = dum[0][1:] 
        boknum = dum[1] 
        expnum = boknight+boknum  # concatenate the two numbers 
    exptime = head.get('exptime')
    if exptime is None:
        raise ValueError('no exptime in header')
        dum = fits.getdata(catfiles[0],1)
        hd1 = dum['field_header_card']
        exptime = hd1['exptime']
    dateobs = head['date-obs']
    airmass = head['airmass']
    t = Time(dateobs, format='fits')
    mjd = t.mjd
    #mjd = date2jd(dateobs,/mjd) 
    logger.info('FILTER = '+filt)
    logger.info('EXPTIME = %.2f sec.' % exptime)
    logger.info('MJD = %.4f' % mjd)
    
    # Set some catalog values 
    meas['filter'] = filt
    meas['mjd'] = mjd 
    meas['measid'] = instrument+'.'+str(expnum)+'.'+np.char.array(meas['ccdnum'].astype(str))+'.'+np.char.array(meas['number'].astype(str))
    
    # Start the exposure-level table
    edt = [('file',(str,300)),('wtfile',(str,300)),('maskfile',(str,300)),('instrument',(str,10)),
           ('base',(str,100)),('expnum',int),('ra',float),('dec',float),('dateobs',(str,50)),('mjd',float),
           ('filter',(str,50)),('exptime',float),('airmass',float),('wcscal',(str,50)),('nsources',int),
           ('ngoodsources',int),('nmeas',int),('fwhm',float),('nchips',int),('rarms',float),('decrms',float),
           ('ebv',float),('ngaiamatch',int),('ngoodgaiamatch',int),('zptype',int),('zpterm',float),('zptermerr',float),
           ('zptermsig',float),('zpspatialvar_rms',float),('zpspatialvar_range',float),('zpspatialvar_nccd',int),
           ('nrefmatch',int),('ngoodrefmatch',int),('depth95',float),('depth10sig',float)]
    expinfo = np.zeros(1,dtype=np.dtype(edt))
    expinfo = Table(expinfo)
    expinfo['file'] = fluxfile
    expinfo['wtfile'] = wtfile
    expinfo['maskfile'] = maskfile
    expinfo['base'] = base
    expinfo['expnum'] = int(expnum)
    expinfo['dateobs'] = str(dateobs)
    expinfo['filter'] = filt
    expinfo['exptime'] = float(exptime)
    expinfo['nsources'] = int(ncat)
    for c in ['zpterm','zptermerr','zptermsig','zpspatialvar_rms','zpspatialvar_range']:
        expinfo[c] = 999999.
    expinfo['depth95'] = 99.99
    expinfo['depth10sig'] = 99.99
    expinfo['instrument'] = instrument 
    expinfo['ra'] = cenra 
    expinfo['dec'] = cendec 
    expinfo['mjd'] = mjd 
    #expinfo['mjd'] = utils.getmjd('','CTIO',dateobs=dateobs) 
    expinfo['nchips'] = nchips 
    expinfo['airmass'] = airmass 
    expinfo['fwhm'] = medfwhm 
    expinfo['ngoodsources'] = ngoodsources 
    wcscal = head.get('wcscal')
    if wcscal is None:
        wcscal = 'NAN' 
    expinfo['wcscal'] = wcscal 
    
    # Step 2. Load the reference catalogs 
    #------------------------------------ 
    logger.info('')
    logger.info('Step 2. Load the reference catalogs')
    logger.info('------------------------------------') 
         
    # Getting reference catalogs 
    if inpref is None: 
        # Search radius 
        radius = 1.1 * np.sqrt( (0.5*rarange)**2 + (0.5*decrange)**2 )
        ref = query.getrefdata(instrument+'-'+filt,cenra,cendec,radius) 
        if len(ref) == 0: 
            logger.info('No Reference Data')
            return 
    # Using input reference catalog 
    else: 
        logger.info('Reference catalogs input')
        if rawrap == False: 
            gdref, = np.where((inpref['ra'] >= np.min(meas[racol])-0.01) & 
                              (inpref['ra'] <= np.max(meas[racol])+0.01) &
                              (inpref['dec'] >= np.min(meas[deccol])-0.01) & 
                              (inpref['dec'] <= np.max(meas[deccol])+0.01))
        else: 
            ra = meas[racol]
            bdra, = np.where(ra > 180) 
            if len(bdra) > 0 : 
                ra[bdra]-=360 
            gdref, = np.where((inpref['ra'] <= np.max(ra)-0.01) & 
                              (inpref['ra'] >= np.min(ra+360)-0.01) &
                              (inpref['dec'] >= np.min(meas[deccol])-0.01) & 
                              (inpref['dec'] <= np.max(meas[deccol])+0.01))
        ref = inpref[gdref] 
        logger.info(str(ngdref)+' reference stars in our region')

    # Step 3. Astrometric calibration 
    #---------------------------------- 
    # At the chip level, linear fits in RA/DEC 
    logger.info('')
    logger.info('Step 3. Astrometric calibration')
    logger.info('--------------------------------') 
    # Get reference catalog with Gaia values 
    gdgaia, = np.where(ref['source'] > 0) 
    gaia = ref[gdgaia]
    # Match everything to Gaia at once, this is much faster!
    gdmeas, = np.where(np.isfinite(meas[racol]) & np.isfinite(meas[deccol]))
    if len(gdmeas)==0:
        logger.info('No catalog stars with good PSF coordinates')
        return
    ind1,ind2,dist = coords.xmatch(gaia['ra'],gaia['dec'],meas[racol][gdmeas],meas[deccol][gdmeas],1.0)
    ngmatch = len(ind1)
    if ngmatch == 0: 
        logger.info('No gaia matches')
        return 
    allgaiaind = np.zeros(ncat,int)-1
    allgaiaind[gdmeas[ind2]] = ind1 
    allgaiadist = np.zeros(ncat,float)+999999. 
    allgaiadist[gdmeas[ind2]] = coords.sphdist(gaia['ra'][ind1],gaia['dec'][ind1],
                                              meas[racol][gdmeas[ind2]],meas[deccol][gdmeas[ind2]])*3600 
    # CCD loop 
    for i in range(nchips): 
        if chinfo['nsources'][i]==0: 
            continue
        # Get chip sources using CCDNUM 
        chind1,chind2 = dln.match(chinfo['ccdnum'][i],meas['ccdnum'])
        nchmatch = len(chind1)
        meas1 = meas[chind2] 
        # Gaia matches for this chip 
        gaiaind1 = allgaiaind[chind2] 
        gaiadist1 = allgaiadist[chind2] 
        gmatch, = np.where((gaiaind1 > -1) & (gaiadist1 <= 0.5))   # get sources with Gaia matches 
        if len(gmatch) == 0: 
            gmatch, = np.where(gaiaind1 > -1 and gaiadist1 <= 1.0) 
        if len(gmatch) < 5: 
            logger.info('Not enough Gaia matches')
            # Add threshold to astrometric errors 
            meas1['raerr'] = np.sqrt(meas1['raerr']**2 + 0.100**2) 
            meas1['decerr'] = np.sqrt(meas1['decerr']**2 + 0.100**2) 
            meas[chind2] = meas1 
            continue
        #gaia1b = gaia[ind1] 
        #meas1b = meas1[ind2] 
        gaia1b = gaia[gaiaind1[gmatch]] 
        meas1b = meas1[gmatch]
        # Apply quality cuts 
        #  no bad CP flags 
        #  no SE truncated or incomplete data flags 
        #  must have good photometry
        if 'pmra' in gaia1b.colnames and 'pmdec' in gaia1b.colnames:  # we have proper motion information 
            qcuts1, = np.where((meas1b['imaflags_iso'] == 0) & (~((meas1b['flags'] & 8) > 0)) & (~((meas1b['flags'] & 16) > 0)) & 
                               (meas1b['mag_auto'] < 50) & np.isfinite(gaia1b['pmra']) & np.isfinite(gaia1b['pmdec']))
        else:
            qcuts1, = np.where((meas1b['imaflags_iso'] == 0) & (~((meas1b['flags'] & 8) > 0)) & (~((meas1b['flags'] & 16) > 0)) & 
                               (meas1b['mag_auto'] < 50))
        if len(qcuts1) == 0: 
            logger.info('Not enough stars after quality cuts')
            # Add threshold to astrometric errors 
            meas1['raerr'] = np.sqrt(meas1['raerr']**2 + 0.100**2) 
            meas1['decerr'] = np.sqrt(meas1['decerr']**2 + 0.100**2) 
            meas[chind2] = meas1 
            continue
        gaia2 = gaia1b[qcuts1] 
        meas2 = meas1b[qcuts1] 
             
        # Precess the Gaia coordinates to the epoch of the observation 
        # The reference epoch for Gaia DR2 is J2015.5 (compared to the 
        # J2015.0 epoch for Gaia DR1). 

        # https://www.cosmos.esa.int/web/gaia/earlydr3#:~:text=The%20reference%20epoch%20for%20Gaia,full%20Gaia%20DR3)%20is%20J2016.
        # The reference epoch for Gaia DR3 (both Gaia EDR3 and the full Gaia DR3) is J2016.0.

        if 'pmra' in gaia2.colnames and 'pmdec' in gaia2.colnames:
            #gaiamjd = 57206.0   # 7/3/2015, J2015.5
            gaiamjd = 57388.0   # 1/1/2016, J2016.0, EDR3, DR3
            delt = (mjd-gaiamjd)/365.242170 # convert to years 
            # convert from mas/yr->deg/yr and convert to angle in RA 
            gra_epoch = gaia2['ra'] + delt*gaia2['pmra']/3600.0/1000.0/np.cos(np.deg2rad(gaia2['dec'])) 
            gdec_epoch = gaia2['dec'] + delt*gaia2['pmdec']/3600.0/1000.0 
        else: 
            gra_epoch = gaia2['ra']
            gdec_epoch = gaia2['dec'] 

        # Rotate to coordinates relative to the center of the field 
        gaialon,gaialat = coords.rotsphcen(gra_epoch,gdec_epoch,chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        lon1,lat1 = coords.rotsphcen(meas2[racol],meas2[deccol],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        # ---- Fit RA as function of RA/DEC ---- 
        londiff = gaialon-lon1 
        err = None
        if 'ra_error' in gaia2.colnames:
            err = np.sqrt(gaia2['ra_error']**2 + meas2['raerr']**2) 
        if 'ra_error' not in gaia2.colnames and 'e_ra_icrs' in gaia2.colnames:
            err = np.sqrt(gaia2['e_ra_icrs']**2 + meas2['raerr']**2) 
        if err is None:
            err = meas2['raerr']
        lonmed = np.median(londiff) 
        lonsig = np.maximum(dln.mad(londiff), 1e-5)  # 0.036" 
        gdlon, = np.where(np.abs(londiff-lonmed) < 3.0*lonsig)# remove outliers 
        if len(gdlon) > 5:   # use constant if not enough stars 
            npars = 4 
        else: 
            npars = 1 
        initpars = np.zeros(npars,float) 
        initpars[0] = np.median([londiff]) 
        racoef,racov = curve_fit(utils.func_poly2d_wrap,[lon1[gdlon],lat1[gdlon]],
                                 londiff[gdlon],sigma=err[gdlon],p0=initpars)
        yfitall = utils.func_poly2d(lon1,lat1,*racoef)
        rarms1 = dln.mad((londiff[gdlon]-yfitall[gdlon])*3600.) 
        rastderr = rarms1/np.sqrt(len(gdlon))
        # Use bright stars to get a better RMS estimate 
        gdstars, = np.where((meas2['fwhm_world']*3600 < 2*medfwhm) & (1.087/meas2['magerr_auto'] > 50))
        if len(gdstars) < 20: 
            gdstars, = np.where((meas2['fwhm_world']*3600 < 2*medfwhm) & (1.087/meas2['magerr_auto'] > 30))
        if len(gdstars) > 5: 
            diff = (londiff-yfitall)*3600. 
            rarms = dln.mad(diff[gdstars]) 
            rastderr = rarms/np.sqrt(len(gdstars))
        else:
            rarms = rarms1 
        # ---- Fit DEC as function of RA/DEC ----- 
        latdiff = gaialat-lat1 
        err = None
        if 'dec_error' in gaia2.colnames: 
            err = np.sqrt(gaia2['dec_error']**2 + meas2['decerr']**2) 
        if 'dec_error' not in gaia2.colnames and 'e_dec_icrs' in gaia2.colnames:
            err = np.sqrt(gaia2['e_de_icrs']**2 + meas2['decerr']**2) 
        if err is None: 
            err = meas2['decerr']
        latmed = np.median(latdiff) 
        latsig = np.maximum(dln.mad(latdiff), 1e-5)   # 0.036" 
        gdlat, = np.where(np.abs(latdiff-latmed) < 3.0*latsig)  # remove outliers 
        if len(gdlat) > 5:     # use constant if not enough stars 
            npars = 4 
        else: 
            npars = 1 
        initpars = np.zeros(npars,float) 
        initpars[0] = np.median(latdiff) 
        deccoef,deccov = curve_fit(utils.func_poly2d_wrap,[lon1[gdlat],lat1[gdlat]],
                                   latdiff[gdlat],sigma=err[gdlat],p0=initpars)
        yfitall = utils.func_poly2d(lon1,lat1,*deccoef)
        decrms1 = dln.mad((latdiff[gdlat]-yfitall[gdlat])*3600.) 
        decstderr = decrms1/np.sqrt(len(gdlat))
        # Use bright stars to get a better RMS estimate 
        if len(gdstars) > 5: 
            diff = (latdiff-yfitall)*3600. 
            decrms = dln.mad(diff[gdstars]) 
            decstderr = decrms/np.sqrt(len(gdstars))
        else:
            decrms = decrms1 
        logger.info('  CCDNUM=%3d  NSOURCES=%5d  %5d/%5d GAIA matches  RMS(RA/DEC)=%7.4f/%7.4f STDERR(RA/DEC)=%7.4f/%7.4f arcsec' % 
                    (chinfo['ccdnum'][i],nchmatch,len(gmatch),len(qcuts1),rarms,decrms,rastderr,decstderr))
        # Apply to all sources 
        lon,lat = coords.rotsphcen(meas1['alpha_j2000'],meas1['delta_j2000'],
                                   chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        lon2 = lon + utils.func_poly2d(lon,lat,*racoef)
        lat2 = lat + utils.func_poly2d(lon,lat,*deccoef) 
        ra2,dec2 = coords.rotsphcen(lon2,lat2,chinfo['cenra'][i],chinfo['cendec'][i],reverse=True,gnomic=True)
        meas1['ra'] = ra2 
        meas1['dec'] = dec2
        if version >= 'v4':
            gdmeas, = np.where(np.isfinite(meas1['rapsf']) & np.isfinite(meas1['decpsf']))
            lon,lat = coords.rotsphcen(meas1['rapsf'][gdmeas],meas1['decpsf'][gdmeas],
                                       chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
            lon2 = lon + utils.func_poly2d(lon,lat,*racoef)
            lat2 = lat + utils.func_poly2d(lon,lat,*deccoef) 
            ra2,dec2 = coords.rotsphcen(lon2,lat2,chinfo['cenra'][i],chinfo['cendec'][i],reverse=True,gnomic=True)
            meas1['rapsf'][gdmeas] = ra2 
            meas1['decpsf'][gdmeas] = dec2
        # Add to astrometric errors 
        meas1['raerr'] = np.sqrt(meas1['raerr']**2 + rarms**2) 
        meas1['decerr'] = np.sqrt(meas1['decerr']**2 + decrms**2) 
        # Stuff back into the main structure 
        meas[chind2] = meas1 
        chinfo['ngaiamatch'][i] = len(gmatch)
        chinfo['ngoodgaiamatch'][i] = len(qcuts1) 
        chinfo['rarms'][i] = rarms 
        chinfo['rastderr'][i] = rastderr 
        chinfo['racoef'][i] = racoef 
        chinfo['decrms'][i] = decrms 
        chinfo['decstderr'][i] = decstderr 
        chinfo['deccoef'][i] = deccoef 
    
    # Get reddening
    coo = SkyCoord(ra=meas['ra'],dec=meas['dec'],unit='deg')
    sfd = SFDQuery()
    ebv = sfd(coo)                
    meas['ebv'] = ebv 
                 
    # Put in exposure-level information 
    expinfo['rarms'] = np.median(chinfo['rarms']) 
    expinfo['decrms'] = np.median(chinfo['decrms']) 
    expinfo['ebv'] = np.median(ebv) 
    #expinfo['gaianmatch'] = median(chinfo['gaianmatch']) 
    expinfo['ngaiamatch'] = np.sum(chinfo['ngaiamatch']) 
    expinfo['ngoodgaiamatch'] = np.sum(chinfo['ngoodgaiamatch']) 
    
    # Step 4. Photometric calibration 
    #-------------------------------- 
    logger.info('')
    logger.info('Step 4. Photometric calibration')
    logger.info('-------------------------------')
    instfilt = instrument+'-'+filt  # instrument-filter combination 
                 
    # Now crossmatch with our catalog 
    dcr = 1.0
    gdmeas, = np.where(np.isfinite(meas['ra']) & np.isfinite(meas['dec']))
    ind1,ind2,dist = coords.xmatch(ref['ra'],ref['dec'],meas['ra'][gdmeas],meas['dec'][gdmeas],dcr)
    nmatch = len(ind1)
    logger.info(str(nmatch)+' matches to reference catalog')
    if nmatch == 0: 
        logger.info('No matches to reference catalog')
        return
    ref1 = ref[ind1]
    meas1 = meas[gdmeas[ind2]]
    # Use Gaia XP synthetic photometry
    # Get the model magnitudes 
    mmags = modelmag.modelmag(ref1,instfilt,cendec,eqnfile) 
    if len(mmags) == 1 and mmags[0] < -1000: 
        print('No good model mags')
        return
    
    # Get the zero-points
    mmexpinfo,mmchinfo = getzpterm(meas1,ref1,mmags,expinfo.copy(),chinfo.copy(),kind='modelmag')
    gexpinfo,gchinfo = getzpterm(meas1,ref1,mmags,expinfo.copy(),chinfo.copy(),kind='gaiaxpsynth')    

    print('Using Gaia for everything now')
    expinfo = gexpinfo.copy()
    chinfo = gchinfo.copy()
    
    # Apply the zero-point to the full catalogs 
    # USE CHIP-LEVEL ZERO-POINTS WHEN POSSIBLE!!! 
    # Create an output catalog for each chip 
    nsrc = np.cumsum(chinfo['nsources'])
    lo = np.append(0,nsrc[0:nchips])
    hi = nsrc-1 
    for i in range(nchips): 
        nmeas1 = hi[i]-lo[i]+1 
        if nmeas1 > 0: 
            meas1 = meas[lo[i]:hi[i]] 
            #; Use chip-level zero-point 
            #if chinfo[i].nrefmatch gt 5 then begin 
            #  chinfo[i].zptype = 1 
            #; Use exposure-level zero-point 
            #endif else begin 
            #  chinfo[i].zpterm = expinfo['zpterm']
            #  chinfo[i].zptermerr = expinfo['zptermerr']
            #  chinfo[i].zptype = 2 
            #endelse 
            # Always use exposure-level zero-points,  they are good enough 
            chinfo['zpterm'][i] = expinfo['zpterm']
            chinfo['zptermerr'][i] = expinfo['zptermerr']
            #chinfo['zptype'][i] = 2 
                     
            #gdmeasmag, = np.where(meas1['magpsf'] < 50) 
            #meas1['cmag'][gdmeasmag] = meas1['magpsf'][gdmeasmag] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            #meas1['cerr'][gdmeasmag] = np.sqrt(meas1['errpsf'][gdmeasmag]**2 + chinfo['zptermerr'][i]**2) # add calibration error in quadrature 
            # Print out the results 
            logger.info('  CCDNUM=%3d  NREFSOURCES=%5d  ZPTYPE=%2d  ZPTERM=%7.4f+/-%7.4f' %
                        (chinfo['ccdnum'][i],chinfo['nrefmatch'][i],chinfo['zptype'][i],chinfo['zpterm'][i],chinfo['zptermerr'][i]))
            #meas[lo[i]:hi[i]] = meas1  # stuff back in 
             
    # Print out the results 
    logger.info('NPHOTREFMATCH=%d' % expinfo['nrefmatch'])
    logger.info('EXPOSURE ZPTERM=%.4f +/- %.4f  SIG=%.4f mag' % (expinfo['zpterm'],expinfo['zptermerr'],expinfo['zptermsig']))
    logger.info('ZPSPATIALVAR:  RMS=%.4f RANGE=%.4f NCCD=%d' % 
                (expinfo['zpspatialvar_rms'],expinfo['zpspatialvar_range'],expinfo['zpspatialvar_nccd']))

             
    # Measure the depth 
    #   need good photometry 
    gdmag, = np.where(meas['magpsf'] < 50) 
    ngdmag = len(gdmag)
    if ngdmag > 0: 
        # Get 95% percentile depth 
        cmag = meas['magpsf'][gdmag] + 2.5*np.log10(exptime) + expinfo['zpterm'][0]
        cerr = np.sqrt(meas['errpsf'][gdmag]**2 + expinfo['zptermerr'][0]**2) # add calibration error in quadrature 
        si = np.argsort(cmag) 
        cmag = cmag[si]
        cerr = cerr[si]
        depth95 = cmag[int(np.round(0.95*ngdmag))-1] 
        expinfo['depth95'] = depth95 
        chinfo['depth95'] = depth95 
        logger.info('95% percentile depth'+' = %.2f mag' % depth95)
        # Get 10 sigma depth 
        #  S/N = 1.087/err 
        #  so S/N=5 is for err=1.087/5=0.2174 
        #  S/N=10 is for err=1.087/10=0.1087 
        depth10sig = 99.99 
        depind, = np.where((cmag < 50) & (cmag > depth95-3.0) & (cerr >= 0.0987) & (cerr <= 0.1187))
        if len(depind) < 5: 
            depind, = np.where((cmag < 50) & (cmag > depth95-3.0) & (cerr >= 0.0787) & (cerr <= 0.1387))
        if len(depind) > 5: 
            depth10sig = np.median(cmag[depind]) 
        else: 
            depind, = np.where(cmag < 50) 
            if len(depind) > 0: 
                depth10sig = np.max(cmag[depind]) 
        logger.info('10sigma depth = %.2f mag' % depth10sig)
        expinfo['depth10sig'] = depth10sig 
        chinfo['depth10sig'] = depth10sig 
        
    # Step 5. Write out the final catalogs and metadata 
    #-------------------------------------------------- 
    #if redo and selfcal and expinfo['ztype']==2:
    #    # Create backup of original versions 
    #    logger.info('Copying meas and meta files to v1 versions')
    #    metafile = expdir+'/'+base+'_meta.fits' 
    #    if os.path.exists(metafile): 
    #        dln.file_copy(metafile,expdir+'/'+base+'_meta.v1.fits',overwrite=True)


    hdu = fits.HDUList()   # main table
    mhdu = fits.HDUList()  # meta-data
    mhdu.append(fits.table_to_hdu(expinfo))

    
    # Create an output catalog for each chip 
    nsrc = np.cumsum(chinfo['nsources'])
    lo = np.append(0,nsrc[0:nchips]) 
    hi = nsrc-1 
    for i in range(nchips): 
        nmeas1 = hi[i]-lo[i]+1 
        if nmeas1 == 0: 
            logger.info('No sources for CCDNUM='+str(chinfo['ccdnum'][i]))
            continue
        meas1 = meas[lo[i]:hi[i]+1] 
        logger.info('CCDNUM={:4d}  NSOURCES={:d}'.format(meas1['ccdnum'][0],len(meas1)))
        
        # Apply QA cuts 
        #---------------- 
        badflag = np.zeros(len(meas1),int)
        
        # Mask bad chip data 
        # Half of chip 31 for MJD>56660 
        #  c4d_131123_025436_ooi_r_v2 with MJD=56619 has problems too 
        #  if the background b/w the left and right half is large then BAd 
        lft31, = np.where((meas1[xcol] < 1024) & (meas1['ccdnum'] == 31))
        rt31, = np.where((meas1[xcol] >= 1024) & (meas1['ccdnum'] == 31))
        if len(lft31) > 10 and len(rt31) > 10: 
            lftback = np.median(meas1['background'][lft31]) 
            rtback = np.median(meas1['background'][rt31]) 
            mnback = 0.5*(lftback+rtback) 
            sigback = dln.mad(meas1['background']) 
            if np.abs(lftback-rtback) > (np.sqrt(mnback)>sigback): 
                jump31 = True
            else: 
                jump31 = False
            #logger.info('  Big jump in CCDNUM 31 background levels' 
        else:
            jump31 = False 
        if expinfo['mjd'] > 56600 or jump31 == 1: 
            badchip31 = True 
            # Remove bad measurements 
            # X: 1-1024 okay 
            # X: 1025-2049 bad 
            # use 1000 as the boundary since sometimes there's a sharp drop 
            # at the boundary that causes problem sources with SExtractor 
            bdind, = np.where((meas1[xcol] > 1000) & (meas1['ccdnum'] == 31))
            ngdind = len(meas1)-len(bdind)
            if len(bdind) > 0:  # some bad ones found
                logger.info('  '+str(len(bdind))+' bad chip 31 measurements')
                badflag[bdind] += 1
            else:
                badchip31 = False  # chip 31 
                         
        # Mask sourcds with bad quality mask flag (IMAFLAGS_ISO) 
        bdmeas, = np.where(meas1['imaflags_iso'] > 0)
        if len(bdmeas) > 0: 
            logger.info('  '+str(len(bdmeas))+' sources with bad CP flags.')
            badflag[bdmeas] += 2 
            
        # Mask sources with bad SE FLAGS 
        #   this masks problematic truncatd sources near chip edges 
        bdseflags, = np.where( ((meas1['flags'] & 8) > 0) |   # object truncated 
                               ((meas1['flags'] & 16) > 0))   # aperture truncate 
        if len(bdseflags) > 0: 
            logger.info('  '+str(len(bdseflags))+' truncated sources')
            badflag[bdseflags] += 4
            
        # Mask low-S/N sources 
        #  snr = 1.087/err 
        snrcut = 5.0 
        bdsnr, = np.where(1.087/meas1['magerr_auto'] < snrcut) 
        if len(bdsnr) > 0: 
            logger.info('  '+str(len(bdsnr))+' sources with S/N<'+str(snrcut))
            badflag[bdsnr] += 8
            
        # Convert to final format 
        if nmeas1 > 0:
            mdt = [('measid',(str,100)),('objectid',(str,100)),('exposure',(str,50)),
                   ('ccdnum',np.int8),('filter',(str,50)),('mjd',float),('x',np.float32),('y',np.float32),
                   ('ra',float),('raerr',np.float32),('dec',float),('decerr',np.float32),('mag_auto',np.float32),
                   ('magerr_auto',np.float32),('mag_aper1',np.float32),('magerr_aper1',np.float32),('mag_aper2',np.float32),
                   ('magerr_aper2',np.float32),('mag_aper4',np.float32),('magerr_aper4',np.float32),
                   ('mag_aper6',np.float32),('magerr_aper6',np.float32),('mag_aper8',np.float32),
                   ('magerr_aper8',np.float32),('mag_iso',np.float32),('magerr_iso',np.float32),
                   ('kron_radius',np.float32),('background',np.float32),('threshold',np.float32),('isoarea_image',np.float32),
                   ('isoarea_world',np.float32),('x2_world',np.float32),('y2_world',np.float32),('xy_world',np.float32),
                   ('asemi',np.float32),('asemierr',np.float32),('bsemi',np.float32),('bsemierr',np.float32),
                   ('theta',np.float32),('thetaerr',np.float32),('ellipticity',np.float32),
                   ('errx2_world',np.float32),('erry2_world',np.float32),('errxy_world',np.float32),
                   ('fwhm',np.float32),('flags',np.int16),('imaflags_iso',np.int32),('nimaflags_iso',np.int32),
                   ('class_star',np.float32),
                   ('ndet_iter',int),('repeat',int),('xpsf',float),('ypsf',float),('magpsf',float),
                   ('errpsf',float),('skypsf',float),('iter',int),('chi',float),('sharp',float),
                   ('rapsf',float),('decpsf',float),('ebv',np.float32),('haspsf',bool),('snr',np.float32),
                   ('badflag',int)]
            final = np.zeros(nmeas1,dtype=np.dtype(mdt))
            final = Table(final)
            # Add units
            for c in ['ra','dec','rapsf','decpsf','theta','thetaerr']:
                final[c].unit = 'degree'
            for c in ['raerr','decerr','fwhm','asemi','asemierr','bsemi','bsemierr','kron_radius']:
                final[c].unit = 'arcsec'
            for c in ['mag_auto','magerr_auto','mag_aper1','magerr_aper1','mag_aper2','magerr_aper2',
                      'mag_aper4','magerr_aper4','mag_aper6','magerr_aper6','mag_aper8','magerr_aper8',
                      'mag_iso','magerr_iso','magpsf','errpsf','ebv']:
                final[c].unit = 'magnitude'
            for c in ['x','y','xpsf','ypsf']:
                final[c].unit = 'pixel'
            for c in ['background','skypsf']:
                final[c].unit = 'ADU'
            final['mjd'].unit = 'day'
            # Copy over all columns in common
            for n in final.colnames:
                if n in meas1.colnames:
                    final[n][:] = meas1[n]
            # Special handling of some columns
            final['exposure'][:] = base
            final['x'][:] = meas1['x_image']
            final['y'][:] = meas1['y_image']
            # Calibrate the photometry
            #  only calibrate good photometric measurements
            zpoff = 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            for c in ['mag_auto','mag_aper1','mag_aper2','mag_aper4','mag_aper6','mag_aper8','magpsf']:
                if c[:8]=='mag_aper':
                    indx = {'1':0,'2':1,'4':2,'6':3,'8':4}[c[-1]]
                    mag = meas1['mag_aper'][:,indx]
                else:
                    mag = meas1[c]
                gdmag = ((mag < 50) & np.isfinite(mag))
                final[c][:] = np.nan   # bad by default
                final[c][gdmag] = mag[gdmag]+zpoff
            final['magerr_auto'][:] = meas1['magerr_auto']
            final['magerr_aper1'][:] = meas1['magerr_aper'][:,0] 
            final['magerr_aper2'][:] = meas1['magerr_aper'][:,1] 
            final['magerr_aper4'][:] = meas1['magerr_aper'][:,2]
            final['magerr_aper6'][:] = meas1['magerr_aper'][:,3] 
            final['magerr_aper8'][:] = meas1['magerr_aper'][:,4]
            final['asemi'][:] = meas1['a_world'] * 3600.         # convert to arcsec 
            final['asemierr'][:] = meas1['erra_world'] * 3600.   # convert to arcsec 
            final['bsemi'][:] = meas1['b_world'] * 3600.         # convert to arcsec 
            final['bsemierr'][:] = meas1['errb_world'] * 3600.   # convert to arcsec 
            final['theta'][:] = 90-meas1['theta_world']          # make CCW E of N 
            final['thetaerr'][:] = meas1['errtheta_world'] 
            final['fwhm'][:] = meas1['fwhm_world'] * 3600.       # convert to arcsec 
            final['skypsf'][:] = meas1['sky']
            final['haspsf'][:] = meas1['magpsf'] < 50
            final['snr'][:] = 1.087/meas1['magerr_auto']
            final['badflag'][:] = badflag
            # Other columns we are keeping in case we need to recreate the original measurement file
            # already copied over in for loop above
            
            #if redo and selfcal and expinfo['zptyper']==2: # Create backup of original versions 
            #    if os.path.exists(outfile) == 1 : 
            #        dln.file_copy(outfile,expdir+'/'+base+'_'+str(chinfo[i].ccdnum,2)+'_meas.v1.fits',overwrite=True)

            # Append to the output HDUList
            hdu1 = fits.table_to_hdu(final)
            hdu1.header['EXTNAME'] = ccdnum
            hdu.append(hdu1)
            mhdu1 = fits.table_to_hdu(chinfo[i:i+1])
            mhdu1.header['EXTNAME'] = ccdnum
            mhdu.append(mhdu1)                    # add metadata for this chip

    # Write to file 
    outfile = expdir+'/'+base+'_meas2.fits'
    logger.info('Writing table to '+outfile)    
    hdu.writeto(outfile,overwrite=True)
    hdu.close()
                     
    # Meta-data file 
    metafile = expdir+'/'+base+'_meta2.fits' 
    logger.info('Writing metadata to '+metafile)
    mhdu.writeto(metafile,overwrite=True)
    mhdu.close()
            
    dt = time.time()-t00 
    logger.info('dt = %.2f sec.' % dt)

def calibrate_healpix(pix,version,nside=64,redo=False):
    """
    This program is a wrapper around NSC_INSTCAL_CALIBRATE
    for all exposures in the same region of the sky.
    It retrieves the reference catalogs ONCE which saves time.

    Parameters
    ----------
    pix       The HEALPix pixel number.
    version   The version name, e.g. 'v3'. 
    =nside    The HEALPix nside to use.  Default is 64.
    /redo     Rerun on exposures that were previously processed.

    Returns
    -------

    Example
    -------

    calibrate_healpix(1045,'v3')

    By D. Nidever 2017
    Translated to Python by D. Nidever, May 2022
    """

    # Main NOAO DECam source catalog
    lsdir,mssdir,localdir = utils.rootdirs()
    fdir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
    tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
    if os.path.exists(fdir)==False:
        os.makedirs(fdir+'logs/')
    if os.path.exists(tmpdir)==False:
        os.makedirs(tmpdir)

    t00 = time.time()

    # Load the list of exposures
    listfile = fdir+'/lists/nsc_calibrate_healpix_list.fits'
    if os.path.exists(listfile)==False:
        print(listfile,' NOT FOUND')
        return
    hplist = Table.read(listfile)

    # Get the exposures for this healpix
    print('Calibrating InstCal SExtractor catalogs for Healpix pixel = '+str(pix))
    ind1,ind2 = dln.match(hplist['pix'],pix)
    nind = len(ind1)
    if nind==0:
        print('No exposures')
        return
    print('NEXPOSURES = '+str(nind))
    hplist1 = hplist[ind]
    hplist1['expdir'] = np.char.array(hplist1['expdir']).strip()
    hplist1['instrument ']= np.char.array(hplist1['instrument']).strip()
    hplist1['filter'] = np.char.array(hplist1['filter']).strip()

    # Central coordinates
    cenra,cendec = hp.pix2ang(nside,pix,lonlat=True)
    print('RA  = %.6f' % cenra)
    print('DEC = %.6f' % cendec)
    cencoo = SkyCoord(ra=cenra,dec=cendec,unit='deg')
    glon = cencoo.galactic.l.degree
    glat = cencoo.galactic.b.degree
    print('L = %.6f' % glon)
    print('B = %.6f' % glat)

    # List of instrument-filters
    filters = np.char.array(hplist1['instrument']).strip()+'-'+np.char.array([f.strip()[0:2] for f in hplist1['filter']]).strip()
    filters = np.unique(filters)

    # Get required radius
    #  DECam      needs 1.1 deg
    #  Mosaic3    needs 0.43 deg
    #  Bok90prime needs 0.75 deg
    nc4d = np.sum(np.char.array(hplist1['instrument']).find('c4d') > -1)
    nksb = np.sum(np.char.array(hplist1['instrument']).find('ksb') > -1)
    minradius = 0.43
    if nksb>0:
        minradius = np.maximum(minradius, 0.75)
    if nc4d>0: 
        minradius = np.maximum(minradius, 1.1)
    # Add extra area for the size of the healpix pixel
    #   nside=128 is roughly 27' across
    #   nside=64 is roughly 55' across
    #   nside=32 is roughly 110' across
    radius = minradius + 0.5

    # Get all of the reference data that we need
    print('')
    ref = getrefdata(filters,cenra,cendec,radius)

    # Loop over the exposures
    for i in range(nind):
        print('')
        print('---- EXPOSURE '+str(i+1)+' OF '+str(nind)+' ----')
        print('')
        expdir = hplist1['expdir'][i]
        lo = expdir.find('/d1')
        expdir = dldir + expdir[lo+5:]
        calibrate(expdir,ref,redo=redo)

    print('')
    print('Total time = %.2f sec' % (time.time()-t00))


 
                
