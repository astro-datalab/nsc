#!/usr/bin/env python

from argparse import ArgumentParser
from astropy.io import fits
from astropy.table import *
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from glob import glob
import healpy as hp
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy import stats
import shutil
import socket
import subprocess
import sys
import time

from nsc import utils,query,modelmag
#from utils_tempest import *
import warnings
warnings.resetwarnings()
warnings.filterwarnings('ignore',category=UserWarning,append=True)

from slurm_funcs import *


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
    #sumfile = dldir+'dnidever/nsc/instcal/'+version+'/lists/nsc_calibrate_summary.fits'
    sumfile = dldir+'/nsc/instcal/'+version+'/lists/nsc_calibrate_summary.fits'
    if os.path.exists(sumfile) == False:
        if silent==False:
            logger.info(sumfile+' NOT FOUND')
        return
    logger.info('Loading '+sumfile)
    suminfo = Table.read(sumfile,1)
    for n in suminfo.colnames: suminfo[n].name = n.lower()  # lowercase names
    allinstfilt = np.char.array(suminfo['instrument'].astype(str))+'-'+np.char.array(suminfo['filter'].astype(str))
    gdfilt, = np.where((allinstfilt == instfilt) & (suminfo['nrefmatch'] > 10) & (suminfo['success'] == 1) & (suminfo['fwhm'] < 2.0))
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
        #catfile = str(refexp['expdir'][i]).strip()+'/'+str(refexp['base'][i]).strip()+'_meas.fits'
        catfile = str(refexp['expdir'][i]).strip()+'/'+str(refexp['base'][i]).strip()+'.fits'
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
            dt = [('sourceid',(str,100)),('ra',float),('dec',float),('ndet',int),('cmag',float),('cerr',float)]
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


def calibrate(expdir,inpref=None,refcatfile=None,logfilename=None,eqnfile=None,redo=False,selfcal=False,saveref=False,ncpu=1,logger=None,gsynthphot=False,psf=False):
    """
    Perform photometry and astrometric calibration of an NSC exposure using
    external catalogs.

    Parameters
    ----------
    expdir : str
       The absolute path to the tarred and zipped NSCexposure directory.
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
    gsynthphot: bool, optional
       If True, perform photometric calibration with
       Gaia synthetic photometry.  Default is False.
    psf : bool, optional
        If True, perform astrometric calibration with
        PSF coordinates instead of SE coordinates.
        Default is False.
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

    #print("psf,gsynthphot = ",psf,gsynthphot,type(psf),type(gsynthphot))
    if type(psf)==str: psf = eval(psf)
    if type(gsynthphot)==str: gsynthphot = eval(gsynthphot)

    # Calibrate catalogs for one exposure
    dldir,mssdir,localdir = utils.rootdirs()
    basedir = "/home/x25h971/nsc/instcal/v4/"
    tardir = basedir+(expdir.split("/")[-3])+"/"+(expdir.split("/")[-2])
    print("tardir = ",tardir)

    # Make sure the directory exists
    if os.path.exists(expdir) == False:
        raise ValueError(expdir+' NOT FOUND')

    t00 = time.time()

    if expdir[-1]=='/':
        expdir = expdir[0:-1]
    base = os.path.basename(expdir) #.split(".tar.gz")[0]
    if logger is None:
        logger = dln.basiclogger("calibrate_"+str(expdir))
    #outfile = expdir.split(".tar.gz")[0]+'/'+base+'_meta.fits'
    outfile = expdir+'/'+base+'_meta.fits'
    # get version number
    lo = expdir.find('nsc/instcal/')
    dum = expdir[lo+12:]
    version = dum[0:dum.find('/')]

    logger.info('Calibrate catalogs for exposure '+base+' in '+expdir)
    logger.info('Version = '+version)
    logger.info("calibrating with gsynth photometry: "+str(gsynthphot))
    logger.info("calibrating with psf coordinates: "+str(psf))

    # Check for output file
    if os.path.exists(outfile) == True and redo==False:
        logger.info(outfile+' already exists and /redo not set.')
        # Re-compress the exposure repo!
        os.chdir(tardir)
        print("tar -czf "+(expdir.split("/")[-1])+".tar.gz "+(expdir.split("/")[-1]))
        os.system("tar -czf "+(expdir.split("/")[-1])+".tar.gz "+(expdir.split("/")[-1]))
        #print("shutil.rmtree("+str(expdir.split("/")[-1])+")")
        #shutil.rmtree(expdir.split("/")[-1])
        #print("os.rmdir("+str(expdir.split("/")[-1])+")")
        #os.rmdir(expdir.split("/")[-1])
        os.chdir(basedir)
        return(True)

    # What instrument is this?
    instrument = 'c4d'# by default
    if expdir.find('/k4m/') > -1:
        instrument = 'k4m'
    if expdir.find('/ksb/') > -1:
        instrument = 'ksb'
    logger.info('This is a '+instrument+' exposure')

    # Model magnitude equation file
    if eqnfile is None:
        #eqnfile = dldir+'dnidever/nsc/instcal/'+version+'/config/modelmag_equations.txt'
        if gsynthphot==False: eqnfile = dldir+'nsc/instcal/'+version+'/config/modelmag_equations.txt'
        else: eqnfile = dldir+'nsc/instcal/'+version+'/config/modelmag_equations_gsynth.txt'
    logger.info('Using model magnitude equation file '+eqnfile)
    if os.path.exists(eqnfile) == False:
        raise ValueError(eqnfile+' NOT FOUND')
    eqnlines = dln.readlines(eqnfile)
    for l in eqnlines: logger.info(l)


    # Get the exposure files
    # Load the logfile and get absolute flux filename
    #loglines = dln.readlines(expdir.split(".tar.gz")[0]+'/'+base+'.log')
    if os.path.exists(expdir+'/'+base+'.log'):
        loglines = dln.readlines(expdir+'/'+base+'.log')
        #ind = dln.grep(loglines,'Step #2: Copying InstCal images from mass store archive',index=True)
        ind = dln.grep(loglines,"Copying InstCal images downloaded from Astro Data Archive",index=True)
        fline = loglines[ind[0]+1]
        #lo = fline.find('/archive') # as in /mss1/archive/pipe/.........
        #fluxfile = mssdir+str(fline[lo+1:])
        fluxfile = fline.split("   ")[-1].strip()
        wline = loglines[ind[0]+2]
        #lo = wline.find('/archive')
        #wtfile = mssdir+str(wline[lo+1:])
        wtfile = wline.split("   ")[-1].strip()
        mline = loglines[ind[0]+3]
        #lo = mline.find('/archive')
        #maskfile = mssdir+str(mline[lo+1:])
        maskfile = mline.split("   ")[-1].strip()
    else:
        explist = Table.read(basedir+"lists/nscdr3_instcal_list.fits.gz")
        expline = explist[explist['base']==expdir.split("/")[-1].split(".tar.gz")[0]]
        fluxfile = expline['fluxfile'][0].split("/")[-1]
        wtfile = expline['wtfile'][0].split("/")[-1]
        maskfile = expline['maskfile'][0].split("/")[-1]
    logger.info(fluxfile+wtfile+maskfile+" = flux,wt,maskfiles")


    # Step 1. Read in the catalogs 
    #----------------------------- 
    logger.info('')
    logger.info('Step 1. Read in the catalogs')
    logger.info('-----------------------------')
    catfiles = []
    catfiles1 = glob(expdir+'/'+base+'_[1-9].fits')
    if len(catfiles1) > 0:
        catfiles += catfiles1
    catfiles2 = glob(expdir+'/'+base+'_[0-9][0-9].fits') 
    if len(catfiles2) > 0:
        catfiles += catfiles2
    ncatfiles = len(catfiles) 
    if ncatfiles == 0:
        raise ValueError('No catalog files found')
    nchips = ncatfiles
    logger.info(str(ncatfiles)+' catalogs found')
     
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
     
    # Figure out the number of sources 
    ncat = 0 
    ncatarr = np.zeros(ncatfiles,int) 
    for i in range(ncatfiles):
        head = fits.getheader(catfiles[i],1) # was extension 2
        ncatarr[i] = head['NAXIS2']
        #print(head['NAXIS2']," = head['NAXIS2']")
    ncat = int(np.sum(ncatarr)) 
    logger.info(str(ncat)+' total sources')
    if ncat == 0: 
        logger.info('No sources')
        return 
    # Create structure, exten=1 has header now 
    ind, = np.where(ncatarr > 0)
    cat1 = fits.getdata(catfiles[ind[0]],1) # was 2
    cdt = cat1.dtype.descr
    #if 'ccdnum' not in cdt:
    cdt += [('ccdnum',int),('ebv',float),('ra',float),('dec',float),('raerr',float),
            ('decerr',float),('cmag',float),('cerr',float),('sourceid',(str,50)),
            ('filter',(str,50)),('mjd',float)]
    cat = np.zeros(ncat,dtype=np.dtype(cdt))
    cat = Table(cat)
    for c in cat.colnames:
        cat[c].name = c.lower()
    # Start the chips summary structure
    dt = [('expdir',(str,300)),('instrument',(str,10)),('filename',(str,300)),('measfile',(str,300)),
          ('ccdnum',int),('nsources',int),('nsources_psf',int),('nmeas',int),('cenra',float),('cendec',float),('ngaiamatch',int),
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

    # prepare headers from fluxfile
    #with fits.open(fluxfile) as hdulist:
    #    nextensions = len(hdulist)-1
    nextensions = len(subprocess.getoutput("ls "+expdir+"/*.fits").split("\n"))
    print("extensions? ",nextensions)
    xtns = []
    ccds = []
    #for nex in range(nextensions):
    #    hd_tmp = fits.getheader(fluxfile,nex+1)
    #    ccds.append(hd_tmp['CCDNUM'])
    #    xtns.append(nex+1)
    for nex in subprocess.getoutput("ls "+expdir+"/*.fits").split("\n"):
        if nex.split("_")[-1]!="meta.fits":
            hd_tmp = fits.getheader(nex)
            ccds.append(hd_tmp['CCDNUM'])
            xtns.append(nex.split("_")[-1].split(".")[0])
    ccds = np.array(ccds)
    xtns = np.array(xtns)
    logger.info("ccds = "+str(ccds)+", extensions = "+str(xtns))
    # Load the files 
    cnt = 0
    for i in range(ncatfiles):
        dum = os.path.splitext(os.path.basename(catfiles[i]))[0].split('_')
        ccdnum = int(dum[-1])
        logger.info("reading chip "+str(ccdnum))
        hd = fits.getheader(catfiles[i],1) # was 2
        cat1 = Table.read(catfiles[i],1) # was ext. 2
        if psf==True:
            cat1 = cat1[[not i for i in MaskedColumn(cat1['RAPSF']).mask]]
            #print("shrunk cat to PSF sources")
        for c in cat1.colnames:
            cat1[c].name = c.lower()
        ncat1 = hd['naxis2']   # handles empty catalogs; will be changed
        #logger.info(str(ncat1)+" = ncat1 for now")
        # Fix negative FWHM values 
        #  use A_WORLD and B_WORLD which are never negative 
        bdfwhm, = np.where(cat1['fwhm_world'] < 0.1)
        if len(bdfwhm) > 0: 
            cat['fwhm_world'][bdfwhm] = np.sqrt(cat1['a_world'][bdfwhm]**2+cat1['b_world'][bdfwhm]**2)*2.35 
         
        #ncat1 = n_elements(cat1) 
        chinfo['filename'][i] = catfiles[i] 
        chinfo['ccdnum'][i] = ccdnum 
        chinfo['nsources'][i] = ncat1 
        # Get the chip corners
        dum = fits.getdata(catfiles[i],1)
        #print("dum[0] = ",dum[0])
        #hd1 = fits.Header.fromstring('\n'.join(dum[0][0]),sep='\n')
        #hdflux = fits.getheader(catfiles[i],ccdnum-1)
        fbool = ccds==int(ccdnum)
        #logger.info("fbool = "+str(fbool))
        logger.info("get header for ccd "+str(ccds[fbool][0])+", extension "+str(xtns[fbool][0]))
        #hdflux = fits.getheader(fluxfile,xtns[fbool][0]) #flux alert!
        #nx = hdflux['NAXIS1']
        #ny = hdflux['NAXIS2']
        loglines = dln.readlines(catfiles[i].split(".fits")[0]+".logs")
        imline = np.array(dln.grep(loglines,"Picture size")[0].split(" "))
        imline = imline[imline!=""]
        nx = int(imline[np.where(imline=="size:")[0][0]+1])
        ny = int(imline[np.where(imline=="size:")[0][0]+2])
        #print("nx,ny = ",nx,ny)
        try:
            hd1 = fits.getheader(catfiles[i])
            wcs1 = WCS(hd1)
        except:
            logger.info('Problem with WCS in header '+catfiles[i])    
            continue
        #extast,hd1,ast,noparams# check the WCS 
        #if noparams <= 0: 
        #    logger.info('Problem with WCS in header '+catfiles[i])
        #    goto,BOMB1 
        #print("wcs1 = ",wcs1)
        vcoo = wcs1.pixel_to_world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
        #print("vcoo = ",vcoo)
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        chinfo['vra'][i] = vra 
        chinfo['vdec'][i] = vdec 
        if ncat1 > 0: 
            if psf==True:
                ncat1 = len(cat1)
                chinfo['nsources_psf'] = ncat1
                logger.info("assigned chinfo['nsources_psf']")
            logger.info(str(ncat1)+" = length of cat at chinfo['nsources_psf'] assignment")
            temp = cat[cnt:cnt+ncat1]
            for n in cat1.colnames:
                temp[n] = cat1[n]
            #STRUCT_ASSIGN,cat1,temp,/nozero 
            temp['ccdnum'] = ccdnum 
            if psf==False:
                temp['ra'] = cat1['alpha_j2000']  # add these here in case the astrometric correction fails later on 
                temp['dec'] = cat1['delta_j2000'] 
            else:
                temp['ra'] = cat1['rapsf']  # add these here in case the astrometric correction fails later on 
                temp['dec'] = cat1['decpsf'] 
            # Add coordinate uncertainties 
            #   sigma = 0.644*FWHM/SNR 
            #   SNR = 1.087/magerr 
            snr = 1.087/temp['magerr_auto']
            bderr, = np.where((temp['magerr_auto'] > 10) & (temp['magerr_iso'] < 10))
            if len(bderr) > 0 : 
                snr[bderr] = 1.087/temp[bderr]['magerr_iso']
            bderr, = np.where((temp['magerr_auto'] > 10) & (temp['magerr_iso'] > 10))
            if len(bderr) > 0 : 
                snr[bderr] = 1 
            coorderr = 0.664*(temp['fwhm_world']*3600)/snr 
            temp['raerr'] = coorderr 
            temp['decerr'] = coorderr 
            # Stuff into main structure 
            cat[cnt:cnt+ncat1] = temp 
            cnt += ncat1 
            # Wrapping around RA=0 
            if psf==False:
                cenra = np.mean(dln.minmax(cat1['alpha_j2000']))
                ra_range=cat1['alpha_j2000']
            else:
                cenra = np.mean(dln.minmax(cat1['rapsf']))
                ra_range=cat1['rapsf']
            #logger.info(str(cenra)+" = cenra")
            #if dln.valrange(cat1['alpha_j2000']) > 100: 
            if dln.valrange(ra_range) > 100:
                #ra = cat1['alpha_j2000']
                ra = ra_range
                #print("ra min max =",np.min(ra),np.max(ra))
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
            if psf==False: chinfo['cendec'][i] = np.mean(dln.minmax(cat1['delta_j2000'])) 
            else: chinfo['cendec'][i] = np.mean(dln.minmax(cat1['decpsf'])) 

    # Exposure level values 
    if psf==True:
        #print(np.mean(chinfo['nsources_psf']),np.mean(chinfo['cenra'])," mean nsources and cenra")
        gdchip, = np.where((chinfo['nsources_psf'] > 0) & (chinfo['cenra'] < 400))
    else:
        #print(np.mean(chinfo['nsources']),np.mean(chinfo['cenra'])," mean nsources and cenra")
        gdchip, = np.where((chinfo['nsources'] > 0) & (chinfo['cenra'] < 400))
    #print(gdchip," = gdchip")
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
    #glactc,cenra,cendec,2000.0,glon,glat,1,/deg 
    logger.info('GLON = %.5f' % glon)
    logger.info('GLAT = %.5f' % glat)
    # Number of good sources
    goodsources, = np.where((cat['imaflags_iso']==0) & (~((cat['flags']& 8)>0)) &
                            (~((cat['flags']&16)>0)) & (cat['mag_auto'] < 50))
    ngoodsources = len(goodsources)
    logger.info('GOOD SRCS = '+str(ngoodsources))

    # Measure median seeing FWHM
    gdcat, = np.where((cat['mag_auto'] < 50) & (cat['magerr_auto'] < 0.05) & (cat['class_star'] > 0.8))
    medfwhm = np.median(cat['fwhm_world'][gdcat]*3600.)
    logger.info('FWHM = %.2f arcsec' % medfwhm)


    # Load the meta-data from the original header
    #READLINE,expdir+'/'+base+'.head',head
    #head = fits.getheader(fluxfile,0) # flux alert!
    head = fits.getheader(catfiles[i])
    #head = hd1
    #print("head = ",head)
    filterlong = head['filter']
    if filterlong is None:
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
    cat['filter'] = filt
    cat['mjd'] = mjd
    cat['sourceid'] = instrument+'.'+str(expnum)+'.'+np.char.array(cat['ccdnum'].astype(str))+'.'+np.char.array(cat['number'].astype(str))

    # Start the exposure-level structure
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
    if inpref is None and refcatfile is None:
        # Search radius
        radius = 1.1 * np.sqrt( (0.5*rarange)**2 + (0.5*decrange)**2 )
        ref = query.getrefdata(instrument+'-'+filt,cenra,cendec,radius)
        if len(ref) == 0:
            logger.info('No Reference Data')
            return

    elif inpref is None and refcatfile is not None:
        ref = Table.read(refcatfile)
        logger.info("Reference catalog read in from "+refcatfile)

    # Using input reference catalog
    else:
        logger.info('Reference catalogs input')
        if rawrap == False:
            #gdref, = np.where((inpref['ra'] >= np.min(cat['alpha_j2000'])-0.01) &
            #                  (inpref['ra'] <= np.max(cat['alpha_j2000'])+0.01) &
            #                  (inpref['dec'] >= np.min(cat['delta_j2000'])-0.01) &
            #                  (inpref['dec'] <= np.max(cat['delta_j2000'])+0.01))
            if psf==False:
                catra = cat['alpha_j2000']
                catdec = cat['delta_j2000']
            else:
                catra = cat['rapsf']
                catdec = cat['decpsf']
            gdref, = np.where((inpref['ra'] >= np.min(catra)-0.01) &
                              (inpref['ra'] <= np.max(catra)+0.01) &
                              (inpref['dec'] >= np.min(catdec)-0.01) &
                              (inpref['dec'] <= np.max(catdec)+0.01))
        else:
            if psf==False:
                ra = cat['alpha_j2000']
                dec = cat['delta_j2000']
            else:
                ra = cat['rapsf']
                dec = cat['decpsf']
            bdra, = np.where(ra > 180)
            if len(bdra) > 0 :
                ra[bdra]-=360
            #gdref, = np.where((inpref['ra'] <= np.max(ra)-0.01) &
            #                  (inpref['ra'] >= np.min(ra+360)-0.01) &
            #                  (inpref['dec'] >= np.min(cat['delta_j2000'])-0.01) &
            #                  (inpref['dec'] <= np.max(cat['delta_j2000'])+0.01))
            gdref, = np.where((inpref['ra'] <= np.max(ra)-0.01) &
                              (inpref['ra'] >= np.min(ra+360)-0.01) &
                              (inpref['dec'] >= np.min(dec)-0.01) &
                              (inpref['dec'] <= np.max(dec)+0.01))
        ref = inpref[gdref]
        ngdref = len(gdref)
        logger.info(str(ngdref)+' reference stars in our region')



    # Step 3. Astrometric calibration
    #----------------------------------
    # At the chip level, linear fits in RA/DEC
    logger.info('')
    logger.info('Step 3. Astrometric calibration')
    logger.info('--------------------------------')
    # Get reference catalog with Gaia values
    #logger.info("ref cat columns = "+str(ref.colnames))
    gdgaia, = np.where(ref['source'] > 0)
    gaia = ref[gdgaia]
    # Match everything to Gaia at once, this is much faster!
    if psf==False: ind1,ind2,dist = coords.xmatch(gaia['ra'],gaia['dec'],cat['alpha_j2000'],cat['delta_j2000'],1.0)
    else: ind1,ind2,dist = coords.xmatch(gaia['ra'],gaia['dec'],cat['rapsf'],cat['decpsf'],1.0)
    ngmatch = len(ind1)
    if ngmatch == 0:
        logger.info('No gaia matches')
        return
    allgaiaind = np.zeros(ncat,int)-1
    allgaiaind[ind2] = ind1
    allgaiadist = np.zeros(ncat,float)+999999. 
    if psf==False: allgaiadist[ind2] = coords.sphdist(gaia['ra'][ind1],gaia['dec'][ind1],cat['alpha_j2000'][ind2],cat['delta_j2000'][ind2])*3600 
    else: allgaiadist[ind2] = coords.sphdist(gaia['ra'][ind1],gaia['dec'][ind1],cat['rapsf'][ind2],cat['decpsf'][ind2])*3600 
    # CCD loop 
    for i in range(nchips): 
        if (psf==True and chinfo['nsources_psf'][i]==0) or (psf==False and chinfo['nsources'][i]==0):
            continue
        # Get chip sources using CCDNUM 
        chind1,chind2 = dln.match(chinfo['ccdnum'][i],cat['ccdnum'])
        nchmatch = len(chind1)
        cat1 = cat[chind2] 
        logger.info("chip sources = "+str(len(cat1))+" "+str(nchmatch))
        # Gaia matches for this chip 
        gaiaind1 = allgaiaind[chind2] 
        gaiadist1 = allgaiadist[chind2] 
        gmatch, = np.where((gaiaind1 > -1) & (gaiadist1 <= 0.5))   # get sources with Gaia matches 
        if len(gmatch) == 0: 
            gmatch, = np.where((gaiaind1 > -1) & (gaiadist1 <= 1.0)) 
        if len(gmatch) < 5: 
            logger.info('Not enough Gaia matches')
            # Add threshold to astrometric errors 
            cat1['raerr'] = np.sqrt(cat1['raerr']**2 + 0.100**2) 
            cat1['decerr'] = np.sqrt(cat1['decerr']**2 + 0.100**2) 
            cat[chind2] = cat1 
            continue
        #gaia1b = gaia[ind1] 
        #cat1b = cat1[ind2] 
        gaia1b = gaia[gaiaind1[gmatch]] 
        cat1b = cat1[gmatch] 
        # Apply quality cuts 
        #  no bad CP flags 
        #  no SE truncated or incomplete data flags 
        #  must have good photometry
        #print("gaia1b = ",gaia1b)
        if 'pmra' in gaia1b.colnames and 'pmdec' in gaia1b.colnames:  # we have proper motion information 
            qcuts1, = np.where((cat1b['imaflags_iso'] == 0) & (~((cat1b['flags'] & 8) > 0)) & (~((cat1b['flags'] & 16) > 0)) & 
                               (cat1b['mag_auto'] < 50) & np.isfinite(gaia1b['pmra']) & np.isfinite(gaia1b['pmdec']))
        else:
            qcuts1, = np.where((cat1b['imaflags_iso'] == 0) & (~((cat1b['flags'] & 8) > 0)) & (~((cat1b['flags'] & 16) > 0)) & 
                               (cat1b['mag_auto'] < 50))
        if len(qcuts1) == 0: 
            logger.info('Not enough stars after quality cuts')
            # Add threshold to astrometric errors 
            cat1['raerr'] = np.sqrt(cat1['raerr']**2 + 0.100**2) 
            cat1['decerr'] = np.sqrt(cat1['decerr']**2 + 0.100**2) 
            cat[chind2] = cat1 
            continue
        gaia2 = gaia1b[qcuts1] 
        cat2 = cat1b[qcuts1] 
             
        # Precess the Gaia coordinates to the epoch of the observation 
        # The reference epoch for Gaia DR2 is J2015.5 (compared to the 
        # J2015.0 epoch for Gaia DR1). 

        # https://www.cosmos.esa.int/web/gaia/earlydr3#:~:text=The%20reference%20epoch%20for%20Gaia,full%20Gaia%20DR3)%20is%20J2016.
        # The reference epoch for Gaia DR3 (both Gaia EDR3 and the full Gaia DR3) is J2016.0.

        if 'PMRA' in gaia2.colnames and 'PMDEC' in gaia2.colnames:
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
        #ROTSPHCEN,gaia2.ra,gaia2.dec,chinfo[i].cenra,chinfo[i].cendec,gaialon,gaialat,/gnomic 
        gaialon,gaialat = coords.rotsphcen(gra_epoch,gdec_epoch,chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        if psf==False: lon1,lat1 = coords.rotsphcen(cat2['alpha_j2000'],cat2['delta_j2000'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        else: lon1,lat1 = coords.rotsphcen(cat2['rapsf'],cat2['decpsf'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        # ---- Fit RA as function of RA/DEC ---- 
        londiff = gaialon-lon1 
        err = None
        if 'RA_ERROR' in gaia2.colnames:
            err = np.sqrt(gaia2['ra_error']**2 + cat2['raerr']**2) 
        if 'RA_ERROR' not in gaia2.colnames and 'E_RA_ICRS' in gaia2.colnames:
            err = np.sqrt(gaia2['e_ra_icrs']**2 + cat2['raerr']**2) 
        if err is None:
            err = cat2['raerr']
        lonmed = np.median(londiff) 
        lonsig = np.maximum(dln.mad(londiff), 1e-5)  # 0.036" 
        gdlon, = np.where(np.abs(londiff-lonmed) < 3.0*lonsig)# remove outliers 
        if len(gdlon) > 5:   # use constant if not enough stars 
            npars = 4 
        else:
            npars = 1 
        initpars = np.zeros(npars,float) 
        initpars[0] = np.median([londiff]) 
        #print("initpars = ",initpars)
        racoef,racov = curve_fit(utils.func_poly2d_wrap,[lon1[gdlon],lat1[gdlon]],londiff[gdlon],sigma=err[gdlon],p0=initpars)
        yfitall = utils.func_poly2d(lon1,lat1,*racoef)
        #print("gdlons:")
        #print(londiff[gdlon])
        #print("yfitall = ",yfitall)
        #print(yfitall[gdlon])
        if len(gdlon)>5: rarms1 = dln.mad((londiff[gdlon]-yfitall[gdlon])*3600.) 
        else: rarms1 = dln.mad((londiff[gdlon]-np.repeat(yfitall,len(gdlon)))*3600.)
        rastderr = rarms1/np.sqrt(len(gdlon))
        # Use bright stars to get a better RMS estimate 
        gdstars, = np.where((cat2['fwhm_world']*3600 < 2*medfwhm) & (1.087/cat2['magerr_auto'] > 50))
        if len(gdstars) < 20: 
            gdstars, = np.where((cat2['fwhm_world']*3600 < 2*medfwhm) & (1.087/cat2['magerr_auto'] > 30))
        if len(gdstars) > 5: 
            diff = (londiff-yfitall)*3600. 
            rarms = dln.mad(diff[gdstars]) 
            rastderr = rarms/np.sqrt(len(gdstars))
        else:
            rarms = rarms1 
        # ---- Fit DEC as function of RA/DEC ----- 
        latdiff = gaialat-lat1 
        err = None
        if 'DEC_ERROR' in gaia2.colnames: 
            err = sqrt(gaia2['dec_error']**2 + cat2['decerr']**2) 
        if 'DEC_ERROR' not in gaia2.colnames and 'E_DEC_ICRS' in gaia2.colnames:
            err = np.sqrt(gaia2['e_de_icrs']**2 + cat2['decerr']**2) 
        if err is None: 
            err = cat2['decerr']
        latmed = np.median(latdiff) 
        latsig = np.maximum(dln.mad(latdiff), 1e-5)   # 0.036" 
        gdlat, = np.where(np.abs(latdiff-latmed) < 3.0*latsig)  # remove outliers 
        if len(gdlat) > 5:     # use constant if not enough stars 
            npars = 4 
        else: 
            npars = 1 
        initpars = np.zeros(npars,float) 
        initpars[0] = np.median(latdiff) 
        deccoef,deccov = curve_fit(utils.func_poly2d_wrap,[lon1[gdlat],lat1[gdlat]],latdiff[gdlat],sigma=err[gdlat],p0=initpars)
        yfitall = utils.func_poly2d(lon1,lat1,*deccoef)
        if len(gdlat)>5: decrms1 = dln.mad((latdiff[gdlat]-yfitall[gdlat])*3600.) 
        else: decrms1 = dln.mad((latdiff[gdlat]-np.repeat(yfitall,len(gdlat)))*3600.) 
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
        if psf==False: lon,lat = coords.rotsphcen(cat1['alpha_j2000'],cat1['delta_j2000'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        else: lon,lat = coords.rotsphcen(cat1['rapsf'],cat1['decpsf'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        lon2 = lon + utils.func_poly2d(lon,lat,*racoef)
        lat2 = lat + utils.func_poly2d(lon,lat,*deccoef) 
        ra2,dec2 = coords.rotsphcen(lon2,lat2,chinfo['cenra'][i],chinfo['cendec'][i],reverse=True,gnomic=True)
        cat1['ra'] = ra2 
        cat1['dec'] = dec2 
        # Add to astrometric errors 
        cat1['raerr'] = np.sqrt(cat1['raerr']**2 + rarms**2) 
        cat1['decerr'] = np.sqrt(cat1['decerr']**2 + decrms**2) 
        # Stuff back into the main structure 
        cat[chind2] = cat1 
        chinfo['ngaiamatch'][i] = len(gmatch)
        chinfo['ngoodgaiamatch'][i] = len(qcuts1) 
        chinfo['rarms'][i] = rarms 
        chinfo['rastderr'][i] = rastderr 
        chinfo['racoef'][i] = racoef 
        chinfo['decrms'][i] = decrms 
        chinfo['decstderr'][i] = decstderr 
        chinfo['deccoef'][i] = deccoef 

                 
    # Get reddening
    coo = SkyCoord(ra=cat['ra'],dec=cat['dec'],unit='deg')
    sfd = SFDQuery()
    ebv = sfd(coo)                
    #glactc,cat['ra'],cat['dec'],2000.0,glon,glat,1,/deg 
    #ebv = dust_getval(glon,glat,/noloop,/interp)
    cat['ebv'] = ebv 
                 
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
    ind1,ind2,dist = coords.xmatch(ref['ra'],ref['dec'],cat['ra'],cat['dec'],dcr)
    nmatch = len(ind1)
    logger.info(str(nmatch)+' matches to reference catalog')
    if nmatch == 0:
        logger.info('No matches to reference catalog')
        goto,ENDBOMB
    ref1 = ref[ind1]
    cat1 = cat[ind2]
    # Get the model magnitude
    mmags = modelmag.modelmag(ref1,instfilt,cendec,eqnfile) 
    if len(mmags) == 1 and mmags[0] < -1000: 
        print('No good model mags')
        return 
    # Get the good sources 
    gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) & (~((cat1['flags'] & 16) > 0)) &
                      (cat1['mag_auto'] < 50) &  (cat1['magerr_auto'] < 0.05) & (cat1['class_star'] > 0.8) &
                      (cat1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
    #  if the seeing is bad then class_star sometimes doens't work well 
    if medfwhm > 1.8 and len(gdcat) < 100:
        gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) & (~((cat1['flags'] & 16) > 0)) &
                          (cat1['mag_auto'] < 50) &  (cat1['magerr_auto'] < 0.05) & 
                          (cat1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
    if len(gdcat) > 0:
        ref2 = ref1[gdcat]
        mmags2 = mmags[gdcat,:]
        cat2 = cat1[gdcat]
        # Matched structure
        mag2 = cat2['mag_auto'] + 2.5*np.log10(exptime)   # correct for the exposure time 
        mstr = {'col':mmags2[:,2],'mag':mag2,'model':mmags2[:,0],'err':mmags2[:,1],'ccdnum':cat2['ccdnum']} 
        #print("min,max,median,mean model mags = ",min(mmags2[:,0]),max(mmags2[:,0]),np.median(mmags2[:,0]),np.mean(mmags2[:,0]))
        #print("median,mean model mag errs = ",np.median(mmags2[:,1]),np.mean(mmags2[:,1]))
        #print("min,max,median,mean nsc mags = ",min(mag2),max(mag2),np.median(mag2),np.mean(mag2))
        # Measure the zero-point
        expinfo,chinfo = fitzpterm(mstr,expinfo,chinfo)
        expinfo['zptype'] = 1 
    else:
        logger.info('No good reference sources')
    
    # Use self-calibration 
    if expinfo['nrefmatch'] <= 5 and selfcal:
        expinfo,chinfo = selfcalzpterm(expdir,cat,expinfo,chinfo)
        expinfo['zptype'] = 3 

    # Apply the zero-point to the full catalogs 
    # USE CHIP-LEVEL ZERO-POINTS WHEN POSSIBLE!!! 
    # Create an output catalog for each chip 
    ## Old way to get source
    ##if psf==True: nsrc = np.cumsum(chinfo['nsources_psf'])
    ##else: nsrc = np.cumsum(chinfo['nsources'])
    ##lo = np.append(0,nsrc[0:nchips])
    ##hi = nsrc-1 
    for i in range(nchips): 
        ##ncat1 = hi[i]-lo[i]+1 
        # Get chip sources using CCDNUM 
        chind1,chind2 = dln.match(chinfo['ccdnum'][i],cat['ccdnum'])
        nchmatch = len(chind1)
        if nchmatch > 0:
        ##if ncat1 > 0: 
            ##cat1 = cat[lo[i]:hi[i]] 
            cat1 = cat[chind2] 
            logger.info("chip sources = "+str(len(cat1))+" "+str(nchmatch))
            ##logger.info("ncat1, length cat1 = "+str(ncat1)+" "+str(len(cat1))+" during phot.cali")
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
            chinfo['zptype'][i] = 2 
                     
            gdcatmag, = np.where(cat1['mag_auto'] < 50)
            cat1['cmag'][gdcatmag] = cat1['mag_auto'][gdcatmag] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            cat1['cerr'][gdcatmag] = np.sqrt(cat1['magerr_auto'][gdcatmag]**2 + chinfo['zptermerr'][i]**2) # add calibration error in quadrature 
            #for j=0,n_elements(cat1.mag_aper)-1 do begin 
            #  gdcatmag = where(cat1.mag_aper[j] lt 50,ngd) 
            #  cat1[gdcatmag].cmag = cat1[gdcatmag].mag_auto + 2.5*alog10(exptime) + zpterm 
            #  cat1[gdcatmag].cerr = sqrt(cat1[gdcatmag].magerr_auto^2 + zptermerr^2)  ; add calibration error in quadrature 
            #endfor 
            # Print out the results 
            logger.info('  CCDNUM=%3d  NREFSOURCES=%5d  ZPTYPE=%2d  ZPTERM=%7.4f+/-%7.4f' %
                        (chinfo['ccdnum'][i],chinfo['nrefmatch'][i],chinfo['zptype'][i],chinfo['zpterm'][i],chinfo['zptermerr'][i]))
            ##cat[lo[i]:hi[i]] = cat1  # stuff back in 
            cat[chind2] = cat1
             
    # Print out the results 
    logger.info('NPHOTREFMATCH=%d' % expinfo['nrefmatch'])
    logger.info('EXPOSURE ZPTERM=%.4f +/- %.4f  SIG=%.4f mag' % (expinfo['zpterm'],expinfo['zptermerr'],expinfo['zptermsig']))
    logger.info('ZPSPATIALVAR:  RMS=%.4f RANGE=%.4f NCCD=%d' % 
                (expinfo['zpspatialvar_rms'],expinfo['zpspatialvar_range'],expinfo['zpspatialvar_nccd']))

             
    # Measure the depth 
    #   need good photometry 
    gdmag, = np.where(cat['cmag'] < 50) 
    ngdmag = len(gdmag)
    if ngdmag > 0: 
        # Get 95% percentile depth 
        cmag = cat['cmag'][gdmag]
        si = np.argsort(cmag) 
        cmag = cmag[si] 
        depth95 = cmag[int(np.round(0.95*ngdmag))-1] 
        expinfo['depth95'] = depth95 
        chinfo['depth95'] = depth95 
        logger.info('95% percentile depth'+' = %.2f mag' % depth95)
        # Get 10 sigma depth 
        #  S/N = 1.087/err 
        #  so S/N=5 is for err=1.087/5=0.2174 
        #  S/N=10 is for err=1.087/10=0.1087 
        depth10sig = 99.99 
        depind, = np.where((cat['cmag'] < 50) & (cat['cmag'] > depth95-3.0) & (cat['cerr'] >= 0.0987) & (cat['cerr'] <= 0.1187))
        if len(depind) < 5 : 
            depind, = np.where((cat['cmag'] < 50) & (cat['cmag'] > depth95-3.0) & (cat['cerr'] >= 0.0787) & (cat['cerr'] <= 0.1387))
        if len(depind) > 5: 
            depth10sig = np.median(cat['cmag'][depind]) 
        else: 
            depind, = np.where(cat['cmag'] < 50) 
            if len(depind) > 0: 
                depth10sig = np.max(cat['cmag'][depind]) 
        logger.info('10sigma depth = %.2f mag' % depth10sig)
        expinfo['depth10sig'] = depth10sig 
        chinfo['depth10sig'] = depth10sig 
             
    # Step 5. Write out the final catalogs and metadata 
    #-------------------------------------------------- 
    if redo and selfcal and expinfo['ztype']==2:
        # Create backup of original versions 
        logger.info('Copying meas and meta files to v1 versions')
        metafile = expdir+'/'+base+'_meta.fits' 
        if os.path.exists(metafile): 
            dln.file_copy(metafile,expdir+'/'+base+'_meta.v1.fits',overwrite=True)
             
    # Create an output catalog for each chip 
    ## Old way of getting chip sources
    ##if psf==True: nsrc = np.cumsum(chinfo['nsources_psf'])
    ##else: nsrc = np.cumsum(chinfo['nsources'])
    ##lo = np.append(0,nsrc[0:nchips]) 
    ##hi = nsrc-1 
    for i in range(nchips): 
        ##ncat1 = hi[i]-lo[i]+1 
        # Get chip sources using CCDNUM 
        chind1,chind2 = dln.match(chinfo['ccdnum'][i],cat['ccdnum'])
        nchmatch = len(chind1)
        if nmatch==0:
        ##if ncat1 == 0: 
            logger.info('No sources for CCDNUM='+str(chinfo['ccdnum'][i]))
            continue
        ##cat1 = cat[lo[i]:hi[i]] 
        ##logger.info("ncat1, len(cat1) = "+str(ncat1)+" "+str(len(cat1))+" before QA cuts")
        # Get chip sources using CCDNUM 
        cat1 = cat[chind2] 
        ncat1 = len(cat1)
        #logger.info("chip sources = "+str(len(cat1))+" "+str(nchmatch)+" before QA cuts")

                 
        # Apply QA cuts 
        #---------------- 
                 
        # Remove bad chip data 
        # Half of chip 31 for MJD>56660 
        #  c4d_131123_025436_ooi_r_v2 with MJD=56619 has problems too 
        #  if the background b/w the left and right half is large then BAd 
        lft31, = np.where((cat1['x_image'] < 1024) & (cat1['ccdnum'] == 31))
        rt31, = np.where((cat1['x_image'] >= 1024) & (cat1['ccdnum'] == 31))
        if len(lft31) > 10 and len(rt31) > 10: 
            lftback = np.median(cat1['background'][lft31]) 
            rtback = np.median(cat1['background'][rt31]) 
            mnback = 0.5*(lftback+rtback) 
            sigback = dln.mad(cat1['background']) 
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
            bdind, = np.where((cat1['x_image'] > 1000) & (cat1['ccdnum'] == 31))
            ngdind = len(cat1)-len(bdind)
            if len(bdind) > 0:  # some bad ones found 
                if ngdind == 0:  # all bad 
                    #logger.info('NO useful measurements in ',list[i].file 
                    cat1 = None
                    ncat1 = 0 
                    continue
                else: 
                    #logger.info('  Removing '+strtrim(nbdind,2)+' bad chip 31 measurements, '+strtrim(ngdind,2)+' left.' 
                    cat1 = np.delete(cat1,bdind)
                    ncat1 = len(cat1)         
            else:
                badchip31 = False  # chip 31 
                         
        # Make a cut on quality mask flag (IMAFLAGS_ISO) 
        bdcat, = np.where(cat1['imaflags_iso'] > 0)
        if len(bdcat) > 0: 
            #logger.info('  Removing '+str(nbdcat)+' sources with bad CP flags.' )
            if len(bdcat) == ncat1: 
                continue
            cat1 = np.delete(cat1,bdcat)
            ncat1 = len(cat1) 
                         
        # Make cuts on SE FLAGS 
        #   this removes problematic truncatd sources near chip edges 
        bdseflags, = np.where( ((cat1['flags'] & 8) > 0) |   # object truncated 
                               ((cat1['flags'] & 16) > 0))   # aperture truncate 
        if len(bdseflags) > 0: 
            #logger.info('  Removing '+str(nbdseflags)+' truncated sources' 
            if len(bdseflags) == ncat1: 
                continue
            cat1 = np.delete(cat1,bdseflags)
            ncat1 = len(cat1) 
                         
        # Removing low-S/N sources 
        #  snr = 1.087/err 
        snrcut = 5.0 
        bdsnr, = np.where(1.087/cat1['magerr_auto'] < snrcut) 
        if len(bdsnr) > 0: 
            #logger.info('  Removing '+str(nbdsnr)+' sources with S/N<',strtrim(snrcut,2) 
            if len(bdsnr) == ncat1: 
                continue
            cat1 = np.delete(cat1,bdsnr)
            ncat1 = len(cat1) 
        #logger.info("ncat1, len(cat1) = "+str(ncat1)+" "+str(len(cat1))+" after QA cuts")

        # Convert to final format
        if ncat1 > 0:
            #print(ncat1,len(cat1)," post-tabling")
            cat1 = Table(cat1)
            logger.info("^^^ for CCD "+str(cat1['ccdnum'][0]))
            #print(ncat1,len(cat1)," post-tabling")
            mdt = [('measid',(str,100)),('objectid',(str,100)),('exposure',(str,50)),
                   ('ccdnum',int),('filter',(str,50)),('mjd',float),('x',float),('y',float),
                   ('ra',float),('raerr',float),('dec',float),('decerr',float),('mag_auto',float),
                   ('magerr_auto',float),('mag_aper1',float),('magerr_aper1',float),('mag_aper2',float),
                   ('magerr_aper2',float),('mag_aper4',float),('magerr_aper4',float),('mag_aper8',float),
                   ('magerr_aper8',float),('kron_radius',float),('asemi',float),('asemierr',float),
                   ('bsemi',float),('bsemierr',float),('theta',float),('thetaerr',float),('fwhm',float),
                   ('flags',int),('class_star',float)] #,
                   #('background',float),('threshold',float),('isoarea_image',float),('isoarea_world',float),
                   #('x2_world',float),('y2_world',float),('xy_world',float),
                   #('errx2_world',float),('erry2_world',float),('errxy_world',float),
                   #('imaflags_iso',float),('nimaflags_iso',float),('ndet_iter',float),('ellipticity',float),
                   #('repeat',float),('xpsf',float),('ypsf',float),('magpsf',float),('errpsf',float),
                   #('sky',float),('iter',float),('chi',float),('sharp',float),
                   #('alpha_j2000),(),(),()]
            meas = np.zeros(ncat1,dtype=np.dtype(mdt))
            meas = Table(meas)
            #logger.info("ncat1,length cat1,length meas = "+str(ncat1)+" "+str(len(cat1))+" "+str(len(meas)))
            #print("cat1 colnames: ",cat1.colnames)
            for n in meas.colnames:
                if n in cat1.colnames:
                    meas[n] = cat1[n]
            meas['measid'] = cat1['sourceid'].astype(str)
            meas['exposure'] = base 
            meas['x'] = cat1['x_image']
            meas['y'] = cat1['y_image']
            meas['mag_auto'] = cat1['cmag']
            meas['magerr_auto'] = cat1['cerr']
            meas['mag_aper1'] = cat1['mag_aper'][:,0] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            meas['magerr_aper1'] = cat1['magerr_aper'][:,0] 
            meas['mag_aper2'] = cat1['mag_aper'][:,1] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            meas['magerr_aper2'] = cat1['magerr_aper'][:,1] 
            meas['mag_aper4'] = cat1['mag_aper'][:,2] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            meas['magerr_aper4'] = cat1['magerr_aper'][:,2] 
            meas['mag_aper8'] = cat1['mag_aper'][:,4] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
            meas['magerr_aper8'] = cat1['magerr_aper'][:,4] 
            meas['asemi'] = cat1['a_world'] * 3600.         # convert to arcsec 
            meas['asemierr'] = cat1['erra_world'] * 3600.   # convert to arcsec 
            meas['bsemi'] = cat1['b_world'] * 3600.         # convert to arcsec 
            meas['bsemierr'] = cat1['errb_world'] * 3600.   # convert to arcsec 
            meas['theta'] = 90-cat1['theta_world']          # make CCW E of N 
            meas['thetaerr'] = cat1['errtheta_world'] 
            meas['fwhm'] = cat1['fwhm_world'] * 3600.       # convert to arcsec 
            meas['class_star'] = cat1['class_star']
                             
            # Write to file 
            #outfile = expdir+'/'+base+'_'+str(chinfo['ccdnum'][i])+'_meas.fits' 
            outfile = expdir+'/'+base+'_'+str(chinfo['ccdnum'][i])+'.fits' 
            chinfo['nmeas'][i] = ncat1   # updating Chinfo 
            chinfo['measfile'][i] = outfile 
                            
            if redo and selfcal and expinfo['zptyper']==2: # Create backup of original versions 
                if os.path.exists(outfile) == 1 : 
                    #dln.file_copy(outfile,expdir+'/'+base+'_'+str(chinfo[i].ccdnum,2)+'_meas.v1.fits',overwrite=True)
                    dln.file_copy(outfile,expdir+'/'+base+'_'+str(chinfo[i]['ccdnum'])+'.v1.fits',overwrite=True)

            hdu = fits.HDUList()
            with fits.open(outfile) as hdul:
                hdu.append(hdul[0])
                hdu.append(hdul[1])
                hdu.append(fits.table_to_hdu(meas))
                hdu.append(fits.table_to_hdu(chinfo[i:i+1]))   # add chip stucture for this chip
                hdu.writeto(outfile,overwrite=True)
            ##hdu.close()
            #fits.update(outfile,data=np.array(meas),ext=2)
            #fits.update(outfile,data=np.array(chinfo[i:i+1]),ext=3)

    # Total number of measurements
    expinfo['nmeas'] = np.sum(chinfo['nmeas'])

    # Meta-data file 
    metafile = expdir+'/'+base+'_meta.fits' 
    logger.info('Writing metadata to '+metafile)
    hdus = fits.HDUList()
    hdus.append(fits.table_to_hdu(expinfo))
    hdus.append(fits.table_to_hdu(chinfo))  # add chip structure to second extension 
    hdus.writeto(metafile,overwrite=True)

    if os.path.exists(metafile):
        logger.info("exposure successfully calibrated!")
        print(logfilename,expdir.split(".tar.gz")[0]+"/"+logfilename.split("/")[-1])
        shutil.copy(logfilename,expdir.split(".tar.gz")[0]+"/"+logfilename.split("/")[-1])
        calbrated = True
    else: calbrated = False
    # Re-compress the exposure repo!
    os.chdir(tardir)
    print("tar -czf "+(expdir.split("/")[-1])+".tar.gz "+expdir.split("/")[-1])
    os.system("tar -czf "+(expdir.split("/")[-1])+".tar.gz "+(expdir.split("/")[-1]))
    #print("shutil.rmtree("+str(expdir.split("/")[-1])+")")
    #shutil.rmtree(expdir.split("/")[-1])
    #print("os.rmdir("+str(expdir.split("/")[-1])+")")
    #os.rmdir(expdir.split("/")[-1])
    os.chdir(basedir)

    dt = time.time()-t00 
    logger.info('dt = %.2f sec.' % dt)
    #return(ref1,ref2,mmags,mmags2)
    
    return(calbrated)


def complete_job(jobname,jstruct,torun,jobinds,jobtype):
    cputime,maxrss,jobid = sacct_cmd(jobname,['cputimeraw','maxrss','jobid'],c=True,m=True)
    jstruct[jobtype+'_cputime'][torun][jobinds] = cputime
    jstruct[jobtype+'_maxrss'][torun][jobinds] = maxrss
    jstruct[jobtype+'_jobid'][torun][jobinds] = jobid
    
def calibrate_healpix(pix,version,nside=32,maxexp=10,logtime="000",gsynth=True,psf=True,redo=False,logger=None):
    """
    This program is a wrapper around calibrate()
    for all exposures in the same region of the sky.
    It retrieves the reference catalogs ONCE which saves time.

    Parameters
    ----------
    pix       The HEALPix pixel number.
    version   The version name, e.g. 'v3'.
    =nside    The HEALPix nside to use.  Default is 32
    =maxexp   The maximum number of exposures to process at once
    /redo     Rerun exposures that were previously processed.

    Returns
    -------

    Example
    -------

    calibrate_healpix(1045,'v3')

    By D. Nidever 2017
    Translated to Python by D. Nidever, May 2022
    """

    # Set up directories
    basedir = "/home/x25h971/nsc/instcal/"+str(version)+"/"
    inddir = basedir+"healpix_indicators/"+str(nside)+"_"+str(int(pix))+"/"
    indfile = inddir+str(pix)+".txt"
    outfiledir = basedir+"outfiles/"
    refcatdir = basedir+"refcats/"
    tmpdir = basedir+'tmp/'
    logdir = basedir+"lists/logs/"
    if os.path.exists(logdir)==False:
        os.mkdir(logdir)
    if os.path.exists(tmpdir)==False:
        os.mkdir(tmpdir)
    t00 = time.time()

    # Do we have a logfile?
    if logger is None:
        logger = dln.basiclogger("calibrate_"+str(pix))

    # Load the list of HEALPix-labeled exposures
    listfile = basedir+'/lists/nsc_calibrate_healpix_list_blackmorerepos1.fits.gz'
    if os.path.exists(listfile)==False:
        logger.info(listfile,' NOT FOUND')
        return
    hplist = Table.read(listfile)

    # Get the exposures for this healpix
    logger.info('Calibrating InstCal SExtractor catalogs for Healpix pixel = '+str(pix))
    #ind1,ind2 = dln.match(hplist['pix'],[pix])
    #imatches,ind1,ind2 = np.intersect1d(hplist['pix'],np.array([pix]),return_indices=True)
    ind1 = np.isin(hplist['pix'],np.array([pix]))
    nind = len(hplist[ind1])
    hplist1 = hplist[ind1]
    if nind==0:
        logger.info('No exposures for pix '+str(pix))
        return
    logger.info('For pix '+str(pix)+', NEXPOSURES = '+str(nind))
    hplist1['expdir'] = [i.strip().split(".fits")[0]+".tar.gz" for i in hplist1['expdir']]
    hplist1['exp_done'] = np.repeat(False,nind) # for previous attempts
    hplist1['exp_calibrated'] = np.repeat(False,nind) # for this attempt
    hplist1['exp_ready'] = np.repeat(False,nind) # has the exposure been measured?  do we know where it is and how to get it?
    hplist1['filter'] = np.char.array(hplist1['filter']).strip()
    hplist1['transfer_group'] = np.repeat(-9,nind)
    hplist1['transfer_in_jobname'] = Column(dtype="U100",length=nind)
    hplist1['transfer_in_cmd'] = Column(dtype="U5000",length=nind)
    hplist1['transfer_in_jobid'] = Column(dtype="U10",length=nind)
    hplist1['transfer_in_jobstatus'] = Column(dtype="U10",length=nind)
    hplist1['transfer_out_jobname'] = Column(dtype="U100",length=nind)
    hplist1['transfer_in_cputime'] = Column(dtype="U10",length=nind)
    hplist1['transfer_in_maxrss'] = Column(dtype="U100",length=nind)
    hplist1['transfer_out_cmd'] = Column(dtype="U5000",length=nind)
    hplist1['transfer_out_jobid'] = Column(dtype="U10",length=nind)
    hplist1['transfer_out_jobstatus'] = Column(dtype="U10",length=nind)
    hplist1['calibrate_jobname'] = Column(dtype="U100",length=nind)
    hplist1['calibrate_cmd'] = Column(dtype="U1000",length=nind)
    hplist1['calibrate_jobid'] = Column(dtype="U10",length=nind)
    hplist1['calibrate_jobstatus'] = Column(dtype="U10",length=nind)
    hplist1['calibrate_cputime'] = Column(dtype="U10",length=nind)
    hplist1['calibrate_maxrss'] = Column(dtype="U100",length=nind)
    hplist1['torun'] = np.repeat(False,nind)
    hplist1['transfer_in_jobstatus'] = "-9"
    hplist1['transfer_out_jobstatus'] = "-9"
    hplist1['transfer_out_cputime'] = Column(dtype="U10",length=nind)
    hplist1['transfer_out_maxrss'] = Column(dtype="U100",length=nind)
    hplist1['calibrate_jobstatus'] = "-9"


    # Get list of already-calibrated exposures from indicator repo
    explist = subprocess.getoutput("ls "+inddir).split("\n")
    print("explist = ",explist)
    if "No such file or directory" in explist[0] or len(explist)<2:
        hplist1['exp_done'][:] = False
        logger.info("no indicator files for this pix "+str(pix))
    else:
        done_exposures = np.array([i.split(".")[0] for i in explist])
        hp_exposures = np.array([i.split("/")[-1].split(".tar.gz")[0]  for i in hplist1['expdir']])
        #print("done exposures = ",done_exposures)
        #print("hp_exposures = ",hp_exposures)
        exp_dones = np.isin(hp_exposures,done_exposures)
        hplist1['exp_done'][exp_dones] = True
        logger.info(str(len(hplist1[exp_dones]))+" exposures completed already for pix "+str(pix))
        if redo:
            logger.info("Re-calibrating those exposures!")
            # Remove all the indicator files
            for ed in explist:
                print(inddir+ed)
                os.remove(str(inddir+ed))
    # Get list of exposures that are ready to calibrate (have been measuresd)
    hplist1['exp_ready'] = [i!="-99" for i in hplist1['blackmore_repo'].filled("-99")]
    # Select exposures to calibrate!
    if redo:
        torun = np.where(hplist1['exp_ready']==True)[0]
    else:
        torun = np.where((hplist1['exp_done']==False) & (hplist1['exp_ready']==True))[0]
    ntorun = len(hplist1[torun])
    hplist1['torun'][torun] = True
    logger.info('NEXPOSURES to run = '+str(ntorun))
    if ntorun==0:
        sys.exit()
    # Get all of the reference data that we need, save to file to delete after
    # Central coordinates
    cenra,cendec = hp.pix2ang(nside,pix,lonlat=True)
    logger.info('RA  = %.6f' % cenra)
    logger.info('DEC = %.6f' % cendec)
    cencoo = SkyCoord(ra=cenra,dec=cendec,unit='deg')
    glon = cencoo.galactic.l.degree
    glat = cencoo.galactic.b.degree
    logger.info('L = %.6f' % glon)
    logger.info('B = %.6f' % glat)
    # List of instrument-filters
    #filters = np.char.array((hplist1['instrument']).strip())+'-'+str&np.char.array([f.strip()[0:2] for f in hplist1['filter']]).strip())
    filters = np.array([i['instrument'].strip()+"-"+i['filter'][0:2] for i in hplist1])
    filters = np.unique(filters)
    for i in range(0,len(filters)):
        if filters[i].split("-")[-1]!="VR" and filters[i].split("-")[-1]!="Y":
            filters[i] = str(str(filters[i].split("-")[0])+"-"+str(np.char.lower(filters[i].split("-")[-1])))
    logger.info("filters = "+", ".join(filters))
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
    radius = minradius + 1 # 0.5
    # Check to see if the file already exists
    refcatfile = refcatdir+"calibrate_healpix32_"+str(pix)+".fits.gz"
    if os.path.exists(refcatfile):
        logger.info("Reference catalog already queried, at "+refcatfile)
        ref = Table.read(refcatfile)
    else:
    # Do the query!
        logger.info('Querying reference catalog for healpix...')
        ref = query.getrefdata(filters,cenra,cendec,radius)
        ref.write(refcatfile,overwrite=True)
        logger.info("refcat saved to "+refcatfile)

    # Split exposures up into transfer batches
    ntransfers = (ntorun//maxexp)+1
    logger.info("We will be submitting "+str(ntransfers)+" transfers for "+str(maxexp)+" exposures per transfer")
    hplist1['transfer_group'][torun] = np.concatenate([np.repeat(i,maxexp) for i in range(ntransfers)])[:ntorun]

    # Save job structure to runfile
    runfile = basedir+"lists/runfiles/calibrate_healpix_"+str(pix)+"_run."+str(logtime)+".fits.gz"
    hplist1.write(runfile,overwrite=True)
    #print("hplist1 colnames = ",hplist1.colnames)
    logger.info("job structure saved to "+runfile)


    # Start submitting exposures for transfer & calibration!
    # ------------------------------------------------------
    if gsynth: gsth = " --gsynth " # Photometric calibration with Gaia Synthetic Photometry?
    else: gsth = ""
    if psf: ps = " --psf"          # Astrometric calibration with PSF coordinates?
    else: ps = ""
    if redo: rdo = " --redo"       # Redo the exposure?
    else: rdo = ""

    sleep_time = 10
    # - Loop through exposure batches (one transfer and one calibration job maintained simultaneously)
    for tns in range(ntransfers):
        logger.info(" ")
        logger.info("Writing & submitting transfer jobs for exposure batch "+str(tns))
        trans = np.where(hplist1['transfer_group'][torun]==tns)[0]
        #print("length of transfer group = ",len(hplist1[torun[trans]]))
        #print(hplist1[torun[trans]])

        # -- Create & submit script for expdir batch transfer from blackmore to tempest
        expdir_srcs = ",".join([i for i in hplist1['blackmore_repo'][torun[trans]]])
        expdir_dests = ",".join([i['expdir'].split("/")[-3]+"/"+i['expdir'].split("/")[-2]+"/"+i['blackmore_repo'].split("/")[-1] for i in hplist1[torun[trans]]])
        cmd_trans = "python "+str(basedir)+"globus-sdk-transfer.py blackmore tempest /phyx-nidever/kfas/nsc_meas/ "+basedir+" "+expdir_srcs+" "+expdir_dests
        print("cmd_trans_in = "+str(cmd_trans))
        jname_trans = "ctransin_"+str(pix)+"_"+str(tns)+"_"+str(logtime)
        hplist1['transfer_in_cmd'][torun[trans]] = cmd_trans
        hplist1['transfer_in_jobname'][torun[trans]] = jname_trans
        # Dependencies!  each transfer in job transiN must wait until transi(N-1) and cal(N-2) are done
        deps = []
        if tns>0: deps.append(hplist1[torun][hplist1[torun]['transfer_group']==(tns-1)]['transfer_in_jobid'][0].strip()) # add jid of last transfer in job
        if tns>1: deps.append(hplist1[torun][hplist1[torun]['transfer_group']==(tns-2)]['calibrate_jobid'][0].strip())   # add jid of last-last calibration job
        # Write job script for exposure transfer to tempest from blackmore
        jfile_trans = write_jscript(jname_trans,tns,"priority",[cmd_trans],outfiledir,"katiefasbender@montana.edu",cpupt=2,mempc="4G",rtime="05:00:00",parallel=False,dependent=deps)
        #os.system("sbatch "+jfile_trans)
        trans_jid = subprocess.getoutput("sbatch "+jfile_trans).split(" ")[-1]
        logger.info("submitted transfer_in job "+jname_trans+" "+str(trans_jid))
        time.sleep(sleep_time) # wait and get the id of the transfer job you just submitted
        #trans_jid = sacct_cmd(jname_trans,['jobid'],c=False,m=False)
        hplist1['transfer_in_jobid'][torun[trans]] = trans_jid

        # -- Create and submit script for calibration of expdir batch
        cmds_cal = []
        for i in range(len(hplist1[torun[trans]])):
            expdir = hplist1['expdir'][torun[trans[i]]]
            #logger.info("expdir = "+expdir)
            #exposure_runs.append(expdir)
            calcmd = "python "+basedir+"calibrate_exposure_tempest.py --version "+str(version)+" --pix "+str(pix)+" --nside "+str(nside)+" --expdir "+str(expdir.split(".tar.gz")[0])+" --refcat "+refcatfile+gsth+ps+rdo
            hplist1['calibrate_cmd'][torun[trans[i]]] = calcmd
            #print("cmd_cal = "+calcmd)
            cmds_cal.append(calcmd) #calibrate(expdir.split(".tar.gz")[0],ref,redo=redo,gsynthphot=True,psf=True)
        jname_cal = "calibrate_"+str(pix)+"_"+str(tns)+"_"+str(logtime)
        hplist1['calibrate_jobname'][torun[trans]] = jname_cal
        # Dependencies!  each calibration job calN must wait until transiN and cal(N-1) are done
        deps = [trans_jid]                                                                                      # add this transfer in job
        if tns>0: deps.append(hplist1[torun][hplist1[torun]['transfer_group']==(tns-1)]['calibrate_jobid'][0].strip()) # add jid of last calibration job
        # Write job script for exposure calibration
        jfile_cal = write_jscript(jname_cal,tns,"priority",cmds_cal,outfiledir,"katiefasbender@montana.edu",cpupt=2,mempc="4G",rtime="05:00:00",parallel=True,dependent=deps)
        #os.system("sbatch "+jfile_cal)
        cal_jid = subprocess.getoutput("sbatch "+jfile_cal).split(" ")[-1]
        hplist1['calibrate_jobid'][torun[trans]] = cal_jid
        logger.info("submitted calibration job "+jname_cal+" "+str(cal_jid))
        time.sleep(sleep_time) # wait and get the id of the transfer job you just submitted

        # -- Save job structure to runfile
        hplist1.write(runfile,overwrite=True)
        logger.info("Saving job structure to "+str(runfile))

    # - Check on exposure jobs every 5 minutes or so
    # Only once an exposure job is done may the exposure file, tarred and zipped with its base_meta.fits file within, be added to the transfer-out
    run_exposures = 0
    while run_exposures<ntorun:
        # loop through batches and check on them, zipping, writing, and transferring the files back when done
        for tns in range(ntransfers):
            logger.info("Checking current jobs of transfer group "+str(tns))
            trans = np.where(hplist1['transfer_group'][torun]==tns)[0]
            for jtype in ["calibrate","transfer_in","transfer_out"]:
                jname = hplist1[torun[trans]][jtype+'_jobname'][0]
                jstat = hplist1[torun[trans]][jtype+'_jobstatus'][0]
                logger.info(jname+" "+jstat)
                if jstat=="RUNNING" or jstat=="PENDING" or jstat=="REQUEUED" or jstat=="-9":
                    if jname.strip()!="":
                        new_jstat = sacct_cmd(jname,['state'],c=False,m=False)
                    else: new_jstat = "-9"
                    if new_jstat.strip()=="": new_jstat = "-9"
                    hplist1[jtype+"_jobstatus"][torun[trans]] = new_jstat
                    if new_jstat=="COMPLETED":
                        complete_job(jname,hplist1,torun,trans,jtype)
                        if jtype=="transfer_out": 
                            calibrated = np.where(hplist1[torun[trans]]['exp_calibrated']==True)[0]
                            for ed in hplist1[torun[trans]]['expdir'][calibrated]:
                                # Check for expdir
                                if os.path.exists(ed):
                                    os.remove(ed)
                                else: print(expdir," does not exist? tns group ",tns)
                        if jtype=="calibrate":
                            edirs = hplist1[torun[trans]]['expdir']
                            for ed,n in zip(edirs,range(len(edirs))):
                                ifile = inddir+ed.split("/")[-1].split(".tar.gz")[0]+".txt"
                                if os.path.exists(ifile): hplist1['exp_calibrated'][torun[trans[n]]] = True
                            calibrated = np.where(hplist1[torun[trans]]['exp_calibrated']==True)[0]
                            # -- Create and submit script for transfer of expdir batch from tempest back to blackmore
                            expdir_srcs = ",".join([i.split("v4/")[-1] for i in hplist1['expdir'][torun[trans[calibrated]]]])
                            expdir_dests = ",".join([i.split("v4/")[-1] for i in hplist1['expdir'][torun[trans[calibrated]]]])
                            cmd_trans = "python "+str(basedir)+"globus-sdk-transfer.py tempest blackmore "+basedir+" /phyx-nidever/kfas/nsc_meas/ "+expdir_srcs+" "+expdir_dests
                            print("cmd_trans_out = "+cmd_trans)
                            jname_trans = "ctransout_"+str(pix)+"_"+str(tns)+"_"+str(logtime)
                            hplist1['transfer_out_cmd'][torun[trans]] = cmd_trans
                            hplist1['transfer_out_jobname'][torun[trans]] = jname_trans
                            # Dependencies!  each transfer out job transoN must wait until transo(N-1) and cal(N) are done
                            #deps = [cal_jid]
                            #if tns>0: deps.append(hplist1[hplist1[torun]['transfer_group']==(tns-1)]['transfer_in_jobid'][0].strip()) # add jid of last transfer out job
                            # Write job script for exposure transfer to tempest from blackmore
                            jfile_trans = write_jscript(jname_trans,tns,"priority",[cmd_trans],outfiledir,"katiefasbender@montana.edu",cpupt=2,mempc="4G",rtime="05:00:00",parallel=False)#,dependent=deps)
                            trans_jid = subprocess.getoutput("sbatch "+jfile_trans).split(" ")[-1]
                            logger.info("submitted tans_out job "+jname_trans+" "+str(trans_jid))
                            #time.sleep(sleep_time) # wait and get the id of the transfer job you just submitted
                            #trans_jid = sacct_cmd(jname_trans,['jobid'],c=False,m=False)
                            hplist1['transfer_out_jobid'][torun[trans]] = trans_jid
                    if jtype=="transfer_out" and new_jstat!="RUNNING" and new_jstat!="PENDING" and new_jstat!="REQUEUED" and new_jstat!="-9": run_exposures = run_exposures+len(hplist1[torun[trans]])
            logger.info(str(run_exposures)+"/"+str(ntorun)+" exposures processed...")
        logger.info("wait a sec...")
        time.sleep(10)
        

    # -  Check for exposures copmleted from this healpix in indicator directory
    # If all exposures have been sucessfully calibrated, save a file indicating healpix completion
    #exps = subprocess.getoutput("ls "+inddir).split("\n")
    #if "No such file or directory" in exps[0]: hplist1['exp_calibrated'][:] = False
    #else:
    #    done_exps = np.array([i.split(".")[0] for i in exps])
    #    hp_exps = np.array([i.split("/")[-2].split(".fits")[0]  for i in hplist1['expdir']])
    #    exp_dones = np.isin(hp_exps,done_exps)
    #    hplist1['exp_calibrated'][exp_dones] = True
    #if len(hplist1['exp_calibrated'][exp_dones])==nind:
    if len(hplist1[hplist1['exp_calibrated']])==nind:
        logger.info("HEALPix Complete!")
        with (indfile,"w") as f:
            f.writelines("done! "+str(logtime)+" to "+str(time.time())+"\n")
            f.close()
    logger.info(str(len(hplist1[hplist1['exp_calibrated']]))+" HEALPix exposures calibrated")

    # - Save job structure to runfile
    hplist1.write(runfile,overwrite=True)
    logger.info("Saving job structure to "+str(runfile))
    logger.info('Total time = %.2f sec' % (time.time()-t00))
