#!/usr/bin/env python

import os
import time
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.table import Table
from astrop.wcs import WCS
from astropy.coordinates import SkyCoord
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from . import utils,query

def calibrate(expdir,inpref=None,eqnfile=None,redo=False,selfcal=False,saveref=False,ncpu=1,logger=None):
    """
    Perform photometry and astrometric calibration of an NSC exposure using
    external catalogs.

    Parameters
    ----------
    expdir : str
       The absolute path to the NSCexposure directory.
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

    calibreate(expdir,inpref,eqnfile)


    Written by D. Nidever in Nov 2017
    Translated to Python by D. Nidever, April 2022
    """

    # Calibrate catalogs for one exposure
    dldir,mssdir,localdir = utils.rootdirs()
     
    # Make sure the directory exists 
    if os.path.exists(expdir) == False:
        raise ValueError(expdir+' NOT FOUND')
     
    t00 = time.time() 
     
    base = os.path.basename(expdir) 
    if logger is None:
        logger = dln.basiclogger()
    outfile = expdir+'/'+base+'_meta.fits' 
    # get version number 
    lo = expdir.find('nsc/instcal/') 
    dum = expdir[lo+12:] 
    version = dum[0:dum.find('/')]
     
    logger.info('Calibrate catalogs for exposure '+base+' in '+expdir)
    logger.info('Version = '+version)
     
    # Check for output file 
    if os.path.exists(outfile) == False and redo==False:
        logger.info(outfile+' already exists and /redo not set.')
     
    # What instrument is this? 
    instrument = 'c4d'# by default
    if expdir.find('/k4m/') > -1:
        instrument = 'k4m'
    if expdir.find('/ksb/') > -1:
        instrument = 'ksb' 
    logger.info('This is a '+instrument+' exposure')
     
    # Model magnitude equation file
    if eqnfile is None:
        eqnfile = dldir+'users/dnidever/nsc/instcal/'+version+'/config/modelmag_equations.txt' 
    logger.info('Using model magnitude equation file '+eqnfile)
    if os.path.exists(eqnfile) == false: 
        raise ValueError(eqnfile+' NOT FOUND')
    eqnlines = dln.readlines(eqnfile)
    for l in eqnlines: logger.info(l)
     
     
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
        head = fits.getheader(catfiles[i],2)
        ncatarr[i] = head['NAXIS2']
    ncat = int(np.sum(ncatarr)) 
    logger.info(str(ncat),' total sources' 
    if ncat == 0: 
        logger.info('No sources')
        return 
    # Create structure, exten=1 has header now 
    ind, = np.where(ncatarr > 0)
    cat1 = fits.getdata(catfiles[ind[0]],2)
    cdt = cat1.dtype
    cdt += [('ccdnum',int),('ebv',float),('ra',float),('dec',float),('raerr',float),
            ('decerr',float),('cmag',float),('cerr',float),('sourceid',(np.str,50)),
            ('filter',(np.str,50)),('mjd',float)]
    cat = np.zeros(ncat,dtype=np.dtype(cdt))
    # Start the chips summary structure
    dt = [('expdir',(np.str,300)),('instrument',(np.str,10)),('filename',(np.str,300)),('measfile',(np.str,300)),
          ('ccdnum',int),('nsources',int),('nmeas',int),('cenra',float),('cendec',float),('ngaiamatch',int),
          ('ngoodgaiamatch',int),('rarms',float),('rastderr',float),('racoef',(float,4)),('decrms',float),
          ('decstderr',float),('decstderr',float),('deccoef',(float,4)),('vra',(float,4)),('vdec',(float,4)),
          ('zptype',int),('zpterm',float),('zptermerr',float),('nrefmatch',int),('depth95',float),('depth10sig',float)]
    chinfo = np.zeros(nchips,dtype=np.dtype(dt))
    for c in ['cenra','cendec','rarms','rastderr','decrms','decstderr','zpterm','zptermerr']:
        chinfo[c] = 999999.
    chinfo['depth95'] = 99.99
    chinfo['depth10sig'] = 99.99
    chinfo['expdir'] = expdir 
    chinfo['instrument'] = instrument 
    # Load the files 
    cnt = 0
    for i in range(ncatfiles):
        dum = os.path.splitext(os.path.basename(catfiles[i]))[0].split('_')
        ccdnum = int(dum[-1])
        hd = fits.getheader(catfiles[i],2)
        cat1 = Table.read(catilfes[i],2)
        ncat1 = hd['naxis2']   # handles empty catalogs 
         
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
        dum = fits.getdata(catfiles[1],1)
        #dum = MRDFITS(catfiles[i],1,/silent) 
        hd1 = dum['field_header_card']
        nx = hd1['NAXIS1']
        ny = hd1['NAXIS2']
        wcs1 = WCS(hd1)
        #extast,hd1,ast,noparams# check the WCS 
        #if noparams <= 0: 
        #    logger.info('Problem with WCS in header '+catfiles[i])
        #    goto,BOMB1 
        vcoo = wcs1.pixel_to_world(,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        chinfo['vra'][i] = vra 
        chinfo['vdec'][i] = vdec 
        if ncat1 > 0: 
            temp = cat[cnt:cnt+ncat1]
            for n in cat1.colnames:
                temp[n] = cat1[n]
            #STRUCT_ASSIGN,cat1,temp,/nozero 
            temp['ccdnum'] = ccdnum 
            temp['ra'] = cat1['alpha_j2000']  # add these here in case the astrometric correction fails later on 
            temp['dec'] = cat1['delta_j2000'] 
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
            cenra = np.mean(minmax(cat1['alpha_j2000'])) 
            # Wrapping around RA=0 
            if dln.valrange(cat1['alpha_j2000)'] > 100: 
                ra = cat1['alpha_j2000']
                bdra, = np.where(ra > 180) 
                if len(bdra) > 0: 
                    ra[bdra]-=360 
                bdra2, = np.where(ra < -180) 
                if len(bdra2) > 0: 
                    ra[bdra2] += 360 
                cenra = np.mean(dln.minmax(ra)) 
                if cenra < 0: 
                    cenra += 360 
            chinfo['cenra'][i] = cenra 
            chinfo['cendec'][i] = np.mean(dln.minmax(cat1['delta_j2000'])) 
        BOMB1: 
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
        rawrap = 1 
     else rawrap=0 
        logger.info('CENRA  = %.5f' % cenra)
        logger.info('CENDEC = %.5f' % cendec)
        glactc,cenra,cendec,2000.0,glon,glat,1,/deg 
        logger.info('GLON = %.5f' % glon)
        logger.info('GLAT = %.5f' % glat)
        # Number of good sources 
        goodsources, = np.where((cat['imaflags_iso']==0) & (~((cat['flags']& 8)>0)) &
                                 (~((cat['flags']&16)>0)) & (cat['mag_auto'] < 50))
        logger.info('GOOD SRCS = '+str(len(goodsources)))
         
        # Measure median seeing FWHM 
        gdcat, = np.where((cat['mag_auto'] < 50) & (cat['magerr_auto'] < 0.05) & (cat['class_star'] > 0.8))
        medfwhm = np.median(cat['fwhm_world'][gdcat]*3600.) 
        logger.info('FWHM = %.2f arcsec' % medfwhm)
         
        # Load the logfile and get absolute flux filename 
        loglines = dln.readlines(expdir+'/'+base+'.log')
        ind = dln.grep(loglines,'Step#2: Copying InstCal images from mass store archive')
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
        #READLINE,expdir+'/'+base+'.head',head
        head = fits.getheader(fluxfile,0)
        filterlong = head.get('filter')
        if filterlong is None:
            dum = mrdfits(catfiles[0],1,/silent) 
            hd1 = dum.field_header_card 
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
            hd1 = dum.field_header_card 
            exptime = hd1['exptime']
        dateobs = head['date-obs']
        airmass = head['airmass']
        mjd = date2jd(dateobs,/mjd) 
        logger.info('FILTER = '+filt)
        logger.info('EXPTIME = %.2f sec.' % exptime)
        logger.info('MJD = %.4f' % mjd)
         
        # Set some catalog values 
        cat['filter'] = filt
        cat['mjd'] = mjd 
        cat['sourceid'] = instrument+'.'+str(expnum)+'.'+str(cat['ccdnum'])+'.'+str(cat['number']) 
         
        # Start the exposure-level structure
        edt = [('file',(np.str,300)),('wtfile',(np.str,300)),('maskfile',(np.str,300)),('instrument',(np.str,10)),
               ('base',(np.str,100)),('exptnum',int),('ra',float),('dec',float),('dateobs',(np.str,50)),('mjd',float),
               ('filter',(np.str,50)),('exptime',float),('airmass',float),('wcscal',(np.str,50)),('nsources',int),
               ('ngoodsources',int),('nmeas',int),('fwhm',float),('nchips',int),('rarms',float),('decrms',float),
               ('ebv',float),('ngaiamatch',int),('ngoodgaiamatch',int),('zptype',int),('zpterm',float),('zptermerr',float),
               ('zptermsig',float),('zpspatialvar_rms',float),('zpspatialvar_range',float),('zpspatialvar_nccd',int),
               ('nrefmatch',int),('ngoodrefmatch',int),('depth95',float),('depth10sig',float)]
        expinfo = np.zeros(1,dtype=np.dtype(edt))
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
        wcscal = sxpar(head,'wcscal',count=nwcscal) 
        if nwcscal == 0 : 
            wcscal='NAN' 
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
            ref = getrefdata(instrument+'-'+filt,cenra,cendec,radius) 
            if count == 0: 
                logger.info('No Reference Data')
                return 
        # Using input reference catalog 
        else: 
            logger.info('Reference catalogs input')
            if rawrap == 0: 
                gdref, = np.where(inpref.ra >= min(cat.alpha_j2000)-0.01 and inpref.ra <= max(cat.alpha_j2000)+0.01 and
                                  inpref.dec >= min(cat.delta_j2000)-0.01 and inpref.dec <= max(cat.delta_j2000)+0.01,ngdref) 
            else: 
                ra = cat.alpha_j2000 
                bdra , = np.where(ra > 180,nbdra) 
                if nbdra > 0 : 
                    ra[bdra]-=360 
                gdref , = np.where((inpref.ra <= max(ra)+0.01 or inpref.ra >= min(ra+360)-0.01) and
                                   inpref.dec >= min(cat.delta_j2000)-0.01 and inpref.dec <= max(cat.delta_j2000)+0.01,ngdref) 
            ref = inpref[gdref] 
            logger.info(str(ngdref)+' reference stars in our region')
         
         
        # Step 3. Astrometric calibration 
        #---------------------------------- 
        # At the chip level, linear fits in RA/DEC 
        logger.info('')
        logger.info('Step 3. Astrometric calibration')
        logger.info('--------------------------------') 
        # Get reference catalog with Gaia values 
        gdgaia , = np.where(ref['source'] > 0,ngdgaia) 
        gaia = ref[gdgaia] 
        # Match everything to Gaia at once, this is much faster! 
        ind1,ind2,dist = coords.xmatch(gaia['ra'],gaia['dec'],cat['alpha_j2000'],cat['delta_j2000'],1.0)
        ngmatch = len(ind1)
        if ngmatch == 0: 
            logger.info('No gaia matches')
            return 
        allgaiaind = np.zeros(ncat,int)-1 
        allgaiaind[ind2] = ind1 
        allgaiadist = np.zeros(ncat,float)+999999. 
        allgaiadist[ind2] = coords.sphdist(gaia['ra'][ind1],gaia['dec'][ind1],cat['alpha_j2000'][ind2],cat['delta_j2000'][ind2])*3600 
        # CCD loop 
        for i in range(nchips): 
            if chinfo['nsources'][i]==0: 
                goto,BOMB 
            # Get chip sources using CCDNUM 
            chind1,chind2 = dln.match(chinfo['ccdnum'][i],cat['ccdnum'])
            cat1 = cat[chind2] 
            # Gaia matches for this chip 
            gaiaind1 = allgaiaind[chind2] 
            gaiadist1 = allgaiadist[chind2] 
            gmatch, = np.where(gaiaind1 > -1 and gaiadist1 <= 0.5)   # get sources with Gaia matches 
            if len(gmatch) == 0: 
                gmatch, = np.where(gaiaind1 > -1 and gaiadist1 <= 1.0) 
            if len(gmatch) < 5: 
                logger.info('Not enough Gaia matches')
                # Add threshold to astrometric errors 
                cat1['raerr'] = np.sqrt(cat1['raerr']**2 + 0.100**2) 
                cat1['decerr'] = np.sqrt(cat1['decerr']**2 + 0.100**2) 
                cat[chind2] = cat1 
                goto,BOMB 
            #gaia1b = gaia[ind1] 
            #cat1b = cat1[ind2] 
            gaia1b = gaia[gaiaind1[gmatch]] 
            cat1b = cat1[gmatch] 
            # Apply quality cuts 
            #  no bad CP flags 
            #  no SE truncated or incomplete data flags 
            #  must have good photometry
            if 'pmra' in gaia1b and 'pmdec' in gaia1b:  # we have proper motion information 
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
                goto,BOMB 
            gaia2 = gaia1b[qcuts1] 
            cat2 = cat1b[qcuts1] 
             
            # Precess the Gaia coordinates to the epoch of the observation 
            # The reference epoch for Gaia DR2 is J2015.5 (compared to the 
            # J2015.0 epoch for Gaia DR1). 
            if tag_exist(gaia2,'PMRA') and tag_exist(gaia2,'PMDEC'): 
                gaiamjd = 57206.0
                delt = (mjd-gaiamjd)/365.242170# convert to years 
                # convert from mas/yr->deg/yr and convert to angle in RA 
                gra_epoch = gaia2['ra'] + delt*gaia2['pmra']/3600.0/1000.0/np.cos(np.deg2rad(gaia2['dec'])) 
                gdec_epoch = gaia2['dec'] + delt*gaia2['pmdec']/3600.0/1000.0 
            else: 
                gra_epoch = gaia2['ra']
                gdec_epoch = gaia2['dec'] 
             
            # Rotate to coordinates relative to the center of the field 
            #ROTSPHCEN,gaia2.ra,gaia2.dec,chinfo[i].cenra,chinfo[i].cendec,gaialon,gaialat,/gnomic 
            gaialon,gaialat = coords.rotsphcen(gra_epoch,gdec_epoch,chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
            lon1,lat1 = coords.rotsphcen(cat2['alpha_j2000'],cat2['delta_j2000'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
            # ---- Fit RA as function of RA/DEC ---- 
            londiff = gaialon-lon1 
            err = None
            if tag_exist(gaia2,'RA_ERROR') : 
                err = np.sqrt(gaia2['ra_error']**2 + cat2['raerr']**2) 
            if tag_exist(gaia2,'RA_ERROR') == 0 and tag_exist(gaia2,'E_RA_ICRS') : 
                err = np.sqrt(gaia2['e_ra_icrs']**2 + cat2['raerr']**2) 
            if err is None:
                err = cat2['raerr']
            lonmed = np.median([londiff]) 
            lonsig = dln.mad([londiff]) > 1e-5# 0.036" 
            gdlon, = np.where(np.abs(londiff-lonmed) < 3.0*lonsig)# remove outliers 
            if len(gdlon) > 5:   # use constant if not enough stars 
                npars = 4 
            else: 
                npars=1 
            initpars = np.zeros(npars,float) 
            initpars[0] = np.median([londiff]) 
            parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars) 
            racoef = MPFIT2DFUN('func_poly2d',lon1[gdlon],lat1[gdlon],londiff[gdlon],err[gdlon],initpars,status=status,dof=dof,
                                bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet) 
            yfitall = FUNC_POLY2D(lon1,lat1,racoef) 
            rarms1 = dln.mad((londiff[gdlon]-yfit)*3600.) 
            rastderr = rarms1/sqrt(ngdlon) 
            # Use bright stars to get a better RMS estimate 
            gdstars, = np.where((cat2['fwhm_world']*3600 < 2*medfwhm) & (1.087/cat2['magerr_auto'] > 50))
            if len(gdstars) < 20: 
                gdstars, = np.where((cat2['fwhm_world']*3600 < 2*medfwhm) & (1.087/cat2['magerr_auto'] > 30))
            if len(gdstars) > 5: 
                diff = (londiff-yfitall)*3600. 
                rarms = dln.mad(diff[gdstars]) 
                rastderr = rarms/np.sqrt(ngdstars) 
             else rarms=rarms1 
                # ---- Fit DEC as function of RA/DEC ----- 
                latdiff = gaialat-lat1 
                undefine,err 
                if tag_exist(gaia2,'DEC_ERROR') : 
                    err = sqrt(gaia2['dec_error']**2 + cat2['decerr']**2) 
                if tag_exist(gaia2,'DEC_ERROR') == 0 and tag_exist(gaia2,'E_DEC_ICRS') : 
                    err = np.sqrt(gaia2['e_de_icrs']**2 + cat2['decerr']**2) 
                if len(err) == 0 : 
                    err = cat2['decerr']
                latmed = np.median(latdiff) 
                latsig = dln.mad(latdiff) > 1e-5# 0.036" 
                gdlat, = np.where(np.abs(latdiff-latmed) < 3.0*latsig)  # remove outliers 
                if len(gdlat) > 5:     # use constant if not enough stars 
                    npars = 4 
                else: 
                    npars = 1 
                initpars = np.zeros(npars,float) 
                initpars[0] = np.median(latdiff) 
                parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars) 
                deccoef = MPFIT2DFUN('func_poly2d',lon1[gdlat],lat1[gdlat],latdiff[gdlat],err[gdlat],initpars,status=status,dof=dof,
                                     bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet) 
                yfitall = FUNC_POLY2D(lon1,lat1,deccoef) 
                decrms1 = dln.mad((latdiff[gdlat]-yfit)*3600.) 
                decstderr = decrms1/sqrt(ngdlat) 
                # Use bright stars to get a better RMS estimate 
                if ngdstars > 5: 
                    diff = (latdiff-yfitall)*3600. 
                    decrms = dln.mad(diff[gdstars]) 
                    decstderr = decrms/sqrt(ngdstars) 
                 else decrms=decrms1 
                    logger.info('  CCDNUM=',string(chinfo['ccdnum'][i],format='(i3)'),'  NSOURCES=',string(nchmatch,format='(i5)'),'  ',
                                string(ngmatch,format='(i5)'),'/',string(nqcuts1,format='(i5)'),' GAIA matches  RMS(RA/DEC)=',
                                string(rarms,format='(f7.4)')+'/'+string(decrms,format='(f7.4)'),' STDERR(RA/DEC)=',
                                string(rastderr,format='(f7.4)')+'/'+string(decstderr,format='(f7.4)'),' arcsec' 
                    # Apply to all sources 
                    lon2,lat = coords.rotsphcen(cat1['alpha_j2000'],cat1['delta_j2000'],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
                    lon2 = lon + FUNC_POLY2D(lon,lat,racoef) 
                    lat2 = lat + FUNC_POLY2D(lon,lat,deccoef) 
                    ra2,dec2 = coords.rotsphcen(lon2,lat2,chinfo['cenra'][i],chinfo['cendec'][i],reverse=True,gnomic=True)
                    cat1['ra'] = ra2 
                    cat1['dec'] = dec2 
                    # Add to astrometric errors 
                    cat1['raerr'] = np.sqrt(cat1['raerr']**2 + rarms**2) 
                    cat1['decerr'] = np.sqrt(cat1['decerr']**2 + decrms**2) 
                    # Stuff back into the main structure 
                    cat[chind2] = cat1 
                    chinfo['ngaiamatch'][i] = ngmatch 
                    chinfo['ngoodgaiamatch'][i] = nqcuts1 
                    chinfo['rarms'][i] = rarms 
                    chinfo['rastderr'][i] = rastderr 
                    chinfo['racoef'][i] = racoef 
                    chinfo['decrms'][i] = decrms 
                    chinfo['decstderr'][i] = decstderr 
                    chinfo['deccoef'][i] = deccoef 
                    BOMB: 
                 
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
                # Get the model magnitudes 
                mmags = getmodelmag(ref1,instfilt,cendec,eqnfile) 
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
                if len(gdcat) == 0: 
                    logger.info('No good reference sources')
                    goto,ENDBOMB 
                ref2 = ref1[gdcat] 
                mmags2 = mmags[gdcat,:] 
                cat2 = cat1[gdcat] 
                # Matched structure 
                mag2 = cat2['mag_auto'] + 2.5*np.log10(exptime)   # correct for the exposure time 
                mstr = {col:float(mmags2[:,2]),mag:float(mag2),model:float(mmags2[:,0]),err:float(mmags2[:,1]),ccdnum:int(cat2.ccdnum)} 
                # Measure the zero-point
                expinfo = nsc_instcal_calibrate_fitzpterm(mstr,expinfo,chinfo)
                expinfo['zptype'] = 1 
                 
 
             
            # Use self-calibration 
            if expinfo['nrefmatch'] <= 5 and selfcal:
                NSC_INSTCAL_CALIBRATE_SELFCALZPTERM,expdir,cat,expinfo,chinfo 
                expinfo['zptype'] = 3 
            # Apply the zero-point to the full catalogs 
            # USE CHIP-LEVEL ZERO-POINTS WHEN POSSIBLE!!! 
            # Create an output catalog for each chip 
            nsrc = long64(np.sum(chinfo.nsources,/cum)) 
            lo = [0L,nsrc[0:nchips-2]] 
            hi = nsrc-1 
            for i in range(nchips): 
                ncat1 = hi[i]-lo[i]+1 
                if ncat1 > 0: 
                    cat1 = cat[lo[i]:hi[i]] 
                    #; Use chip-level zero-point 
                    #if chinfo[i].nrefmatch gt 5 then begin 
                    #  chinfo[i].zptype = 1 
                    #; Use exposure-level zero-point 
                    #endif else begin 
                    #  chinfo[i].zpterm = expinfo.zpterm 
                    #  chinfo[i].zptermerr = expinfo.zptermerr 
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
                    logger.info('  CCDNUM=',string(chinfo[i].ccdnum,format='(i3)'),'  NREFSOURCES=',string(chinfo[i].nrefmatch,format='(i5)'),'  ZPTYPE=',                  string(chinfo[i].zptype,format='(i2)'),'  ZPTERM=',string(chinfo[i].zpterm,format='(f7.4)'),'+/-',string(chinfo[i].zptermerr,format='(f7.4)') 
                    cat[lo[i]:hi[i]] = cat1# stuff back in 
             
            # Print out the results 
            logger.info('NPHOTREFMATCH=',str(expinfo['nrefmatch']))
            logger.info('EXPOSURE ZPTERM=',stringize(expinfo.zpterm,ndec=4),'+/-',stringize(expinfo.zptermerr,ndec=4),'  SIG=',stringize(expinfo.zptermsig,ndec=4),'mag' 
            logger.info('ZPSPATIALVAR:  RMS=',stringize(expinfo.zpspatialvar_rms,ndec=4),' ',         'RANGE=',stringize(expinfo.zpspatialvar_range,ndec=4),' NCCD=',str(expinfo.zpspatialvar_nccd,2) 
             
            # Measure the depth 
            #   need good photometry 
            gdmag, = np.where(cat['cmag'] < 50) 
            if len(gdmag) > 0: 
                # Get 95% percentile depth 
                cmag = cat['cmag'][gdmag]
                si = np.argsort(cmag) 
                cmag = cmag[si] 
                depth95 = cmag[int(np.round(0.95*ngdmag)-1] 
                expinfo['depth95'] = depth95 
                chinfo['depth95'] = depth95 
                logger.info('95% percentile depth = %.2f mag' % depth95)
                # Get 10 sigma depth 
                #  S/N = 1.087/err 
                #  so S/N=5 is for err=1.087/5=0.2174 
                #  S/N=10 is for err=1.087/10=0.1087 
                depth10sig = 99.99 
                depind, = np.where(cat.cmag < 50 and cat.cmag > depth95-3.0 and cat.cerr >= 0.0987 and cat.cerr <= 0.1187,ndepind) 
                if len(depind) < 5 : 
                    depind, = np.where(cat.cmag < 50 and cat.cmag > depth95-3.0 and cat.cerr >= 0.0787 and cat.cerr <= 0.1387,ndepind) 
                if len(depind) > 5: 
                    depth10sig = np.median([cat[depind].cmag]) 
                else: 
                    depind, = np.where(cat.cmag < 50,ndepind) 
                    if len(depind) > 0: 
                        depth10sig = np.max(cat['cmag'][depind]) 
                logger.info('10sigma depth = %.2f mag' % depth10sig)
                expinfo['depth10sig'] = depth10sig 
                chinfo['depth10sig'] = depth10sig 
             
            # Step 5. Write out the final catalogs and metadata 
            #-------------------------------------------------- 
            if keyword_set(redo) and keyword_set(selfcal) and expinfo['zptype'] == 2: 
                # Create backup of original versions 
                logger.info('Copying meas and meta files to v1 versions')
                metafile = expdir+'/'+base+'_meta.fits' 
                if os.path.exists(metafile): 
                    FILE_COPY,metafile,expdir+'/'+base+'_meta.v1.fits',/overwrite 
             
             
            # Create an output catalog for each chip 
            nsrc = long64(np.sum(chinfo.nsources,/cum)) 
            lo = [0L,nsrc[0:nchips-2]] 
            hi = nsrc-1 
            for i in range(nchips): 
                ncat1 = hi[i]-lo[i]+1 
                if ncat1 == 0: 
                    logger.info('No sources for CCDNUM='+str(chinfo['ccdnum'][i]))
                    goto,CHIPBOMB 
                cat1 = cat[lo[i]:hi[i]] 
                 
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
                        jump31 = 1 
                    else: 
                        jump31 = 0 
                    #logger.info('  Big jump in CCDNUM 31 background levels' 
                 else jump31=0 
                    if expinfo.mjd > 56600 or jump31 == 1: 
                        badchip31 = 1 
                        # Remove bad measurements 
                        # X: 1-1024 okay 
                        # X: 1025-2049 bad 
                        # use 1000 as the boundary since sometimes there's a sharp drop 
                        # at the boundary that causes problem sources with SExtractor 
                        bdind, = np.where((cat1['x_image'] > 1000) & (cat1['ccdnum'] == 31))
                        if len(bdind) > 0:# some bad ones found 
                            if ngdind == 0:# all bad 
                                #logger.info('NO useful measurements in ',list[i].file 
                                undefine,cat1 
                                ncat1 = 0 
                                goto,CHIPBOMB 
                            else: 
                                #logger.info('  Removing '+strtrim(nbdind,2)+' bad chip 31 measurements, '+strtrim(ngdind,2)+' left.' 
                                REMOVE,bdind,cat1 
                                ncat1 = len(cat1) 
                    # some bad ones to remove 
                     else badchip31=0# chip 31 
                         
                        # Make a cut on quality mask flag (IMAFLAGS_ISO) 
                        bdcat, = np.where(cat1['imaflags_iso'] > 0)
                        if len(bdcat) > 0: 
                            #logger.info('  Removing ',strtrim(nbdcat,2),' sources with bad CP flags.' 
                            if nbdcat == ncat1 : 
                                goto,CHIPBOMB 
                            REMOVE,bdcat,cat1 
                            ncat1 = len(cat1) 
                         
                        # Make cuts on SE FLAGS 
                        #   this removes problematic truncatd sources near chip edges 
                        bdseflags, = np.where( ((cat1['flags'] & 8) > 0) |   # object truncated 
                                               ((cat1['flags'] & 16) > 0))   # aperture truncate 
                        if len(bdseflags) > 0: 
                            #logger.info('  Removing ',strtrim(nbdseflags,2),' truncated sources' 
                            if nbdseflags == ncat1 : 
                                goto,CHIPBOMB 
                            REMOVE,bdseflags,cat1 
                            ncat1 = len(cat1) 
                         
                        # Removing low-S/N sources 
                        #  snr = 1.087/err 
                        snrcut = 5.0 
                        bdsnr, = np.where(1.087/cat1['magerr_auto'] < snrcut) 
                        if len(bdsnr) > 0: 
                            #logger.info('  Removing ',strtrim(nbdsnr,2),' sources with S/N<',strtrim(snrcut,2) 
                            if nbdsnr == ncat1 : 
                                goto,CHIPBOMB 
                            REMOVE,bdsnr,cat1 
                            ncat1 = len(cat1) 
                         
                        # Convert to final format 
                        if ncat1 > 0:
                            mdt = [('measid',(np.str,100)),('objectid',(np.str,100)),('exposure',(np.str,50)),
                                   ('ccdnum',int),('filter',(np.str,50)),('mjd',float),('x',float),('y',float),
                                   ('ra',float),('raerr',float),('dec',float),('decerr',float),('mag_auto',float),
                                   ('magerr_auto',float),('mag_aper1',float),('magerr_aper1',float),('mag_aper2',float),
                                   ('magerr_aper2',float),('mag_aper4',float),('magerr_aper4',float),('mag_aper8',float),
                                   ('magerr_aper8',float),('kron_radius',float),('asemi',float),('asemierr',float),
                                   ('bsemi',float),('bsemierr',float),('theta',float),('thetaerr',float),('fwhm',float),
                                   ('flags',int),('class_star',float)]
                            meas = np.zeros(ncat1,dtype=np.dtype(mdt))
                            meas = Table(meas)
                            for n in meas.colnames:
                                meas[n] = cat1[n]
                            meas['measid'] = str(cat1['sourceid']) 
                            meas['exposure'] = base 
                            meas['x'] = cat1['x_image']
                            meas['y'] = cat1['y_image']
                            meas['mag_auto'] = cat1['cmag']
                            meas['magerr_auto'] = cat1['cerr']
                            meas['mag_aper1'] = cat1['mag_aper'][0] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
                            meas['magerr_aper1'] = cat1['magerr_aper'][0] 
                            meas['mag_aper2'] = cat1['mag_aper'][1] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
                            meas['magerr_aper2'] = cat1['magerr_aper'][1] 
                            meas['mag_aper4'] = cat1['mag_aper'][2] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
                            meas['magerr_aper4'] = cat1['magerr_aper'][2] 
                            meas['mag_aper8'] = cat1['mag_aper'][4] + 2.5*np.log10(exptime) + chinfo['zpterm'][i]
                            meas['magerr_aper8'] = cat1['magerr_aper'][4] 
                            meas['asemi'] = cat1['a_world'] * 3600.         # convert to arcsec 
                            meas['asemierr'] = cat1['erra_world'] * 3600.   # convert to arcsec 
                            meas['bsemi'] = cat1['b_world'] * 3600.         # convert to arcsec 
                            meas['bsemierr'] = cat1['errb_world'] * 3600.   # convert to arcsec 
                            meas['theta'] = 90-cat1['theta_world']          # make CCW E of N 
                            meas['thetaerr'] = cat1['errtheta_world'] 
                            meas['fwhm'] = cat1['fwhm_world'] * 3600.       # convert to arcsec 
                            meas['class_star'] = cat1['class_star']
                             
                            # Write to file 
                            outfile = expdir+'/'+base+'_'+str(chinfo['ccdnum'][i])+'_meas.fits' 
                            chinfo['nmeas'][i] = ncat1   # updating Chinfo 
                            chinfo['measfile'][i] = outfile 
                             
                            if keyword_set(redo) and keyword_set(selfcal) and expinfo.zptype == 2 :# Create backup of original versions 
                                $ 
                            if os.path.exists(outfile) == 1 : 
                                FILE_COPY,outfile,expdir+'/'+base+'_'+str(chinfo[i].ccdnum,2)+'_meas.v1.fits',/overwrite 
                             
                            MWRFITS,meas,outfile,/create 
                            MWRFITS,chinfo[i],outfile,/silent# add chip stucture for this chip 
                        CHIPBOMB: 
                     
                    # Total number of measurements 
                    expinfo['nmeas'] = np.sum(chinfo['nmeas']) 
                     
                    # Meta-data file 
                    metafile = expdir+'/'+base+'_meta.fits' 
                    logger.info('Writing metadata to '+metafile)
                    hdus = fits.HDUList()
                    hdu.append(fits.table_to_hdu(expinfo))
                    hdu.append(fits.table_to_hdu(chinfo))  # add chip structure to second extension 
                    hdu.writeto(metafile,overwrite=True)
                     
                    dt = time.time()-t00 
                    logger.info('dt = %.2f sec.' % dt)
                     
 
                
