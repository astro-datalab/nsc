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
from . import utils,query,modelmag

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

def loadheader(headfile):
    """
    Load header file
    """

    headlines = dln.readlines(headfile)
    lo = dln.grep(headlines,'^SIMPLE  =',index=True)
    begind = dln.grep(headlines,'^XTENSION',index=True)
    begind = lo+begind
    #endind = dln.grep(headlines,'^END',index=True)
    endind = np.concatenate((np.array(begind[1:])-1,np.array([len(headlines)-1])))
    headdict = {}
    # Loop over the extendions
    for i in range(len(begind)):
        lines = headlines[begind[i]:endind[i]+1]
        if i==0:
            headdict['main'] = lines
        else:
            numind = dln.grep(lines,'^CCDNUM',index=True)
            line1 = lines[numind[0]]
            ccdnum = int(line1.split('/')[0].split('=')[1])
            extind = dln.grep(lines,'^EXTNUM',index=True)
            headdict[ccdnum] = lines
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

    outdir = expdir.replace('instcal/v4/','instcal/v4/gaiaxpsynthtest/')
    outfile = outdir+'/'+base+'_meta_modelmags.fits' 
    if os.path.exists(outdir)==False:
        os.makedirs(outdir)
    
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
    if expdir.endswith('.tar.gz'):
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
    else:
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
        headfile = os.path.join(expdir,base+'.hdr')
        if os.path.exists(headfile)==False:
            headfile = os.path.join(dldir,'dnidever','nsc','instcal',version,
                                    'header',instrument,night,base+'.hdr')
            if os.path.exists(headfile)==False:
                print(headfile+' not found')
                return
                #raise ValueError(headfile+' not found')
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
    ncatarr = np.zeros(ncatfiles,int) 
    for i in range(ncatfiles):
        head = fits.getheader(catfiles[i],1)
        ncatarr[i] = head['NAXIS2']
    ncat = int(np.sum(ncatarr)) 
    logger.info(str(ncat)+' total sources')
    if ncat == 0: 
        logger.info('No sources')
        return 
    # Create structure, exten=1 has header now 
    ind, = np.where(ncatarr > 0)
    cat1 = fits.getdata(catfiles[ind[0]],1)
    #cat1 = fits.getdata(catfiles[ind[0]],2)    
    cdt = cat1.dtype.descr
    cdt += [('ccdnum',int),('ebv',float),('ra',float),('dec',float),('raerr',float),
            ('decerr',float),('cmag',float),('cerr',float),('sourceid',(str,50)),
            ('filter',(str,50)),('mjd',float)]
    cat = np.zeros(ncat,dtype=np.dtype(cdt))
    cat = Table(cat)
    for c in cat.colnames:
        cat[c].name = c.lower()
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
    for i in range(ncatfiles):
        dum = os.path.splitext(os.path.basename(catfiles[i]))[0].split('_')
        ccdnum = int(dum[-1])
        hd = fits.getheader(catfiles[i],1)
        cat1 = Table.read(catfiles[i],1)
        #hd = fits.getheader(catfiles[i],2)        
        #cat1 = Table.read(catfiles[i],2)
        for c in cat1.colnames:
            cat1[c].name = c.lower()
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
        if version < 'v4':
            dum = fits.getdata(catfiles[i],1)
            hd1 = fits.Header.fromstring('\n'.join(dum[0][0]),sep='\n')
        else:
            dum = headdict[int(ccdnum)]
            hd1 = fits.Header.fromstring('\n'.join(dum),sep='\n')
        nx = hd1['NAXIS1']
        ny = hd1['NAXIS2']
        try:
            wcs1 = WCS(hd1)
        except:
            logger.info('Problem with WCS in header '+catfiles[i])    
            continue
        #extast,hd1,ast,noparams# check the WCS 
        #if noparams <= 0: 
        #    logger.info('Problem with WCS in header '+catfiles[i])
        #    goto,BOMB1 
        vcoo = wcs1.pixel_to_world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
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
            if len(bderr) > 0: 
                snr[bderr] = 1.087/temp[bderr]['magerr_iso']
            bderr, = np.where((temp['magerr_auto'] > 10) & (temp['magerr_iso'] > 10))
            if len(bderr) > 0: 
                snr[bderr] = 1 
            coorderr = 0.664*(temp['fwhm_world']*3600)/snr 
            temp['raerr'] = coorderr 
            temp['decerr'] = coorderr 
            # Stuff into main structure 
            cat[cnt:cnt+ncat1] = temp
            cnt += ncat1
            cenra = np.mean(dln.minmax(cat1['alpha_j2000'])) 
            # Wrapping around RA=0 
            if dln.valrange(cat1['alpha_j2000']) > 100: 
                ra = cat1['alpha_j2000']
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
            chinfo['cendec'][i] = np.mean(dln.minmax(cat1['delta_j2000'])) 

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
        #READLINE,expdir+'/'+base+'.head',head
        head = fits.getheader(fluxfile,0)
    else:
        dum = headdict['main']
        head = fits.Header.fromstring('\n'.join(dum),sep='\n')        
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
    cat['filter'] = filt
    cat['mjd'] = mjd 
    cat['sourceid'] = instrument+'.'+str(expnum)+'.'+np.char.array(cat['ccdnum'].astype(str))+'.'+np.char.array(cat['number'].astype(str))
         
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
            gdref, = np.where((inpref['ra'] >= np.min(cat[racol])-0.01) & 
                              (inpref['ra'] <= np.max(cat[racol])+0.01) &
                              (inpref['dec'] >= np.min(cat[deccol])-0.01) & 
                              (inpref['dec'] <= np.max(cat[deccol])+0.01))
        else: 
            ra = cat[racol]
            bdra, = np.where(ra > 180) 
            if len(bdra) > 0 : 
                ra[bdra]-=360 
            gdref, = np.where((inpref['ra'] <= np.max(ra)-0.01) & 
                              (inpref['ra'] >= np.min(ra+360)-0.01) &
                              (inpref['dec'] >= np.min(cat[deccol])-0.01) & 
                              (inpref['dec'] <= np.max(cat[deccol])+0.01))
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
    gdcat, = np.where(np.isfinite(cat[racol]) & np.isfinite(cat[deccol]))
    if len(gdcat)==0:
        logger.info('No catalog stars with good PSF coordinates')
        return
    ind1,ind2,dist = coords.xmatch(gaia['ra'],gaia['dec'],cat[racol][gdcat],cat[deccol][gdcat],1.0)
    ngmatch = len(ind1)
    if ngmatch == 0: 
        logger.info('No gaia matches')
        return 
    allgaiaind = np.zeros(ncat,int)-1
    allgaiaind[gdcat[ind2]] = ind1 
    allgaiadist = np.zeros(ncat,float)+999999. 
    allgaiadist[gdcat[ind2]] = coords.sphdist(gaia['ra'][ind1],gaia['dec'][ind1],
                                              cat[racol][gdcat[ind2]],cat[deccol][gdcat[ind2]])*3600 
    # CCD loop 
    for i in range(nchips): 
        if chinfo['nsources'][i]==0: 
            continue
        # Get chip sources using CCDNUM 
        chind1,chind2 = dln.match(chinfo['ccdnum'][i],cat['ccdnum'])
        nchmatch = len(chind1)
        cat1 = cat[chind2] 
        # Gaia matches for this chip 
        gaiaind1 = allgaiaind[chind2] 
        gaiadist1 = allgaiadist[chind2] 
        gmatch, = np.where((gaiaind1 > -1) & (gaiadist1 <= 0.5))   # get sources with Gaia matches 
        if len(gmatch) == 0: 
            gmatch, = np.where(gaiaind1 > -1 and gaiadist1 <= 1.0) 
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
        #ROTSPHCEN,gaia2.ra,gaia2.dec,chinfo[i].cenra,chinfo[i].cendec,gaialon,gaialat,/gnomic 
        gaialon,gaialat = coords.rotsphcen(gra_epoch,gdec_epoch,chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        lon1,lat1 = coords.rotsphcen(cat2[racol],cat2[deccol],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
        # ---- Fit RA as function of RA/DEC ---- 
        londiff = gaialon-lon1 
        err = None
        if 'ra_error' in gaia2.colnames:
            err = np.sqrt(gaia2['ra_error']**2 + cat2['raerr']**2) 
        if 'ra_error' not in gaia2.colnames and 'e_ra_icrs' in gaia2.colnames:
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
        racoef,racov = curve_fit(utils.func_poly2d_wrap,[lon1[gdlon],lat1[gdlon]],
                                 londiff[gdlon],sigma=err[gdlon],p0=initpars)
        yfitall = utils.func_poly2d(lon1,lat1,*racoef)
        rarms1 = dln.mad((londiff[gdlon]-yfitall[gdlon])*3600.) 
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
        if 'dec_error' in gaia2.colnames: 
            err = np.sqrt(gaia2['dec_error']**2 + cat2['decerr']**2) 
        if 'dec_error' not in gaia2.colnames and 'e_dec_icrs' in gaia2.colnames:
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
        lon,lat = coords.rotsphcen(cat1[racol],cat1[deccol],chinfo['cenra'][i],chinfo['cendec'][i],gnomic=True)
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

    # DO WE WANT TO USE THE SOURCE EXTRACTOR OF DAOPHOT COORDINATES
    # AS THE DEFAULT RA/DEC???
    # don't have daophot coordinates for all the sources, why?
        
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
    gdcat, = np.where(np.isfinite(cat['ra']) & np.isfinite(cat['dec']))
    ind1,ind2,dist = coords.xmatch(ref['ra'],ref['dec'],cat['ra'][gdcat],cat['dec'][gdcat],dcr)
    nmatch = len(ind1)
    logger.info(str(nmatch)+' matches to reference catalog')
    if nmatch == 0: 
        logger.info('No matches to reference catalog')
        return
    ref1 = ref[ind1]
    cat1 = cat[gdcat[ind2]]
    # Use Gaia XP synthetic photometry
    # Get the model magnitudes 
    mmags = modelmag.modelmag(ref1,instfilt,cendec,eqnfile) 
    if len(mmags) == 1 and mmags[0] < -1000: 
        print('No good model mags')
        return
    # Using model magnitudes
    logger.info('Getting model magnitudes zeropoint')
    # Get the good sources 
    gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) & (~((cat1['flags'] & 16) > 0)) &
                      (cat1['mag_auto'] < 50) & (cat1['magerr_auto'] < 0.05) & (cat1['class_star'] > 0.8) &
                      (cat1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
    #  if the seeing is bad then class_star sometimes doesn't work well 
    if medfwhm > 1.8 and len(gdcat) < 100:
        gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) & (~((cat1['flags'] & 16) > 0)) &
                          (cat1['mag_auto'] < 50) &  (cat1['magerr_auto'] < 0.05) & 
                          (cat1['fwhm_world']*3600 < 2*medfwhm) & (mmags[:,0] < 50) & (mmags[:,1] < 5))
    if len(gdcat) == 0:
        logger.info('No good reference sources')
        return
    ref2 = ref1[gdcat] 
    mmags2 = mmags[gdcat,:] 
    cat2 = cat1[gdcat]
    # Matched structure
    #mag2 = cat2['mag_auto'] + 2.5*np.log10(exptime)   # correct for the exposure time
    mag2 = cat2['magpsf'] + 2.5*np.log10(exptime)   # correct for the exposure time    
    mstr = {'col':mmags2[:,2],'mag':mag2,'model':mmags2[:,0],
            'err':mmags2[:,1],'ccdnum':cat2['ccdnum']}
    mmatchedcatfile = outdir+'/'+base+'_modelmags_matchedcat.fits'
    logger.info('Saving modelmags matched catalogs to '+mmatchedcatfile)
    temp = Table(mstr)
    temp['source'] = ref2['source']
    temp['ra'] = ref2['ra']
    temp['dec'] = ref2['dec']    
    temp.write(mmatchedcatfile,overwrite=True)
    # Measure the zero-point
    mexpinfo = expinfo.copy()
    mchinfo = chinfo.copy()
    mexpinfo,mchinfo = fitzpterm(mstr,mexpinfo,mchinfo)

    # Print out the results 
    logger.info('NPHOTREFMATCH=%d' % mexpinfo['nrefmatch'])
    logger.info('EXPOSURE ZPTERM=%.4f +/- %.4f  SIG=%.4f mag' % (mexpinfo['zpterm'],mexpinfo['zptermerr'],mexpinfo['zptermsig']))
    logger.info('ZPSPATIALVAR:  RMS=%.4f RANGE=%.4f NCCD=%d' % 
                (mexpinfo['zpspatialvar_rms'],mexpinfo['zpspatialvar_range'],mexpinfo['zpspatialvar_nccd']))
        
    # Using Gaia XP synthetic photometry
    logger.info('Getting Gaia XP synthetic photometry zeropoint')
    refmagcol = 'gsynth_'+filt.lower()+'mag'
    referrcol = 'e_gsynth_'+filt.lower()+'mag'
    # Get the good sources 
    gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) &
                      (~((cat1['flags'] & 16) > 0)) & (cat1['magpsf'] < 50) &
                      (cat1['errpsf'] < 0.05) & (cat1['class_star'] > 0.8) &
                      (cat1['fwhm_world']*3600 < 2*medfwhm) &
                      (np.abs(cat1['sharp']) < 1) & (cat1['chi']<5) &
                      (ref1['bp'] > 0) & (ref1['bp'] < 50) &
                      (ref1['rp'] > 0) & (ref1['rp'] < 50) &
                      (ref1[refmagcol] > 0) & (ref1[refmagcol] < 50))
    #  if the seeing is bad then class_star sometimes doesn't work well 
    if medfwhm > 1.8 and len(gdcat) < 100:
        gdcat, = np.where((cat1['imaflags_iso'] == 0) & (~((cat1['flags'] & 8) > 0)) &
                          (~((cat1['flags'] & 16) > 0)) &
                          (cat1['magpsf'] < 50) & (cat1['errpsf'] < 0.05) & 
                          (cat1['fwhm_world']*3600 < 2*medfwhm) &
                          (ref1['bp'] > 0) & (ref1['bp'] < 50) &
                          (ref1['rp'] > 0) & (ref1['rp'] < 50) &
                          (ref1[refmagcol] > 0) & (ref1[refmagcol] < 50))
    if len(gdcat) == 0:
        logger.info('No good reference sources')
        return
    ref2 = ref1[gdcat] 
    mmags2 = mmags[gdcat,:]
    cat2 = cat1[gdcat]
    # Matched structure
    mag2 = cat2['magpsf'] + 2.5*np.log10(exptime)   # correct for the exposure time
    mstr = {'col':ref2['bp']-ref2['rp'],'mag':mag2,'model':ref2[refmagcol],
            'err':ref2[referrcol],'ccdnum':cat2['ccdnum']}
    gmatchedcatfile = outdir+'/'+base+'_gaiaxpsynth_matchedcat.fits'
    logger.info('Saving gaia matched catalogs to '+gmatchedcatfile)
    temp = Table(mstr)
    temp['source'] = ref2['source']
    temp['ra'] = ref2['ra']
    temp['dec'] = ref2['dec']
    temp.write(gmatchedcatfile,overwrite=True)    
    # Measure the zero-point
    gexpinfo = expinfo.copy()
    gchinfo = chinfo.copy()
    gexpinfo,gchinfo = fitzpterm(mstr,gexpinfo,gchinfo)

    # Print out the results 
    logger.info('NPHOTREFMATCH=%d' % gexpinfo['nrefmatch'])
    logger.info('EXPOSURE ZPTERM=%.4f +/- %.4f  SIG=%.4f mag' % (gexpinfo['zpterm'],gexpinfo['zptermerr'],gexpinfo['zptermsig']))
    logger.info('ZPSPATIALVAR:  RMS=%.4f RANGE=%.4f NCCD=%d' % 
                (gexpinfo['zpspatialvar_rms'],gexpinfo['zpspatialvar_range'],gexpinfo['zpspatialvar_nccd']))
                     
    # Meta-data file 
    mmetafile = outdir+'/'+base+'_meta_modelmags.fits' 
    logger.info('Writing model magnitudes metadata to '+mmetafile)
    hdu = fits.HDUList()
    hdu.append(fits.table_to_hdu(mexpinfo))
    hdu.append(fits.table_to_hdu(mchinfo))  # add chip structure to second extension 
    hdu.writeto(mmetafile,overwrite=True)
    hdu.close()
    gmetafile = outdir+'/'+base+'_meta_gaiaxpsynth.fits' 
    logger.info('Writing gaia xp synthetic photometry metadata to '+gmetafile)
    hdu = fits.HDUList()
    hdu.append(fits.table_to_hdu(gexpinfo))
    hdu.append(fits.table_to_hdu(gchinfo))  # add chip structure to second extension 
    hdu.writeto(gmetafile,overwrite=True)    
    hdu.close()
    
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


 
                
