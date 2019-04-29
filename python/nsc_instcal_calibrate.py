#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from astropy.table import Table, Column
import time
import shutil
import re
import subprocess
import glob
import logging
import socket
#from scipy.signal import convolve2d
from scipy.ndimage.filters import convolve

if __name__ == "__main__":

# Calibrate exposures for one NSC exposure

    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # on thing/hulk use
    if (host == "thing") | (host == "hulk"):
        dir = "/dl1/users/dnidever/nsc/instcal/"+verdir
        mssdir = "/mss1/"
        tmproot = "/d0/dnidever/nsc/instcal/"+verdir+"tmp/"
    # on gp09 use
    if (host == "gp09") | (host == "gp08") | (host == "gp07") | (host == "gp06") | (host == "gp05"):
        dir = "/net/dl1/users/dnidever/nsc/instcal/"+verdir
        mssdir = "/net/mss1/"
        tmproot = "/data0/dnidever/nsc/instcal/"+verdir+"tmp/"

    t0 = time.time()

    #print sys.argv

    # Not enough inputs
    n = len(sys.argv)
    if n < 1:
        print "Syntax - nsc_instcal_calibrate.py expdir"
        sys.exit()

    # Inputs
    expdir = sys.argv[1]
    # Check that the directory exist
    if os.path.exists(expdir) == False:
        print(expdir+" NOT FOUND")
        sys.exit()

    # Get exposure base name
    dum = expdir.split('/')
    base = (dum[-1] if dum[-1]!="" else dum[-2])

    # Get version number
    lo = expdir.find('nsc/instcal/')
    dum = expdir[lo+12:]
    version = dum[0:dum.find('/')]


    # Set up logging to screen and logfile
    #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    logfile = tmpdir+"/"+base+".log"
    #fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, fileName))
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    rootLogger.setLevel(logging.NOTSET)

    rootLogger.info("Calibrate catalogs for exposure "+base+" in "+expdir+" on host="+host)
    rootLogger.info("Version = "+version)

    # Check for output file
    if file_test(outfile) eq 1 and not keyword_set(redo) then begin
         rootLogger.info("outfile,' already exists and /redo not set.")
        return

    # What instrument is this?
    instrument = "c4d"  # by default
    if expdir.find("/k4m/") != -1: instrument = "k4m"
    if expdir.find("/ksb/") != -1: instrument = "ksb"
    rootLogger.info("This is a "+instrument+" exposure")


    # Model magnitude equation file
    eqnfile = dldir+"users/dnidever/nsc/instcal/"+version+"/config/modelmag_equations.txt"
    rootLogger.info("Using model magnitude equation file "+eqnfile)
    if os.path.exists(eqnfile) is False:
        rootLogger.info(eqnfile+" NOT FOUND")
        sys.exit()
    READLINE,eqnfile,eqnlines
     rootLogger.info(eqnlines)


    # Step 1. Read in the catalogs 
    #-----------------------------
    rootLogger.info("")
    rootLogger.info("Step 1. Read in the catalogs")
    rootLogger.info("-----------------------------")
    catfiles = []
    catfiles1 = glob.glob(expdir+"/"+base+"_[1-9].fits")
    if len(catfiles1) > 0: catfiles.append(catfiles1)
    catfiles2 = glob.glob(expdir+"/"+base+"_[1-9][0-9].fits")
    if len(catfiles2) > 0: catfiles.append(catfiles2)    
    ncatfiles = len(catfiles)
    if ncatfiles == 0:
        rootLogger.info("No catalog files found")        
        sys.exit()
    nchips = ncatfiles
    rootLogger.info(str(ncatfiles)+" catalogs found")


    # Check that this isn't a problematic Mosaic3 exposure
    if stregex(expdir,'/k4m/',/boolean) eq 1:
        dum = fits.getdata(catfiles[0],1)
        head0 = dum["field_header_card"]
        pixcnt = head0["PIXCNT*"] #,count=npixcnt)
        if pixcnt > 0:
            rootLogger.info("This is a Mosaic3 exposure with pixel shift problems")
            sys.exit()
        wcscal = head0["WCSCAL"]  #,count=nwcscal)
        if (nwcscal > 0) and strtrim(wcscal,2) eq 'Failed':
            rootLogger.info("This is a Mosaic3 exposure with failed WCS calibration")
            sys.exit()

    # Figure out the number of sources
    ncat = 0
    for fil in catfiles:
        head = fits.getheader(fil,2)
        ncat += head["NAXIS2"]
    rootLogger.info(str(ncat)+" total sources")

    # Create structure, exten=1 has header now
    cat1 = fits.getdata(catfiles[0],2,header=False)
    if type(cat) != "astropy.io.fits.fitsrec.FITS_rec":
        rootLogger.info("Chip 1 catalog is empty.")
        sys.exit()
    cat1 = Table(cat1)
    #catdtype = cat1.dtype
    newcat = cat1.copy()[0:1]
    newnames = ['CCDNUM','EBV','RA','DEC','RAERR','DECERR','CMAG','CERR','SOURCEID','FILTER','MJD']
    newtypes = ['int','float','float64','float64','float','float','float','float','S100','S50','float64']
    newvals = [0L,0.0,0.0,0.0,0.0,0.0,99.99,9.99,' ',' ',0.0]
    newcols = []
    for n,t,v in zip(newnames,newtypes,newvals):
        col = Column(name=n,length=1,dtype=t)
        col[:] = v
        newcols.append(col)
    newcat.add_columns(newcols)
    cat = np.zeros(ncat,dtype=newcat)
    # Start the chips summary structure
    chdtype = np.dtype([('expdir',np.str,100),('instrument',np.str,10),('filename',np.str,200),
                        ('measfile',np.str,200),('ccdnum',int),('nsources',np.long),('nmeas',np.long),
                        ('cenra',np.float64),('cendec',np.float64),('ngaiamatch',np.long),
                        ('ngoodgaiamatch',np.long),('rarms',float),('rastderr',float),('racoef',np.float64,4),
                        ('decrms',float),('decstderr',float),('deccoef',np.float64,4),('vra',np.flat64,4),
                        ('vdec',np.float64,4),('zpterm',float),('zptermerr',float),('nrefmatch',np.long),
                        ('depth95',float),('depth10sig',float)])
    chstr = np.zeros(nchips,dtype=chdtype)
    chstr['expdir'] = expdir
    chstr['instrument'] = instrument
    chstr['cenra'] = 999999
    chstr['cendec'] = 999999
    chstr['rarms'] = 999999
    chstr['rastderr'] = 999999
    chstr['decrms'] = 999999
    chstr['decstderr'] = 999999
    chstr['zpterm'] = 999999
    chstr['zptermerr'] = 999999
    chstr['depth95'] = 99.99
    chstr['depth95sig'] = 99.99    
    # Load the files
    cnt = 0L
    for i in range(ncatfiles):
        dum = os.path.splitext(os.path.basename(catfiles[i]))[0]
        ccdnum = int(dum.split('_')[-1])
        cat1, hd = fits.getdata(catfiles[i],2,header=True)
        ncat1 = hd['naxis2'] # handles empty catalogs
        chstr['filename'][i] = catfiles[i]
        chstr['ccdnum'][i] = ccdnum
        chstr['nsources'][i] = ncat1
        # Get the chip corners
        dum = fits.getdata(catfiles[i],1)
        hd1 = dum['FIELD_HEADER_CARD']
        nx = hd1['NAXIS1']
        ny = hd1['NAXIS2']
        w = WCS(hd1)
        if w.is_celestial is True:
            vra,vdec = w.wcs_pix2world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1],0)
            #foot = w.calc_footprint()   # different order
            #vra = foot[:,0]
            #vdec = foot[:,1]
            chstr['vra'][i] = vra
            chstr['vdec'][i] = vdec
            if ncat1 > 0:
                temp = cat[cnt:cnt+ncat1-1]
                STRUCT_ASSIGN,cat1,temp,/nozero
                temp.ccdnum = ccdnum
                temp.ra = cat1.alpha_j2000  # add these here in case the astrometric correction fails later on
                temp.dec = cat1.delta_j2000
                # Add coordinate uncertainties
                #   sigma = 0.644*FWHM/SNR
                #   SNR = 1.087/magerr
                snr = 1.087/temp['magerr_auto']
                bderr = where(temp.magerr_auto gt 10 and temp.magerr_iso lt 10,nbderr)
                if nbderr gt 0 then snr[bderr]=1.087/temp[bderr].magerr_iso
                bderr = where(temp.magerr_auto gt 10 and temp.magerr_iso gt 10,nbderr)
                if nbderr gt 0 then snr[bderr] = 1
                coorderr = 0.664*(temp.fwhm_world*3600)/snr
                temp['raerr'] = coorderr
                temp['decerr'] = coorderr
                # Stuff into main structure
                cat[cnt:cnt+ncat1-1] = temp
                cnt += ncat1
                cenra = mean(minmax(cat1.alpha_j2000))
                # Wrapping around RA=0
                if range(cat1['alpha_j2000']) > 100:
                    ra = cat1['alpha_j2000']
                    bdra = where(ra gt 180,nbdra)
                    if nbdra gt 0 then ra[bdra]-=360
                    bdra2 = where(ra lt -180,nbdra2)
                    if nbdra2 gt 0 then ra[bdra2]+=360
                    cenra = mean(minmax(ra))
                    if cenra lt 0 then cenra+=360
                chstr['cenra'][i] = cenra
                chstr['cendec'][i] = mean(minmax(cat1.delta_j2000))
        # no good WCS
        rootLogger.info("Problem with WCS in header "+catfiles[i])


    # Exposure level values
    gdchip = where(chstr.nsources gt 0 and chstr.cenra lt 400,ngdchip)
    if ngdchip eq 0:
        rootLogger.info("No good chip catalogs with good WCS.")
        sys.exit()
    # Central coordinates of the entire field
    cendec = mean(minmax(chstr[gdchip].cendec))
    decrange = range(chstr[gdchip].vdec)
    cenra = mean(minmax(chstr[gdchip].cenra))
    rarange = range(chstr[gdchip].vra)*cos(cendec/!radeg)
    # Wrapping around RA=0
    if range(minmax(chstr[gdchip].cenra)) > 100:
        ra = chstr[gdchip].cenra
        bdra = where(ra gt 180,nbdra)
        if nbdra gt 0 then ra[bdra]-=360
        bdra2 = where(ra lt -180,nbdra2)
        if nbdra2 gt 0 then ra[bdra2]+=360
        cenra = mean(minmax(ra))
        if cenra lt 0 then cenra+=360
        # use chip VRA to get RA range
        vra = chstr[gdchip].vra
        bdra = where(vra gt 180,nbdra)
        if nbdra gt 0 then vra[bdra]-=360
        bdra2 = where(vra lt -180,nbdra2)
        if nbdra2 gt 0 then vra[bdra2]+=360
        rarange = range(vra)*cos(cendec/!radeg)
        rawrap = 1
    else:
        rawrap=0
    rootLogger.info("CENRA  = "+str(cenra))
    rootLogger.info("CENDEC = "+str(cendec))
    glactc,cenra,cendec,2000.0,glon,glat,1,/deg
    rootLogger.info("GLON = "+str(glon))
    rootLogger.info("GLAT = "+str(glat))
    # Number of good sources
    goodsources = where(cat.imaflags_iso eq 0 and not ((cat.flags and 8) eq 8) and not ((cat.flags and 16) eq 16) and $
                    cat.mag_auto lt 50,ngoodsources)
    rootLogger.info("GOOD SRCS = "+str(ngoodsources))

    # Measure median seeing FWHM
    gdcat = where(cat.mag_auto lt 50 and cat.magerr_auto lt 0.05 and cat.class_star gt 0.8,ngdcat)
    medfwhm = median(cat[gdcat].fwhm_world*3600.)
    rootLogger.info("FWHM = "+str(medfwhm)+" arcsec")

    # Load the logfile and get absolute flux filename
    READLINE,expdir+'/'+base+'.log',loglines
    ind = where(stregex(loglines,'Step #2: Copying InstCal images from mass store archive',/boolean) eq 1,nind)
    fline = loglines[ind[0]+1]
    lo = strpos(fline,'/archive')
    # make sure the mss1 directory is correct for this server
    fluxfile = mssdir+strtrim(strmid(fline,lo+1),2)
    wline = loglines[ind[0]+2]
    lo = strpos(wline,'/archive')
    wtfile = mssdir+strtrim(strmid(wline,lo+1),2)
    mline = loglines[ind[0]+3]
    lo = strpos(mline,'/archive')
    maskfile = mssdir+strtrim(strmid(mline,lo+1),2)
    
    # Load the meta-data from the original header
    #READLINE,expdir+'/'+base+'.head',head
    head = headfits(fluxfile,exten=0)
    filterlong = strtrim(sxpar(head,'filter',count=nfilter),2)
    if nfilter == 0:
        dum = mrdfits(catfiles[0],1,/silent)
        hd1 = dum.field_header_card
        filterlong = strtrim(sxpar(hd1,'filter'),2)
    if strmid(filterlong,0,2) eq 'VR' then filter='VR' else filter=strmid(filterlong,0,1)
    if filterlong eq 'bokr' then filter='r'
    expnum = sxpar(head,'expnum')
    if instrument == 'ksb':  # Bok doesn't have expnum
        #DTACQNAM= '/data1/batch/bok/20160102/d7390.0049.fits.fz'   
        dtacqnam = sxpar(head,'DTACQNAM',count=ndtacqnam)
        if ndtacqnam == 0:
            rootLogger.info("I cannot create an EXPNUM for this Bok exposure")
            sys.exit()
        bokbase = file_basename(dtacqnam)
        dum = strsplit(bokbase,'.',/extract)
        boknight = strmid(dum[0],1)
        boknum = dum[1]
        expnum = boknight+boknum  # concatenate the two numbers
    exptime = sxpar(head,'exptime',count=nexptime)
    if nexptime == 0:
        dum = mrdfits(catfiles[0],1,/silent)
        hd1 = dum.field_header_card
        exptime = sxpar(hd1,'exptime')
    dateobs = sxpar(head,'date-obs')
    airmass = sxpar(head,'airmass')
    mjd = date2jd(dateobs,/mjd)
    rootLogger.info("FILTER = "+filter)
    rootLogger.info("EXPTIME = "+str(exptime+" sec.")
    rootLogger.info("MJD = "+str(mjd))

    # Set some catalog values
    cat.filter = filter
    cat.mjd = mjd
    cat.sourceid = instrument+'.'+strtrim(expnum,2)+'.'+strtrim(cat.ccdnum,2)+'.'+strtrim(cat.number,2)

    # Start the exposure-level structure
    expstr = {file:fluxfile,wtfile:wtfile,maskfile:maskfile,instrument:'',base:base,expnum:long(expnum),ra:0.0d0,dec:0.0d0,dateobs:string(dateobs),$
              mjd:0.0d,filter:filter,exptime:float(exptime),airmass:0.0,nsources:long(ncat),ngoodsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,ebv:0.0,ngaiamatch:0L,$
              ngoodgaiamatch:0L,zptype:0,zpterm:999999.0,zptermerr:99999.0,zptermsig:999999.0,zpspatialvar_rms:999999.0,zpspatialvar_range:999999.0,$
              zpspatialvar_nccd:0,nrefmatch:0L,ngoodrefmatch:0L,depth95:99.99,depth10sig:99.99}
    expstr['instrument'] = instrument
    expstr['ra'] = cenra
    expstr['dec'] = cendec
    expstr['mjd'] = mjd
    #expstr['mjd'] = photred_getmjd('','CTIO',dateobs=dateobs)
    expstr['nchips'] = nchips
    expstr['airmass'] = airmass
    expstr['fwhm'] = medfwhm
    expstr['ngoodsources'] = ngoodsources

    # Step 2. Load the reference catalogs
    #------------------------------------
    rootLogger.info("")
    rootLogger.info("Step 2. Load the reference catalogs")
    rootLogger.info("------------------------------------")

    # Getting reference catalogs
    if n_elements(inpref) == 0:
        # Search radius
        radius = 1.1 * sqrt( (0.5*rarange)^2 + (0.5*decrange)^2 ) 
        #ref = GETREFDATA_V3(filter,cenra,cendec,radius,count=count)
        ref = GETREFDATA(instrument+'-'+filter,cenra,cendec,radius,count=count)
        if count == 0:
            printlog,logfi,'No Reference Data'
            sys.exit()
    # Using input reference catalog
    else:
        rootLogger.info("Reference catalogs input")
        if rawrap == 0:
            gdref = where(inpref.ra ge min(cat.alpha_j2000)-0.01 and inpref.ra le max(cat.alpha_j2000)+0.01 and $
                          inpref.dec ge min(cat.delta_j2000)-0.01 and inpref.dec le max(cat.delta_j2000)+0.01,ngdref)
        else:
            ra = cat.alpha_j2000
            bdra = where(ra gt 180,nbdra)
            if nbdra gt 0 then ra[bdra]-=360
            gdref = where((inpref.ra le max(ra)+0.01 or inpref.ra ge min(ra+360)-0.01) and $
                          inpref.dec ge min(cat.delta_j2000)-0.01 and inpref.dec le max(cat.delta_j2000)+0.01,ngdref)
        ref = inpref[gdref]
        rootLogger.info(str(ngdref)+" reference stars in our region")

    # Step 3. Astrometric calibration
    #----------------------------------
    # At the chip level, linear fits in RA/DEC
    rootLogger.info("")
    rootLogger.info("Step 3. Astrometric calibration")
    rootLogger.info("--------------------------------")
    # Get reference catalog with Gaia values
    gdgaia = where(ref.source gt 0,ngdgaia)
    gaia = ref[gdgaia]
    # Match everything to Gaia at once, this is much faster!
    SRCMATCH,gaia.ra,gaia.dec,cat.alpha_j2000,cat.delta_j2000,1.0,ind1,ind2,/sph,count=ngmatch
    if ngmatch == 0:
        rootLogger.info("No gaia matches")
        sys.exit()
    allgaiaind = lonarr(ncat)-1
    allgaiaind[ind2] = ind1
    allgaiadist = fltarr(ncat)+999999.
    allgaiadist[ind2] = sphdist(gaia[ind1].ra,gaia[ind1].dec,cat[ind2].alpha_j2000,cat[ind2].delta_j2000,/deg)*3600
    # CCD loop
    for i in range(nchips):
        if chstr[i].nsources eq 0 then goto,BOMB
        # Get chip sources using CCDNUM
        MATCH,chstr[i].ccdnum,cat.ccdnum,chind1,chind2,/sort,count=nchmatch
        cat1 = cat[chind2]
        # Gaia matches for this chip
        gaiaind1 = allgaiaind[chind2]
        gaiadist1 = allgaiadist[chind2]
        gmatch = where(gaiaind1 gt -1 and gaiadist1 le 0.5,ngmatch)  # get sources with Gaia matches
        if ngmatch eq 0 then gmatch = where(gaiaind1 gt -1 and gaiadist1 le 1.0,ngmatch)
        if ngmatch < 5:
            rootLogger.info("Not enough Gaia matches")
            # Add threshold to astrometric errors
            cat1.raerr = sqrt(cat1.raerr^2 + 0.100^2)
            cat1.decerr = sqrt(cat1.decerr^2 + 0.100^2)
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
        if tag_exist(gaia1b,'PMRA') and tag_exist(gaia1b,'PMDEC'):   # we have proper motion information
            qcuts1 = where(cat1b.imaflags_iso eq 0 and not ((cat1b.flags and 8) eq 8) and not ((cat1b.flags and 16) eq 16) and $
                           cat1b.mag_auto lt 50 and finite(gaia1b.pmra) eq 1 and finite(gaia1b.pmdec) eq 1,nqcuts1)
        else:
            qcuts1 = where(cat1b.imaflags_iso eq 0 and not ((cat1b.flags and 8) eq 8) and not ((cat1b.flags and 16) eq 16) and $
                   cat1b.mag_auto lt 50,nqcuts1)
        if nqcuts1 == 0:
            rootLogger.info("Not enough stars after quality cuts")
            # Add threshold to astrometric errors
            cat1.raerr = sqrt(cat1.raerr^2 + 0.100^2)
            cat1.decerr = sqrt(cat1.decerr^2 + 0.100^2)
            cat[chind2] = cat1
            goto,BOMB
        gaia2 = gaia1b[qcuts1]
        cat2 = cat1b[qcuts1]

        # Precess the Gaia coordinates to the epoch of the observation
        # The reference epoch for Gaia DR2 is J2015.5 (compared to the
        # J2015.0 epoch for Gaia DR1).
        if tag_exist(gaia2,'PMRA') and tag_exist(gaia2,'PMDEC'):
            gaiamjd = 57206.0d0
            delt = (mjd-gaiamjd)/365.242170d0   # convert to years
            # convert from mas/yr->deg/yr and convert to angle in RA
            gra_epoch = gaia2.ra + delt*gaia2.pmra/3600.0d0/1000.0d0/cos(gaia2.dec*d2r)
            gdec_epoch = gaia2.dec + delt*gaia2.pmdec/3600.0d0/1000.0d0
        else:
            gra_epoch = gaia2.ra
            gdec_epoch = gaia2.dec

        # Rotate to coordinates relative to the center of the field
        #ROTSPHCEN,gaia2.ra,gaia2.dec,chstr[i].cenra,chstr[i].cendec,gaialon,gaialat,/gnomic
        ROTSPHCEN,gra_epoch,gdec_epoch,chstr[i].cenra,chstr[i].cendec,gaialon,gaialat,/gnomic
        ROTSPHCEN,cat2.alpha_j2000,cat2.delta_j2000,chstr[i].cenra,chstr[i].cendec,lon1,lat1,/gnomic
        # ---- Fit RA as function of RA/DEC ----
        londiff = gaialon-lon1
        undefine,err
        if tag_exist(gaia2,'RA_ERROR') then err = sqrt(gaia2.ra_error^2 + cat2.raerr^2)
        if tag_exist(gaia2,'RA_ERROR') eq 0 and tag_exist(gaia2,'E_RA_ICRS') then err = sqrt(gaia2.e_ra_icrs^2 + cat2.raerr^2)
        if n_elements(err) eq 0 then err=cat2.raerr
        lonmed = median([londiff])
        lonsig = mad([londiff]) > 1e-5   # 0.036"
        gdlon = where(abs(londiff-lonmed) lt 3.0*lonsig,ngdlon)  # remove outliers
        if ngdlon gt 5 then npars = 4 else npars=1  # use constant if not enough stars
        initpars = dblarr(npars)
        initpars[0] = median([londiff])
        parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
        racoef = MPFIT2DFUN('func_poly2d',lon1[gdlon],lat1[gdlon],londiff[gdlon],err[gdlon],initpars,status=status,dof=dof,$
                            bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)
        yfitall = FUNC_POLY2D(lon1,lat1,racoef)
        rarms1 = MAD((londiff[gdlon]-yfit)*3600.)
        rastderr = rarms1/sqrt(ngdlon)
        # Use bright stars to get a better RMS estimate
        gdstars = where(cat2.fwhm_world*3600 lt 2*medfwhm and 1.087/cat2.magerr_auto gt 50,ngdstars)
        if ngdstars lt 20 then gdstars = where(cat2.fwhm_world*3600 lt 2*medfwhm and 1.087/cat2.magerr_auto gt 30,ngdstars)
        if ngdstars > 5:
            diff = (londiff-yfitall)*3600.
            rarms = MAD(diff[gdstars])
            rastderr = rarms/sqrt(ngdstars)
        else:
            rarms=rarms1
        # ---- Fit DEC as function of RA/DEC -----
        latdiff = gaialat-lat1
        undefine,err
        if tag_exist(gaia2,'DEC_ERROR') then err = sqrt(gaia2.dec_error^2 + cat2.decerr^2)
        if tag_exist(gaia2,'DEC_ERROR') eq 0 and tag_exist(gaia2,'E_DEC_ICRS') then err = sqrt(gaia2.e_de_icrs^2 + cat2.decerr^2)
        if n_elements(err) eq 0 then err=cat2.decerr
        latmed = median([latdiff])
        latsig = mad([latdiff]) > 1e-5  # 0.036"
        gdlat = where(abs(latdiff-latmed) lt 3.0*latsig,ngdlat)  # remove outliers
        if ngdlat gt 5 then npars = 4 else npars=1  # use constant if not enough stars
        initpars = dblarr(npars)
        initpars[0] = median([latdiff])
        parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
        deccoef = MPFIT2DFUN('func_poly2d',lon1[gdlat],lat1[gdlat],latdiff[gdlat],err[gdlat],initpars,status=status,dof=dof,$
                             bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)
        yfitall = FUNC_POLY2D(lon1,lat1,deccoef)
        decrms1 = MAD((latdiff[gdlat]-yfit)*3600.)
        decstderr = decrms1/sqrt(ngdlat)
        # Use bright stars to get a better RMS estimate
        if ngdstars > 5:
            diff = (latdiff-yfitall)*3600.
            decrms = MAD(diff[gdstars])
            decstderr = decrms/sqrt(ngdstars)
        else:
            decrms=decrms1
        rootLogger.info("  CCDNUM="+str(chstr[i].ccdnum)+"  NSOURCES="+str(nchmatch)+"  "+str(ngmatch)+"/"+str(nqcuts1)+" GAIA matches  RMS(RA/DEC)="+str(rarms)+"/"+str(decrms)+" STDERR(RA/DEC)="+str(rastderr)+"/"+str(decstderr)+" arcsec"
        # Apply to all sources
        ROTSPHCEN,cat1.alpha_j2000,cat1.delta_j2000,chstr[i].cenra,chstr[i].cendec,lon,lat,/gnomic
        lon2 = lon + FUNC_POLY2D(lon,lat,racoef)
        lat2 = lat + FUNC_POLY2D(lon,lat,deccoef)
        ROTSPHCEN,lon2,lat2,chstr[i].cenra,chstr[i].cendec,ra2,dec2,/reverse,/gnomic
        cat1.ra = ra2
        cat1.dec = dec2
        # Add to astrometric errors
        cat1.raerr = sqrt(cat1.raerr^2 + rarms^2)
        cat1.decerr = sqrt(cat1.decerr^2 + decrms^2)
        # Stuff back into the main structure
        cat[chind2] = cat1
        chstr[i].ngaiamatch = ngmatch
        chstr[i].ngoodgaiamatch = nqcuts1
        chstr[i].rarms = rarms
        chstr[i].rastderr = rastderr
        chstr[i].racoef = racoef
        chstr[i].decrms = decrms
        chstr[i].decstderr = decstderr
        chstr[i].deccoef = deccoef
        BOMB:

    # Get reddening
    glactc,cat.ra,cat.dec,2000.0,glon,glat,1,/deg
    ebv = dust_getval(glon,glat,/noloop,/interp)
    cat.ebv = ebv
        
    # Put in exposure-level information
    expstr.rarms = median(chstr.rarms)
    expstr.decrms = median(chstr.decrms)
    expstr.ebv = median(ebv)
    #expstr.gaianmatch = median(chstr.gaianmatch)
    expstr.ngaiamatch = total(chstr.ngaiamatch)
    expstr.ngoodgaiamatch = total(chstr.ngoodgaiamatch)


    # Step 4. Photometric calibration
    #--------------------------------
    rootLogger.info("")
    rootLogger.info("Step 4. Photometric calibration")
    rootLogger.info("-------------------------------")
    instfilt = instrument+'-'+filter    # instrument-filter combination

    # Now crossmatch with our catalog
    dcr = 1.0
    SRCMATCH,ref.ra,ref.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
    rootLogger.info(str(nmatch)+" matches to reference catalog")
    if nmatch == 0:
         rootLogger.info("No matches to reference catalog'
        goto,ENDBOMB
    ref1 = ref[ind1]
    cat1 = cat[ind2]
    # Get the model magnitudes
    mmags = GETMODELMAG(ref1,instfilt,cendec,eqnfile)
    # Get the good sources
    gdcat = where(cat1.imaflags_iso eq 0 and not ((cat1.flags and 8) eq 8) and not ((cat1.flags and 16) eq 16) and $
                  cat1.mag_auto lt 50 and cat1.magerr_auto lt 0.05 and cat1.class_star gt 0.8 and $
                  cat1.fwhm_world*3600 lt 2*medfwhm and mmags[*,0] lt 50 and mmags[*,1] lt 5,ngdcat)
    #  if the seeing is bad then class_star sometimes doens't work well
    if medfwhm > 1.8 and ngdcat < 100:
        gdcat = where(cat1.imaflags_iso eq 0 and not ((cat1.flags and 8) eq 8) and not ((cat1.flags and 16) eq 16) and $
                      cat1.mag_auto lt 50 and cat1.magerr_auto lt 0.05 and $
                      cat1.fwhm_world*3600 lt 2*medfwhm and mmags[*,0] lt 50 and mmags[*,1] lt 5,ngdcat)
    if ngdcat == 0:
         rootLogger.info("No good reference sources'
        goto,ENDBOMB
    ref2 = ref1[gdcat]
    mmags2 = mmags[gdcat,*]
    cat2 = cat1[gdcat]
    # Matched structure
    mag2 = cat2.mag_auto + 2.5*alog10(exptime)  # correct for the exposure time
    mstr = {col:float(mmags2[*,2]),mag:float(mag2),model:float(mmags2[*,0]),err:float(mmags2[*,1]),ccdnum:long(cat2.ccdnum)}
    # Measure the zero-point
    NSC_INSTCAL_CALIBRATE_FITZPTERM,mstr,expstr,chstr
    expstr.zptype = 1

    ENDBOMB:

    # Use self-calibration
    if expstr.nrefmatch le 5 and keyword_set(selfcal):
        NSC_INSTCAL_CALIBRATE_SELFCALZPTERM,expdir,cat,expstr,chstr
        expstr.zptype = 3
    # Apply the zero-point to the full catalogs
    # USE CHIP-LEVEL ZERO-POINTS WHEN POSSIBLE!!!
    print,'USE CHIP-LEVEL ZERO-POINTS WHEN POSSIBLE!!!'
    # Create an output catalog for each chip
    nsrc = long64(total(chstr.nsources,/cum))
    lo = [0L,nsrc[0:nchips-2]]
    hi = nsrc-1
    for i in range(nchips):
        ncat1 = hi[i]-lo[i]+1
        if ncat1 > 0:
            cat1 = cat[lo[i]:hi[i]]
            if chstr[i].nrefmatch gt 5 then begin
            chstr[i].zptype = 1
        else:
            chstr[i].zpterm = expstr.zperm
            chstr[i].zptermerr = expstr.zptermerr
            chstr[i].zptype = 2

        gdcatmag = where(cat1.mag_auto lt 50,ngd)
        cat1[gdcatmag].cmag = cat1[gdcatmag].mag_auto + 2.5*alog10(exptime) + chstr[i].zpterm
        cat1[gdcatmag].cerr = sqrt(cat1[gdcatmag].magerr_auto^2 + chstr[i].zptermerr^2)  # add calibration error in quadrature
        #for j=0,n_elements(cat1.mag_aper)-1 do begin
        #  gdcatmag = where(cat1.mag_aper[j] lt 50,ngd)
        #  cat1[gdcatmag].cmag = cat1[gdcatmag].mag_auto + 2.5*alog10(exptime) + zpterm
        #  cat1[gdcatmag].cerr = sqrt(cat1[gdcatmag].magerr_auto^2 + zptermerr^2)  # add calibration error in quadrature
        #endfor
        # Print out the results
         rootLogger.info("  CCDNUM=',strtrim(chstr[i].ccdnum,2),'  NREFSOURCES=',strtrim(chstr[i].nrefmatch,2),'  ZPTYPE=',strtrim(chstr[i].zptype,2),$
        '  ZPTERM=',stringize(zpterm,ndec=4),'+/-',stringize(zptermerr,ndec=4)
        cat[lo[i]:hi[i]] = cat1  # stuff back in
        #stop
    #stop
    #gdcatmag = where(cat.mag_auto lt 50,ngd)
    #cat[gdcatmag].cmag = cat[gdcatmag].mag_auto + 2.5*alog10(exptime) + expstr.zpterm
    #cat[gdcatmag].cerr = sqrt(cat[gdcatmag].magerr_auto^2 + expstr.zptermerr^2)  # add calibration error in quadrature
    # Print out the results
     rootLogger.info("NPHOTREFMATCH=',strtrim(expstr.nrefmatch,2)
     rootLogger.info("EXPOSURE ZPTERM=',stringize(expstr.zpterm,ndec=4),'+/-',stringize(expstr.zptermerr,ndec=4),'  SIG=',stringize(expstr.zptermsig,ndec=4),'mag'
     rootLogger.info("ZPSPATIALVAR:  RMS=',stringize(expstr.zpspatialvar_rms,ndec=3),' ',$
    'RANGE=',stringize(expstr.zpspatialvar_range,ndec=3),' NCCD=',strtrim(expstr.zpspatialvar_nccd,2)

    # Measure the depth
    #   need good photometry
    gdmag = where(cat.cmag lt 50,ngdmag)
    if ngdmag > 0:
        # Get 95% percentile depth
        cmag = cat[gdmag].cmag
        si = sort(cmag)
        cmag = cmag[si]
        depth95 = cmag[round(0.95*ngdmag)-1]
        expstr.depth95 = depth95
        chstr.depth95 = depth95
         rootLogger.info("95% percentile depth = '+stringize(depth95,ndec=2)+' mag'
        # Get 10 sigma depth
        #  S/N = 1.087/err
        #  so S/N=5 is for err=1.087/5=0.2174
        #  S/N=10 is for err=1.087/10=0.1087
        depth10sig = 99.99
        depind = where(cat.cmag lt 50 and cat.cmag gt depth95-3.0 and cat.cerr ge 0.0987 and cat.cerr le 0.1187,ndepind)
        if ndepind lt 5 then depind = where(cat.cmag lt 50 and cat.cmag gt depth95-3.0 and cat.cerr ge 0.0787 and cat.cerr le 0.1387,ndepind)
        if ndepind > 5:
            depth10sig = median([cat[depind].cmag])
        else:
            depind = where(cat.cmag lt 50,ndepind)
            if ndepind gt 0 then depth10sig=max([cat[depind].cmag])
         rootLogger.info("10sigma depth = '+stringize(depth10sig,ndec=2)+' mag'
        expstr.depth10sig = depth10sig
        chstr.depth10sig = depth10sig

    # Step 5. Write out the final catalogs and metadata
    #--------------------------------------------------
    if keyword_set(redo) and keyword_set(selfcal) and expstr.zptype == 2:
        # Create backup of original versions
         rootLogger.info("Copying meas and meta files to v1 versions'
        metafile = expdir+'/'+base+'_meta.fits'
        if file_test(metafile) eq 1 then FILE_COPY,metafile,expdir+'/'+base+'_meta.v1.fits',/overwrite

    # Create an output catalog for each chip
    nsrc = long64(total(chstr.nsources,/cum))
    lo = [0L,nsrc[0:nchips-2]]
    hi = nsrc-1
    for i in range(nchips):
        ncat1 = hi[i]-lo[i]+1
        if ncat1 == 0:
             rootLogger.info("No sources for CCDNUM='+strtrim(chstr[i].ccdnum,2)
            goto,CHIPBOMB
        cat1 = cat[lo[i]:hi[i]]

        # Apply QA cuts
        #----------------

        # Remove bad chip data
        # Half of chip 31 for MJD>56660
        #  c4d_131123_025436_ooi_r_v2 with MJD=56619 has problems too
        #  if the background b/w the left and right half is large then BAd
        lft31 = where(cat1.x_image lt 1024 and cat1.ccdnum eq 31,nlft31)
        rt31 = where(cat1.x_image ge 1024 and cat1.ccdnum eq 31,nrt31)
        if nlft31 > 10 and nrt31 > 10:
            lftback = median([cat1[lft31].background])
            rtback = median([cat1[rt31].background])
            mnback = 0.5*(lftback+rtback)
            sigback = mad(cat1.background)
            if abs(lftback-rtback) gt (sqrt(mnback)>sigback) then jump31=1 else jump31=0
            # rootLogger.info("  Big jump in CCDNUM 31 background levels'
        else:
            jump31=0

        badchip31 = 0
        if (expstr.mjd > 56600) | (jump31 == 1):
            badchip31 = 1
            # Remove bad measurements
            # X: 1-1024 okay
            # X: 1025-2049 bad
            # use 1000 as the boundary since sometimes there's a sharp drop
            # at the boundary that causes problem sources with SExtractor
            bdind = where(cat1.x_image gt 1000 and cat1.ccdnum eq 31,nbdind,comp=gdind,ncomp=ngdind)
            if nbdind > 0:   # some bad ones found
                if ngdind == 0:   # all bad
                    # rootLogger.info("NO useful measurements in ',list[i].file
                    undefine,cat1
                    ncat1 = 0
                    goto,CHIPBOMB
                else:
                    # rootLogger.info("  Removing '+strtrim(nbdind,2)+' bad chip 31 measurements, '+strtrim(ngdind,2)+' left.'
                    REMOVE,bdind,cat1
                    ncat1 = n_elements(cat1)

        # Make a cut on quality mask flag (IMAFLAGS_ISO)
        bdcat = where(cat1.imaflags_iso gt 0,nbdcat)
        if nbdcat > 0:
            # rootLogger.info("  Removing ',strtrim(nbdcat,2),' sources with bad CP flags.'
            if nbdcat eq ncat1 then goto,CHIPBOMB
            REMOVE,bdcat,cat1
            ncat1 = n_elements(cat1)

        # Make cuts on SE FLAGS
        #   this removes problematic truncatd sources near chip edges
        bdseflags = where( ((cat1.flags and 8) eq 8) or $             # object truncated
                           ((cat1.flags and 16) eq 16),nbdseflags)    # aperture truncate
        if nbdseflags > 0:
            # rootLogger.info("  Removing ',strtrim(nbdseflags,2),' truncated sources'
            if nbdseflags eq ncat1 then goto,CHIPBOMB
            REMOVE,bdseflags,cat1
            ncat1 = n_elements(cat1)

        # Removing low-S/N sources
        #  snr = 1.087/err
        snrcut = 5.0
        bdsnr = where(1.087/cat1.magerr_auto lt snrcut,nbdsnr)
        if nbdsnr > 0:
            # rootLogger.info("  Removing ',strtrim(nbdsnr,2),' sources with S/N<',strtrim(snrcut,2)
            if nbdsnr eq ncat1 then goto,CHIPBOMB
            REMOVE,bdsnr,cat1
            ncat1 = n_elements(cat1)

        # Convert to final format
        if ncat1 > 0:
            schema = {measid:'',objectid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,raerr:0.0,dec:0.0d0,decerr:0.0,$
                      mag_auto:0.0,magerr_auto:0.0,mag_aper1:0.0,magerr_aper1:0.0,mag_aper2:0.0,magerr_aper2:0.0,$
                      mag_aper4:0.0,magerr_aper4:0.0,mag_aper8:0.0,magerr_aper8:0.0,kron_radius:0.0,$
                      asemi:0.0,asemierr:0.0,bsemi:0.0,bsemierr:0.0,theta:0.0,thetaerr:0.0,fwhm:0.0,flags:0,class_star:0.0}
            meas = replicate(schema,ncat1)
            STRUCT_ASSIGN,cat1,meas,/nozero
            meas.measid = strtrim(cat1.sourceid,2)
            meas.exposure = base
            meas.x = cat1.x_image
            meas.y = cat1.y_image
            meas.mag_auto = cat1.cmag
            meas.magerr_auto = cat1.cerr
            meas.mag_aper1 = cat1.mag_aper[0] + 2.5*alog10(exptime) + chstr[i].zpterm
            meas.magerr_aper1 = cat1.magerr_aper[0]
            meas.mag_aper2 = cat1.mag_aper[1] + 2.5*alog10(exptime) + chstr[i].zpterm
            meas.magerr_aper2 = cat1.magerr_aper[1]
            meas.mag_aper4 = cat1.mag_aper[2] + 2.5*alog10(exptime) + chstr[i].zpterm
            meas.magerr_aper4 = cat1.magerr_aper[2]
            meas.mag_aper8 = cat1.mag_aper[4] + 2.5*alog10(exptime) + chstr[i].zpterm
            meas.magerr_aper8 = cat1.magerr_aper[4]
            meas.asemi = cat1.a_world * 3600.            # convert to arcsec
            meas.asemierr = cat1.erra_world * 3600.      # convert to arcsec
            meas.bsemi = cat1.b_world * 3600.            # convert to arcsec
            meas.bsemierr = cat1.errb_world * 3600.      # convert to arcsec
            meas.theta = 90-cat1.theta_world             # make CCW E of N
            meas.thetaerr = cat1.errtheta_world
            meas.fwhm = cat1.fwhm_world * 3600.          # convert to arcsec
            meas.class_star = cat1.class_star

            # Write to file
            outfile = expdir+'/'+base+'_'+strtrim(chstr[i].ccdnum,2)+'_meas.fits'
            chstr[i].nmeas = ncat1  # updating CHSTr
            chstr[i].measfile = outfile

            if keyword_set(redo) and keyword_set(selfcal) and expstr.zptype == 2:   # Create backup of original versions 
            if file_test(outfile) eq 1 then FILE_COPY,outfile,expdir+'/'+base+'_'+strtrim(chstr[i].ccdnum,2)+'_meas.v1.fits',/overwrite

            MWRFITS,meas,outfile,/create
            MWRFITS,chstr[i],outfile,/silent   # add chip stucture for this chip
        CHIPBOMB:


    # Meta-data file
    metafile = expdir+'/'+base+'_meta.fits'
     rootLogger.info("Writing metadata to ',metafile
    MWRFITS,expstr,metafile,/create
    MWRFITS,chstr,metafile,/silent  # add chip structure to second extension

    dt = systime(1)-t00
     rootLogger.info("dt = ',stringize(dt,ndec=2),' sec.'
