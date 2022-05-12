#!/usr/bin/env python

import os
import time
import numpy as np
import socket
from astropy.io import fits,ascii
from astropy.table import Table,vstack
import subprocess
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from astroquery.vizier import Vizier
from astropy.coordinates import Angle,SkyCoord
import healpy as hp
import astropy.units as u
from . import utils

Vizier.TIMEOUT = 600
Vizier.ROW_LIMIT = -1
Vizier.cache_location = None

def tempest_query(cenra,cendec,radius,refcat,nside=32,silent=False,logger=None):
    """
    Get reference catalog information from FITS files on tempest.
 
    Parameters
    ----------
    cenra : float
       Central RA for the search. 
    cendec : float
       Central DEC for the search. 
    radius : float
       Search radius in degrees. 
    refcat : table
       Reference catalog name (e.g. 2MASS, Gaia, etc.)
    silent : bool, optional
       Don't print anything to the screen. 
    logger : logging object, optional
       Logging object to use for printing messages.
 
    Returns
    -------
    ref : astropy table
        Search results from the reference catalog. 
 
    Example
    -------

    cat = tempest_query(cenra,cendec,radius,refcat) 
 
    """

    t0 = time.time() 

    if logger is None:
        logger = dln.basiclogger()


    # What healpix do we need to load
    upix = hp.ang2pix(nside,cenra,cendec,lonlat=True)
    cenvec = hp.ang2vec(cenra,cendec,lonlat=True)
    allpix = hp.query_disc(nside,cenvec,np.deg2rad(radius),inclusive=True)

    # Loop over healpix
    ref = None
    for p in allpix:
        reffile = '/home/x51j468/catalogs/'+refcat+'/ring'+str(nside)+'/'+str(p//1000)+'/'+str(p)+'.fits'
        if os.path.exists(reffile)==False:
            print(reffile,' NOT FOUND')
            continue
        tab = Table.read(reffile)
        # Sometimes the catalog files have a 2nd dimension
        for n in tab.colnames:
            if tab[n].ndim==2 and tab[n].shape[1]==1:
                tab[n] = tab[n].flatten()
        ntab = len(tab)
        if ref is None:
            ref = tab
        else:
            ref = vstack((ref,tab))

    if ref is None:
        return []
            
    # Do radius cut
    dist = coords.sphdist(cenra,cendec,ref['ra'],ref['dec'])
    if dist.ndim==2:
        dist = dist.flatten()
    gd, = np.where(dist <= radius)
    ref = ref[gd]
    
    return ref


def getrefcat(cenra,cendec,radius,refcat,version=None,saveref=False,
              savefile=None,nside=32,silent=False,logger=None):
    """
    Get reference catalog information from DL database 
 
    Parameters
    ----------
    cenra : float
       Central RA for the search. 
    cendec : float
       Central DEC for the search. 
    radius : float
       Search radius in degrees. 
    refcat : table
       Reference catalog name (e.g. 2MASS, Gaia, etc.)
    version : str
       Version of NSC reduction.. 
    saveref : bool, optional
       Save the output to SAVEFILE or a default filename. Default is False.
    savefile : str, optional
       The file to save to or search for existing catalog. 
    silent : bool, optional
       Don't print anything to the screen. 
    logger : logging object, optional
       Logging object to use for printing messages.
 
    Returns
    -------
    ref : astropy table
        Search results from the reference catalog. 
 
    Example
    -------

    cat = getrefcat(cenra,cendec,radius,refcat,savefile=savefile,saveref=saveref) 
 
    By D. Nidever  Sep 2017 
    Translated to Python by D. Nidever, April 2022
    """
     
    count = 0 
    t0 = time.time() 

    host = socket.gethostname()
    host = host.split('.')[0]
    if 'tempest' in host.lower():
        tempest = True
    else:
        tempest = False

    # Run query on tempest
    if tempest:
        return tempest_query(cenra,cendec,radius,refcat,nside=nside,silent=silent,logger=logger)

    if logger is None:
        logger = dln.basiclogger()
    
    # Check that we have psql installed 
    out = subprocess.check_output(['which','psql'],shell=False)
    if type(out) is bytes:
        out = out.decode()
    out = out.strip()
    if dln.size(out)>1:
        out = out[0]
    if os.path.exists(out) == 0: 
        raise ValueError('No PSQL found on this sytem.')
     
    # Temporary directory 
    # /tmp is often small and get get fille dup
    dldir,mssdir,localdir = utils.rootdirs()
    if version is None:
        version = 'v3' 
    fdir = dldir+'users/dnidever/nsc/instcal/'+version+'/' 
    tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/' 
     
    # FLIP THIS AROUND, INPUT SHOULD BE THE "EASY" VERSION!!! 
    refname = str(refcat).upper()
    if refname == 'II/312/AIS': 
        refname = 'GALEX' 
    elif refname == '2MASS-PSC': 
        refname = 'TMASS' 
    elif refname == '2MASS': 
        refname = 'TMASS' 
    elif refname == 'GAIA/GAIA': 
        refname = 'GAIA' 
    elif refname == 'Skymapper': 
        refname = 'SKYMAPPER' 
    elif refname == 'GLIMPSE': 
        refname = 'II/293/glimpse' 
    elif refname == 'SAGE': 
        refname = 'II/305/archive' 
    elif refname == 'ATLASREFCAT2': 
        refname = 'ATLAS' 
     
    if savefile is None:
        savefile = tmpdir+'ref_%.5f_%.5f_%.3f_%s.fits' % (cenra,cendec,radius,refname)
    if os.path.exists(os.path.abspath(os.path.dirname(savefile)))==False:
        os.makedirs(os.path.abspath(os.path.dirname(savefile)))

    if silent==False:
        logger.info('Querying %s: RA=%.5f DEC=%.5f Radius=%.3f' % (refname,cenra,cendec,radius))
     
    # Loading previously loaded file 
    if os.path.exists(savefile): 
        if silent==False:
            logger.info('Loading previously-saved file '+savefile)
        ref = Table.read(savefile) 
         
    # Do the Query 
    #-------------- 
    else: 
         
        # Use DataLab database search 
        #---------------------------- 
        if refname in ['TMASS','GAIA','GAIADR2','GAIAEDR3','PS','SKYMAPPER','SKYMAPPERDR2','ALLWISE','ATLAS']:
            if refname == 'TMASS': 
                tablename = 'twomass.psc' 
                cols = 'designation,ra as raj2000,dec as dej2000,j_m as jmag,j_cmsig as e_jmag,h_m as hmag,h_cmsig as e_hmag,k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg' 
                ##server = 'gp04.datalab.noao.edu' 
                #server = 'gp01.datalab.noirlab.edu' 
                ##server = 'dldb1.sdm.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
            racol = 'ra' 
            deccol = 'dec' 
            if refname == 'GAIA': 
                tablename = 'gaia_dr1.gaia_source' 
                cols = 'source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,dec as de_icrs,dec_error as e_de_icrs,'
                cols += 'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,phot_g_mean_mag as gmag' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                ##server = 'dldb1.sdm.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
            if refname == 'GAIADR2': 
                tablename = 'gaia_dr2.gaia_source' 
                cols = 'source_id as source,ra,ra_error,dec,dec_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'
                cols += 'phot_g_mean_mag as gmag,phot_bp_mean_mag as bp,phot_bp_mean_flux as fbp,phot_bp_mean_flux_error as e_fbp,'
                cols += 'phot_rp_mean_mag as rp,phot_rp_mean_flux as frp,phot_rp_mean_flux_error as e_frp' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
            if refname == 'GAIAEDR3': 
                tablename = 'gaia_edr3.gaia_source' 
                cols = 'source_id as source,ra,ra_error,dec,dec_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'
                cols += 'phot_g_mean_mag as gmag,phot_bp_mean_mag as bp,phot_bp_mean_flux as fbp,phot_bp_mean_flux_error as e_fbp,'
                cols += 'phot_rp_mean_mag as rp,phot_rp_mean_flux as frp,phot_rp_mean_flux_error as e_frp' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
            if refname == 'PS': 
                #tablename = 'cp_calib.ps1' 
                tablename = 'public.ps1' 
                cols = 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag' 
                ##server = 'gp02.datalab.noirlab.edu' 
                #server = 'gp01.datalab.noirlab.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
            if refname == 'SKYMAPPER': 
                tablename = 'skymapper_dr1.master' 
                cols = 'raj2000, dej2000, u_psf as sm_umag, e_u_psf as e_sm_umag, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag,'
                cols += 'e_r_psf as e_sm_rmag, i_psf as sm_imag,e_i_psf as e_sm_imag, z_psf as sm_zmag, e_z_psf as e_sm_zmag' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
                racol = 'raj2000' 
                deccol = 'dej2000' 
            if refname == 'SKYMAPPERDR2': 
                tablename = 'skymapper_dr2.master' 
                cols = 'raj2000, dej2000, u_psf as sm_umag, e_u_psf as e_sm_umag, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag,'
                cols += 'e_r_psf as e_sm_rmag, i_psf as sm_imag,e_i_psf as e_sm_imag, z_psf as sm_zmag, e_z_psf as e_sm_zmag' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
                racol = 'raj2000' 
                deccol = 'dej2000' 
            if refname == 'ALLWISE': 
                tablename = 'allwise.source' 
                #cols = 'ra, dec, w1mdef as w1mag, w1sigmdef as e_w1mag, w2mdef as w2mag, w2sigmdef as e_w2mag' 
                cols = 'ra, dec, w1mpro as w1mag, w1sigmpro as e_w1mag, w2mpro as w2mag, w2sigmpro as e_w2mag' 
                #server = 'gp04.datalab.noao.edu' 
                #server = 'gp01.datalab.noirlab.edu' 
                server = 'db02.datalab.noirlab.edu' 
                user = 'dlquery'
            if refname == 'ATLAS': 
                tablename = 'atlasrefcat2' 
                cols = 'objid,ra,dec,plx as parallax,dplx as parallax_error,pmra,dpmra as pmra_error,pmdec,dpmdec as pmdec_error,gaia,dgaia as gaiaerr,'
                cols += 'bp,dbp as bperr,rp,drp as rperr,teff,agaia,dupvar,ag,rp1,r1,r10,g as gmag,dg as gerr,gchi,gcontrib,'
                cols += 'r as rmag, dr as rerr,rchi,rcontrib,i as imag,di as ierr,ichi,icontrib,z as zmag,dz as zerr,zchi,zcontrib,nstat,'
                cols += 'j as jmag,dj as jerr,h as hmag,dh as herr,k as kmag,dk as kerr' 
                server = 'gp10.datalab.noirlab.edu' 
                user = 'datalab'

            # Use Postgres command with q3c cone search 
            refcattemp = savefile.replace('.fits','.txt') 
            cmd = "psql -h "+server+" -U "+user+" -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename
            cmd += " WHERE q3c_radial_query(%s,%s,%.5f,%.5f,%.3f)'" % (racol,deccol,cenra,cendec,radius)
            cmd += " > "+refcattemp
            dln.remove(refcattemp,allow=True)
            dln.remove(savefile,allow=True) 
            out = subprocess.check_output(cmd,shell=True)
            # Check for empty query
            tlines = dln.readlines(refcattemp,nreadline=4)
            if len(tlines) < 4: 
                if silent==False:
                    logger.info('No Results')
                return []
            #  Load ASCII file and create the FITS file 
            ref = ascii.read(refcattemp,data_start=3,delimiter='|')
            #ref = importascii(refcattemp,/header,delim='|',skipline=2,/silent) 
            if saveref:
                logger.info('Saving catalog to file '+savefile)
                ref.write(savefile,overwrite=True)
            dln.remove(refcattemp,allow=True)
             
            # Fix 0.0 mags/errs in ATLAS 
            if refname == 'ATLAS': 
                magcols = ['gaia','bp','rp','gmag','rmag','imag','zmag','jmag','hmag','kmag'] 
                errcols = ['gaiaerr','bperr','rperr','gerr','rerr','ierr','zerr','jerr','herr','kerr'] 
                cols = ref.colnames
                # Set mags with 0.0 to 99.99 
                for j in range(len(magcols)):
                    if magcols[j] in ref.colnames:
                        bdmag = (ref[magcols[j]] <= 0.0)
                        if np.sum(bdmag)>0:
                            ref[magcols[j]][bdmag] = 99.99
                # Set errors with 0.0 to 9.99 
                for j in range(len(errcols)):
                    if errcols[j] in ref.colnames:
                        bderr = (ref[errcols[j]] <= 0.0)
                        if np.sum(bderr)>0:
                            ref[errcols[j]][bderr] = 9.99
             
        # Use astroquery vizier
        #---------------- 
        #   for low density with 2MASS/GAIA and always for GALEX and APASS 
        else: 
             
            # Use QUERYVIZIER for GALEX (python code has problems) 
            #if refname == 'II/312/ais' or refname == 'GALEX': 
            # if refcat eq 'APASS' then cfa=0 else cfa=1  ; cfa doesn't have APASS 
            #cfa = 1   # problems with CDS VizieR and cfa has APASS now 
            #if refcat == 'SAGE': 
            #    cfa = 0 
            result = Vizier.query_region(SkyCoord(ra=cenra, dec=cendec, unit='deg'),
                                         radius=Angle(radius, "deg"), catalog=refname)
            # Check for failure 
            if len(result)==0:
                if silent==False : 
                    logger.info('Failure or No Results')
                return []                
            ref = result[0]
            ref.meta['description'] = ref.meta['description'][0:50]
            #ref = QUERYVIZIER(refname,[cenra,cendec],radius*60,cfa=cfa,timeout=600,/silent) 
               
            # Fix/homogenize the GAIA tags 
            if refname == 'GAIA': 
                nref = len(ref) 
                orig = ref.copy()
                dt = [('source',int),('ra_icrs',float),('e_ra_icrs',float),('de_icrs',float),('e_de_icrs',float),
                      ('fg',float),('e_fg',float),('gmag',float)]
                ref = np.zeros(nref,dtype=np.dtype(dt))
                ref = Table(ref)
                for n in orig.colnames:
                    ref[n] = orig[n]
                ref['fg'] = orig['_fg_']
                ref['e_fg'] = orig['e__fg_']
                ref['gmag'] = orig['_gmag_']
                del orig
            # Fix/homogenize the 2MASS tags 
            elif refname == 'TMASS': 
                nref = len(ref) 
                orig = ref.copy()
                dt = [('designation',(np.str,50)),('raj2000',float),('dej2000',float),('jmag',float),('e_jmag',float),
                      ('hmag',float),('e_hmag',float),('kmag',float),('e_kmag',float),('qflg',(np.str,20))]
                ref = np.zeros(nref,dtype=np.dtype(dt))
                for n in orig.colnames:
                    ref[n] = orig[n]
                ref['designation'] = orig['_2mass']
                del orig
            # Fix NANs in ALLWISE 
            elif refname == 'ALLWISE': 
                bd, = np.where(np.isfinite(ref['_3_6_']) == False)
                if len(bd) > 0: 
                    ref['_3_6_'][bd] = 99.99 
                    ref['e__3_6_'][bd] = 9.99 
                bd, = np.where(np.isfinite(ref['_4_5_']) == False)
                if len(bd) > 0: 
                    ref['_4_5_'][bd] = 99.99 
                    ref['e__4_5_'][bd] = 9.99 
             
            # Save the file 
            if saveref:
                if silent==False: 
                    logger.info('Saving catalog to file '+savefile)
                ref.write(savefile,overwrite=True)  # only save if necessary
                          
    if silent==False:
        logger.info('%d sources found   dt=%.1f sec.' % (len(ref),time.time()-t0))
     
    return ref 


def getexttype(cenra,cendec,radius):
    """
    Get the extinction type for the NSC. 
 
    Parameters
    ----------
    cenra : float
       Right Ascension at the center of the image. 
    cendec : float
       Declination at the center of the image. 
    radius : float
       Radius of the image. 
 
    Returns
    -------
    exttype : int
       Extinction type: 
            1 - SFD 
            2 - RJCE ALLWISE 
            3 - RJCE GlIMPSE 
            4 - RJCE SAGE 
 
    Example
    -------

    ext_type = getexttype(cenra,cendec,radius) 
    
    By D. Nidever  Feb 2019 
    Translated to Python by D. Nidever, April 2022
    """
     
    # Figure out what reddening method we are using 
    #---------------------------------------------- 
    # Extinction types: 
    # 1 - SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2 
    # 2 - RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE 
    #     data available 
    # 3 - RJCE GLIMPSE, GLIMPSE data available 
    # 4 - RJCE SAGE, SAGE data available 
    ext_type = 0 
    cencoo = SkyCoord(ra=cenra,dec=cendec,unit='deg')
    cengl = cencoo.galactic.l.deg
    cengb = cencoo.galactic.b.deg
                          
    # Get grid of SFD EBV values across the field to see 
    #  what regime we are working in 
    x = np.linspace(-radius,radius,100).reshape(-1,1) + np.zeros(100,float).reshape(1,-1)
    y = np.zeros(100,float).reshape(-1,1) + np.linspace(-radius,radius,100)
    rr,dd = coords.rotsphcen(x,y,cenra,cendec,gnomic=True,reverse=True)
    coo = SkyCoord(ra=rr,dec=dd,unit='deg')
    sfd = SFDQuery()
    ebv_grid = sfd(coo)
    maxebv = np.max(ebv_grid) 

    # Check if there is any GLIMPSE data available 
    #  do 3x3 grid with small matching radius 
    if np.abs(cengb) < 5 and (cengl < 65 or cengl > 290): 
        x = np.linspace(-radius,radius,3).reshape(-1,1) +np.zeros(3,float).reshape(1,-1)
        y = np.zeros(3,float).reshape(-1,1) + np.linspace(-radius,radius,3)
        rr,dd = coords.rotsphcen(x,y,cenra,cendec,gnomic=True,reverse=True)
        rr = rr.flatten()
        dd = dd.flatten()
        nresult = 0 
        cnt = 0
        while (nresult == 0) and (cnt < 9): 
            result = Vizier.query_region(SkyCoord(ra=rr[cnt],dec=dd[cnt],unit='deg'),
                                         radius=Angle(1.0,"arcmin"),catalog='II/293/glimpse')
            nresult = len(result)
            cnt += 1 
        if nresult > 0: 
            ext_type = 3 
     
    # Check if there is any SAGE data available 
    #  do 3x3 grid with small matching radius 
    lmcrad = coords.sphdist(81.9,-69.867,cenra,cendec) 
    smcrad = coords.sphdist(13.183,-72.8283,cenra,cendec) 
    if lmcrad < 5.0 or smcrad < 4.0: 
        x = np.linspace(-radius,radius,3).reshape(-1,1) +np.zeros(3,float).reshape(1,-1)
        y = np.zeros(3,float).reshape(-1,1) + np.linspace(-radius,radius,3)
        rr,dd = coords.rotsphcen(x,y,cenra,cendec,gnomic=True,reverse=True)
        rr = rr.flatten() 
        dd = dd.flatten() 
        nresult = 0 
        cnt = 0
        while (nresult == 0) and (cnt < 9): 
            result = Vizier.query_region(SkyCoord(ra=rr[cnt],dec=dd[cnt],unit='deg'),
                                         radius=Angle(1.0,"arcmin"),catalog='II/305/archive')
            nresult = len(result)
            cnt += 1 
        if nresult > 0: 
            ext_type = 4 
     
    # Use RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE 
    #     data available 
    if ext_type == 0 and (np.abs(cengb) < 16 or maxebv > 0.2) : 
        ext_type = 2 
     
    # SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2 
    if ext_type == 0: 
        ext_type = 1 
     
    return ext_type 

                          
def getrefdata(filt,cenra,cendec,radius,saveref=False,silent=False,
               dcr=0.5,modelmags=False,logger=None):
    """
    Get reference catalog information needed for a given filter. 
 
    Parameters
    ----------
    filt : str
       Filter and instrument name, e.g. "c4d-u". 
         This can be an array of multiple filters and 
         then all the reference needed for all the filters 
         will be included. 
    cenra : float
       Central RA for the search. 
    cendec : float
       Central DEC for the search. 
    radius : float
       Search radius in degrees. 
    dcr : float, optional
       The cross-matching radius in arcsec.  Default is 0.5". 
    saveref : bool, optional
       Save the output to FILE. 
    silent : bool, optional
       Don't print anything to the screen. 
    modelmags : bool, optional
       Return the model magnitudes as well. 
    logger : logging object
       Logging object to use for printing ot the screen.
 
    Returns
    -------
    ref : astropy table
       Search results from the reference catalog all in one 
         catalog
 
    Example
    -------

    cat = getrefdata('g',cenra,cendec,radius,saveref=saveref) 
 
    By D. Nidever  Sep 2017 
    Translated by D. Nidever, April 2022
    """
                          
    t0 = time.time() 

    if dln.size(filt)==1 and type(filt) is str:
        filt = [filt]
     
    # Check that we have psql installed 
    out = subprocess.check_output(['which','psql'],shell=False)
    if type(out) is bytes:
        out = out.decode()
    out = out.strip()
    if dln.size(out)>1:
        out = out[0]
    if os.path.exists(out) == 0: 
        raise ValueError('No PSQL found on this sytem.')

    if logger is None:
        logger = dln.basiclogger()
     
    # Figure out what reddening method we are using 
    #---------------------------------------------- 
    # Extinction types: 
    # 1 - SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2 
    # 2 - RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE 
    #     data available 
    # 3 - RJCE GLIMPSE, GLIMPSE data available 
    # 4 - RJCE SAGE, SAGE data available 
    ext_type = getexttype(cenra,cendec,radius) 
     
    # If there is GLIMPSE or SAGE data, also get ALLWISE data 
    #  in case there is only partial coverage 
     
     
    if silent==False: 
        logger.info('Getting reference catalogs for:')
        logger.info('FILTER(S) = '+', '.join(filt))
        logger.info('CENRA  = '+str(cenra)) 
        logger.info('CENDEC = '+str(cendec)) 
        logger.info('RADIUS = '+str(radius)+' deg')
        if ext_type==1:
            logger.info('Extinction Type: 1 - SFD')
        elif ext_type==2:
            logger.info('Extinction Type: 2 - RJCE ALLWISE')
        elif ext_type==3:
            logger.info('Extinction Type: 3 - RJCE GLIMPSE')
        elif ext_type==4:
            logger.info('Extinction Type: 4 - RJCE SAGE')
 
     
    # Load the reference catalogs 
    #----------------------------- 
    # Figure out the reference catalogs that we need based on 
    #  filter-instrument combination
    refcat = []
    for i in range(len(filt)): 
        instfilt = filt[i] 
        # If no instrument information then assume it's DECam 
        if instfilt.find('-') == -1: 
            instfilt = 'c4d-'+instfilt 
        # Check the cases
        # DECam u-band                           
        if instfilt=='c4d-u':
            # Use GAIA, 2MASS and GALEX to calibrate 
            refcat += ['2MASS-PSC','II/312/ais'] 
            if cendec <= 0: 
                refcat += ['Skymapper']
        # DECam g-band
        elif instfilt=='c4d-g':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                # Use 2MASS and Skymapper to calibrate 
                #push,refcat,['2MASS-PSC','Skymapper'] 
                refcat += ['2MASS-PSC','ATLAS'] 
        # DECam r-band
        elif instfilt=='c4d-r':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                # Use 2MASS and Skymapper to calibrate 
                #push,refcat,['2MASS-PSC','Skymapper'] 
                refcat += ['2MASS-PSC','ATLAS'] 
        # DECam i-band
        elif instfilt=='c4d-i':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                # Use Skymapper and 2MASS to calibrate 
                #push,refcat,['2MASS-PSC','Skymapper'] 
                refcat += ['2MASS-PSC','ATLAS'] 
        # DECam z-band 
        elif instfilt=='c4d-z':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                # Use Skymapper and 2MASS to calibrate 
                #push,refcat,['2MASS-PSC','Skymapper'] 
                refcat += ['2MASS-PSC','ATLAS'] 
        # DECam Y-band 
        elif instfilt=='c4d-Y':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                # Use 2MASS to calibrate 
                refcat += ['2MASS-PSC'] 
        # DECam VR-band 
        elif instfilt=='c4d-VR':
            # Use PS1 if possible 
            if cendec > -29: 
                refcat += ['2MASS-PSC','PS'] 
            else: 
                refcat += ['2MASS-PSC','ATLAS'] 
        # Bok+90Prime g-band 
        elif instfilt=='ksb-g':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Bok+90Prime r-band 
        elif instfilt=='ksb-r':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Mosaic3 g-band 
        elif instfilt=='k4m-g':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Mosaic3 r-band 
        elif instfilt=='k4m-r':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Mosaic3 i-band 
        elif instfilt=='k4m-i':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Mosaic3 z-band 
        elif instfilt=='k4m-z':
            # Use PS1 
            refcat += ['2MASS-PSC','PS'] 
        # Mosaic3 VR-band 
        elif instfilt=='k4m-VR':
            # Use GAIA G-band to calibrate 
            refcat += ['2MASS-PSC','PS'] 
            #push,refcat,['2MASS-PSC'] 
        else:
            logger.info(filt+' not currently supported')

 
    # Extinction catalogs 
    if ext_type >= 2: 
        refcat += ['ALLWISE']
    elif ext_type == 3: 
        refcat += ['GLIMPSE']
    elif ext_type == 4: 
        refcat += ['SAGE']
    # If near DEC=-30 then load BOTH PS and ATLAS 
    if (cendec > -34 and cendec < -26) and (('PS' in refcat) or ('ATLAS' in refcat)):
        refcat += ['PS','ATLAS']  # code below removes duplicates 
    # Some catalogs 
    if len(refcat) > 0: 
        # Get the unique ones
        refcat,ui = np.unique(refcat,return_index=True)
        # Always add Gaia at the beginning 
        refcat = ['GAIAEDR3']+list(refcat)
    else:
        refcat = ['GAIAEDR3']
    nrefcat = len(refcat) 
 
    # Figure out the new columns that need to be added 
    newcols = []
    for i in range(nrefcat):
        if refcat[i]=='GAIADR2' or refcat[i]=='GAIAEDR3':
            newcols +=['source','ra','ra_error','dec','dec_error','pmra','pmra_error','pmdec','pmdec_error','gmag','e_gmag','bp','e_bp','rp','e_rp'] 
        elif refcat[i]=='2MASS-PSC':
            newcols += ['jmag','e_jmag','hmag','e_hmag','kmag','e_kmag','qflg'] 
        elif refcat[i]=='PS':
            newcols += ['ps_gmag','ps_rmag','ps_imag','ps_zmag','ps_ymag'] 
        elif refcat[i]=='APASS':
            newcols += ['apass_gmag','e_apass_gmag','apass_rmag','e_apass_rmag'] 
        elif refcat[i]=='II/312/ais':
            newcols += ['nuv','e_nuv']  # Galex 
        elif refcat[i]=='Skymapper':
            newcols += ['sm_umag','e_sm_umag','sm_gmag','e_sm_gmag','sm_rmag','e_sm_rmag','sm_imag','e_sm_imag','sm_zmag','e_sm_zmag']  # Skymapper DR1 
        elif refcat[i]=='ALLWISE':
            newcols += ['w1mag','e_w1mag','w2mag','e_w2mag'] 
        elif refcat[i]=='GLIMPSE':
            newcols += ['gl_36mag','e_gl_36mag','gl_45mag','e_gl_45mag'] 
        elif refcat[i]=='SAGE':
            newcols += ['sage_36mag','e_sage_36mag','sage_45mag','e_sage_45mag'] 
        elif refcat[i]=='ATLAS':
            newcols += ['atlas_gmag','e_atlas_gmag','atlas_gcontrib','atlas_rmag','e_atlas_rmag','atlas_rcontrib',
                        'atlas_imag','e_atlas_imag','atlas_icontrib','atlas_zmag','e_atlas_zmag','atlas_zcontrib'] 
        else:
            raise ValueError(refcat[i]+' NOT SUPPORTED')
    newcols += ['ebv_sfd','ejk','e_ejk','ext_type'] 
    if modelmags:
        newcols += ['model_mag']
    nnewcols = len(newcols) 
 
    # Load the necessary catalogs 
    nrefcat = len(refcat) 
    if silent==False: 
        logger.info(str(nrefcat)+' reference catalogs to load: '+', '.join(refcat))
    ref = None
    for i in range(nrefcat): 
        t0 = time.time() 
        if silent==False: 
                logger.info('Loading '+refcat[i]+' reference catalog')
        # Load the catalog 
        ref1 = getrefcat(cenra,cendec,radius,refcat[i],silent=silent,logger=logger) 
        nref1 = len(ref1)
        if nref1 == 0: 
            continue
 
        # Initialize the catalog 
        dt = []
        for j in range(nnewcols): 
            dtype1 = float
            if newcols[j] == 'qflg': 
                dtype1 = (np.str,20)
            if newcols[j] == 'ext_type': 
                dtype1 = int
            if newcols[j] == 'source': 
                dtype1 = int
            if newcols[j] == 'ra': 
                dtype1 = float 
            if newcols[j] == 'dec': 
                dtype1 = float
            if newcols[j].find('contrib') > -1:
                dtype1 = int
            dt += [(newcols[j],dtype1)]
 
        # First successful one, initialize the catalog 
        if ref is None:
            ref = np.zeros(nref1,dtype=np.dtype(dt))
            ref = Table(ref)
            for n in ref.colnames:
                if n in ref1.colnames:
                    ref[n] = ref1[n]
            #ref = replicate(schema,nref1) 
            #struct_assign,ref1,ref,/nozero
            if 'ra' not in ref1.colnames:
                ref['ra'] = ref['ra_icrs']
                ref['dec'] = ref['de_icrs']
            ind1 = np.arange(nref1).astype(int) 
            ind2 = np.arange(nref1).astype(int)
            nmatch = nref1 
 
        # Second and later 
        else: 
 
            # Get RA/DEC columns 
            # 2MASS, Galex, APASS use RAJ2000 
            # PS uses RA/DEC 
            if (refcat[i] != 'PS' and refcat[i] != 'ALLWISE' and refcat[i] != 'ATLAS'):
                racol = 'raj2000'
            else:
                racol = 'ra'
            if (refcat[i] != 'PS' and refcat[i] != 'ALLWISE' and refcat[i] != 'ATLAS'):
                deccol = 'dej2000'
            else:
                deccol = 'dec'
                        
            # Crossmatch 
            ind1,ind2,dist = coords.xmatch(ref['ra'],ref['dec'],ref1[racol],ref1[deccol],dcr)
            if silent==False: 
                logger.info(str(nmatch)+' matches')
 
        # Add magnitude columns 
        if nmatch > 0:
            if refcat[i]=='GAIADR2' or refcat[i]=='GAIAEDR3':
                temp = ref[ind1]
                for n in temp.colnames:
                    if n in ref1.colnames:
                        temp[n] = ref1[n][ind2]
                #struct_assign,ref1[ind2],temp,/nozero 
                temp['e_gmag'] = 2.5*np.log10(1.0+ref1['e_fg'][ind2]/ref1['fg'][ind2]) 
                temp['e_bp'] = 2.5*np.log10(1.0+ref1['e_fbp'][ind2]/ref1['fbp'][ind2]) 
                temp['e_rp'] = 2.5*np.log10(1.0+ref1['e_frp'][ind2]/ref1['frp'][ind2]) 
                ref[ind1] = temp 
            elif refcat[i]=='2MASS-PSC':
                ref['jmag'][ind1] = ref1['jmag'][ind2]
                ref['e_jmag'][ind1] = ref1['e_jmag'][ind2]
                ref['hmag'][ind1] = ref1['hmag'][ind2]
                ref['e_hmag'][ind1] = ref1['e_hmag'][ind2]
                ref['kmag'][ind1] = ref1['kmag'][ind2]
                ref['e_kmag'][ind1] = ref1['e_kmag'][ind2]
                ref['qflg'][ind1] = ref1['qflg'][ind2]
            elif refcat[i]=='PS':
                ref['ps_gmag'][ind1] = ref1['gmag'][ind2]
                ref['ps_rmag'][ind1] = ref1['rmag'][ind2]
                ref['ps_imag'][ind1] = ref1['imag'][ind2]
                ref['ps_zmag'][ind1] = ref1['zmag'][ind2]
                ref['ps_ymag'][ind1] = ref1['ymag'][ind2]
            elif refcat[i]=='APASS':
                ref['apass_gmag'][ind1] = ref1['g_mag'][ind2]
                ref['e_apass_gmag'][ind1] = ref1['e_g_mag'][ind2]
                ref['apass_rmag'][ind1] = ref1['r_mag'][ind2]
                ref['e_apass_rmag'][ind1] = ref1['e_r_mag'][ind2]
            elif refcat[i]=='II/312/ais':
                if tag_exist(ref1,'NUV'): 
                    ref['nuv'][ind1] = ref1['nuv'][ind2]
                    ref['e_nuv'][ind1] = ref1['e_nuv'][ind2]
                else: 
                    ref['nuv'][ind1] = ref1['nuvmag'][ind2]
                    ref['e_nuv'][ind1] = ref1['e_nuvmag'][ind2]
            elif refcat[i]=='Skymapper':
                ref['sm_umag'][ind1] = ref1['sm_umag'][ind2]
                ref['e_sm_umag'][ind1] = ref1['e_sm_umag'][ind2]
                ref['sm_gmag'][ind1] = ref1['sm_gmag'][ind2]
                ref['e_sm_gmag'][ind1] = ref1['e_sm_gmag'][ind2]
                ref['sm_rmag'][ind1] = ref1['sm_rmag'][ind2]
                ref['e_sm_rmag'][ind1] = ref1['e_sm_rmag'][ind2]
                ref['sm_imag'][ind1] = ref1['sm_imag'][ind2]
                ref['e_sm_imag'][ind1] = ref1['e_sm_imag'][ind2]
                ref['sm_zmag'][ind1] = ref1['sm_zmag'][ind2]
                ref['e_sm_zmag'][ind1] = ref1['e_sm_zmag'][ind2]
            elif refcat[i]=='ATLAS':
                ref['atlas_gmag'][ind1] = ref1['gmag'][ind2]
                ref['e_atlas_gmag'][ind1] = ref1['gerr'][ind2]
                ref['atlas_gcontrib'][ind1] = ref1['gcontrib'][ind2]
                ref['atlas_rmag'][ind1] = ref1['rmag'][ind2]
                ref['e_atlas_rmag'][ind1] = ref1['rerr'][ind2]
                ref['atlas_rcontrib'][ind1] = ref1['rcontrib'][ind2]
                ref['atlas_imag'][ind1] = ref1['imag'][ind2]
                ref['e_atlas_imag'][ind1] = ref1['ierr'][ind2]
                ref['atlas_icontrib'][ind1] = ref1['icontrib'][ind2]
                ref['atlas_zmag'][ind1] = ref1['zmag'][ind2]
                ref['e_atlas_zmag'][ind1] = ref1['zerr'][ind2]
                ref['atlas_zcontrib'][ind1] = ref1['zcontrib'][ind2]
            elif refcat[i]=='ALLWISE':
                ref['w1mag'][ind1] = ref1['w1mag'][ind2]
                ref['e_w1mag'][ind1] = ref1['e_w1mag'][ind2]
                ref['w2mag'][ind1] = ref1['w2mag'][ind2]
                ref['e_w2mag'][ind1] = ref1['e_w2mag'][ind2]
            elif refcat[i]=='GLIMPSE':
                ref['gl_36mag'][ind1] = ref1['_3_6mag'][ind2]
                ref['e_gl_36mag'][ind1] = ref1['e_3_6mag'][ind2]
                ref['gl_45mag'][ind1] = ref1['_4_5mag'][ind2]
                ref['e_gl_45mag'][ind1] = ref1['e_4_5mag'][ind2]
            elif refcat[i]=='SAGE':
                ref['sage_36mag'][ind1] = ref1['__3_6_'][ind2]
                ref['e_sage_36mag'][ind1] = ref1['e__3_6_'][ind2]
                ref['sage_45mag'][ind1] = ref1['__4_5_'][ind2]
                ref['e_sage_45mag'][ind1] = ref1['e__4_5_'][ind2]
            else:
                raise ValueError(catname+' NOT SUPPORTED')
 
 
        # Add leftover ones 
        if nmatch < len(ref1): 
            left1 = ref1.copy()
            if nmatch > 0:
                left1 = np.delete(left1,ind2)
                #remove,ind2,left1 
            nleft1 = len(left1)
            new = np.zeros(nleft1,dtype=np.dtype(dt))
            new = Table(new)
            new['ra'] = left1[racol]
            new['dec'] = left1[deccol]
            if refcat[i]=='GAIADR2' or refcat[i]=='GAIAEDR3':
                temp = ref[ind1]
                for n in left1.colnames:
                     new[n] = left1[n]
                new['e_gmag'] = 2.5*np.log10(1.0+left1['e_fg']/left1['fg']) 
                new['e_bp'] = 2.5*np.log10(1.0+left1['e_fbp']/left1['fbp']) 
                new['e_rp'] = 2.5*np.log10(1.0+left1['e_frp']/left1['frp']) 
            elif refcat[i]=='2MASS-PSC':
                new['jmag'] = left1['jmag']
                new['e_jmag'] = left1['e_jmag']
                new['hmag'] = left1['hmag']
                new['e_hmag'] = left1['e_hmag']
                new['kmag'] = left1['kmag']
                new['e_kmag'] = left1['e_kmag'] 
                new['qflg'] = left1['qflg']
            elif refcat[i]=='PS':
                new['ps_gmag'] = left1['gmag'] 
                new['ps_rmag'] = left1['rmag']
                new['ps_imag'] = left1['imag']
                new['ps_zmag'] = left1['zmag']
                new['ps_ymag'] = left1['ymag']
            elif refcat[i]=='APASS':
                new['apass_gmag'] = left1['g_mag'] 
                new['e_apass_gmag'] = left1['e_g_mag'] 
                new['apass_rmag'] = left1['r_mag']
                new['e_apass_rmag'] = left1['e_r_mag']
            elif refcat[i]=='II/312/ais':
                if 'NUV' in left1.colnames:
                    new['nuv'] = left1['nuv']
                    new['e_nuv'] = left1['e_nuv'] 
                else: 
                    new['nuv'] = left1['nuvmag'] 
                    new['e_nuv'] = left1['e_nuvmag'] 
            elif refcat[i]=='Skymapper':
                new['sm_umag'] = left1['sm_umag'] 
                new['e_sm_umag'] = left1['e_sm_umag'] 
                new['sm_gmag'] = left1['sm_gmag']
                new['e_sm_gmag'] = left1['e_sm_gmag'] 
                new['sm_rmag'] = left1['sm_rmag']
                new['e_sm_rmag'] = left1['e_sm_rmag'] 
                new['sm_imag'] = left1['sm_imag']
                new['e_sm_imag'] = left1['e_sm_imag'] 
                new['sm_zmag'] = left1['sm_zmag']
                new['e_sm_zmag'] = left1['e_sm_zmag'] 
            elif refcat[i]=='ATLAS':
                new['atlas_gmag'] = left1['gmag']
                new['e_atlas_gmag'] = left1['gerr'] 
                new['atlas_gcontrib'] = left1['gcontrib'] 
                new['atlas_rmag'] = left1['rmag']
                new['e_atlas_rmag'] = left1['rerr'] 
                new['atlas_rcontrib'] = left1['rcontrib'] 
                new['atlas_imag'] = left1['imag']
                new['e_atlas_imag'] = left1['ierr'] 
                new['atlas_icontrib'] = left1['icontrib'] 
                new['atlas_zmag'] = left1['zmag']
                new['e_atlas_zmag'] = left1['zerr'] 
                new['atlas_zcontrib'] = left1['zcontrib'] 
            elif refcat[i]=='ALLWISE':
                new['w1mag'] = left1['w1mag']
                new['e_w1mag'] = left1['e_w1mag'] 
                new['w2mag'] = left1['w2mag']
                new['e_w2mag'] = left1['e_w2mag'] 
            elif refcat[i]=='GLIMPSE':
                new['gl_36mag'] = left1['_3_6mag']
                new['e_gl_36mag'] = left1['e_3_6mag'] 
                new['gl_45mag'] = left1['_4_5mag']
                new['e_gl_45mag'] = left1['e_4_5mag'] 
            elif refcat[i]=='SAGE':
                new['sage_36mag'] = left1['__3_6_']
                new['e_sage_36mag'] = left1['e__3_6_'] 
                new['sage_45mag'] = left1['__4_5_']
                new['e_sage_45mag'] = left1['e__4_5_']
            else:
                raise ValueError(catname+' NOT SUPPORTED')
 
 
            # Combine the two 
            old = ref.copy()
            ref = np.zeros(len(old)+nleft1,dtype=np.dtype(dt))
            ref = Table(ref)
            ref[0:len(old)] = old 
            ref[len(old):] = new
            del old
            del new
            del left1
  
    # Get extinction 
    #---------------- 
    ref = getreddening(ref,ext_type)
 
    count = len(ref) 
 
    logger.info('dt=%.1f sec' % (time.time()-t0))

    return ref 
 

def getreddening(ref,ext_type):
    """
    This calculates E(J-Ks) reddening using reference catalog data for 
    the NSC. 
 
    Parameters
    ----------
    ref : catalog
       The catalog with the reference catalog data. 
    exttype : int
       Extinction type: 
             1 - SFD 
             2 - RJCE ALLWISE 
             3 - RJCE GlIMPSE 
             4 - RJCE SAGE 
 
    Returns
    -------
    The REF catalog EJK column and EXT_TYPE columns are updated. 
 
    Example
    -------

    ref = getreddening(ref,ext_type)
 
    By D. Nidever  Feb 2019 
    Translated to Python by D. Nidever, April 2022
    """
     
    # Add SFD reddening
    coo = SkyCoord(ra=ref['ra'],dec=ref['dec'],unit='deg')
    sfd = SFDQuery()
    ebv = sfd(coo)
    ref['ebv_sfd'] = ebv 
     
    # Start with SFD extinction for all 
    ejk_sfd = 1.5*0.302*ref['ebv_sfd']
    ref['ejk'] = ejk_sfd 
    ref['e_ejk'] = 0.1 # not sure what to put here 
    bd, = np.where(ref['ebv_sfd'] > 0.3)   # E(B-V)>0.3 is roughly where things "break down" 
    if len(bd) > 0 : 
        ref['e_ejk'][bd] = 1.0 
    ref['ext_type'] = 1 
     
    # RJCE extinction 
    if ext_type >= 2: 
         
        # RJCE GLIMPSE, type=3 
        if ext_type == 3: 
            gdglimpse, = np.where((ref['jmag'] < 50) & (ref['hmag'] < 50) & (ref['kmag'] < 50) & (ref['gl_45mag'] < 50))
            if len(gdglimpse) > 0: 
                ejk = 1.5*0.918*(ref['hmag'][gdglimpse]-ref['gl_45mag'][gdglimpse]-0.08) 
                e_ejk = 1.5*0.918*np.sqrt(ref['e_hmag'][gdglimpse]**2+ref['e_gl_45mag'][gdglimpse]**2) 
                #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdglimpse] and e_ejk lt 0.2,ngdejk) 
                gdejk, = np.where(ejk < ejk_sfd[gdglimpse],ngdejk) 
                if len(gdejk) > 0: 
                    ref['ejk'][gdglimpse[gdejk]] = np.maximum(ejk[gdejk],0)
                    ref['e_ejk'][gdglimpse[gdejk]] = e_ejk[gdejk] 
                    ref['ext_type'][gdglimpse[gdejk]] = 3 
         
        # RJCE SAGE, type=4 
        if ext_type == 4: 
            gdsage, = np.where((ref['jmag'] < 50) & (ref['hmag'] < 50) & (ref['kmag'] < 50) & (ref['sage_45mag'] < 50))
            if len(gdsage) > 0: 
                ejk = 1.5*0.918*(ref['hmag'][gdsage]-ref['sage_45mag'][gdsage]-0.08) 
                e_ejk = 1.5*0.918*np.sqrt(ref['e_hmag'][gdsage]**2+ref['e_sage_45mag'][gdsage]**2) 
                #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdsage] and e_ejk lt 0.2,ngdejk) 
                gdejk, = np.where(ejk < ejk_sfd[gdsage]) 
                if len(gdejk) > 0: 
                    ref['ejk'][gdsage[gdejk]] = np.maximum(ejk[gdejk],0)
                    ref['e_ejk'][gdsage[gdejk]] = e_ejk[gdejk] 
                    ref['ext_type'][gdsage[gdejk]] = 4 
         
        # RJCE ALLWISE, type=2 
        #   Use ALLWISE for those that don't have IRAC data 
        gdwise, = np.where((ref['jmag'] < 50) & (ref['hmag'] < 50) & (ref['kmag'] < 50) &
                           (ref['w2mag'] < 50) & (ref['ext_type'] <= 1))
        if len(gdwise) > 0: 
            ejk = 1.5*0.918*(ref['hmag'][gdwise]-ref['w2mag'][gdwise]-0.05) 
            e_ejk = 1.5*0.918*np.sqrt(ref['e_hmag'][gdwise]**2+ref['e_w2mag'][gdwise]**2) 
            #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdwise] and e_ejk lt 0.2,ngdejk) 
            gdejk, = np.where(ejk < ejk_sfd[gdwise])
            if len(gdejk) > 0: 
                ref['ejk'][gdwise[gdejk]] = np.maximum(ejk[gdejk],0)
                ref['e_ejk'][gdwise[gdejk]] = e_ejk[gdejk] 
                ref['ext_type'][gdwise[gdejk]] = 2 
     
    # No good reddening 
    bd, = np.where((ref['ext_type'] == 1) & (ref['ebv_sfd'] > 0.3))   # E(B-V)>0.3 is roughly where things "break down" 
    if len(bd) > 0: 
        ref['ejk'][bd] = 999999.0 
        ref['e_ejk'][bd] = 999999.0 
        ref['ext_type'][bd] = 0 
     
    # Fix NANs in E_EJK 
    bd, = np.where(np.isfinite(ref['e_ejk']) == False)
    if len(bd) > 0 : 
        ref['e_ejk'][bd] = 9.99 

    return ref
 
 
