#!/usr/bin/env python

import os
import time
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
import subprocess
from dlnpyutils import utils as dln,coords
from dustmaps.sfd import SFDQuery
from astroquery.vizier import Vizier
from astropy.coordinates import Angle,SkyCoord
import astropy.units as u
from . import utils

Vizier.TIMEOUT = 600
Vizier.ROW_LIMIT = -1
Vizier.cache_location = None


def getdata(refcat,minra,redo=False,silent=False,logger=None):
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

    cat = getrefcat(cenra,cendec,radius,refcat,file=file,saveref=saveref) 
 
    By D. Nidever  Sep 2017 
    Translated to Python by D. Nidever, April 2022
    """
     
    count = 0 
    t0 = time.time() 

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
        catname = 'II/293/glimpse' 
    elif refname == 'SAGE': 
        catname = 'II/305/archive' 
    elif refname == 'ATLASREFCAT2': 
        refname = 'ATLAS' 


    ra0 = float(minra)
    ra1 = minra+1.0


    outdir = '/net/dl1/users/dnidever/nsc/refcatalogs/'+refname+'/'
    if os.path.exists(outdir)==False:
        os.makedirs(outdir)
    savefile = outdir+'ref_%.6f_%.6f_%s.fits' % (ra0,ra1,refname)
    if os.path.exists(os.path.abspath(os.path.dirname(savefile)))==False:
        os.makedirs(os.path.abspath(os.path.dirname(savefile)))

    if silent==False:
        logger.info('Querying %s: %.6f <= RA < %.6f' % (refname,ra0,ra1))

    
     
    # Loading previously loaded file 
    if os.path.exists(savefile) and redo==False:
        logger.info(savefile+' already exists and redo==False')
        return None
         
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
                cols = 'raj2000 as ra, dej2000 as dec, u_psf as sm_umag, e_u_psf as e_sm_umag, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag,'
                cols += 'e_r_psf as e_sm_rmag, i_psf as sm_imag,e_i_psf as e_sm_imag, z_psf as sm_zmag, e_z_psf as e_sm_zmag' 
                #server = 'gp04.datalab.noirlab.edu' 
                ##server = 'gp01.datalab.noao.edu' 
                server = 'db02.datalab.noirlab.edu'
                user = 'dlquery'
                racol = 'raj2000' 
                deccol = 'dej2000' 
            if refname == 'SKYMAPPERDR2': 
                tablename = 'skymapper_dr2.master' 
                cols = 'raj2000 as ra, dej2000 as dec, u_psf as sm_umag, e_u_psf as e_sm_umag, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag,'
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
            if refname=='SKYMAPPER' or refname=='SKYMAPPERDR2':
                cmd += " WHERE raj2000 >= %.6f and raj2000 < %.6f'" % (ra0,ra1)
            else:
                cmd += " WHERE ra >= %.6f and ra < %.6f'" % (ra0,ra1)
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

            if refname.upper()=='GALEX':
                cols = ['RAJ2000','DEJ2000','FUVmag','e_FUVmag','NUVmag','e_NUVmag']
                catname = 'II/335/galex_ais'
            elif refname.upper()=='GLIMPSE':
                # Only includes GLIMPSE I,II,3D
                cols = ['RAJ2000','DEJ2000','_2MASS','_3.6mag','e_3.6mag','_4.5mag','e_4.5mag']
                catname = 'II/293/glimpse'
            elif refname.upper()=='SAGE':
                cols = ['RAJ2000','DEJ2000','__3.6_','e__3.6_','__4.5_','e__4.5_']
                catname = 'II/305/catalog'

            Vizier.ROW_LIMIT = -1
            Vizier.TIMEOUT = 1000000
            Vizier.columns = cols
            result = Vizier.query_constraints(catalog=catname,
                                              RA='>='+str(ra0)+' & <'+str(ra1))

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
             
        # Convert all mags and errmags to float32
        for n in ref.colnames:
            if n.find('mag') > -1:
                ref[n] = ref[n].astype(np.float32)
            if n.find('e_') > -1 and n.find('mag') > -1:
                ref[n] = ref[n].astype(np.float32)
        # Lowercase column names
        for n in ref.colnames:
            ref[n].name = n.lower()
        # Convert raj2000/dej2000 to ra/dec
        if 'raj2000' in ref.colnames:
            ref['raj2000'].name = 'ra'
        if 'dej2000' in ref.colnames:
            ref['dej2000'].name = 'dec'  

        # Save the file 
        logger.info('Saving catalog to file '+savefile)
        ref.write(savefile,overwrite=True)
                
    if silent==False:
        logger.info('%d sources found   dt=%.1f sec.' % (len(ref),time.time()-t0))
     
    return ref 

def download(refcat,redo=False):
    """ Download data."""

    for i in range(360):
        dum = getdata(refcat,i,redo=redo)
