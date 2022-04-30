#!/usr/bin/env python

import os
import time
import numpy as np
 
def getrefcat(cenra,cendec,radius,refcat,count=count,file=file,saveref=saveref,silent=silent,logfile=logfile):
    """
    Get reference catalog information from DL database 
 
    INPUTS: 
    cenra     Central RA for the search. 
    cendec    Central DEC for the search. 
    radius    Search radius in degrees. 
    refcat    Reference catalog name (e.g. 2MASS, Gaia, etc.). 
    =file     The file to save to or search for existing catalog. 
    /saveref  Save the output to FILE. 
    /silent   Don't print anything to the screen. 
    =logfile  Filename for logging output. 
 
    OUTPUTS: 
    ref       Search results from the reference catalog. 
    =count    Number of elements in REF. 
 
    USAGE: 
    IDL>cat = getrefcat(cenra,cendec,radius,refcat,file=file,saveref=saveref) 
 
    By D. Nidever  Sep 2017 
    """
     
    count = 0 
    t0 = time.time() 
     
    # Check that we have psql installed 
    spawn,['which','psql'],out,errout,/noshell 
    if os.path.exists(out[0]) == 0: 
        print('No PSQL found on this sytem.' 
        return -1 
     
    # Defaults 
    if len(logfile) == 0: 
        logf=-1 
    else: 
        logf=logfile 
     
    # Temporary directory 
    # /tmp is often small and get get fille dup 
    NSC_ROOTDIRS,dldir,mssdir,localdir 
    if len(version) == 0 : 
        version='v2' 
    dir = dldir+'users/dnidever/nsc/instcal/'+version+'/' 
    tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/' 
     
    # FLIP THIS AROUND, INPUT SHOULD BE THE "EASY" VERSION!!! 
    refname = strupcase(refcat) 
    if refname == 'II/312/AIS' : 
        refname='GALEX' 
    if refname == '2MASS-PSC' : 
        refname='TMASS' 
    if refname == '2MASS' : 
        refname='TMASS' 
    if refname == 'GAIA/GAIA' : 
        refname='GAIA' 
    if refname == 'Skymapper' : 
        refname='SKYMAPPER' 
    if refname == 'GLIMPSE' : 
        refname='II/293/glimpse' 
    if refname == 'SAGE' : 
        refname='II/305/archive' 
    if refname == 'ATLASREFCAT2' : 
        refname='ATLAS' 
     
    if len(file) == 0 : 
        file=tmpdir+'ref_'+stringize(cenra,ndec=5)+'_'+stringize(cendec,ndec=5)+'_'+stringize(radius,ndec=3)+'_'+refname+'.fits' 
     
    if not keyword_set(silent) : 
        printlog,logf,'Querying '+refname+': RA='+stringize(cenra,ndec=5)+' DEC='+stringize(cendec,ndec=5)+' Radius='+stringize(radius,ndec=3) 
     
    # Loading previously loaded file 
    If os.path.exists(file) == 1: 
        if not keyword_set(silent) : 
            printlog,logf,'Loading previously-saved file ',file 
        ref = MRDFITS(file,1,/silent) 
         
        # Do the Query 
        #-------------- 
    else: 
         
        # Use DataLab database search 
        #---------------------------- 
        If (refname == 'TMASS' or refname == 'GAIA' or refname == 'GAIADR2' or refname == 'PS' or refname == 'SKYMAPPER' or       refname == 'ALLWISE' or refname == 'ATLAS'): 
            if refname == 'TMASS': 
                tablename = 'twomass.psc' 
                cols = 'designation,ra as raj2000,dec as dej2000,j_m as jmag,j_cmsig as e_jmag,h_m as hmag,h_cmsig as e_hmag,k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg' 
                #server = 'gp04.datalab.noao.edu' 
                server = 'gp01.datalab.noao.edu' 
                #server = 'dldb1.sdm.noao.edu' 
            racol = 'ra' 
            deccol = 'dec' 
            if refname == 'GAIA': 
                tablename = 'gaia_dr1.gaia_source' 
                cols = 'source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,dec as de_icrs,dec_error as e_de_icrs,'+             'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,phot_g_mean_mag as gmag' 
                server = 'gp04.datalab.noao.edu' 
                #server = 'gp01.datalab.noao.edu' 
                #server = 'dldb1.sdm.noao.edu' 
            if refname == 'GAIADR2': 
                tablename = 'gaia_dr2.gaia_source' 
                cols = 'source_id as source,ra,ra_error,dec,dec_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'+             'phot_g_mean_mag as gmag,phot_bp_mean_mag as bp,phot_bp_mean_flux as fbp,phot_bp_mean_flux_error as e_fbp,'+             'phot_rp_mean_mag as rp,phot_rp_mean_flux as frp,phot_rp_mean_flux_error as e_frp' 
                server = 'gp04.datalab.noao.edu' 
                #server = 'gp01.datalab.noao.edu' 
            if refname == 'PS': 
                #tablename = 'cp_calib.ps1' 
                tablename = 'public.ps1' 
                cols = 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag' 
                server = 'gp02.datalab.noao.edu' 
            if refname == 'SKYMAPPER': 
                tablename = 'skymapper_dr1.master' 
                cols = 'raj2000, dej2000, u_psf as sm_umag, e_u_psf as e_sm_umag, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag, e_r_psf as e_sm_rmag, i_psf as sm_imag, '+             'e_i_psf as e_sm_imag, z_psf as sm_zmag, e_z_psf as e_sm_zmag' 
                server = 'gp04.datalab.noao.edu' 
                #server = 'gp01.datalab.noao.edu' 
                racol = 'raj2000' 
                deccol = 'dej2000' 
            if refname == 'ALLWISE': 
                tablename = 'allwise.source' 
                cols = 'ra, dec, w1mdef as w1mag, w1sigmdef as e_w1mag, w2mdef as w2mag, w2sigmdef as e_w2mag' 
                #server = 'gp04.datalab.noao.edu' 
                server = 'gp01.datalab.noao.edu' 
            if refname == 'ATLAS': 
                tablename = 'atlasrefcat2' 
                cols = 'objid,ra,dec,plx as parallax,dplx as parallax_error,pmra,dpmra as pmra_error,pmdec,dpmdec as pmdec_error,gaia,dgaia as gaiaerr,'+             'bp,dbp as bperr,rp,drp as rperr,teff,agaia,dupvar,ag,rp1,r1,r10,g as gmag,dg as gerr,gchi,gcontrib,'+             'r as rmag, dr as rerr,rchi,rcontrib,i as imag,di as ierr,ichi,icontrib,z as zmag,dz as zerr,zchi,zcontrib,nstat,'+             'j as jmag,dj as jerr,h as hmag,dh as herr,k as kmag,dk as kerr' 
                server = 'gp10.datalab.noao.edu' 
             
            # Use Postgres command with q3c cone search 
            refcattemp = repstr(file,'.fits','.txt') 
            cmd = "psql -h "+server+" -U datalab -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename+          " WHERE q3c_radial_query("+racol+","+deccol+","+stringize(cenra,ndec=4,/nocomma)+","+stringize(cendec,ndec=4,/nocomma)+          ","+stringize(radius,ndec=3)+")' > "+refcattemp 
            os.remove(refcattemp,/allow 
            os.remove(file,/allow 
            spawn,cmd,out,outerr 
            # Check for empty query 
            READLINE,refcattemp,tlines,nlineread=4 
            if len(tlines) < 4: 
                if not keyword_set(silent) : 
                    printlog,logf,'No Results' 
                ref = -1 
                nref = 0 
                return ref 
            #  Load ASCII file and create the FITS file 
            ref = importascii(refcattemp,/header,delim='|',skipline=2,/silent) 
            if keyword_set(saveref): 
                printlog,logf,'Saving catalog to file '+file 
                MWRFITS,ref,file,/create 
            os.remove(refcattemp,/allow 
             
            # Fix 0.0 mags/errs in ATLAS 
            if refname == 'ATLAS': 
                magcols = ['gaia','bp','rp','gmag','rmag','imag','zmag','jmag','hmag','kmag'] 
                errcols = ['gaiaerr','bperr','rperr','gerr','rerr','ierr','zerr','jerr','herr','kerr'] 
                tags = tag_names(ref) 
                # Set mags with 0.0 to 99.99 
                for j in range(len(magcols)): 
                    colind , = np.where(strupcase(tags) == strupcase(magcols[j]),ncolind) 
                    if colind > 0: 
                        bdmag , = np.where(ref.(colind[0]) <= 0.0,nbdmag) 
                        if nbdmag > 0 : 
                            ref[bdmag].(colind[0])=99.99 
                # Set errors with 0.0 to 9.99 
                for j in range(len(errcols)): 
                    colind , = np.where(strupcase(tags) == strupcase(errcols[j]),ncolind) 
                    if colind > 0: 
                        bderr , = np.where(ref.(colind[0]) <= 0.0,nbderr) 
                        if nbderr > 0 : 
                            ref[bderr].(colind[0])=9.99 
             
             
            # Use QUERYVIZIER 
            #---------------- 
            #   for low density with 2MASS/GAIA and always for GALEX and APASS 
        else: 
             
            # Use QUERYVIZIER for GALEX (python code has problems) 
            if refname == 'II/312/ais' or refname == 'GALEX': 
                #if refcat eq 'APASS' then cfa=0 else cfa=1  ; cfa doesn't have APASS 
                cfa = 1# problems with CDS VizieR and cfa has APASS now 
                if refcat == 'SAGE' : 
                    cfa=0 
                ref = QUERYVIZIER(refname,[cenra,cendec],radius*60,cfa=cfa,timeout=600,/silent) 
                # Check for failure 
                if size(ref,/type) != 8: 
                    if not keyword_set(silent) : 
                        printlog,logf,'Failure or No Results' 
                    ref = -1 
                    nref = 0 
                    return ref 
                 
                # Use Python code 
            else: 
                # Use python code, it's much faster, ~18x 
                tempfile = MKTEMP('vzr') 
                os.remove(tempfile+'.fits',/allow 
                pylines = 'python -c "from astroquery.vizier import Vizier#'+                'import astropy.units as u;'+                'import astropy.coordinates as coord;'+                'Vizier.TIMEOUT = 600;'+                'Vizier.ROW_LIMIT = -1;'+                'Vizier.cache_location = None;'+                'result = Vizier.query_region(coord.SkyCoord(ra='+strtrim(cenra,2)+', dec='+strtrim(cendec,2)+                ",unit=(u.deg,u.deg),frame='icrs'),width='"+strtrim(radius*60,2)+"m',catalog='"+refname+"');"+                "df=result[0];"+                "df.meta['description']=df.meta['description'][0:50];"+                "df.write('"+tempfile+".fits')"+'"' 
                spawn,pylines,out,errout 
                if os.path.exists(tempfile+'.fits') == 0: 
                    if not keyword_set(silent) : 
                        printlog,logf,'No Results' 
                    ref = -1 
                    nref = 0 
                    os.remove([tempfile,tempfile+'.fits'],/allow 
                    return ref 
                ref = MRDFITS(tempfile+'.fits',1,/silent) 
                os.remove([tempfile,tempfile+'.fits'],/allow 
             
            # Fix/homogenize the GAIA tags 
            if refname == 'GAIA': 
                nref = len(ref) 
                orig = ref 
                ref = replicate({source:0LL,ra_icrs:0.0d0,e_ra_icrs:0.0d0,de_icrs:0.0d0,e_de_icrs:0.0d0,fg:0.0d0,e_fg:0.0d0,gmag:0.0d0},nref) 
                struct_assign,orig,ref 
                ref.fg = orig._fg_ 
                ref.e_fg = orig.e__fg_ 
                ref.gmag = orig._gmag_ 
                undefine,orig 
            # Fix/homogenize the 2MASS tags 
            if refname == 'TMASS': 
                nref = len(ref) 
                orig = ref 
                ref = replicate({designation:'',raj2000:0.0d0,dej2000:0.0d0,jmag:0.0,e_jmag:0.0,hmag:0.0,e_hmag:0.0,kmag:0.0,e_kmag:0.0,qflg:''},nref) 
                struct_assign,orig,ref 
                ref.designation = orig._2mass 
                undefine,orig 
            # Fix NANs in ALLWISE 
            if refname == 'ALLWISE': 
                bd , = np.where(finite(ref._3_6_) == 0,nbd) 
                if nbd > 0: 
                    ref[bd]._3_6_ = 99.99 
                    ref[bd].e__3_6_ = 9.99 
                bd , = np.where(finite(ref._4_5_) == 0,nbd) 
                if nbd > 0: 
                    ref[bd]._4_5_ = 99.99 
                    ref[bd].e__4_5_ = 9.99 
             
            # Save the file 
            if keyword_set(saveref): 
                if not keyword_set(silent) : 
                    printlog,logf,'Saving catalog to file '+file 
                MWRFITS,ref,file,/create# only save if necessary 
    # use queryvizier.pro 
# do the query 
     
    if silent==False:
        printlog,logf,str(len(ref),2)+' sources found   dt=',stringize(time.time()-t0,ndec=1),' sec.' 
     
    count = len(ref) 
     
    return ref 


 
def getrefdata(filt,cenra,cendec,radius,count=count,saveref=saveref,silent=silent,
               dcr=dcr,modelmags=modelmags,logfile=logfile):


#+ 
# 
# GETREFDATA 
# 
# Get reference catalog information needed for a given filter. 
# 
# INPUTS: 
#  filter    Filter and instrument name, e.g. "c4d-u". 
#              This can be an array of multiple filters and 
#              then all the reference needed for all the filters 
#              will be included. 
#  cenra     Central RA for the search. 
#  cendec    Central DEC for the search. 
#  radius    Search radius in degrees. 
#  =dcr      The cross-matching radius in arcsec.  Default is 0.5". 
#  /saveref  Save the output to FILE. 
#  /silent   Don't print anything to the screen. 
#  /modelmags  Return the model magnitudes as well. 
#  =logfile  Filename to write output to. 
# 
# OUTPUTS: 
#  ref       Search results from the reference catalog all in one 
#              structure. 
#  =count    Number of elements in REF. 
# 
# USAGE: 
#  IDL>cat = getrefdata('g',cenra,cendec,radius,saveref=saveref) 
# 
# By D. Nidever  Sep 2017 
#- 
                          
    t0 = time.time() 
    undefine,ref 
    count = 0 
     
    # Not enough inputs 
    if len(filter) == 0 or len(cenra) == 0 or len(cendec) == 0 or    len(radius) == 0: 
        print('Syntax - cat = getrefdata(filter,cenra,cendec,radius,saveref=saveref,dcr=dcr,' 
        print('                          modelmags=modelmags,logfile=logfile)' 
        return -1 
     
    # Defaults 
    if len(dcr) == 0 :# arcsec 
        dcr = 0.5 
    if len(logfile) == 0: 
        logf = -1 
    else: 
        logf=logfile 
     
    # Check that we have psql installed 
    spawn,['which','psql'],out,errout,/noshell 
    if os.path.exists(out[0]) == 0: 
        print('No PSQL found on this sytem.' 
        return -1 
     
    # Figure out what reddening method we are using 
    #---------------------------------------------- 
    # Extinction types: 
    # 1 - SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2 
    # 2 - RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE 
    #     data available 
    # 3 - RJCE GLIMPSE, GLIMPSE data available 
    # 4 - RJCE SAGE, SAGE data available 
    ext_type = GETEXTTYPE(cenra,cendec,radius) 
     
    # If there is GLIMPSE or SAGE data, also get ALLWISE data 
    #  in case there is only partial coverage 
     
     
    if not keyword_set(silent): 
        printlog,logf,'Getting reference catalogs for:' 
        printlog,logf,'FILTER(S) = ',strjoin(filter,', ') 
        printlog,logf,'CENRA  = ',str(cenra,2) 
        printlog,logf,'CENDEC = ',str(cendec,2) 
        printlog,logf,'RADIUS = ',str(radius,2),' deg' 
        case ext_type of 
            1: printlog,logf,'Extinction Type: 1 - SFD' 
            2: printlog,logf,'Extinction Type: 2 - RJCE ALLWISE' 
            3: printlog,logf,'Extinction Type: 3 - RJCE GLIMPSE' 
            4: printlog,logf,'Extinction Type: 4 - RJCE SAGE' 
            else: 
 
     
    # Load the reference catalogs 
    #----------------------------- 
    # Figure out the reference catalogs that we need based on 
    #  filter-instrument combination 
    for i in range(len(filter)): 
        instfilt = filter[i] 
        # If no instrument information then assume it's DECam 
        if strpos(instfilt,'-') == -1 : 
            instfilt='c4d-'+instfilt 
        # Check the cases 
        CASE instfilt of 
            # DECam u-band 
            'c4d-u': begin 
            # Use GAIA, 2MASS and GALEX to calibrate 
            push,refcat,['2MASS-PSC','II/312/ais'] 
            if cendec <= 0 : 
                push,refcat,'Skymapper' 
 
        # DECam g-band 
        'c4d-g': begin 
        # Use PS1 if possible 
        if cendec > -29: 
            push,refcat,['2MASS-PSC','PS'] 
        else: 
            # Use 2MASS and Skymapper to calibrate 
            #push,refcat,['2MASS-PSC','Skymapper'] 
            push,refcat,['2MASS-PSC','ATLAS'] 
 
    # DECam r-band 
    'c4d-r': begin 
    # Use PS1 if possible 
    if cendec > -29: 
        push,refcat,['2MASS-PSC','PS'] 
    else: 
        # Use 2MASS and Skymapper to calibrate 
        #push,refcat,['2MASS-PSC','Skymapper'] 
        push,refcat,['2MASS-PSC','ATLAS'] 
 
# DECam i-band 
'c4d-i': begin 
# Use PS1 if possible 
if cendec > -29: 
    push,refcat,['2MASS-PSC','PS'] 
else: 
    # Use Skymapper and 2MASS to calibrate 
    #push,refcat,['2MASS-PSC','Skymapper'] 
    push,refcat,['2MASS-PSC','ATLAS'] 
 
# DECam z-band 
'c4d-z': begin 
# Use PS1 if possible 
if cendec > -29: 
push,refcat,['2MASS-PSC','PS'] 
else: 
# Use Skymapper and 2MASS to calibrate 
#push,refcat,['2MASS-PSC','Skymapper'] 
push,refcat,['2MASS-PSC','ATLAS'] 
 
# DECam Y-band 
'c4d-Y': begin 
# Use PS1 if possible 
if cendec > -29: 
push,refcat,['2MASS-PSC','PS'] 
else: 
# Use 2MASS to calibrate 
push,refcat,['2MASS-PSC'] 
 
# DECam VR-band 
'c4d-VR': begin 
# Use PS1 if possible 
if cendec > -29: 
push,refcat,['2MASS-PSC','PS'] 
else: 
push,refcat,['2MASS-PSC','ATLAS'] 
 
# Bok+90Prime g-band 
'ksb-g': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Bok+90Prime r-band 
'ksb-r': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Mosaic3 g-band 
'k4m-g': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Mosaic3 r-band 
'k4m-r': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Mosaic3 i-band 
'k4m-i': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Mosaic3 z-band 
'k4m-z': begin 
# Use PS1 
push,refcat,['2MASS-PSC','PS'] 
 
# Mosaic3 VR-band 
'k4m-VR': begin 
# Use GAIA G-band to calibrate 
push,refcat,['2MASS-PSC','PS'] 
#push,refcat,['2MASS-PSC'] 
 
else: begin 
printlog,logf,filter,' not currently supported' 
#return,-1 
 
 
# filter loop 
# Extinction catalogs 
if ext_type >= 2 : 
push,refcat,'ALLWISE' 
if ext_type == 3 : 
push,refcat,'GLIMPSE' 
if ext_type == 4 : 
push,refcat,'SAGE' 
# If near DEC=-30 then load BOTH PS and ATLAS 
if (cendec > -34 and cendec < -26) and    ((where(refcat == 'PS'))[0] != -1 or (where(refcat == 'ATLAS'))[0] != -1): 
push,refcat,['PS','ATLAS']# code below removes duplicates 
# Some catalogs 
if len(refcat) > 0: 
# Get the unique ones 
ui = np.uniq(refcat,np.argsort(refcat)) 
refcat = refcat[ui] 
# Always add Gaia at the beginning 
refcat = ['GAIADR2',refcat] 
 else refcat='GAIADR2'#'GAIA/GAIA' 
nrefcat = len(refcat) 
 
# Figure out the new columns that need to be added 
undefine,newtags 
for i in range(nrefcat): 
case refcat[i] of 
'GAIADR2': push,newtags,['source','ra','ra_error','dec','dec_error','pmra','pmra_error','pmdec','pmdecerr','gmag','e_gmag','bp','e_bp','rp','e_rp'] 
'2MASS-PSC': push,newtags,['jmag','e_jmag','hmag','e_hmag','kmag','e_kmag','qflg'] 
'PS': push,newtags,['ps_gmag','ps_rmag','ps_imag','ps_zmag','ps_ymag'] 
'APASS': push,newtags,['apass_gmag','e_apass_gmag','apass_rmag','e_apass_rmag'] 
'II/312/ais': push,newtags,['nuv','e_nuv']# Galex 
'Skymapper': push,newtags,['sm_umag','e_sm_umag','sm_gmag','e_sm_gmag','sm_rmag','e_sm_rmag','sm_imag','e_sm_imag','sm_zmag','e_sm_zmag']# Skymapper DR1 
'ALLWISE': push,newtags,['w1mag','e_w1mag','w2mag','e_w2mag'] 
'GLIMPSE': push,newtags,['gl_36mag','e_gl_36mag','gl_45mag','e_gl_45mag'] 
'SAGE': push,newtags,['sage_36mag','e_sage_36mag','sage_45mag','e_sage_45mag'] 
'ATLAS': push,newtags,['atlas_gmag','e_atlas_gmag','atlas_gcontrib','atlas_rmag','e_atlas_rmag','atlas_rcontrib',                         'atlas_imag','e_atlas_imag','atlas_icontrib','atlas_zmag','e_atlas_zmag','atlas_zcontrib'] 
else: import pdb; pdb.set_trace(),refcat[i]+' NOT SUPPORTED' 
 
push,newtags,['ebv_sfd','ejk','e_ejk','ext_type'] 
if keyword_set(modelmags) : 
push,newtags,'model_mag' 
nnewtags = len(newtags) 
 
# Load the necessary catalogs 
nrefcat = len(refcat) 
if not keyword_set(silent) : 
printlog,logf,str(nrefcat,2),' reference catalogs to load: '+strjoin(refcat,', ') 
for i in range(nrefcat): 
t0 = time.time() 
if not keyword_set(silent) : 
printlog,logf,'Loading ',refcat[i],' reference catalog' 
 
# Load the catalog 
ref1 = GETREFCAT(cenra,cendec,radius,refcat[i],count=nref1,silent=silent,logfile=logf) 
if nref1 == 0 : 
goto,BOMB 
tags1 = tag_names(ref1) 
 
# Initialize the catalog 
undefine,schema 
for j in range(nnewtags): 
val0 = 99.99 
if newtags[j] == 'qflg' : 
val0='' 
if newtags[j] == 'ext_type' : 
val0=0 
if newtags[j] == 'source' : 
val0=0LL 
if newtags[j] == 'ra' : 
val0=0.0d0 
if newtags[j] == 'dec' : 
val0=0.0d0 
if stregex(newtags[j],'contrib',/boolean) == 1 : 
val0=-1L 
if len(schema) == 0: 
schema=create_struct(newtags[j],val0) 
else: 
schema = create_struct(schema,newtags[j],val0) 
 
 
# First successful one, initialize the catalog 
if len(ref) == 0: 
ref = replicate(schema,nref1) 
struct_assign,ref1,ref,/nozero 
if tag_exist(ref1,'RA') == 0: 
ref.ra = ref.ra_icrs 
ref.dec = ref.de_icrs 
ind1 = lindgen(nref1) 
ind2 = lindgen(nref1) 
nmatch = nref1 
 
# Second and later 
else: 
 
# Get RA/DEC columns 
# 2MASS, Galex, APASS use RAJ2000 
# PS uses RA/DEC 
if (refcat[i] != 'PS' and refcat[i] != 'ALLWISE' and refcat[i] != 'ATLAS'): 
raind, = np.where(tags1 == 'RAJ2000',nraind) 
else: 
raind, = np.where(tags1 == 'RA',nraind) 
if (refcat[i] != 'PS' and refcat[i] != 'ALLWISE' and refcat[i] != 'ATLAS'): 
decind, = np.where(tags1 == 'DEJ2000',ndecind) 
else: 
decind, = np.where(tags1 == 'DEC',ndecind) 
 
# Crossmatch 
SRCMATCH,ref.ra,ref.dec,ref1.(raind),ref1.(decind),dcr,ind1,ind2,/sph,count=nmatch 
if not keyword_set(silent) : 
printlog,logf,str(nmatch,2)+' matches' 
 
# Add magnitude columns 
if nmatch > 0: 
case refcat[i] of 
'GAIADR2': begin 
temp = ref[ind1] 
struct_assign,ref1[ind2],temp,/nozero 
temp.e_gmag = 2.5*alog10(1.0+ref1[ind2].e_fg/ref1[ind2].fg) 
temp.e_bp = 2.5*alog10(1.0+ref1[ind2].e_fbp/ref1[ind2].fbp) 
temp.e_rp = 2.5*alog10(1.0+ref1[ind2].e_frp/ref1[ind2].frp) 
ref[ind1] = temp 
 
'2MASS-PSC': begin 
ref[ind1].jmag = ref1[ind2].jmag 
ref[ind1].e_jmag = ref1[ind2].e_jmag 
ref[ind1].hmag = ref1[ind2].hmag 
ref[ind1].e_hmag = ref1[ind2].e_hmag 
ref[ind1].kmag = ref1[ind2].kmag 
ref[ind1].e_kmag = ref1[ind2].e_kmag 
ref[ind1].qflg = ref1[ind2].qflg 
 
'PS': begin 
ref[ind1].ps_gmag = ref1[ind2].gmag 
ref[ind1].ps_rmag = ref1[ind2].rmag 
ref[ind1].ps_imag = ref1[ind2].imag 
ref[ind1].ps_zmag = ref1[ind2].zmag 
ref[ind1].ps_ymag = ref1[ind2].ymag 
 
'APASS': begin 
ref[ind1].apass_gmag = ref1[ind2].g_mag 
ref[ind1].e_apass_gmag = ref1[ind2].e_g_mag 
ref[ind1].apass_rmag = ref1[ind2].r_mag 
ref[ind1].e_apass_rmag = ref1[ind2].e_r_mag 
 
'II/312/ais': begin 
if tag_exist(ref1,'NUV'): 
ref[ind1].nuv = ref1[ind2].nuv 
ref[ind1].e_nuv = ref1[ind2].e_nuv 
else: 
ref[ind1].nuv = ref1[ind2].nuvmag 
ref[ind1].e_nuv = ref1[ind2].e_nuvmag 
 
'Skymapper': begin 
ref[ind1].sm_umag = ref1[ind2].sm_umag 
ref[ind1].e_sm_umag = ref1[ind2].e_sm_umag 
ref[ind1].sm_gmag = ref1[ind2].sm_gmag 
ref[ind1].e_sm_gmag = ref1[ind2].e_sm_gmag 
ref[ind1].sm_rmag = ref1[ind2].sm_rmag 
ref[ind1].e_sm_rmag = ref1[ind2].e_sm_rmag 
ref[ind1].sm_imag = ref1[ind2].sm_imag 
ref[ind1].e_sm_imag = ref1[ind2].e_sm_imag 
ref[ind1].sm_zmag = ref1[ind2].sm_zmag 
ref[ind1].e_sm_zmag = ref1[ind2].e_sm_zmag 
 
'ATLAS': begin 
ref[ind1].atlas_gmag = ref1[ind2].gmag 
ref[ind1].e_atlas_gmag = ref1[ind2].gerr 
ref[ind1].atlas_gcontrib = ref1[ind2].gcontrib 
ref[ind1].atlas_rmag = ref1[ind2].rmag 
ref[ind1].e_atlas_rmag = ref1[ind2].rerr 
ref[ind1].atlas_rcontrib = ref1[ind2].rcontrib 
ref[ind1].atlas_imag = ref1[ind2].imag 
ref[ind1].e_atlas_imag = ref1[ind2].ierr 
ref[ind1].atlas_icontrib = ref1[ind2].icontrib 
ref[ind1].atlas_zmag = ref1[ind2].zmag 
ref[ind1].e_atlas_zmag = ref1[ind2].zerr 
ref[ind1].atlas_zcontrib = ref1[ind2].zcontrib 
 
'ALLWISE': begin 
ref[ind1].w1mag = ref1[ind2].w1mag 
ref[ind1].e_w1mag = ref1[ind2].e_w1mag 
ref[ind1].w2mag = ref1[ind2].w2mag 
ref[ind1].e_w2mag = ref1[ind2].e_w2mag 
 
'GLIMPSE': begin 
ref[ind1].gl_36mag = ref1[ind2]._3_6mag 
ref[ind1].e_gl_36mag = ref1[ind2].e_3_6mag 
ref[ind1].gl_45mag = ref1[ind2]._4_5mag 
ref[ind1].e_gl_45mag = ref1[ind2].e_4_5mag 
 
'SAGE': begin 
ref[ind1].sage_36mag = ref1[ind2].__3_6_ 
ref[ind1].e_sage_36mag = ref1[ind2].e__3_6_ 
ref[ind1].sage_45mag = ref1[ind2].__4_5_ 
ref[ind1].e_sage_45mag = ref1[ind2].e__4_5_ 
 
else: import pdb; pdb.set_trace(),catname+' NOT SUPPORTED' 
 
 
# Add leftover ones 
if nmatch < len(ref1): 
left1 = ref1 
if nmatch > 0 : 
remove,ind2,left1 
nleft1 = len(left1) 
new = replicate(schema,nleft1) 
new.ra = left1.(raind) 
new.dec = left1.(decind) 
 
case refcat[i] of 
'GAIADR2': begin 
temp = ref[ind1] 
struct_assign,left1,new 
new.e_gmag = 2.5*alog10(1.0+left1.e_fg/left1.fg) 
new.e_bp = 2.5*alog10(1.0+left1.e_fbp/left1.fbp) 
new.e_rp = 2.5*alog10(1.0+left1.e_frp/left1.frp) 
 
'2MASS-PSC': begin 
new.jmag = left1.jmag 
new.e_jmag = left1.e_jmag 
new.hmag = left1.hmag 
new.e_hmag = left1.e_hmag 
new.kmag = left1.kmag 
new.e_kmag = left1.e_kmag 
new.qflg = left1.qflg 
 
'PS': begin 
new.ps_gmag = left1.gmag 
new.ps_rmag = left1.rmag 
new.ps_imag = left1.imag 
new.ps_zmag = left1.zmag 
new.ps_ymag = left1.ymag 
 
'APASS': begin 
new.apass_gmag = left1.g_mag 
new.e_apass_gmag = left1.e_g_mag 
new.apass_rmag = left1.r_mag 
new.e_apass_rmag = left1.e_r_mag 
 
'II/312/ais': begin 
if tag_exist(left1,'NUV'): 
new.nuv = left1.nuv 
new.e_nuv = left1.e_nuv 
else: 
new.nuv = left1.nuvmag 
new.e_nuv = left1.e_nuvmag 
 
'Skymapper': begin 
new.sm_umag = left1.sm_umag 
new.e_sm_umag = left1.e_sm_umag 
new.sm_gmag = left1.sm_gmag 
new.e_sm_gmag = left1.e_sm_gmag 
new.sm_rmag = left1.sm_rmag 
new.e_sm_rmag = left1.e_sm_rmag 
new.sm_imag = left1.sm_imag 
new.e_sm_imag = left1.e_sm_imag 
new.sm_zmag = left1.sm_zmag 
new.e_sm_zmag = left1.e_sm_zmag 
 
'ATLAS': begin 
new.atlas_gmag = left1.gmag 
new.e_atlas_gmag = left1.gerr 
new.atlas_gcontrib = left1.gcontrib 
new.atlas_rmag = left1.rmag 
new.e_atlas_rmag = left1.rerr 
new.atlas_rcontrib = left1.rcontrib 
new.atlas_imag = left1.imag 
new.e_atlas_imag = left1.ierr 
new.atlas_icontrib = left1.icontrib 
new.atlas_zmag = left1.zmag 
new.e_atlas_zmag = left1.zerr 
new.atlas_zcontrib = left1.zcontrib 
 
'ALLWISE': begin 
new.w1mag = left1.w1mag 
new.e_w1mag = left1.e_w1mag 
new.w2mag = left1.w2mag 
new.e_w2mag = left1.e_w2mag 
 
'GLIMPSE': begin 
new.gl_36mag = left1._3_6mag 
new.e_gl_36mag = left1.e_3_6mag 
new.gl_45mag = left1._4_5mag 
new.e_gl_45mag = left1.e_4_5mag 
 
'SAGE': begin 
new.sage_36mag = left1.__3_6_ 
new.e_sage_36mag = left1.e__3_6_ 
new.sage_45mag = left1.__4_5_ 
new.e_sage_45mag = left1.e__4_5_ 
 
else: import pdb; pdb.set_trace(),catname+' NOT SUPPORTED' 
 
 
# Combine the two 
old = ref 
ref = replicate(schema,len(old)+nleft1) 
ref[0:len(old)-1] = old 
ref[len(old)::] = new 
undefine,old,new,left1 
 
BOMB: 
 
# Get extinction 
#---------------- 
GETREDDENING,ref,ext_type 
 
count = len(ref) 
 
print('dt=',stringize(time.time()-t0,ndec=1),' sec' 
 
return ref 
 

def getreddening(ref,ext_type):
    """
#+ 
# 
# GETREDDENING 
# 
# This calculates E(J-Ks) reddening using reference catalog data for 
# the NSC. 
# 
# INPUTS: 
#  ref      The structure with the reference catalog data. 
#  exttype  Extinction type: 
#             1 - SFD 
#             2 - RJCE ALLWISE 
#             3 - RJCE GlIMPSE 
#             4 - RJCE SAGE 
# 
# OUTPUTS: 
#  The REF structure EJK column and EXT_TYPE columns are updated. 
# 
# USAGE: 
#  IDL>getreddening,ref,ext_type 
# 
# By D. Nidever  Feb 2019 
#- 
    """
      
    # Not enough inputs 
    if len(ref) == 0 or ext_type == 0: 
        print('Syntax - getreddening,ref,ext_type' 
        return 
     
    # Add SFD reddening 
    GLACTC,ref.ra,ref.dec,2000.0,glon,glat,1,/deg 
    ebv = dust_getval(glon,glat,/noloop,/interp) 
    ref.ebv_sfd = ebv 
     
    # Start with SFD extinction for all 
    ejk_sfd = 1.5*0.302*ref.ebv_sfd 
    ref.ejk = ejk_sfd 
    ref.e_ejk = 0.1# not sure what to put here 
    bd , = np.where(ref.ebv_sfd > 0.3,nbd)# E(B-V)>0.3 is roughly where things "break down" 
    if nbd > 0 : 
        ref[bd].e_ejk=1.0 
    ref.ext_type = 1 
     
    # RJCE extinction 
    if ext_type >= 2: 
         
        # RJCE GLIMPSE, type=3 
        if ext_type == 3: 
            gdglimpse , = np.where(ref.jmag < 50 and ref.hmag < 50 and ref.kmag < 50 and                       ref.gl_45mag < 50,ngdglimpse) 
            if ngdglimpse > 0: 
                ejk = 1.5*0.918*(ref[gdglimpse].hmag-ref[gdglimpse].gl_45mag-0.08) 
                e_ejk = 1.5*0.918*sqrt(ref[gdglimpse].e_hmag**2+ref[gdglimpse].e_gl_45mag**2) 
                #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdglimpse] and e_ejk lt 0.2,ngdejk) 
                gdejk , = np.where(ejk < ejk_sfd[gdglimpse],ngdejk) 
                if ngdejk > 0: 
                    ref[gdglimpse[gdejk]].ejk = ejk[gdejk] > 0.0 
                    ref[gdglimpse[gdejk]].e_ejk = e_ejk[gdejk] 
                    ref[gdglimpse[gdejk]].ext_type = 3 
         
        # RJCE SAGE, type=4 
        if ext_type == 4: 
            gdsage , = np.where(ref.jmag < 50 and ref.hmag < 50 and ref.kmag < 50 and                    ref.sage_45mag < 50,ngdsage) 
            if ngdsage > 0: 
                ejk = 1.5*0.918*(ref[gdsage].hmag-ref[gdsage].sage_45mag-0.08) 
                e_ejk = 1.5*0.918*sqrt(ref[gdsage].e_hmag**2+ref[gdsage].e_sage_45mag**2) 
                #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdsage] and e_ejk lt 0.2,ngdejk) 
                gdejk , = np.where(ejk < ejk_sfd[gdsage],ngdejk) 
                if ngdejk > 0: 
                    ref[gdsage[gdejk]].ejk = ejk[gdejk] > 0.0 
                    ref[gdsage[gdejk]].e_ejk = e_ejk[gdejk] 
                    ref[gdsage[gdejk]].ext_type = 4 
         
        # RJCE ALLWISE, type=2 
        #   Use ALLWISE for those that don't have IRAC data 
        gdwise , = np.where(ref.jmag < 50 and ref.hmag < 50 and ref.kmag < 50 and                  ref.w2mag < 50 and ref.ext_type <= 1,ngdwise) 
        if ngdwise > 0: 
            ejk = 1.5*0.918*(ref[gdwise].hmag-ref[gdwise].w2mag-0.05) 
            e_ejk = 1.5*0.918*sqrt(ref[gdwise].e_hmag**2+ref[gdwise].e_w2mag**2) 
            #gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdwise] and e_ejk lt 0.2,ngdejk) 
            gdejk , = np.where(ejk < ejk_sfd[gdwise],ngdejk) 
            if ngdejk > 0: 
                ref[gdwise[gdejk]].ejk = ejk[gdejk] > 0.0 
                ref[gdwise[gdejk]].e_ejk = e_ejk[gdejk] 
                ref[gdwise[gdejk]].ext_type = 2 
     
    # No good reddening 
    bd , = np.where(ref.ext_type == 1 and ref.ebv_sfd > 0.3,nbd)# E(B-V)>0.3 is roughly where things "break down" 
    if nbd > 0: 
        ref[bd].ejk = 999999.0 
        ref[bd].e_ejk = 999999.0 
        ref[bd].ext_type = 0 
     
    # Fix NANs in E_EJK 
    bd , = np.where(finite(ref.e_ejk) == 0,nbd) 
    if nbd > 0 : 
        ref[bd].e_ejk = 9.99 
     
 
 
