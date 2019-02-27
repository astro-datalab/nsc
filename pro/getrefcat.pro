;+
;
; GETREFCAT
;
; Get reference catalog information from DL database
;
; INPUTS:
;  cenra     Central RA for the search.
;  cendec    Central DEC for the search.
;  radius    Search radius in degrees.
;  refcat    Reference catalog name (e.g. 2MASS, Gaia, etc.).
;  =file     The file to save to or search for existing catalog.
;  /saveref  Save the output to FILE.
;  /silent   Don't print anything to the screen.
;  =logfile  Filename for logging output.
;
; OUTPUTS:
;  ref       Search results from the reference catalog.
;  =count    Number of elements in REF.
;
; USAGE:
;  IDL>cat = getrefcat(cenra,cendec,radius,refcat,file=file,saveref=saveref)
;
; By D. Nidever  Sep 2017
;-

function getrefcat,cenra,cendec,radius,refcat,count=count,file=file,saveref=saveref,silent=silent,logfile=logfile

undefine,ref
count = 0
t0 = systime(1)

; Not enough inputs
if n_elements(cenra) eq 0 or n_elements(cendec) eq 0 or n_elements(radius) eq 0 or $
   n_elements(refcat) eq 0 then begin
  print,'Syntax - cat = getrefcat(cenra,cendec,radius,refcat,file=file,saveref=saveref)'
  return,-1
endif

; Check that we have psql installed
spawn,['which','psql'],out,errout,/noshell
if file_test(out[0]) eq 0 then begin
  print,'No PSQL found on this sytem.'
  return,-1
endif

; Defaults
if n_elements(logfile) eq 0 then logf=-1 else logf=logfile

; Temporary directory
; /tmp is often small and get get fille dup
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'

; FLIP THIS AROUND, INPUT SHOULD BE THE "EASY" VERSION!!!
refname = strupcase(refcat)
if refname eq 'II/312/AIS' then refname='GALEX'
if refname eq '2MASS-PSC' then refname='TMASS'
if refname eq '2MASS' then refname='TMASS'
if refname eq 'GAIA/GAIA' then refname='GAIA'
if refname eq 'Skymapper' then refname='SKYMAPPER'
if refname eq 'GLIMPSE' then refname='II/293/glimpse'
if refname eq 'SAGE' then refname='II/305/archive'

if n_elements(file) eq 0 then file=tmpdir+'ref_'+stringize(cenra,ndec=5)+'_'+stringize(cendec,ndec=5)+'_'+stringize(radius,ndec=3)+'_'+refname+'.fits'

if not keyword_set(silent) then $
  printlog,logf,'Querying '+refname+': RA='+stringize(cenra,ndec=5)+' DEC='+stringize(cendec,ndec=5)+' Radius='+stringize(radius,ndec=3)

; Loading previously loaded file
if file_test(file) eq 1 then begin
  if not keyword_set(silent) then $
    printlog,logf,'Loading previously-saved file ',file
  ref = MRDFITS(file,1,/silent)

; Do the Query
;--------------
endif else begin

  ; Use DataLab database search for Gaia and 2MASS if density is high
  if (refname eq 'TMASS' or refname eq 'GAIA' or refname eq 'GAIADR2' or refname eq 'PS' or refname eq 'SKYMAPPER' or refname eq 'ALLWISE') then begin
    if refname eq 'TMASS' then begin
      tablename = 'twomass.psc'
      cols = 'designation,ra as raj2000,dec as dej2000,j_m as jmag,j_cmsig as e_jmag,h_m as hmag,h_cmsig as e_hmag,k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg'
      server = 'gp01.datalab.noao.edu'
      ;server = 'dldb1.sdm.noao.edu'
    endif
    racol = 'ra'
    deccol = 'dec'
    if refname eq 'GAIA' then begin
      tablename = 'gaia_dr1.gaia_source'
      cols = 'source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,dec as de_icrs,dec_error as e_de_icrs,'+$
             'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,phot_g_mean_mag as gmag'
      server = 'gp01.datalab.noao.edu'
      ;server = 'dldb1.sdm.noao.edu'
    endif
    if refname eq 'GAIADR2' then begin
      tablename = 'gaia_dr2.gaia_source'
      cols = 'source_id as source,ra,ra_error,dec,dec_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'+$
             'phot_g_mean_mag as gmag,phot_bp_mean_mag as bp,phot_rp_mean_mag as rp'
      server = 'gp01.datalab.noao.edu'
    endif
    if refname eq 'PS' then begin
      ;tablename = 'cp_calib.ps1'
      tablename = 'public.ps1'
      cols = 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag'
      server = 'gp02.datalab.noao.edu'
    endif
    if refname eq 'SKYMAPPER' then begin
      tablename = 'skymapper_dr1.master'
      cols = 'raj2000, dej2000, g_psf as sm_gmag, e_g_psf as e_sm_gmag, r_psf as sm_rmag, e_r_psf as e_sm_rmag, i_psf as sm_imag, '+$
             'e_i_psf as e_sm_imag, z_psf as sm_zmag, e_z_psf as e_sm_zmag'
      server = 'gp01.datalab.noao.edu'
      racol = 'raj2000'
      deccol = 'dej2000'
    endif
    if refname eq 'ALLWISE' then begin
      tablename = 'allwise.source'
      cols = 'ra, dec, w1mpro as w1mag, w1sigmpro as e_w1mag, w2mpro as w2mag, w2sigmpro as e_w2mag'
      server = 'gp01.datalab.noao.edu'
    endif
    
    ; Use Postgres command with q3c cone search
    refcattemp = repstr(file,'.fits','.txt')
    cmd = "psql -h "+server+" -U datalab -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename+$
          " WHERE q3c_radial_query("+racol+","+deccol+","+stringize(cenra,ndec=4,/nocomma)+","+stringize(cendec,ndec=4,/nocomma)+$
          ","+stringize(radius,ndec=3)+")' > "+refcattemp
    file_delete,refcattemp,/allow
    file_delete,file,/allow
    spawn,cmd,out,outerr
    ; Check for empty query
    READLINE,refcattemp,tlines,nlineread=4
    if n_elements(tlines) lt 4 then begin
      if not keyword_set(silent) then printlog,logf,'No Results'
      ref = -1
      nref = 0
      return,ref
    endif
    ;  Load ASCII file and create the FITS file
    ref = importascii(refcattemp,/header,delim='|',skipline=2,/silent)
    if keyword_set(saveref) then begin
      printlog,logf,'Saving catalog to file '+file
      MWRFITS,ref,file,/create
    endif
    file_delete,refcattemp,/allow

  ; Use QUERYVIZIER
  ;   for low density with 2MASS/GAIA and always for GALEX and APASS
  endif else begin

    ; Use python code, it's much faster, ~18x
    tempfile = MKTEMP('vzr')
    file_delete,tempfile+'.fits',/allow
    pylines = 'python -c "from astroquery.vizier import Vizier;'+$
              'import astropy.units as u;'+$
              'import astropy.coordinates as coord;'+$
              'Vizier.TIMEOUT = 600;'+$
              'Vizier.ROW_LIMIT = -1;'+$
              'result = Vizier.query_region(coord.SkyCoord(ra='+strtrim(cenra,2)+', dec='+strtrim(cendec,2)+$
              ",unit=(u.deg,u.deg),frame='icrs'),width='"+strtrim(radius*60,2)+"m',catalog='"+refname+"');"+$
              "df=result[0];df.write('"+tempfile+".fits')"+'"'
    spawn,pylines,out,errout
    if file_test(tempfile+'.fits') eq 0 then begin
      if not keyword_set(silent) then printlog,logf,'No Results'
      ref = -1
      nref = 0
      return,ref
    endif
    ref = MRDFITS(tempfile+'.fits',1,/silent)
    file_delete,tempfile+'.fits'

    ;;if refcat eq 'APASS' then cfa=0 else cfa=1  ; cfa doesn't have APASS
    ;cfa = 1  ; problems with CDS VizieR and cfa has APASS now
    ;if refcat eq 'SAGE' then cfa=0
    ;ref = QUERYVIZIER(refname,[cenra,cendec],radius*60,cfa=cfa,timeout=600,/silent)
    ;; Check for failure
    ;if size(ref,/type) ne 8 then begin
    ;  if not keyword_set(silent) then printlog,logf,'Failure or No Results'
    ;  ref = -1
    ;  nref = 0
    ;  return,ref
    ;endif

    ; Fix/homogenize the GAIA tags
    if refname eq 'GAIA' then begin
      nref = n_elements(ref)
      orig = ref
      ref = replicate({source:0LL,ra_icrs:0.0d0,e_ra_icrs:0.0d0,de_icrs:0.0d0,e_de_icrs:0.0d0,fg:0.0d0,e_fg:0.0d0,gmag:0.0d0},nref)
      struct_assign,orig,ref
      ref.fg = orig._fg_
      ref.e_fg = orig.e__fg_
      ref.gmag = orig._gmag_
      undefine,orig
    endif
    ; Fix/homogenize the 2MASS tags
    if refname eq 'TMASS' then begin
      nref = n_elements(ref)
      orig = ref
      ref = replicate({designation:'',raj2000:0.0d0,dej2000:0.0d0,jmag:0.0,e_jmag:0.0,hmag:0.0,e_hmag:0.0,kmag:0.0,e_kmag:0.0,qflg:''},nref)
      struct_assign,orig,ref
      ref.designation = orig._2mass
      undefine,orig
    endif
    ;; Fix NANs in ALLWISE
    if refname eq 'ALLWISE' then begin
      bd = where(finite(ref._3_6_) eq 0,nbd)
      if nbd gt 0 then begin
        ref[bd]._3_6_ = 99.99
        ref[bd].e__3_6_ = 9.99
      endif
      bd = where(finite(ref._4_5_) eq 0,nbd)
      if nbd gt 0 then begin
        ref[bd]._4_5_ = 99.99
        ref[bd].e__4_5_ = 9.99
      endif
    endif

    ; Save the file
    if keyword_set(saveref) then begin
      if not keyword_set(silent) then $
        printlog,logf,'Saving catalog to file '+file
      MWRFITS,ref,file,/create  ; only save if necessary
    endif
  endelse
endelse

if not keyword_set(silent) then $
  printlog,logf,strtrim(n_elements(ref),2)+' sources found   dt=',stringize(systime(1)-t0,ndec=1),' sec.'

count = n_elements(ref)

return,ref

stop

end
