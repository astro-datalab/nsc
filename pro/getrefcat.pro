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

function getrefcat,cenra,cendec,radius,refcat,count=count,file=file,saveref=saveref,silent=silent

undefine,ref
count = 0

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

logf = -1

; FLIP THIS AROUND, INPUT SHOULD BE THE "EASY" VERSION!!!
refname = strupcase(refcat)
if refname eq 'II/312/AIS' then refname='GALEX'
if refname eq '2MASS-PSC' then refname='TMASS'
if refname eq '2MASS' then refname='TMASS'
if refname eq 'GAIA/GAIA' then refname='GAIA'

if n_elements(file) eq 0 then file='/tmp/ref_'+stringize(cenra,ndec=5)+'_'+stringize(cendec,ndec=5)+'_'+stringize(radius,ndec=3)+'_'+refname+'.fits'

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
  if (refname eq 'TMASS' or refname eq 'GAIA' or refname eq 'PS') then begin
    if refname eq 'TMASS' then begin
      tablename = 'twomass.psc'
      cols = 'designation,ra as raj2000,dec as dej2000,j_m as jmag,j_cmsig as e_jmag,h_m as hmag,h_cmsig as e_hmag,k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg'
      server = 'dldb1.sdm.noao.edu'
    endif
    if refname eq 'GAIA' then begin
      tablename = 'gaia_dr1.gaia_source'
      cols = 'source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,dec as de_icrs,dec_error as e_de_icrs,'+$
             'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,phot_g_mean_mag as gmag'
      server = 'dldb1.sdm.noao.edu'
    endif
    if refname eq 'PS' then begin
      tablename = 'cp_calib.ps1'
      cols = 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag'
      server = 'gp02.datalab.noao.edu'
    endif

    ; Use Postgres command with q3c cone search                                                                                                                                    
    refcattemp = repstr(file,'.fits','.txt')
    cmd = "psql -h "+server+" -U datalab -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename+$
          " WHERE q3c_radial_query(ra,dec,"+stringize(cenra,ndec=4,/nocomma)+","+stringize(cendec,ndec=4,/nocomma)+$
          ","+stringize(radius,ndec=3)+")' > "+refcattemp
    file_delete,refcattemp,/allow
    file_delete,file,/allow
    spawn,cmd,out,outerr
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
    if refcat eq 'APASS' then cfa=0 else cfa=1  ; cfa doesn't have APASS                                                                                                        
    ref = QUERYVIZIER(refcat,[cenra,cendec],radius*60,cfa=cfa)

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

    ; Save the file
    if keyword_set(saveref) then begin
      if not keyword_set(silent) then $
        printlog,logf,'Saving catalog to file '+file
      MWRFITS,ref,file,/create  ; only save if necessary
    endif
  endelse
endelse

if not keyword_set(silent) then $
  printlog,logf,strtrim(n_elements(ref),2)+' sources found'

count = n_elements(ref)

return,ref

stop

end
