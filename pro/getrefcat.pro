function getrefcat,cenra,cendec,refcat,radius,refcatfile=refcatfile,saveref=saveref


varname = refcat
if varname eq 'II/312/ais' then varname='GALEX'
if varname eq '2MASS-PSC' then varname='TMASS'
if varname eq 'GAIA/GAIA' then varname='GAIA'
;refcatfile = expdir+'/'+base+'_'+varname+'.fits'

;if file_test(refcatfile) eq 1 and not keyword_set(redo) then begin                                                                                                                
; Loading previously loaded file                                                                                                                                                   
if file_test(refcatfile) eq 1 then begin
  printlog,logf,'  Loading previously-saved file ',refcatfile
  ;(SCOPE_VARFETCH(varname,/enter)) = MRDFITS(refcatfile,1,/silent)
  ;ref = SCOPE_VARFETCH(varname)
  ref = MRDFITS(refcatfile,1,/silent)

; Do the Query                                                                                                                                                                     
;--------------                                                                                                                                                                    
endif else begin

  ; Use DataLab database search for Gaia and 2MASS if density is high                                                                                                              
  if (varname eq 'TMASS' or varname eq 'GAIA' or varname eq 'PS') then begin
    if varname eq 'TMASS' then begin
      tablename = 'twomass.psc'
      cols = 'designation,ra as raj2000,dec as dej2000,j_m as jmag,j_cmsig as e_jmag,h_m as hmag,h_cmsig as e_hmag,k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg'
      server = 'dldb1.sdm.noao.edu'
    endif
    if varname eq 'GAIA' then begin
      tablename = 'gaia_dr1.gaia_source'
      cols = 'source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,dec as de_icrs,dec_error as e_de_icrs,'+$
             'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,phot_g_mean_mag as gmag'
      server = 'dldb1.sdm.noao.edu'
    endif
    if varname eq 'PS' then begin
      tablename = 'cp_calib.ps1'
      cols = 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag'
      server = 'gp02.datalab.noao.edu'
    endif

    ; Use Postgres command with q3c cone search                                                                                                                                    
    refcattemp = repstr(refcatfile,'.fits','.txt')
    cmd = "psql -h "+server+" -U datalab -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename+$
          " WHERE q3c_radial_query(ra,dec,"+stringize(cenra,ndec=4,/nocomma)+","+stringize(cendec,ndec=4,/nocomma)+$
          ","+stringize(radius,ndec=3)+")' > "+refcattemp
    file_delete,refcattemp,/allow
    file_delete,refcatfile,/allow
    spawn,cmd,out,outerr
    ;  Load ASCII file and create the FITS file                                                                                                                                     
    ref = importascii(refcattemp,/header,delim='|',skipline=2,/silent)
    if keyword_set(saveref) then MWRFITS,ref,refcatfile,/create      ; only save if necessary                                                                                      
    file_delete,refcattemp,/allow
    ;(SCOPE_VARFETCH(varname,/enter)) = ref

  ; Use QUERYVIZIER                                                                                                                                                                
  ;   for low density with 2MASS/GAIA and always for GALEX and APASS                                                                                                               
  endif else begin
    ;ref = QUERYVIZIER(refcat[i],[cenra,cendec],[rarange*1.1*60,decrange*1.1*60],/cfa)                                                                                             
    if refcat[i] eq 'APASS' then cfa=0 else cfa=1  ; cfa doesn't have APASS                                                                                                        
    ref = QUERYVIZIER(refcat[i],[cenra,cendec],radius*60,cfa=cfa)

    ; Fix/homogenize the GAIA tags                                                                                                                                                 
    if varname eq 'GAIA' then begin
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
    if varname eq 'TMASS' then begin
      nref = n_elements(ref)
      orig = ref
      ref = replicate({designation:'',raj2000:0.0d0,dej2000:0.0d0,jmag:0.0,e_jmag:0.0,hmag:0.0,e_hmag:0.0,kmag:0.0,e_kmag:0.0,qflg:''},nref)
      struct_assign,orig,ref
      ref.designation = orig._2mass
      undefine,orig
    endif

    ;(SCOPE_VARFETCH(varname,/enter)) = ref
    ; Save the file                                                                                                                                                                
    if keyword_set(saveref) then MWRFITS,ref,refcatfile,/create  ; only save if necessary                                                                                          
  endelse
endelse

return,ref

stop

end
