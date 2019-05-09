;+
;
; GETREFDATA
;
; Get reference catalog information needed for a given filter.
;
; INPUTS:
;  filter    Filter and instrument name, e.g. "c4d-u".
;              This can be an array of multiple filters and
;              then all the reference needed for all the filters
;              will be included.
;  cenra     Central RA for the search.
;  cendec    Central DEC for the search.
;  radius    Search radius in degrees.
;  =dcr      The cross-matching radius in arcsec.  Default is 0.5".
;  /saveref  Save the output to FILE.
;  /silent   Don't print anything to the screen.
;  /modelmags  Return the model magnitudes as well.
;  =logfile  Filename to write output to.
;
; OUTPUTS:
;  ref       Search results from the reference catalog all in one
;              structure.
;  =count    Number of elements in REF.
;
; USAGE:
;  IDL>cat = getrefdata('g',cenra,cendec,radius,saveref=saveref)
;
; By D. Nidever  Sep 2017
;-

function getrefdata,filter,cenra,cendec,radius,count=count,saveref=saveref,silent=silent,$
                    dcr=dcr,modelmags=modelmags,logfile=logfile

t0 = systime(1)
undefine,ref
count = 0

; Not enough inputs
if n_elements(filter) eq 0 or n_elements(cenra) eq 0 or n_elements(cendec) eq 0 or $
   n_elements(radius) eq 0 then begin
  print,'Syntax - cat = getrefdata(filter,cenra,cendec,radius,saveref=saveref,dcr=dcr,'
  print,'                          modelmags=modelmags,logfile=logfile)'
  return,-1
endif

; Defaults
if n_elements(dcr) eq 0 then dcr = 0.5  ; arcsec
if n_elements(logfile) eq 0 then logf = -1 else logf=logfile

; Check that we have psql installed
spawn,['which','psql'],out,errout,/noshell
if file_test(out[0]) eq 0 then begin
  print,'No PSQL found on this sytem.'
  return,-1
endif

;; Figure out what reddening method we are using
;;----------------------------------------------
;; Extinction types:
;; 1 - SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2
;; 2 - RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE
;;     data available
;; 3 - RJCE GLIMPSE, GLIMPSE data available
;; 4 - RJCE SAGE, SAGE data available
ext_type = GETEXTTYPE(cenra,cendec,radius)

;; If there is GLIMPSE or SAGE data, also get ALLWISE data
;;  in case there is only partial coverage


if not keyword_set(silent) then begin
  printlog,logf,'Getting reference catalogs for:'
  printlog,logf,'FILTER(S) = ',strjoin(filter,', ')
  printlog,logf,'CENRA  = ',strtrim(cenra,2)
  printlog,logf,'CENDEC = ',strtrim(cendec,2)
  printlog,logf,'RADIUS = ',strtrim(radius,2),' deg'
  case ext_type of
  1: printlog,logf,'Extinction Type: 1 - SFD'
  2: printlog,logf,'Extinction Type: 2 - RJCE ALLWISE'
  3: printlog,logf,'Extinction Type: 3 - RJCE GLIMPSE'
  4: printlog,logf,'Extinction Type: 4 - RJCE SAGE'
  else:
  endcase
endif

; Load the reference catalogs
;-----------------------------
; Figure out the reference catalogs that we need based on
;  filter-instrument combination
For i=0,n_elements(filter)-1 do begin
  instfilt = filter[i]
  ; If no instrument information then assume it's DECam
  if strpos(instfilt,'-') eq -1 then instfilt='c4d-'+instfilt
  ; Check the cases
  CASE instfilt of
  ; DECam u-band
  'c4d-u': begin
    ; Use GAIA, 2MASS and GALEX to calibrate
    push,refcat,['2MASS-PSC','II/312/ais']
  end
  ; DECam g-band
  'c4d-g': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ;; Use 2MASS and Skymapper to calibrate
      ;push,refcat,['2MASS-PSC','Skymapper']
      push,refcat,['2MASS-PSC','ATLAS']
    endelse
  end
  ; DECam r-band
  'c4d-r': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ;; Use 2MASS and Skymapper to calibrate
      ;push,refcat,['2MASS-PSC','Skymapper']
      push,refcat,['2MASS-PSC','ATLAS']
    endelse
  end
  ; DECam i-band
  'c4d-i': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ;; Use Skymapper and 2MASS to calibrate
      ;push,refcat,['2MASS-PSC','Skymapper']
      push,refcat,['2MASS-PSC','ATLAS']
    endelse
  end
  ; DECam z-band
  'c4d-z': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ;; Use Skymapper and 2MASS to calibrate
      ;push,refcat,['2MASS-PSC','Skymapper']
      push,refcat,['2MASS-PSC','ATLAS']
    endelse
  end
  ; DECam Y-band
  'c4d-Y': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ; Use 2MASS to calibrate
      push,refcat,['2MASS-PSC']
    endelse
  end
  ; DECam VR-band
  'c4d-VR': begin
    ; Use GAIA G-band to calibrate
    push,refcat,['2MASS-PSC']
  end
  ; Bok+90Prime g-band
  'ksb-g': begin
    ; Use PS1
    push,refcat,'PS'
  end
  ; Bok+90Prime r-band
  'ksb-r': begin
    ; Use PS1
    push,refcat,'PS'
  end
  ; Mosaic3 z-band
  'k4m-z': begin
    ; Use PS1
    push,refcat,'PS'
  end
  else: begin
    printlog,logf,filter,' not currently supported'
    ;return,-1
  end
  ENDCASE
Endfor ; filter loop
;; Extinction catalogs
if ext_type ge 2 then push,refcat,'ALLWISE'
if ext_type eq 3 then push,refcat,'GLIMPSE'
if ext_type eq 4 then push,refcat,'SAGE'
; Some catalogs
if n_elements(refcat) gt 0 then begin
  ; Get the unique ones
  ui = uniq(refcat,sort(refcat))
  refcat = refcat[ui]
  ; Always add Gaia at the beginning
  refcat = ['GAIADR2',refcat]
endif else refcat='GAIADR2'   ;'GAIA/GAIA'
nrefcat = n_elements(refcat)

; Figure out the new columns that need to be added
undefine,newtags
for i=0,nrefcat-1 do begin
  case refcat[i] of
  'GAIADR2': push,newtags,['source','ra','ra_error','dec','dec_error','pmra','pmra_error','pmdec','pmdecerr','gmag','e_gmag','bp','e_bp','rp','e_rp']
  '2MASS-PSC': push,newtags,['jmag','e_jmag','hmag','e_hmag','kmag','e_kmag','qflg']
  'PS': push,newtags,['ps_gmag','ps_rmag','ps_imag','ps_zmag','ps_ymag']
  'APASS': push,newtags,['apass_gmag','e_apass_gmag','apass_rmag','e_apass_rmag']
  'II/312/ais': push,newtags,['nuv','e_nuv']  ; Galex
  'Skymapper': push,newtags,['sm_gmag','e_sm_gmag','sm_rmag','e_sm_rmag','sm_imag','e_sm_imag','sm_zmag','e_sm_zmag']  ; Skymapper DR1
  'ALLWISE': push,newtags,['w1mag','e_w1mag','w2mag','e_w2mag']
  'GLIMPSE': push,newtags,['gl_36mag','e_gl_36mag','gl_45mag','e_gl_45mag']
  'SAGE': push,newtags,['sage_36mag','e_sage_36mag','sage_45mag','e_sage_45mag']
  'ATLAS': push,newtags,['atlas_gmag','e_atlas_gmag','atlas_rmag','e_atlas_rmag','atlas_imag','e_atlas_imag','atlas_zmag','e_atlas_zmag']
  else: stop,refcat[i]+' NOT SUPPORTED'
  endcase
endfor
push,newtags,['ebv_sfd','ejk','e_ejk','ext_type']
if keyword_set(modelmags) then push,newtags,'model_mag'
nnewtags = n_elements(newtags)

; Load the necessary catalogs
nrefcat = n_elements(refcat)
if not keyword_set(silent) then $
  printlog,logf,strtrim(nrefcat,2),' reference catalogs to load: '+strjoin(refcat,', ')
for i=0,nrefcat-1 do begin
  t0 = systime(1)
  if not keyword_set(silent) then $
    printlog,logf,'Loading ',refcat[i],' reference catalog'

  ; Load the catalog
  ref1 = GETREFCAT(cenra,cendec,radius,refcat[i],count=nref1,silent=silent,logfile=logf)
  if nref1 eq 0 then goto,BOMB
  tags1 = tag_names(ref1)


  ;; Initialize the catalog
  undefine,schema
  for j=0,nnewtags-1 do begin
    val0 = 99.99
    if newtags[j] eq 'qflg' then val0=''
    if newtags[j] eq 'ext_type' then val0=0
    if newtags[j] eq 'source' then val0=0LL
    if newtags[j] eq 'ra' then val0=0.0d0
    if newtags[j] eq 'dec' then val0=0.0d0
    if n_elements(schema) eq 0 then schema=create_struct(newtags[j],val0) else $
      schema = create_struct(schema,newtags[j],val0)
  endfor


  ;; First successful one, initialize the catalog
  if n_elements(ref) eq 0 then begin
    ref = replicate(schema,nref1)
    struct_assign,ref1,ref,/nozero
    if tag_exist(ref1,'RA') eq 0 then begin
      ref.ra = ref.ra_icrs
      ref.dec = ref.de_icrs
    endif
    ind1 = lindgen(nref1)
    ind2 = lindgen(nref1)
    nmatch = nref1

  ;; Second and later
 endif else begin

    ; Get RA/DEC columns
    ; 2MASS, Galex, APASS use RAJ2000
    ; PS uses RA/DEC
    if (refcat[i] ne 'PS' and refcat[i] ne 'ALLWISE' and refcat[i] ne 'ATLAS') then raind=where(tags1 eq 'RAJ2000',nraind) else $
       raind=where(tags1 eq 'RA',nraind)
    if (refcat[i] ne 'PS' and refcat[i] ne 'ALLWISE' and refcat[i] ne 'ATLAS') then decind=where(tags1 eq 'DEJ2000',ndecind) else $
       decind=where(tags1 eq 'DEC',ndecind)

    ; Crossmatch
    SRCMATCH,ref.ra,ref.dec,ref1.(raind),ref1.(decind),dcr,ind1,ind2,/sph,count=nmatch
    if not keyword_set(silent) then $
      printlog,logf,strtrim(nmatch,2)+' matches'
 endelse

  ; Add magnitude columns
  if nmatch gt 0 then begin
    case refcat[i] of
    'GAIADR2': begin
      temp = ref[ind1]
      struct_assign,ref1[ind2],temp,/nozero
      temp.e_gmag = 2.5*alog10(1.0+ref1[ind2].e_fg/ref1[ind2].fg)
      temp.e_bp = 2.5*alog10(1.0+ref1[ind2].e_fbp/ref1[ind2].fbp)
      temp.e_rp = 2.5*alog10(1.0+ref1[ind2].e_frp/ref1[ind2].frp)
      ref[ind1] = temp
    end
    '2MASS-PSC': begin
       ref[ind1].jmag = ref1[ind2].jmag
       ref[ind1].e_jmag = ref1[ind2].e_jmag
       ref[ind1].hmag = ref1[ind2].hmag
       ref[ind1].e_hmag = ref1[ind2].e_hmag
       ref[ind1].kmag = ref1[ind2].kmag
       ref[ind1].e_kmag = ref1[ind2].e_kmag
       ref[ind1].qflg = ref1[ind2].qflg
    end
    'PS': begin
       ref[ind1].ps_gmag = ref1[ind2].gmag
       ref[ind1].ps_rmag = ref1[ind2].rmag
       ref[ind1].ps_imag = ref1[ind2].imag
       ref[ind1].ps_zmag = ref1[ind2].zmag
       ref[ind1].ps_ymag = ref1[ind2].ymag
    end
    'APASS': begin
       ref[ind1].apass_gmag = ref1[ind2].g_mag
       ref[ind1].e_apass_gmag = ref1[ind2].e_g_mag
       ref[ind1].apass_rmag = ref1[ind2].r_mag
       ref[ind1].e_apass_rmag = ref1[ind2].e_r_mag
    end
    'II/312/ais': begin
       if tag_exist(ref1,'NUV') then begin
         ref[ind1].nuv = ref1[ind2].nuv
         ref[ind1].e_nuv = ref1[ind2].e_nuv
       endif else begin
         ref[ind1].nuv = ref1[ind2].nuvmag
         ref[ind1].e_nuv = ref1[ind2].e_nuvmag
       endelse
    end
    'Skymapper': begin
       ref[ind1].sm_gmag = ref1[ind2].sm_gmag
       ref[ind1].e_sm_gmag = ref1[ind2].e_sm_gmag
       ref[ind1].sm_rmag = ref1[ind2].sm_rmag
       ref[ind1].e_sm_rmag = ref1[ind2].e_sm_rmag
       ref[ind1].sm_imag = ref1[ind2].sm_imag
       ref[ind1].e_sm_imag = ref1[ind2].e_sm_imag
       ref[ind1].sm_zmag = ref1[ind2].sm_zmag
       ref[ind1].e_sm_zmag = ref1[ind2].e_sm_zmag
    end
    'ATLAS': begin
       ref[ind1].atlas_gmag = ref1[ind2].gmag
       ref[ind1].e_atlas_gmag = ref1[ind2].gerr
       ref[ind1].atlas_rmag = ref1[ind2].rmag
       ref[ind1].e_atlas_rmag = ref1[ind2].rerr
       ref[ind1].atlas_imag = ref1[ind2].imag
       ref[ind1].e_atlas_imag = ref1[ind2].ierr
       ref[ind1].atlas_zmag = ref1[ind2].zmag
       ref[ind1].e_atlas_zmag = ref1[ind2].zerr
    end
    'ALLWISE': begin
       ref[ind1].w1mag = ref1[ind2].w1mag
       ref[ind1].e_w1mag = ref1[ind2].e_w1mag
       ref[ind1].w2mag = ref1[ind2].w2mag
       ref[ind1].e_w2mag = ref1[ind2].e_w2mag
    end
    'GLIMPSE': begin
       ref[ind1].gl_36mag = ref1[ind2]._3_6mag
       ref[ind1].e_gl_36mag = ref1[ind2].e_3_6mag
       ref[ind1].gl_45mag = ref1[ind2]._4_5mag
       ref[ind1].e_gl_45mag = ref1[ind2].e_4_5mag
    end
    'SAGE': begin
       ref[ind1].sage_36mag = ref1[ind2].__3_6_
       ref[ind1].e_sage_36mag = ref1[ind2].e__3_6_
       ref[ind1].sage_45mag = ref1[ind2].__4_5_
       ref[ind1].e_sage_45mag = ref1[ind2].e__4_5_
    end
    else: stop,catname+' NOT SUPPORTED'
    endcase
  endif

  ; Add leftover ones
  if nmatch lt n_elements(ref1) then begin
    left1 = ref1
    remove,ind2,left1
    nleft1 = n_elements(left1)
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
    end
    '2MASS-PSC': begin
       new.jmag = left1.jmag
       new.e_jmag = left1.e_jmag
       new.hmag = left1.hmag
       new.e_hmag = left1.e_hmag
       new.kmag = left1.kmag
       new.e_kmag = left1.e_kmag
       new.qflg = left1.qflg
    end
    'PS': begin
       new.ps_gmag = left1.gmag
       new.ps_rmag = left1.rmag
       new.ps_imag = left1.imag
       new.ps_zmag = left1.zmag
       new.ps_ymag = left1.ymag
    end
    'APASS': begin
       new.apass_gmag = left1.g_mag
       new.e_apass_gmag = left1.e_g_mag
       new.apass_rmag = left1.r_mag
       new.e_apass_rmag = left1.e_r_mag
    end
    'II/312/ais': begin
       if tag_exist(left1,'NUV') then begin
         new.nuv = left1.nuv
         new.e_nuv = left1.e_nuv
       endif else begin
         new.nuv = left1.nuvmag
         new.e_nuv = left1.e_nuvmag
       endelse
    end
    'Skymapper': begin
       new.sm_gmag = left1.sm_gmag
       new.e_sm_gmag = left1.e_sm_gmag
       new.sm_rmag = left1.sm_rmag
       new.e_sm_rmag = left1.e_sm_rmag
       new.sm_imag = left1.sm_imag
       new.e_sm_imag = left1.e_sm_imag
       new.sm_zmag = left1.sm_zmag
       new.e_sm_zmag = left1.e_sm_zmag
    end
    'ATLAS': begin
       ref.atlas_gmag = left1.gmag
       ref.e_atlas_gmag = left1.gerr
       ref.atlas_rmag = left1.rmag
       ref.e_atlas_rmag = left1.rerr
       ref.atlas_imag = left1.imag
       ref.e_atlas_imag = left1.ierr
       ref.atlas_zmag = left1.zmag
       ref.e_atlas_zmag = left1.zerr
    end
    'ALLWISE': begin
       new.w1mag = left1.w1mag
       new.e_w1mag = left1.e_w1mag
       new.w2mag = left1.w2mag
       new.e_w2mag = left1.e_w2mag
    end
    'GLIMPSE': begin
       new.gl_36mag = left1._3_6mag
       new.e_gl_36mag = left1.e_3_6mag
       new.gl_45mag = left1._4_5mag
       new.e_gl_45mag = left1.e_4_5mag
    end
    'SAGE': begin
       new.sage_36mag = left1.__3_6_
       new.e_sage_36mag = left1.e__3_6_
       new.sage_45mag = left1.__4_5_
       new.e_sage_45mag = left1.e__4_5_
    end      
    else: stop,catname+' NOT SUPPORTED'
    endcase
      
    ; Combine the two
    old = ref
    ref = replicate(schema,n_elements(old)+nleft1)
    ref[0:n_elements(old)-1] = old
    ref[n_elements(old):*] = new
    undefine,old,new,left1
  endif

  BOMB:
endfor

;; Get extinction
;;----------------
GETREDDENING,ref,ext_type

count = n_elements(ref)

print,'dt=',stringize(systime(1)-t0,ndec=1),' sec'

return,ref

end
