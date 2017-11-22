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
;  IDL>cat = getrefcat('g',cenra,cendec,radius,saveref=saveref)
;
; By D. Nidever  Sep 2017
;-

function getrefdata,filter,cenra,cendec,radius,count=count,saveref=saveref,silent=silent,$
                    dcr=dcr,modelmags=modelmags,logfile=logfile

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

if not keyword_set(silent) then begin
  printlog,logf,'Getting reference catalogs for:'
  printlog,logf,'FILTER(S) = ',strjoin(filter,', ')
  printlog,logf,'CENRA  = ',strtrim(cenra,2)
  printlog,logf,'CENDEC = ',strtrim(cendec,2)
  printlog,logf,'RADIUS = ',strtrim(radius,2),' deg'
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
      ; Use 2MASS and APASS to calibrate
      push,refcat,['2MASS-PSC','APASS']
    endelse
  end
  ; DECam r-band
  'c4d-r': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ; Use 2MASS and APASS to calibrate
      push,refcat,['2MASS-PSC','APASS']
    endelse
  end
  ; DECam i-band
  'c4d-i': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ; Use GAIA and 2MASS to calibrate
      push,refcat,['2MASS-PSC']
    endelse
  end
  ; DECam z-band
  'c4d-z': begin
    ; Use PS1 if possible
    if cendec gt -29 then begin
      push,refcat,['2MASS-PSC','PS']
    endif else begin
      ; Use GAIA and 2MASS to calibrate  
      push,refcat,['2MASS-PSC']
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
; Some catalogs
if n_elements(refcat) gt 0 then begin
  ; Get the unique ones
  ui = uniq(refcat,sort(refcat))
  refcat = refcat[ui]
  ; Always add Gaia at the beginning
  refcat = ['GAIA/GAIA',refcat]
endif else refcat='GAIA/GAIA'
nrefcat = n_elements(refcat)

; Figure out the new columns that need to be added
undefine,newtags
for i=0,nrefcat-1 do begin
  case refcat[i] of
  'GAIA/GAIA':   ; do nothing
  '2MASS-PSC': push,newtags,['jmag','e_jmag','kmag','e_kmag','qflg']
  'PS': push,newtags,['ps_gmag','ps_rmag','ps_zmag']
  'APASS': push,newtags,['apass_gmag','e_apass_gmag','apass_rmag','e_apass_rmag']
  'II/312/ais': push,newtags,['nuv','e_nuv']  ; Galex
  else: stop,refcat[i]+' NOT SUPPORTED'
  endcase
endfor
push,newtags,'ebv'
if keyword_set(modelmags) then push,newtags,'model_mag'
nnewtags = n_elements(newtags)

; Load the necessary catalogs
nrefcat = n_elements(refcat)
if not keyword_set(silent) then $
  printlog,logf,strtrim(nrefcat,2),' reference catalogs to load'
for i=0,nrefcat-1 do begin
  t0 = systime(1)
  if not keyword_set(silent) then $
    printlog,logf,'Loading ',refcat[i],' reference catalog'

  ; Load the catalog
  ref1 = GETREFCAT(cenra,cendec,radius,refcat[i],count=nref1,silent=silent,logfile=logf)
  if nref1 eq 0 then goto,BOMB
  tags1 = tag_names(ref1)


  ; First successful one, initialize the catalog
  if n_elements(ref) eq 0 then begin
    ; New format, add tags
    schema = ref1[0]
    struct_assign,{dum:''},schema
    schema = create_struct(schema,'ra',999999d0,'dec',999999.0d0)
    for j=0,nnewtags-1 do begin
      val0 = 99.99
      if newtags[j] eq 'qflg' then val0=''
      schema = create_struct(schema,newtags[j],val0)
    endfor
    ref = replicate(schema,n_elements(ref1))
    struct_assign,ref1,ref,/nozero
    ref.ra = ref.ra_icrs
    ref.dec = ref.de_icrs

  ; Crossmatch and add magnitudes
  endif else begin

    ; Get RA/DEC columns
    ; 2MASS, Galex, APASS use RAJ2000
    ; PS uses RA/DEC
    if refcat[i] ne 'PS' then raind=where(tags1 eq 'RAJ2000',nraind) else $
       raind=where(tags1 eq 'RA',nraind)
    if refcat[i] ne 'PS' then decind=where(tags1 eq 'DEJ2000',ndecind) else $
       decind=where(tags1 eq 'DEC',ndecind)

    ; Crossmatch
    SRCMATCH,ref.ra,ref.dec,ref1.(raind),ref1.(decind),dcr,ind1,ind2,/sph,count=nmatch
    if not keyword_set(silent) then $
      printlog,logf,strtrim(nmatch,2)+' matches'

    ; Add magnitude columns
    if nmatch gt 0 then begin
      case refcat[i] of
      '2MASS-PSC': begin
         ref[ind1].jmag = ref1[ind2].jmag
         ref[ind1].e_jmag = ref1[ind2].e_jmag
         ref[ind1].kmag = ref1[ind2].kmag
         ref[ind1].e_kmag = ref1[ind2].e_kmag
         ref[ind1].qflg = ref1[ind2].qflg
      end
      'PS': begin
         ref[ind1].ps_gmag = ref1[ind2].gmag
         ref[ind1].ps_rmag = ref1[ind2].rmag
         ref[ind1].ps_zmag = ref1[ind2].zmag
      end
      'APASS': begin
         ref[ind1].apass_gmag = ref1[ind2].g_mag
         ref[ind1].e_apass_gmag = ref1[ind2].e_g_mag
         ref[ind1].apass_rmag = ref1[ind2].r_mag
         ref[ind1].e_apass_rmag = ref1[ind2].e_r_mag
      end
      'II/312/ais': begin
         ref[ind1].nuv = ref1[ind2].nuv
         ref[ind1].e_nuv = ref1[ind2].e_nuv
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
      '2MASS-PSC': begin
         new.jmag = left1.jmag
         new.kmag = left1.kmag
      end
      'PS': begin
         new.ps_gmag = left1.gmag
         new.ps_rmag = left1.rmag
         new.ps_zmag = left1.zmag
      end
      'APASS': begin
         new.apass_gmag = left1.g_mag
         new.apass_rmag = left1.r_mag
      end
      'II/312/ais': new.nuv = left1.nuv
      else: stop,catname+' NOT SUPPORTED'
      endcase
      
      ; Combine the two
      old = ref
      ref = replicate(schema,n_elements(old)+nleft1)
      ref[0:n_elements(old)-1] = old
      ref[n_elements(old):*] = new
      undefine,old,new,left1
    endif

  endelse  ; second reference catalog
  BOMB:
endfor

; Add reddening
glactc,ref.ra,ref.dec,2000.0,glon,glat,1,/deg
ebv = dust_getval(glon,glat,/noloop,/interp)
ref.ebv = ebv

; Get the model magnitudes
if keyword_set(modelmags) then begin
  model_mag = GETMODELMAG(ref,filter)
  ref.model_mag = model_mag
  gmodel = where(ref.model_mag lt 50,ngmodel)
  if not keyword_set(silent) then $
    printlog,logf,strtrim(ngmodel,2)+' stars with good model magnitudes'
endif

count = n_elements(ref)

return,ref

end
