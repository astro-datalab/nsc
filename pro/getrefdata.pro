;+
;
; GETREFDATA
;
; Get reference catalog information needed for a given filter.
;
; INPUTS:
;  filter    Filter short name, e.g. "u".
;  cenra     Central RA for the search.
;  cendec    Central DEC for the search.
;  radius    Search radius in degrees.
;  /saveref  Save the output to FILE.
;  /silent   Don't print anything to the screen.
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

function getrefdata,filter,cenra,cendec,radius,count=count,saveref=saveref,silent=silent

undefine,ref
count = 0

; Not enough inputs
if n_elements(filter) eq 0 or n_elements(cenra) eq 0 or n_elements(cendec) eq 0 or $
   n_elements(radius) eq 0 then begin
  print,'Syntax - cat = getrefdata(filter,cenra,cendec,radius,saveref=saveref)'
  return,-1
endif

; Check that we have psql installed
spawn,['which','psql'],out,errout,/noshell
if file_test(out[0]) eq 0 then begin
  print,'No PSQL found on this sytem.'
  return,-1
endif

logf = -1

; If no instrument information then assume it's DECam
instfilt = filter
if strpos(filter,'-') eq -1 then instfilt='c4d-'+filter

; Load the reference catalogs
;-----------------------------
refcat = 'GAIA/GAIA'                ; always need this one
CASE instfilt of
; DECam u-band
'c4d-u': begin
  ; Use GAIA, 2MASS and GALEX to calibrate
  push,refcat,['2MASS-PSC','II/312/ais']
end
; DECam g-band
'c4d-g': begin
  ; Use 2MASS and APASS to calibrate
  push,refcat,['2MASS-PSC','APASS']
end
; DECam r-band
'c4d-r': begin
  ; Use 2MASS and APASS to calibrate
  push,refcat,['2MASS-PSC','APASS']
end
; DECam i-band
'c4d-i': begin
  ; Use GAIA and 2MASS to calibrate
  push,refcat,['2MASS-PSC']
end
; DECam z-band
'c4d-z': begin
  ; Use GAIA and 2MASS to calibrate  
  push,refcat,['2MASS-PSC']
end
; DECam Y-band
'c4d-Y': begin
  ; Use 2MASS to calibrate
  push,refcat,['2MASS-PSC']
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
  return,-1
end
ENDCASE
nrefcat = n_elements(refcat)

; Figure out the new columns that need to be added
undefine,newtags
for i=0,nrefcat-1 do begin
  case refcat[i] of
  'GAIA/GAIA':   ; do nothing
  '2MASS-PSC': push,newtags,['jmag','kmag']
  'PS': push,newtags,['ps_gmag','ps_rmag','ps_zmag']
  'APASS': push,newtags,['apass_gmag','apass_rmag']
  'II/312/ais': push,newtags,'nuv'  ; Galex
  else: stop,refcat[i]+' NOT SUPPORTED'
  endcase
endfor
push,newtags,'ebv'
push,newtags,'model_mag'
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
  ref1 = GETREFCAT(cenra,cendec,radius,refcat[i],count=nref1,silent=silent)
  tags1 = tag_names(ref1)



  ; First one, initialize the catalog
  if i eq 0 then begin
    ; New format, add tags
    schema = ref1[0]
    struct_assign,{dum:''},schema
    schema = create_struct(schema,'ra',0.0d0,'dec',0.0d0)
    for j=0,nnewtags-1 do schema=create_struct(schema,newtags[j],99.99)
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
       raind=where(tags1 eq 'DEC',ndecind)

    ; Crossmatch
    dcr = 0.5  ; arcsec
    SRCMATCH,ref.ra,ref.dec,ref1.(raind),ref1.(decind),dcr,ind1,ind2,/sph,count=nmatch
    if not keyword_set(silent) then $
      printlog,logf,strtrim(nmatch,2)+' matches'

    ; Add magnitude columns
    if nmatch gt 0 then begin
      case refcat[i] of
      '2MASS-PSC': begin
         ref[ind1].jmag = ref1[ind2].jmag
         ref[ind1].kmag = ref1[ind2].kmag
      end
      'PS': begin
         ref[ind1].ps_gmag = ref1[ind2].gmag
         ref[ind1].ps_rmag = ref1[ind2].rmag
         ref[ind1].ps_zmag = ref1[ind2].zmag
      end
      'APASS': begin
         ref[ind1].apass_gmag = ref1[ind2].g_mag
         ref[ind1].apass_rmag = ref1[ind2].r_mag
      end
      'II/312/ais': ref[ind1].nuv = ref1[ind2].nuv
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
endfor

; Add reddening
glactc,ref.ra,ref.dec,2000.0,glon,glat,1,/deg
ebv = dust_getval(glon,glat,/noloop,/interp)
ref.ebv = ebv

; Get the model magnitudes
model_mag = GETMODELMAG(ref,filter)
ref.model_mag = model_mag
gmodel = where(ref.model_mag lt 50,ngmodel)
if not keyword_set(silent) then $
  printlog,logf,strtrim(ngmodel,2)+' stars with good model magnitudes'

count = n_elements(ref)

return,ref

end
