;+
;
; GETMODELMAG
;
; This calculates the model magnitudes for the NSC catalog
; given a catalog with the appropriate information
;
; INPUTS:
;  cat        Structure of source with appropriate magnitude
;               columns.
;  filter     Short filter name, e.g. 'g' or 'c4d-g'.
;
; OUTPUTS:
;  model_mag  The model magnitudes for each source.
;
; USAGE:
;  IDL>model_mag = getmodelmag(cat,'g')
;
; By D. Nidever  Sep 2017
;-

function getmodelmag,cat,filter

; This calculates the model magnitude for stars given the
; the magnitudes in reference catalogs
; NUV - Galex NUV magnitude
; GMAG - Gaia G magnitude
; JMAG - 2MASS J magnitude
; KMAG - 2MASS Ks magnitude
; APASS_GMAG - APASS g magnitue
; APASS_RMAG - APASS r magnitude
; EBV  - E(B-V) reddening

; Not enough inputs
if n_elements(cat) eq 0 or n_elements(filter) eq 0 then begin
  print,'Syntax - model_mag = getmodelmag(cat,filter)'
  return,-1
endif

; If no instrument is given then assume DECam
instfilt = filter
if strpos(filter,'-') eq -1 then instfilt='c4d-'+filter

tags = tag_names(cat)
CASE instfilt of
; ---- DECam u-band ----
'c4d-u': begin
  ; Double-check that we have the right columns
  needtags = ['GMAG','JMAG','NUV','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.gmag lt 50 and cat.jmag lt 50 and cat.nuv lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (G-J)o = G-J-1.12*EBV
    col = cat[gd].gmag - cat[gd].jmag - 1.12*cat[gd].ebv
    ; u = 0.2469*NUV + 0.7501*G + 0.5462*GJ0 + 0.6809*EBV + 0.0052  ; for 0.8<GJ0<1.1
    model_mag[gd] = 0.2469*cat[gd].nuv + 0.7501*cat[gd].gmag + 0.5462*col + 0.6809*cat[gd].ebv + 0.0052
  endif
end

;---- DECam g-band ----
'c4d-g': begin
  ; Double-check that we have the right columns
  needtags = ['JMAG','KMAG','APASS_GMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.jmag lt 50 and cat.kmag lt 50 and cat.apass_gmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (J-Ks)o = J-Ks-0.17*EBV
    col = cat[gd].jmag-cat[gd].kmag-0.17*cat[gd].ebv
    ; g = APASS_G - 0.0421*JK0 - 0.05*EBV - 0.0620
    model_mag[gd] = cat[gd].apass_gmag - 0.0421*col - 0.05*cat[gd].ebv - 0.0620
  endif
end

; ---- DECam r-band ----
'c4d-r': begin
  ; Double-check that we have the right columns
  needtags = ['JMAG','KMAG','APASS_RMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.jmag lt 50 and cat.kmag lt 50 and cat.apass_rmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (J-Ks)o = J-Ks-0.17*EBV
    col = cat[gd].jmag-cat[gd].kmag-0.17*cat[gd].ebv
    ; r = APASS_r - 0.0861884*JK0 + 0.0*EBV + 0.0548607
    model_mag[gd] = cat[gd].apass_rmag - 0.0861884*col + 0.0548607
  endif
end

; ---- DECam i-band ----
'c4d-i': begin
  ; Double-check that we have the right columns
  needtags = ['JMAG','KMAG','GMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.jmag lt 50 and cat.kmag lt 50 and cat.gmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (J-Ks)o = J-Ks-0.17*EBV
    col = cat[gd].jmag-cat[gd].kmag-0.17*cat[gd].ebv
    ; i = G - 0.4587*JK0 - 0.276*EBV + 0.0967721
    model_mag[gd] = cat[gd].gmag - 0.4587*col - 0.276*cat[gd].ebv + 0.0967721
  endif
end

; ---- DECam z-band ----
'c4d-z': begin
  ; Double-check that we have the right columns
  needtags = ['JMAG','KMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.jmag lt 50 and cat.kmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (J-Ks)o = J-Ks-0.17*EBV
    col = cat[gd].jmag-cat[gd].kmag-0.17*cat[gd].ebv
    ; z = J + 0.765720*JK0 + 0.40*EBV +  0.605658
    model_mag[gd] = cat[gd].jmag + 0.765720*col + 0.40*cat[gd].ebv +  0.605658
  endif
end

; ---- DECam Y-band ----
'c4d-Y': begin
  ; Double-check that we have the right columns
  needtags = ['JMAG','KMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.jmag lt 50 and cat.kmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then begin
    ; (J-Ks)o = J-Ks-0.17*EBV
    col = cat[gd].jmag-cat[gd].kmag-0.17*cat[gd].ebv
    ; Y = J + 0.54482*JK0 + 0.20*EBV + 0.663380
    model_mag[gd] = cat[gd].jmag + 0.54482*col + 0.20*cat[gd].ebv + 0.663380
  endif
end

; ---- DECam VR-band ----
'c4d-VR': begin
  ; Double-check that we have the right columns
  needtags = ['GMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.gmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then model_mag[gd] = cat[gd].gmag
end


; ---- Bok+90Prime g-band ----
'ksb-g': begin
  ; Double-check that we have the right columns
  needtags = ['PS_GMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.ps_gmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then model_mag[gd] = cat[gd].ps_gmag
end

; ---- Bok+90Prime r-band ----
'ksb-r': begin  ; Double-check that we have the right columns
  needtags = ['PS_RMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.ps_rmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then model_mag[gd] = cat[gd].ps_rmag
end

; ---- Mosaic3 z-band ----
'k4m-z': begin
  ; Double-check that we have the right columns
  needtags = ['PS_ZMAG','EBV']
  MATCH,tags,needtags,ind1,ind2,/sort,count=ntagmatch  
  if ntagmatch lt n_elements(needtags) then begin
    leftind = indgen(n_elements(needtags))
    if ntagmatch gt 0 then remove,ind2,leftind
    print,'Needed columns missing. '+strjoin(needtags[leftind],' ')
  endif
  ; Selected sources with the needed information
  gd = where(cat.ps_zmag lt 50,ngd)
  model_mag = fltarr(n_elements(cat))+99.99
  if ngd gt 0 then model_mag[gd] = cat[gd].ps_zmag
end
else: stop,filter+' not supported'
ENDCASE

return,model_mag

end
