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
;  filter     Short filter name, e.g. 'g'.
;  dec        The declination of the exposure.
;  eqnfile    File with the model magnitude equations.
;
; OUTPUTS:
;  model_mag  An [Nsource,3] array with model magnitudes, errors and color.
;
; USAGE:
;  IDL>model_mag = getmodelmag(cat,'g',-50.0,'modelmag_equations.txt')
;
; By D. Nidever  Feb 2019
;-

function getmodelmag,cat,filter,dec,eqnfile

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
ncat = n_elements(cat)
if ncat eq 0 or n_elements(filter) eq 0 or n_elements(dec) eq 0 or n_elements(eqnfile) eq 0 then begin
  print,'Syntax - model_mag = getmodelmag(cat,filter,dec,eqnfile)'
  return,-1
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
  error = !ERROR_STATE.MSG
  PHOTRED_ERRORMSG,logfile=logf
  CATCH, /CANCEL 
  return,-999999.
endif

; If no instrument is given then assume DECam
instfilt = filter

tags = tag_names(cat)

;; Load the model magnitude equation information
;; band, dec range, color equation, color min/max range, quality cuts, model mag equation
if file_test(eqnfile) eq 0 then begin
  print,eqnfile,' NOT FOUND'
  return,-999999.
endif
eqnstr = IMPORTASCII(eqnfile,/header,/silent)
neqn = n_elements(eqnstr)
;; Get COLOR range
dum = strsplitter(eqnstr.colorange,',',/extract)
dum = repstr(repstr(dum,'[',''),']','')
add_tag,eqnstr,'colorlim',[0.0,0.0],eqnstr
eqnstr.colorlim[0] = float(reform(dum[0,*]))
eqnstr.colorlim[1] = float(reform(dum[1,*]))
;; Get DEC range
dum = strsplitter(eqnstr.decrange,',',/extract)
dum = repstr(repstr(dum,'[',''),']','')
add_tag,eqnstr,'declim',[0.0,0.0],eqnstr
eqnstr.declim[0] = float(reform(dum[0,*]))
eqnstr.declim[1] = float(reform(dum[1,*]))



;; Get the band for this FILTER and DEC.
gd = where(strtrim(eqnstr.instrument,2)+'-'+strtrim(eqnstr.band,2) eq filter and dec ge eqnstr.declim[0] and dec le eqnstr.declim[1],ngd)
if ngd eq 0 then begin
  print,'No model magnitude equation for FILTER='+filter+' and DEC='+stringize(dec,ndec=2)
  return,-999999.
endif
if ngd gt 1 then begin
  print,'Found multiple magnitude equation for FILTER='+filter+' and DEC='+stringize(dec,ndec=2)+'.  Using the first one.'
  gd = gd[0]
endif
eqnstr1 = eqnstr[gd]

;; No parentheses allowed
if strpos(eqnstr1.coloreqn,'(') ne -1 or strpos(eqnstr1.coloreqn,')') ne -1 or $
   strpos(eqnstr1.modelmageqn,'(') ne -1 or strpos(eqnstr1.modelmageqn,')') ne -1 then begin
  print,'No parentheses allowed in the model magnitude equations'
  return,-999999.
endif

;; Are we using color?
if eqnstr1.colorlim[0] lt -10 and eqnstr1.colorlim[1] gt 10 and $
   strpos(strupcase(eqnstr1.modelmageqn),'COLOR') eq -1 then usecolor=0 else usecolor=1

;; Get all columns that we need
coloreqn = eqnstr1.coloreqn
qualitycuts = eqnstr1.qualitycuts
modelmageqn = eqnstr1.modelmageqn
coloreqn_cols = strsplit(coloreqn,'-+*',/extract)
modelmageqn_cols = strsplit(modelmageqn,'-+*',/extract)
if keyword_set(usecolor) then cols = strupcase([coloreqn_cols, modelmageqn_cols]) else cols=strupcase(modelmageqn_cols)
;; Remove numbers and "COLOR"
bd = where(valid_num(cols) eq 1 or strupcase(cols) eq 'COLOR' or cols eq '??',nbd)
if nbd gt 0 then begin
  if nbd lt n_elements(cols) then REMOVE,bd,cols else undefine,cols
endif
ncols = n_elements(cols)
;; No columns left
if ncols eq 0 then begin
  print,'No columns to use.'
  return,-999999.
endif
;; Only unique columns
ui = uniq(cols,sort(cols))
cols = cols[ui]
ncols = n_elements(cols)
MATCH,tags,cols,ind1,ind2,/sort,count=ntagmatch  
if ntagmatch lt ncols then begin
  leftind = indgen(ncols)
  if ntagmatch gt 0 then remove,ind2,leftind
  print,'Needed columns missing. '+strjoin(cols[leftind],' ')
endif

;; Make the color
;;  replace the columns by CAT[GD].COLUMN
if keyword_set(usecolor) then begin
  coloreqn_cols = strupcase(strsplit(coloreqn,'-+*',/extract))
  coloreqn_cols = coloreqn_cols[uniq(coloreqn_cols,sort(coloreqn_cols))]  ; unique ones
  bd = where(valid_num(coloreqn_cols) eq 1,nbd)  ;; Remove numbers and "COLOR"
  if nbd gt 0 then REMOVE,bd,coloreqn_cols
  colcmd = strupcase(coloreqn)
  for i=0,n_elements(coloreqn_cols)-1 do colcmd=repstr(colcmd,coloreqn_cols[i],'cat.'+coloreqn_cols[i])
  dum = EXECUTE('color='+colcmd)
endif else color=fltarr(ncat)

;; Make quality cuts
magcolsind = where((stregex(cols,'mag$',/boolean,/fold_case) eq 1 and stregex(cols,'^e_',/boolean,/fold_case) eq 0) or $
                   cols eq 'NUV',nmagcolsind)
;; make sure all magnitudes are good (<50) and finite
goodmask = bytarr(ncat)+1
for i=0,nmagcolsind-1 do begin
  magind = where(strupcase(tags) eq strupcase(cols[magcolsind[i]]),nmagind)
  goodmask AND= (cat.(magind[0]) lt 50 and cat.(magind[0]) gt 0 and finite(cat.(magind[0])) eq 1)
endfor
;; input quality cuts
;;  replace <=, <, >, >=, =, &, |
qualitycuts = repstr(qualitycuts,'<=',' le ')
qualitycuts = repstr(qualitycuts,'>=',' ge ')
qualitycuts = repstr(qualitycuts,'>',' gt ')
qualitycuts = repstr(qualitycuts,'<',' lt ')
qualitycuts = repstr(qualitycuts,'=',' eq ')
qualitycuts = repstr(qualitycuts,'&',' and ')
qualitycuts = repstr(qualitycuts,'|',' or ')
;; fix columns names
qualitycuts_cols = strsplit(qualitycuts,' ',/extract)
for i=0,n_elements(qualitycuts_cols)-1 do begin
  col = qualitycuts_cols[i]
  colind = where(tags eq strupcase(col),ncolind)
  if colind gt 0 then qualitycuts = repstr(qualitycuts,col,'cat.'+col)
endfor
dum = EXECUTE('goodmask AND= '+qualitycuts)
;; Apply the color range
if keyword_set(usecolor) then goodmask AND= (color ge eqnstr1.colorlim[0] and color le eqnstr1.colorlim[1])
;; Get the sources that pass all cuts
gd = where(goodmask eq 1,ngd)
if ngd eq 0 then begin
  print,'No good sources left'
  return,-999999.
endif

; Make the model magnitude
;;  replace the columns by CAT[GD].COLUMN
modelmageqn_cols = strupcase(strsplit(modelmageqn,'-+*',/extract))
bd = where(valid_num(modelmageqn_cols) eq 1 or strupcase(modelmageqn_cols) eq 'COLOR',nbd)  ;; Remove numbers and "COLOR"
if nbd gt 0 then REMOVE,bd,modelmageqn_cols
modelmageqn_cols = modelmageqn_cols[uniq(modelmageqn_cols,sort(modelmageqn_cols))]  ; unique ones
magcmd = strupcase(modelmageqn)
for i=0,n_elements(modelmageqn_cols)-1 do magcmd=repstr(magcmd,modelmageqn_cols[i],'cat[gd].'+modelmageqn_cols[i])
magcmd = repstr(magcmd,'COLOR','COLOR[gd]')
dum = EXECUTE('modelmag_gd='+magcmd)
modelmag = fltarr(ncat)+99.99
modelmag[gd] = modelmag_gd

;; Make the error structure
;;  Each magnitude has an E_MAG error except for PS and Gaia GMAG
;; If we are using PS or GMAG then add the errors for the
undefine,adderrtags
if (where(cols eq 'GMAG'))[0] ne -1 then push,adderrtags,'E_GMAG'
psmagind = where(stregex(cols,'^PS_',/boolean) eq 1 and stregex(cols,'MAG$',/boolean) eq 1,npsmagind)
for i=0,npsmagind-1 do push,adderrtags,'E_'+cols[psmagind]
nadderrtags = n_elements(adderrtags)
;; Making error structure
errtagind = where(stregex(tags,'^E_',/boolean) eq 1,nerrtags)
errtags = tags[errtagind]
errschema = create_struct(errtags[0],0.001)
if nerrtags gt 1 then for i=1,nerrtags-1 do errschema=create_struct(errschema,errtags[i],0.001)
if nadderrtags gt 0 then for i=0,nadderrtags-1 do errschema=create_struct(errschema,adderrtags[i],0.001)
err = replicate(errschema,ncat)
struct_assign,cat,err,/nozero
if (where(cols eq 'GMAG'))[0] ne -1 then err.e_gmag = 2.5*alog10(1.0+cat.e_fg/cat.fg)
;; leave the PS errors at 0.001
;; convert NAN or 99.99 to 9.99 to be consistent
for i=0,n_tags(err)-1 do begin
  bd = where(err.(i) gt 10.0 or finite(err.(i)) eq 0,nbd)
  if nbd gt 0 then err[bd].(i) = 9.99
endfor

;; Calculate the color errors
;; get the columns
if keyword_set(usecolor) then begin
  colorerr_cols = strupcase(strsplit(coloreqn,'-+*',/extract))
  colorerr_cols = colorerr_cols[uniq(colorerr_cols,sort(colorerr_cols))]  ; unique ones
  bd = where(valid_num(colorerr_cols) eq 1 or strupcase(colorerr_cols) eq 'EBV',nbd)  ;; Remove numbers and "EBV"
  if nbd gt 0 then REMOVE,bd,colorerr_cols
  ;; use - and + signs to break apart the components that need to be  squared
  coloreqn_terms = strupcase(strsplit(coloreqn,'-+',/extract))
  ;; remove any terms that don't have a COLORERR_COLS in them
  okay = bytarr(n_elements(coloreqn_terms))
  for i=0,n_elements(coloreqn_terms)-1 do $
    for j=0,n_elements(colorerr_cols)-1 do okay[i] OR= stregex(coloreqn_terms[i],colorerr_cols[j],/boolean)
  bd = where(okay eq 0,nbd)
  if nbd gt 0 then REMOVE,bd,coloreqn_terms
  ;; Now create the equation, add in quadrature
  colorerrcmd = 'sqrt( '+strjoin('('+coloreqn_terms+')^2','+')+' )'
  colorerrcmd = strupcase(colorerrcmd)
  for i=0,n_elements(colorerr_cols)-1 do colorerrcmd=repstr(colorerrcmd,colorerr_cols[i],'err[gd].e_'+colorerr_cols[i])
  dum = EXECUTE('colorerr_gd='+colorerrcmd)
  colorerr = fltarr(ncat)+9.99
  colorerr[gd] = colorerr_gd
endif else colorerr=fltarr(ncat)

;; The modelmag errors
;; get the columns
modelmagerr_cols = strupcase(strsplit(modelmageqn,'-+*',/extract))
modelmagerr_cols = modelmagerr_cols[uniq(modelmagerr_cols,sort(modelmagerr_cols))]  ; unique ones
bd = where(valid_num(modelmagerr_cols) eq 1 or strupcase(modelmagerr_cols) eq 'EBV',nbd)  ;; Remove numbers and "EBV"
if nbd gt 0 then REMOVE,bd,modelmagerr_cols
;;   use - and + signs to break apart the components that need to be  squared
modelmageqn_terms = strupcase(strsplit(modelmageqn,'-+',/extract))
;; remove any terms that don't have a COLORERR_COLS in them
okay = bytarr(n_elements(modelmageqn_terms))
for i=0,n_elements(modelmageqn_terms)-1 do $
  for j=0,n_elements(modelmagerr_cols)-1 do okay[i] OR= stregex(modelmageqn_terms[i],modelmagerr_cols[j],/boolean)
bd = where(okay eq 0,nbd)
if nbd gt 0 then REMOVE,bd,modelmageqn_terms
;; Now create the equation, add in quadrature 
modelmagerrcmd = 'sqrt( '+strjoin('('+modelmageqn_terms+')^2','+')+' )'
modelmagerrcmd = strupcase(modelmagerrcmd)
for i=0,n_elements(modelmageqn_cols)-1 do modelmagerrcmd=repstr(modelmagerrcmd,modelmageqn_cols[i],'err[gd].e_'+modelmageqn_cols[i])
modelmagerrcmd = repstr(modelmagerrcmd,'COLOR','COLORERR[gd]')
dum = EXECUTE('modelmagerr_gd='+modelmagerrcmd)
modelmagerr = fltarr(ncat)+9.99
modelmagerr[gd] = modelmagerr_gd

;; combine modelmag and modelmagerr
;; Combine mags and errors
mags = fltarr(n_elements(cat),3)
mags[*,0] = modelmag
mags[*,1] = modelmagerr
mags[*,2] = color

return,mags

end
