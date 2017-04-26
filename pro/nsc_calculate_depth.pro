pro nsc_calculate_depth

; Calculcate the 5 sigma depth for each exposure

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
  
; Restore the calibration summary file
print,'Loading nsc_instcal_calibrate.fits'
str = MRDFITS(dir+'nsc_instcal_calibrate.fits',1,/silent)
str.expdir = strtrim(str.expdir,2)
str.instrument = strtrim(str.instrument,2)
str.metafile = strtrim(str.metafile,2)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)
str.filter = strtrim(str.filter,2)
nstr = n_elements(str)

; Add depth column
print,'Adding depth column'
schema_str = str[0]
struct_assign,{dum:''},schema_str
schema_str = create_struct(schema_str,'depth95',99.99,'depth10sig',99.99)
old = str
str = replicate(schema_str,nstr)
struct_assign,old,str,/nozero
undefine,old

For i=0,nstr-1 do begin
  if (i+1) mod 1000 eq 0 then print,i+1

  catfile = str[i].expdir+'/'+str[i].base+'_cat.fits'
  ;if file_test(catfile) and str[i].success eq 1 then begin
  if str[i].success eq 1 then begin
    cat = MRDFITS(catfile,1,/silent,status=status)  
    if status lt 0 then goto,BOMB

    ; Need good photometry
    gdmag = where(cat.cmag lt 50,ngdmag)
    if ngdmag eq 0 then goto,BOMB

    ; Measure the depth
    ; Get 95% percentile depth
    cmag = cat[gdmag].cmag
    si = sort(cmag)
    cmag = cmag[si]
    depth95 = cmag[round(0.95*ngdmag)]
    str[i].depth95 = depth95
    ; Get 10 sigma depth
    ;  S/N = 1.087/err
    ;  so S/N=5 is for err=1.087/5=0.2174
    ;  S/N=10 is for err=1.087/10=0.1087
    depth10sig = 99.99
    depind = where(cat.cmag lt 50 and cat.cerr ge 0.0987 and cat.cerr le 0.1187,ndepind)
    if ndepind lt 5 then depind = where(cat.cmag lt 50 and cat.cerr ge 0.0787 and cat.cerr le 0.1387,ndepind)
    if ndepind gt 5 then begin
      depth10sig = median([cat[depind].cmag])
    endif else begin
      depind = where(cat.cmag lt 50,ndepind)
      if ndepind gt 0 then depth10sig=max([cat[depind].cmag])
    endelse
    str[i].depth10sig = depth10sig
  endif

  BOMB:

Endfor

;print,'NEED TO SAVE THE FINAL FILE!!!'
MWRFITS,str,dir+'nsc_instcal_calibrate_depth.fits',/create

stop

end
