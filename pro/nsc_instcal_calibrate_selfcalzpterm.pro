pro nsc_instcal_calibrate_selfcalzpterm,expdir,cat,expstr,chstr,logfile=logfile,silent=silent

; Measure the zero-point for this exposure using other
; successfully-calibrated exposures as secondary references 

NSC_ROOTDIRS,dldir,mssdir,localdir
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(expdir) eq 0 or n_elements(cat) eq 0 or n_elements(expstr) eq 0 or n_elements(chstr) eq 0 then begin
  print,'Syntax - nsc_instcal_calibrate_selfcatzpterm,expdir,cat,expstr,chstr'
  return
endif

t0 = systime(1)

if n_elements(logfile) eq 0 then logf=-1 else logf=logfile

; Set zero-point column to bad by default
expstr.zpterm = 999999.
expstr.zptermerr = 999999.
expstr.zptermsig = 999999.
expstr.zpspatialvar_rms = 999999.
expstr.zpspatialvar_range = 999999.
expstr.zpspatialvar_nccd = 0
expstr.nrefmatch = 0
expstr.nrefgdmatch = 0
chstr.zpterm = 999999.
chstr.zptermerr = 999999.
chstr.nrefmatch = 0

; What instrument is this?
instrument = 'c4d'  ; by default
if stregex(expdir,'/k4m/',/boolean) eq 1 then instrument='k4m'
if stregex(expdir,'/ksb/',/boolean) eq 1 then instrument='ksb'
instfilt = instrument+'-'+strtrim(expstr[0].filter,2)
; version
lo = strpos(expdir,'dl1')
dum = strsplit(strmid(expdir,lo+3),'/',/extract)
version = dum[4]
if version eq 'c4d' or version eq 'k4m' or version eq 'ksb' then version=''  ; version1

; Central coordinates
cenra = expstr.ra
cendec = expstr.dec

; Read in the summary structure
sumfile = dldir+'users/dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_calibrate.fits'
if file_test(sumfile) eq 0 then begin
  if not keyword_set(silent) then printlog,logf,sumfile,' NOT FOUND'
  return  
endif
sum = MRDFITS(sumfile,1,/silent)
allinstfilt = strtrim(sum.instrument,2)+'-'+strtrim(sum.filter,2)
gdfilt = where(allinstfilt eq instfilt and sum.nrefmatch gt 10 and sum.success eq 1 and sum.fwhm lt 2.0,ngdfilt)
if ngdfilt eq 0 then begin
  if not keyword_set(silent) then printlog,logf,'No good ',instfilt,' exposures'
  return
endif
sumfilt = sum[gdfilt]

; Minimum matching radius
minradius = 0.43
if instrument eq 'ksb' then minradius=minradius > 0.75
if instrument eq 'c4d' then minradius=minradius > 1.1

; Calculate healpix number
nside = 32L
theta = (90-sumfilt.dec)/radeg
phi = sumfilt.ra/radeg
ANG2PIX_RING,nside,theta,phi,allpix
ANG2PIX_RING,nside,(90-cendec)/radeg,cenra/radeg,cenpix
MATCH,allpix,cenpix,ind1,ind2,/sort,count=nmatch
if nmatch eq 0 then begin
  if not keyword_set(silent) then printlog,logf,'No good ',instfilt,' exposure at this position'
  return
endif
sumfiltpos = sumfilt[ind1]
dist = SPHDIST(sumfiltpos.ra,sumfiltpos.dec,cenra,cendec,/deg)
gddist = where(dist lt minradius,ngddist)
if ngddist eq 0 then begin
  if not keyword_set(silent) then printlog,logf,'No good ',instfilt,' exposures within ',stringize(minradius,ndec=2),' deg'
  return
endif
refexp = sumfiltpos[gddist]
nrefexp = n_elements(refexp)
if not keyword_set(silent) then printlog,logf,strtrim(ngddist,2),' good exposures within ',stringize(minradius,ndec=2),' deg'

; Pick the best ones
si = reverse(sort(refexp.zptermerr))
refexp = refexp[si]

; Loop over the reference exposures
;----------------------------------
refcnt = 0LL
nrefexpused = 0L
FOR i=0,nrefexp-1 do begin

  ; Load the catalog
  catfile = strtrim(refexp[i].expdir,2)+'/'+strtrim(refexp[i].base,2)+'_cat.fits'
  cat1 = MRDFITS(catfile,1,/silent)
  ; Only want good ones
  gdcat1 = where(cat1.imaflags_iso eq 0 and ((cat1.flags and 8) ne 8) and ((cat1.flags and 16) ne 16) and $
                 cat1.cmag lt 50 and 1.087/cat1.cerr gt 10 and cat1.fwhm_world*3600 lt 2.5 and $
                 ; Remove bad measurements
                 ; X: 1-1024 okay
                 ; X: 1025-2049 bad
                 ; use 1000 as the boundary since sometimes there's a sharp drop
                 ; at the boundary that causes problem sources with SExtractor
                 not (cat1.x_image gt 1000 and cat1.ccdnum eq 31),ngdcat1)

  ; Only want sources inside the radius of the exposure

  if not keyword_set(silent) then printlog,logf,i+1,ngdcat1,'   '+refexp[i].base,format='(I5,I10,A-50)'

  ; Some good sources
  if ngdcat1 gt 0 then begin
    nrefexpused++

    ; Convert to small catalog version
    nnew = ngdcat1
    schema = {souceid:'',ra:0.0d0,dec:0.0d0,ndet:1L,cmag:0.0,cerr:0.0}
    newcat = REPLICATE(schema,nnew)
    STRUCT_ASSIGN,cat1[gdcat1],newcat,/nozero

    ; First one, start REF structure
    if n_elements(ref) eq 0 then begin
      ; Start ref structure
      ref = replicate(schema,1e5>nnew)
      nref = n_elements(ref)
      ; Convert CMAG/CMAGERR to cumulative quantities
      cmag = 2.5118864d^newcat.cmag * (1.0d0/newcat.cerr^2)
      cerr = 1.0d0/newcat.cerr^2
      newcat.cmag = cmag
      newcat.cerr = cerr
      ; Add to REF
      ref[0:nnew-1] = newcat
      refcnt += nnew

    ; Second or later, add to them
    endif else begin

      ; Crossmatch
      SRCMATCH,ref[0:refcnt-1].ra,ref[0:refcnt-1].dec,newcat.ra,newcat.dec,0.5,ind1,ind2,/sph,count=nmatch
      if not keyword_set(silent) then printlog,logf,'   ',strtrim(nmatch,2),' crossmatches'

      ; Combine crossmatches
      if nmatch gt 0 then begin
        ; Cumulative sum of weighted average
        ref[ind1].cmag += 2.5118864d^newcat[ind2].cmag * (1.0d0/newcat[ind2].cerr^2)
        ref[ind1].cerr += 1.0d0/newcat[ind2].cerr^2
        ref[ind1].ndet++
        ; Remove
        if nmatch lt n_elements(newcat) then REMOVE,ind2,newcat else undefine,newcat
      endif

      ; Add leftover ones
      if n_elements(newcat) gt 0 then begin
        if not keyword_set(silent) then printlog,logf,'   Adding ',strtrim(n_elements(newcat),2),' leftovers'
        ; Convert CMAG/CMAGERR to cumulative quantities
        cmag = 2.5118864d^newcat.cmag * (1.0d0/newcat.cerr^2)
        cerr = 1.0d0/newcat.cerr^2
        newcat.cmag = cmag
        newcat.cerr = cerr
        ; Add more elements
        if refcnt+n_elements(newcat) gt nref then begin
          oldref = ref
          nnewref = 1e5 > n_elements(newcat)
          ref = REPLICATE(schema,nref+nnewref)
          ref[0:nref-1] = oldref
          nref = n_elements(ref)
          undefine,oldref
        endif
        ; Add new elements to REF
        ref[refcnt:refcnt+n_elements(newcat)-1] = newcat
        refcnt += n_elements(newcat)     
      endif  ; add leftovers
    endelse  ; 2nd or later exposure
  endif  ; some good sources

  ; Do we have enough exposures
  if nrefexpused ge (nrefexp < 10) then goto,DONE

  ;stop

ENDFOR  ; reference exposures loop
DONE:
; No good sources
if refcnt eq 0 then begin
  if not keyword_set(silent) then printlog,logf,'No good sources'
  return
endif
; Remove extra elements
ref = ref[0:refcnt-1]
; Take average for CMAG and CERR
newflux = ref.cmag / ref.cerr
newmag = 2.50*alog10(newflux)
newerr = sqrt(1.0/ref.cerr)
ref.cmag = newmag
ref.cerr = newerr

; Cross-match with the observed data
SRCMATCH,ref.ra,ref.dec,cat.ra,cat.dec,0.5,ind1,ind2,/sph,count=nmatch
if nmatch eq 0 then begin
  if not keyword_set(silent) then printlog,logf,'No reference sources match to the data'
  return
endif
cat2 = cat[ind2]
ref2 = ref[ind1]
; Take a robust mean relative to CMAG
mag = cat2.mag_auto + 2.5*alog10(expstr.exptime)  ; correct for the exposure time
err = cat2.magerr_auto
model_mag = ref2.cmag
col2 = fltarr(nmatch)
mstr = {col:col2,mag:float(mag),model:float(model_mag),err:float(err),ccdnum:long(cat2.ccdnum)}


; MEASURE THE ZERO-POINT
;-----------------------
; Same as from nsc_instcal_calibrate_fitzpterm.pro

; Fit the global and ccd zero-points
n = n_elements(mstr.col)
diff = mstr.model - mstr.mag
err = mstr.err
; Make a sigma cut
med = median([diff])
sig = mad([diff])
gd = where(abs(diff-med) lt 3*sig,ngd)
x = fltarr(n)
zpterm = dln_poly_fit(x[gd],diff[gd],0,measure_errors=err[gd],sigma=zptermerr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
zpterm = zpterm[0] & zptermerr=zptermerr[0]
; Save in exposure structure
expstr.zpterm = zpterm
expstr.zptermerr = zptermerr
expstr.zptermsig = sig
expstr.nrefmatch = n
expstr.nrefgdmatch = ngd

; Measure chip-level zpterm and variations
nchips = n_elements(chstr)
for i=0,nchips-1 do begin
  MATCH,chstr[i].ccdnum,mstr.ccdnum,ind1,ind2,/sort,count=nmatch
  if nmatch ge 5 then begin
    gd1 = where(abs(diff[ind2]-med) lt 3*sig,ngd1)
    if ngd1 ge 5 then begin
      zpterm1 = dln_poly_fit(x[ind2[gd1]],diff[ind2[gd1]],0,measure_errors=err[ind2[gd1]],sigma=zptermerr1,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
      chstr[i].zpterm = zpterm1
      chstr[i].zptermerr = zptermerr1
      chstr[i].nrefmatch = nmatch
   endif
 endif
endfor

; Measure spatial variations
gdchip = where(chstr.zpterm lt 1000,ngdchip)
if ngdchip gt 1 then begin
  expstr.zpspatialvar_rms = stddev([chstr[gdchip].zpterm])
  expstr.zpspatialvar_range = max(chstr[gdchip].zpterm)-min(chstr[gdchip].zpterm)
  expstr.zpspatialvar_nccd = ngdchip
endif

; Print out the results
if not keyword_set(silent) then begin
  printlog,logf,'NMATCH = ',strtrim(nmatch,2)
  printlog,logf,'ZPTERM=',stringize(expstr.zpterm,ndec=4),'+/-',stringize(expstr.zptermerr,ndec=4),'  SIG=',stringize(expstr.zptermsig,ndec=4),'mag'
  printlog,logf,'ZPSPATIALVAR:  RMS=',stringize(expstr.zpspatialvar_rms,ndec=3),' ',$
           'RANGE=',stringize(expstr.zpspatialvar_range,ndec=3),' NCCD=',strtrim(expstr.zpspatialvar_nccd,2)
endif

dt = systime(1)-t0
if not keyword_set(silent) then printlog,logf,'dt = ',stringize(dt,ndec=2),' sec.'

if keyword_set(stp) then stop

end
