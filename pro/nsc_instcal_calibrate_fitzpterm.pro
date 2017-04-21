pro nsc_instcal_calibrate_fitzpterm,mstr,expstr,chstr

; Fit the global and ccd zero-points

n = n_elements(mstr.col)
diff = mstr.model - mstr.mag
err = mstr.err
; Make a sigma cut
med = median(diff)
sig = mad(diff)
gd = where(abs(diff-med) lt 3*sig,ngd)
x = fltarr(n)
zpterm = dln_poly_fit(x[gd],diff[gd],0,measure_errors=err[gd],sigma=zptermerr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
zpterm = zpterm[0] & zptermerr=zptermerr[0]
; Save in exposure structure
expstr.zpterm = zpterm
expstr.zptermerr = zptermerr
expstr.zptermsig = sig
expstr.nrefmatch = n

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

end
