pro qacheck_photscatter,pix,redo=redo,version=version

; Measure photometric scatter for all band in healpix fields

t0 = systime(1)

; Combine all the exposures that fall in this healpix pixel
if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
dldir = '/net/dl2/'
if n_elements(version) eq 0 then version='v3'
dir = dldir+'dnidever/nsc/instcal/'+version+'/'
outdir = dir+'combine/'
radeg = 180.0d0 / !dpi

subdir = strtrim(long(pix)/1000,2)    ; use the thousands to create subdirectory grouping
outfile = outdir+'/qaphotrms/'+subdir+'/'+strtrim(pix,2)+'.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' already exists and /redo not set'
  return
endif

; Check obj file
objfile = outdir+'/'+subdir+'/'+strtrim(pix,2)+'.fits.gz'
if file_test(objfile) eq 0 then begin
  print,objfile,' NOT FOUND'
  return
endif

  str = {pix:long(pix),nexp:lonarr(7),minscatter:fltarr(7)+999999.,maglow:fltarr(7)+999999.,minfitscatter:fltarr(7)+999999.}

  meta = mrdfits(objfile,1,/silent)
  allobj = mrdfits(objfile,2,/silent)
  gdobj = where(allobj.fwhm lt 1.5 and allobj.class_star gt 0.5,ngdobj)
  allobj1 = allobj[gdobj]
  tags = tag_names(allobj)

  filters = ['u','g','r','i','z','y','vr']
  nfilters = n_elements(filters)

  ; Filter loop
  for j=0,nfilters-1 do begin
    minmed = 99.99
    maglow = 99.99  
    magind = where(tags eq strupcase(filters[j])+'MAG',nmagind)
    scatind = where(tags eq strupcase(filters[j])+'RMS',nscatind)
    detind = where(tags eq 'NPHOT'+strupcase(filters[j]),ndetind)
    ;if str.nexp[j] lt 3 then goto,BOMB
    ; Need at least 3 detections to get a good RMS
    gdndet = where(allobj1.(detind) gt 2,ngdndet)
    if ngdndet lt 5 then goto,BOmB
    allobj2 = allobj1[gdndet]
    str.nexp[j] = median([allobj2.(detind)])
    mag = allobj2.(magind)
    scatter = allobj2.(scatind)
    bindata,mag,scatter,bin=0.2,xbin,med,min=15,max=20,/med,gdind=gdind
    bindata,mag,scatter,bin=0.2,xbin,num,min=15,max=20,/hist
    xbin = xbin[gdind]
    med = med[gdind]
    ; sometimes the values are still 99.99
    gdind2 = where(med lt 10,ngdind2)
    xbin = xbin[gdind2]
    med = med[gdind2]
    ; smooth
    smmed = gsmooth(med,2)
    ; bright star, low scatter region
    minmed0 = min(med)
    minind0 = first_el(minloc(med))
    gdlow = where(med lt 2*minmed0 and abs(xbin-xbin[minind0]) lt 2.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 2*minmed0 and abs(xbin-xbin[minind0]) lt 3.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 4*minmed0 and abs(xbin-xbin[minind0]) lt 4.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 5*minmed0 and abs(xbin-xbin[minind0]) lt 5.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 6*minmed0 and abs(xbin-xbin[minind0]) lt 6.0,ngdlow)
    if ngdlow lt 2 then begin
      gdlow1 = where(abs(xbin-xbin[minind0]) lt 3.0,ngdlow1)
      si = sort(med[gdlow1])
      gdlow = gdlow1[si[0:2]]  ; take the lowest 3 points
    endif
    undefine,coef0,coef1,coef2
    coef0 = median([med[gdlow]])
    coef1 = robust_poly_fitq(xbin[gdlow],med[gdlow],1)
    if ngdlow gt 2 then begin
      coef2 = robust_poly_fitq(xbin[gdlow],med[gdlow],2)
      coef = coef2
    endif else coef=coef1
    minmed = min(med[gdlow])
    minind = gdlow[first_el(minloc(med[gdlow]))]
    maglow = xbin[minind]
    minfitmed = min(poly(xbin[gdlow],coef))
    str.minscatter[j] = minmed
    str.maglow[j] = maglow
    str.minfitscatter[j] = minfitmed
    BOMB:
    print,' ',filters[j],' ',str.nexp[j],'  ',minmed,' ',maglow
  endfor

; Make output directories
if file_test(outdir+'/qaphotrms/',/directory) eq 0 then file_mkdir,outdir+'/qaphotrms/'
if file_test(outdir+'/qaphotrms/'+subdir,/directory) eq 0 then file_mkdir,outdir+'/qaphotrms/'+subdir

; Writing output file
MWRFITS,str,outfile,/create

;stop

end
