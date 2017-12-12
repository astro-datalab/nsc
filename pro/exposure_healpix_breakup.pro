pro exposure_healpix_breakup,expdir,version=version,nside=nside,redo=redo

; Break up the exposure catalogs into separate files
; for the healpix regions

nexpdir = n_elements(expdir)
if nexpdir eq 0 then begin
  print,'Syntax - exposure_healpix_breakup,expdir,version=version,nside=nside,redo=redo'
endif

; Defaults
if n_elements(nside) eq 0 then nside = 128
if n_elements(version) eq 0 then version='v2'
radeg = 180.0d0 / !dpi

NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'

; Loop through the exposures
For i=0,nexpdir-1 do begin
  expdir1 = expdir[i]
  base1 = file_basename(expdir1)
  catfile = expdir1+'/'+base1+'_cat.fits'
  print,strtrim(i+1,2),'/'+strtrim(nexpdir,2),' ',catfile
  if file_test(catfile) eq 1 then begin
    cat = MRDFITS(catfile,1,/silent,status=status,error_action=2)
    if status ne 0 then goto,BOMB
    ; Get healpix values for all sources
    theta = (90-cat.dec)/radeg
    phi = cat.ra/radeg
    ANG2PIX_RING,nside,theta,phi,allpix
    ui = uniq(allpix,sort(allpix))
    upix = allpix[ui]
    nupix = n_elements(upix)
    print,strtrim(nupix,2),' unique healpix pixels'
    ; Break them up
    si = sort(allpix)
    allpix = allpix[si]
    cat = cat[si]
    brklo = where(allpix ne shift(allpix,1),nbrk)
    brkhi = [brklo[1:nbrk-1]-1,n_elements(allpix)-1]
    npixsrc = brkhi-brklo+1
    ; Loop through the pixels
    For j=0,nupix-1 do begin
      lo = brklo[j]
      hi = brkhi[j]
      cat1 = cat[lo:hi]
      pix1 = allpix[lo]
      outfile = expdir1+'/'+base1+'_cat.'+strtrim(pix1,2)+'.fits'
      if file_test(outfile) eq 0 or keyword_set(redo) then begin
        print,'  Writing ',outfile
        MWRFITS,cat1,outfile,/create
     endif else print,outfile,' EXISTS and /redo NOT set'
    Endfor  ; healpix pixel loop
  endif else print,catfile,' NOT FOUND'
  BOMB:
Endfor  ; exposure loop

end
