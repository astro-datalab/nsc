pro create_chip_table

; Create the NSC exposure table

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
if host eq 'nidevermacbookpro' then listdir=dldir else listdir = dir+'lists/'
cmbdir = dir +'combine/'

; exposures, EXPOSURE_mets.fits[1], nsc_instcal_calibrate.fits
;     -use the nsc_healpix_list.fits file to figure out which exposures passed the cut
;     -need to figure out how many sources were actually used, I think I need to actually load the object.fits.gz[1] meta-data info to get this
;
; -chips, EXPOSURE_meta.fits[2], nsc_instcal_calibrate.fits
;      -same as exposure

; Load the calibration exposure summary file
expsum = MRDFITS(listdir+'nsc_instcal_calibrate.fits',1)
expsum.expdir = strtrim(expsum.expdir,2)
expsum.base = strtrim(expsum.base,2)
expsum.dateobs = strtrim(expsum.dateobs,2)
expsum.filter = strtrim(expsum.filter,2)
b = where(strmid(expsum.expdir,0,4) ne '/net',nb)
if nb gt 0 then expsum[b].expdir='/net'+expsum[b].expdir
nexpnum = n_elements(expnum)
; load the chip summary file
chsum = MRDFITS(listdir+'nsc_instcal_calibrate.fits',2)
chsum.expdir = strtrim(chsum.expdir,2)
chsum.filename = strtrim(chsum.filename,2)
chsum.dateobs = strtrim(chsum.dateobs,2)
b = where(strmid(chsum.expdir,0,4) ne '/net',nb)
if nb gt 0 then chsum[b].expdir='/net'+chsum[b].expdir
b = where(strmid(chsum.filename,0,4) ne '/net',nb)
if nb gt 0 then chsum[b].filename='/net'+chsum[b].filename
nchsum = n_elements(chsum)

;; Load the final list of combine exposures
;healstr = MRDFITS(listdir+'nsc_healpix_list.fits',1)
;healstr.file = strtrim(healstr.file,2)
;healstr.base = strtrim(healstr.base,2)
;index = MRDFITS(listdir+'nsc_healpix_list.fits',2)
;; Get unique exposure
;ui = uniq(healstr.file,sort(healstr.file))
;uhstr = healstr[ui]
;hexpdir = file_dirname(uhstr.file)+'/'
;b = where(strmid(hexpdir,0,4) ne '/net',nb)
;if nb gt 0 then hexpdir[b]='/net'+hexpdir[b]
;
;; Match HEALSTR exposures to SUM to get the combined exposures list
;MATCH,sexpdir,hexpdir,ind1,ind2,/sort,count=nmatch
;fsum = sum[ind1]  ; final list of exposures
;nexp = nmatch

; I need to add the EXPSTR index number to the HEALSTR structure

; Final columns
tags = ['instrument', 'exposure','expnum','ccdnum','ra','dec','dateobs','mjd','filter','exptime','airmass','nsources','fwhm','rarms','rastderr',$
        'ra_coef1','ra_coef2','ra_coef3','ra_coef4','decrms','decstderr','dec_coef1','dec_coef2','dec_coef3','dec_coef4','ebv','vertex_ra1','vertex_ra2',$
        'vertex_ra3','vertex_ra4','vertex_dec1','vertex_dec2','vertex_dec3','vertex_dec4','ngaiamatch','zpterm','zptermerr','nrefmatch','depth95','depth10sig']
types = [7,7,3,2,5,5,7,5,7,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,3,4,4]
; ngoodgaiamatch, ngoodrefmatch
; zpspatialvar_XXX
schema = create_struct(tags[0],fix('',type=types[0]))
for i=1,n_elements(tags)-1 do schema=create_struct(schema,tags[i],fix('',type=types[i]))
chipstr = REPLICATE(schema,nchsum)
STRUCT_ASSIGN,chsum,chipstr       ; copy over information
chipstr.ra = chsum.cenra
chipstr.dec = chsum.cendec
chipstr.exposure = file_basename(chsum.expdir)
chipstr.ra_coef1 = chsum.racoef[0]
chipstr.ra_coef2 = chsum.racoef[1]
chipstr.ra_coef3 = chsum.racoef[2]
chipstr.ra_coef4 = chsum.racoef[3]
chipstr.dec_coef1 = chsum.deccoef[0]
chipstr.dec_coef2 = chsum.deccoef[1]
chipstr.dec_coef3 = chsum.deccoef[2]
chipstr.dec_coef4 = chsum.deccoef[3]
chipstr.vertex_ra1 = chsum.vra[0]
chipstr.vertex_ra2 = chsum.vra[1]
chipstr.vertex_ra3 = chsum.vra[2]
chipstr.vertex_ra4 = chsum.vra[3]
chipstr.vertex_dec1 = chsum.vdec[0]
chipstr.vertex_dec2 = chsum.vdec[1]
chipstr.vertex_dec3 = chsum.vdec[2]
chipstr.vertex_dec4 = chsum.vdec[3]
chipstr.nsources = 0

; need exposure, expnum, filter, mjd, dateobs, ebv, fwhm, exptime, airmass


; Make index of CHSTR to EXPSUM
si = sort(chipstr.expdir)
;chipstr = chipstr[si]
; find the breaks
brklo = where(chipstr[si].expdir ne shift(chipstr[si].expdir,1),nbrk)
brkhi = [brklo[1:nbrk-1]-1,n_elements(chipstr)-1]
npixexp = brkhi-brklo+1

; Loop over the exposures
for i=0,nexpnum-1 do begin
  chind = si[brklo[i]:brkhi[i]]
  chipstr[chind].expnum = expsum[i].expnum
  chipstr[chind].filter = expsum[i].filter
  chipstr[chind].mjd = expsum[i].mjd
  chipstr[chind].dateobs = expsum[i].dateobs
  chipstr[chind].ebv = expsum[i].ebv
  chipstr[chind].fwhm = expsum[i].fwhm
  chipstr[chind].exptime = expsum[i].exptime
  chipstr[chind].airmass = expsum[i].airmass
endfor

stop

; Get list of pix for each exposure
indhealstr = l64indgen(n_elements(healstr))
allhexpdir = file_dirname(healstr.file)+'/'
b = where(strmid(allhexpdir,0,4) ne '/net',nb)
if nb gt 0 then allhexpdir[b]='/net'+allhexpdir[b]
allpix = healstr.pix
si = sort(allhexpdir)
allhexpdir = allhexpdir[si]
allpix = allpix[si]
indhealstr = indhealstr[si]
; find the breaks
brklo = where(allhexpdir ne shift(allhexpdir,1),nbrk)
brkhi = [brklo[1:nbrk-1]-1,n_elements(allhexpdir)-1]
npixexp = brkhi-brklo+1
; Loop over the exposures
healstr_expind = lonarr(n_elements(healstr))
for i=0L,nexp-1 do begin
  lo = brklo[i]
  hi = brkhi[i]
  indhealstr1 = indhealstr[lo:hi]  ; indices into HEALSTR
  healstr_expind[indhealstr1] = i
endfor

; Loop over the pixels and load the meta-data
npix = n_elements(index)
undefine,badpix,duppix
for i=0,npix-1 do begin
  if i mod 100 eq 0 then print,i
  pixfile = cmbdir+strtrim(long(index[i].pix)/1000,2)+'/'+strtrim(index[i].pix,2)+'.fits.gz'
  if file_test(pixfile) eq 1 then begin
    meta = MRDFITS(pixfile,1,/silent)
    nmeta = n_elements(meta)
    mbase = strtrim(meta.base,2)
    ; Check for duplicate exposures
    uimbase = uniq(mbase,sort(mbase))
    if n_elements(uimbase) ne nmeta then begin
      print,index[i].pix,' has duplicates'
      push,duppix,index[i].pix
      mbaseorig = mbase
      mbase = mbase[uimbase]  ; get unique ones
      nmeta = n_elements(mbase)
    endif
    ; Load the IDs
    idstr = MRDFITS(pixfile,3,/silent)
    idstr.sourceid = strtrim(idstr.sourceid,2)
    nidstr = n_elements(idstr)
    uiidstr = uniq(idstr,sort(idstr))
    ; Get the HEALSTR elements for this pixel
    plo = index[i].lo
    phi = index[i].hi
    file1 = healstr[plo:phi].file
    base1 = file_basename(file1,'_cat.fits')
    expdir1 = file_basename(file1)+'/'
    b = where(strmid(expdir1,0,4) ne '/net',nb)
    if nb gt 0 then expdir1[b]='/net'+expdir1[b]
    healstr_expind1 = healstr_expind[plo:phi]
    ; Match up the exposures, some might be missing in META if they 
    ;  had no good sources
    MATCH,base1,mbase,ind1,ind2,/sort,count=nmatch
    if nmatch lt nmeta then begin
      print,index[i].pix,'  not all exposures matched'
      push,badpix,index[i].pix
      ; Rematch with full exposure structure
      MATCH,chipstr.exposure,mbase,ind1,ind2,/sort,count=nmatch2
      if nmatch2 ne nmeta then stop,'exposures missing in CHIPSTR'
    endif
    if nmatch gt 0 then begin
      meta2 = meta[ind2]
      healstr_expind2 = healstr_expind1[ind1]
      chipstr[healstr_expind2].nsources += meta2.nsources  ; add the nsources
    endif
  endif else print,pixfile,' NOT FOUND'
endfor

; Cut out any exposures with NO sources
bdchip = where(chipstr.nsources eq 0,nbdchip)
print,strtrim(nbchip,2),' chips with NO sources'
;REMOVE,bdexp,chipstr


stop

end
