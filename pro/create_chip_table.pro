pro create_chip_table,version

; Create the NSC exposure table

if n_elements(version) eq 0 then begin
  print,'Syntax - create_chip_table,version'
  return
endif

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
if host eq 'nidevermacbookpro' then listdir=dldir else listdir = dir+'lists/'
cmbdir = dir +'combine/'

; exposures, EXPOSURE_mets.fits[1], nsc_calibrate_summary.fits
;     -use the nsc_healpix_list.fits file to figure out which exposures passed the cut
;     -need to figure out how many sources were actually used, I think I need to actually load the object.fits.gz[1] meta-data info to get this
;
; -chips, EXPOSURE_meta.fits[2], nsc_calibrate_summary.fits
;      -same as exposure

; Load the calibration exposure summary file
expsum = MRDFITS(listdir+'nsc_calibrate_summary.fits.gz',1)
expsum.expdir = strtrim(expsum.expdir,2)
expsum.base = strtrim(expsum.base,2)
expsum.base = file_basename(expsum.expdir)
expsum.dateobs = strtrim(expsum.dateobs,2)
expsum.filter = strtrim(expsum.filter,2)
b = where(strmid(expsum.expdir,0,4) ne '/net',nb)
if nb gt 0 then expsum[b].expdir='/net'+expsum[b].expdir
; remove exposures that failed
bd = where(expsum.success eq 0,nbd)
if nbd gt 0 then REMOVE,bd,expsum
nexpsum = n_elements(expsum)
; load the chip summary file
chsum = MRDFITS(listdir+'nsc_calibrate_summary.fits.gz',2)
chsum.expdir = strtrim(chsum.expdir,2)
chsum.filename = strtrim(chsum.filename,2)
b = where(strmid(chsum.expdir,0,4) ne '/net',nb)
if nb gt 0 then chsum[b].expdir='/net'+chsum[b].expdir
b = where(strmid(chsum.filename,0,4) ne '/net',nb)
if nb gt 0 then chsum[b].filename='/net'+chsum[b].filename
nchsum = n_elements(chsum)

;; Load the final list of combine exposures
;healstr = MRDFITS(listdir+'nsc_healpix_list.fits',1)
;healstr.file = strtrim(healstr.file,2)
;healstr.base = strtrim(healstr.base,2)
index = MRDFITS(listdir+'nsc_instcal_combine_healpix_list.fits.gz',2)
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
        'vertex_ra3','vertex_ra4','vertex_dec1','vertex_dec2','vertex_dec3','vertex_dec4','ngaiamatch','zpterm','zptermerr','nrefmatch','chipzpterm',$
        'chipzptermerr','chipnrefmatch','depth95','depth10sig']
types = [7,7,3,2,5,5,7,5,7,4,4,3,4,4,4,5,5,5,5,4,4,5,5,5,5,4,5,5,5,5,5,5,5,5,3,4,4,3,4,4,3,4,4]
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
chipstr.chipzpterm = chsum.zpterm
chipstr.chipzptermerr = chsum.zptermerr
chipstr.chipnrefmatch = chsum.nrefmatch
chipstr.nsources = 0

; need exposure, expnum, filter, mjd, dateobs, ebv, fwhm, exptime, airmass

; Sort EXPSUM by exposure
si = sort(expsum.base)
expsum = expsum[si]

; Make index of CHSTR to EXPSUM
si = sort(chipstr.exposure)
;chipstr = chipstr[si]
; find the breaks
brklo = where(chipstr[si].exposure ne shift(chipstr[si].exposure,1),nbrk)
brkhi = [brklo[1:nbrk-1]-1,n_elements(chipstr)-1]
nchexp = brkhi-brklo+1
chexpindex = {exposure:expsum.base,lo:brklo,hi:brkhi,num:nchexp,index:si}

; Loop over the exposures and stuff the exposure-level
;  information into CHIPSTR structure
for i=0,nexpsum-1 do begin
  if i mod 5000 eq 0 then print,i
  chind = si[brklo[i]:brkhi[i]]
  if expsum[i].base ne chipstr[chind[0]].exposure then stop,'exposure does not match'
  chipstr[chind].expnum = expsum[i].expnum
  chipstr[chind].filter = expsum[i].filter
  chipstr[chind].mjd = expsum[i].mjd
  chipstr[chind].dateobs = expsum[i].dateobs
  chipstr[chind].ebv = expsum[i].ebv
  chipstr[chind].fwhm = expsum[i].fwhm
  chipstr[chind].exptime = expsum[i].exptime
  chipstr[chind].airmass = expsum[i].airmass
  chipstr[chind].zpterm = expsum[i].zpterm    ; exposure-level zpterm/err/nrefmatch
  chipstr[chind].zptermerr = expsum[i].zptermerr
  chipstr[chind].nrefmatch = expsum[i].nrefmatch
endfor

; Load all of the CHIPSUM files
chipsumfiles = FILE_SEARCH(cmbdir+'chipsum/*/*_chipsum.fits',count=nchipsumfiles)
print,strtrim(nchipsumfiles,2),' chipsum files'
; figure out how many elements we need, TAKES TOO LONG
;count = 0LL
;for i=0,nchipsumfiles-1 do count+=sxpar(headfits(chipsumfiles[i],exten=1),'naxis2')
chipsum_schema = {pix:0L,exposure:'',ccdnum:0L,nsources:0L}
chipsum = replicate(chipsum_schema,1e7)
nchipsum = n_elements(chipsum)
cnt = 0LL
for i=0,nchipsumfiles-1 do begin
  if i mod 1000 eq 0 then print,i
  chipsum1 = MRDFITS(chipsumfiles[i],1,/silent)
  nchipsum1 = n_elements(chipsum1)
  ; add more elements
  if cnt+nchipsum1 gt nchipsum then begin
    print,'Ading more elements'
    old = chipsum
    nnew = 5e6
    chipsum = replicate(chipsum_schema,nchipsum+nnew)
    chipsum[0:nchipsum-1] = old
    nchipsum = n_elements(chipsum)
  endif
  chipsum[cnt:cnt+nchipsum1-1] = chipsum1
  cnt += nchipsum1
endfor
; trim off the extra elements
chipsum = chipsum[0:cnt-1]
nchipsum = n_elements(chipsum)
chipsum.exposure = strtrim(chipsum.exposure,2)
print,strtrim(nchipsum,2),' CHIPSUM elements'


; Find the breaks
chids = chipsum.exposure+'-'+strtrim(chipsum.ccdnum,2)
uichids = uniq(chids,sort(chids))
uchids = chids[uichids]
nuchids = n_elements(uchids)
si = sort(chids)
; find the breaks
brklo = where(chids[si] ne shift(chids[si],1),nbrk)
if nbrk gt 0 then begin
  brkhi = [brklo[1:nbrk-1]-1,n_elements(chids)-1]
  nchsrc = brkhi-brklo+1
endif else begin
  brklo = 0
  brkhi = n_elements(chids)-1
  nchsrc = n_elements(chids)
endelse
;uidstrindex = {id:uchids,exposure:strarr(nuchids),ccdnum:lonarr(nuchids),lo:brklo,hi:brkhi,num:nchsrc,index:si}
; Get unique exposures
dum2 = strsplitter(uchids,'-',/extract)
uchids_exposure_all = reform(dum2[0,*])  ; all
uidstrindex_exposure = reform(dum2[0,*])
uidstrindex_ccdnum = long(reform(dum2[1,*]))

; Now match to the CHIPSTR structure
chipstr_chid = chipstr.exposure+'-'+strtrim(chipstr.ccdnum,2)
MATCH,chipstr_chid,uchids,ind1,ind2,/sort,count=nmatch
chipstr[ind1].nsources = nchsrc[ind2]   ;uidstrindex[ind2].num

; Cut out any chips with NO sources
bdchip = where(chipstr.nsources eq 0,nbdchip)
print,strtrim(nbdchip,2),' chips with NO sources'
;REMOVE,bdexp,chipstr

; save the final file
;MWRFITS,chipstr,dir+'lists/nsc_instcal_chip.fits',/create

stop

end
