pro create_exposure_table,version

; Create the NSC exposure table

if n_elements(version) eq 0 then begin
  print,'Syntax - create_exposure_table,version'
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

; Load the calibration summary file
sum = MRDFITS(listdir+'nsc_calibrate_summary.fits.gz',1)
sum.expdir = strtrim(sum.expdir,2)
sum.base = strtrim(sum.base,2)
sum.dateobs = strtrim(sum.dateobs,2)
sum.filter = strtrim(sum.filter,2)
sexpdir = sum.expdir
b = where(strmid(sexpdir,0,4) ne '/net',nb)
if nb gt 0 then sexpdir[b]='/net'+sexpdir[b]

; Load the final list of combine exposures
healstr = MRDFITS(listdir+'nsc_instcal_combine_healpix_list.fits',1)
healstr.file = strtrim(healstr.file,2)
healstr.base = strtrim(healstr.base,2)
index = MRDFITS(listdir+'nsc_instcal_combine_healpix_list.fits',2)
; Get unique exposure
ui = uniq(healstr.file,sort(healstr.file))
uhstr = healstr[ui]
hexpdir = file_dirname(uhstr.file)+'/'
b = where(strmid(hexpdir,0,4) ne '/net',nb)
if nb gt 0 then hexpdir[b]='/net'+hexpdir[b]

; Match HEALSTR exposures to SUM to get the combined exposures list
MATCH,sexpdir,hexpdir,ind1,ind2,/sort,count=nmatch
fsum = sum[ind1]  ; final list of exposures
nexp = nmatch

; I need to add the EXPSTR index number to the HEALSTR structure


; Final columns
tags = ['instrument', 'exposure','expnum','ra','dec','dateobs','mjd','filter','exptime','airmass','nsources','fwhm','nchips','rarms','rastderr',$
        'decrms','decstderr','ebv','ngaiamatch','zptype','zpterm','zptermerr','zptermsig','nrefmatch','depth95','depth10sig']
types = [7,7,3,5,5,7,5,7,4,4,3,4,2,4,4,4,4,4,3,2,4,4,4,3,4,4]
; ngoodgaiamatch, ngoodrefmatch
; zpspatialvar_XXX
schema = create_struct(tags[0],fix('',type=types[0]))
for i=1,n_elements(tags)-1 do schema=create_struct(schema,tags[i],fix('',type=types[i]))
expstr = REPLICATE(schema,nexp)
STRUCT_ASSIGN,fsum,expstr       ; copy over information
expstr.exposure = fsum.base
expstr.nsources = 0

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

stop

; Loop over the pixels and load the meta-data
npix = n_elements(index)
for i=0,npix-1 do begin
  pixfile = cmbdir+strtrim(long(index[i].pix)/1000,2)+'/'+strtrim(index[i].pix,2)+'.fits.gz'
  if file_test(pixfile) eq 1 then begin
    meta = MRDFITS(pixfile,1,/silent)
    nmeta = n_elements(meta)
    mbase = strtrim(meta.base,2)
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
    if nmatch lt nmeta then stop,'not all exposures matched'
    if nmatch gt 0 then begin
      meta2 = meta[ind2]
      healstr_expind2 = healstr_expind1[ind1]
      expstr[healstr_expind2].nsources += meta.nsources  ; add the nsources
    endif
    
    stop

  endif else print,pixfile,' NOT FOUND'
endfor

; Cut out any exposures with NO sources
bdexp = where(expstr.nsources eq 0,nbdexp)
print,strtrim(nbdexp,2),' exposures with NO sources'
;REMOVE,bdexp,expstr


stop

end
