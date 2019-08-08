pro nsc_instcal_measure_update,expdir,redo=redo

;; Update the measurement catalog with the objectID

; Not enough inputs
if n_elements(expdir) eq 0 then begin
  print,'Syntax - nsc_instcal_measure_update,expdir,redo=redo'
  return
endif

t0 = systime(1)

;; Get version number from exposure directory
lo = strpos(expdir,'nsc/instcal/')
dum = strmid(expdir,lo+12)
version = strmid(dum,0,strpos(dum,'/'))
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
nside = 128
radeg = 180.0d0 / !dpi

; Check if output file already exists
base = file_basename(expdir)
outfile = expdir+'/'+base+'_1_meas.fits'
if file_test(outfile) eq 0 then begin
  print,outfile,' NOT FOUND'
  return
endif

print,'Adding objectID for measurement catalogs for exposure = ',base

;;  Load the exposure and metadata files
metafile = expdir+'/'+base+'_meta.fits'
meta = MRDFITS(metafile,1,/silent)
meta.file = strtrim(meta.file,2)
meta.base = strtrim(meta.base,2)
meta.dateobs = strtrim(meta.dateobs,2)
nmeta = n_elements(meta)
chstr = MRDFITS(metafile,2,/silent)
chstr.measfile = strtrim(chstr.measfile,2)
nchips = n_elements(chstr)

;; Load the catalogs
undefine,meas,nmeas
add_tag,chstr,'meas_index',0L,chstr
count = 0LL
for i=0,nchips-1 do begin
  cat1 = MRDFITS(chstr[i].measfile,1,/silent)
  ncat1 = n_elements(cat1)
  if n_elements(meas) eq 0 then begin
    schema = cat1[0]
    struct_assign,{dum:''},schema
    meas = replicate(schema,long(total(chstr.nmeas)))
  endif
  meas[count:count+ncat1-1] = cat1
  chstr[i].meas_index = count
  count += ncat1
endfor
meas.measid = strtrim(meas.measid,2)

;; Get the OBJECTID from the combined healpix file IDSTR structure
;;  remove any sources that weren't used

;; Figure out which healpix this figure overlaps

theta = (90-meas.dec)/radeg
phi = meas.ra/radeg
ANG2PIX_RING,nside,theta,phi,pix
ui = uniq(pix,sort(pix))
upix = pix[ui]
npix = n_elements(upix)

;; Load the healpix list
;listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_instcal_combine_healpix_list.fits.gz'
;if file_test(listfile) eq 0 then begin
;  print,listfile,' NOT FOUND'
;  return
;endif
;healstr = MRDFITS(listfile,1,/silent)  ; takes ~20s to load
;healstr.file = strtrim(healstr.file,2)
;healstr.base = strtrim(healstr.base,2)
;ind = where(healstr.base eq base,nind)  ; takes ~2s
;upix = healstr[ind].pix

;; Loop over the pixels
for i=0,npix-1 do begin
  objfile = dir+'combine/'+strtrim(long(upix[i])/1000,2)+'/'+strtrim(upix[i],2)+'.fits.gz'
  if file_test(objfile) eq 1 then begin
    idstr = MRDFITS(objfile,3,/silent)
    idstr.measid = strtrim(idstr.measid,2)
    nidstr = n_elements(idstr)
    MATCH,idstr.measid,meas.measid,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then meas[ind2].objectid=strtrim(idstr[ind1].objectid,2)
    print,i+1,upix[i],nmatch,format='(I5,I10,I7)'
  endif else print,objfile,' NOT FOUND'
endfor


;; Only keep sources with an objectid
ind = where(meas.objectid eq '',nind)
if nind gt 0 then STOP,'some objects are missing IDs'
;if nind gt 0 then begin
;  print,'Keeping ',strtrim(nind,2),' of ',strtrim(ncat,2),' sources'
;  meas = meas[ind]
;endif else begin
;  print,'No sources to keep'
;  return
;endelse

;; Output
print,'Updating measurement catalogs'
for i=0,nchips-1 do begin
  lo = chstr[i].meas_index
  hi = lo+chstr[i].nmeas-1
  meas1 = meas[lo:hi]
  meta1 = MRDFITS(chstr[i].measfile,2,/silent)  ;; load the meta extensions
  MWRFITS,meas1,chstr[i].measfile,/create
  MWRFITS,meta1,chstr[i].measfile,/silent
endfor

print,'dt = ',strtrim(systime(1)-t0,2),' sec.'

if keyword_set(stp) then stop

end
