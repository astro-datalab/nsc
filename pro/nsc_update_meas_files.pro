pro nsc_update_meas_files,expdir
  
;; Fix the aperture photometry in the meas files for exptime

NSC_ROOTDIRS,dldir,mssdir,localdir

print,'Updating aperture photometry in measurement catalog for exptime for ',expdir

base = file_basename(expdir)
catfile = expdir+'/'+base+'_cat.fits'
metafile = expdir+'/'+base+'_meta.fits'
measfile = expdir+'/'+base+'_meas.fits'
if file_test(catfile) eq 0 or file_test(metafile) eq 0 or file_test(measfile) eq 0 then begin
  print,'Catalog/metadata/meas file NOT FOUND'
  return
endif

; What instrument is this?
instrument = 'c4d'  ; by default
if stregex(expdir,'/k4m/',/boolean) eq 1 then instrument='k4m'
if stregex(expdir,'/ksb/',/boolean) eq 1 then instrument='ksb'
print,'This is a '+instrument+' exposure'

; Load catalog
cat = MRDFITS(catfile,1,/silent,status=status)  
if status lt 0 then goto,BOMB
ncat = n_elements(cat)

; Load metadata
meta = MRDFITS(metafile,1,/silent,status=status1)
if status1 lt 0 then goto,BOMB

; Load meas
meas = MRDFITS(measfile,1,/silent,status=status2)
if status2 lt 0 then goto,BOMB

; Check that the aperture photometry has not been corrected yet
;  Compare meas and cat mag_aper[0]
; meas.MEASID  c4d.386720.1.1
; cat.SOURCEID
; meas has zpterm applied
ind = where(strtrim(cat.sourceid,2) eq strtrim(meas[0].measid,2),nind)
if abs(cat[ind[0]].mag_aper[0]-meas[0].mag_aper1+meta.zpterm) gt 0.001 then begin
  print,'These aperture magnitudes have been corrected already!'
  return
endif

;; Correct the MEAS aperture photometry
meas.mag_aper1 += 2.5*alog10(meta.exptime)
meas.mag_aper2 += 2.5*alog10(meta.exptime)
meas.mag_aper4 += 2.5*alog10(meta.exptime)
meas.mag_aper8 += 2.5*alog10(meta.exptime)


; Temporarily move the old metadata and catalogs files aside for safe
; keeping
file_move,measfile,measfile+'.bak',/allow,/over

; Write new versions of files
print,'Writing new catalog to ',measfile
MWRFITS,meas,measfile,/create

; Delete old files
file_delete,measfile+'.bak',/allow

BOMB:

end
