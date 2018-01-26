pro make_measurement_catalog,expdir,version=version,redo=redo,stp=stp

;; Make the measurement catalog for a single exposure

t0 = systime(1)

nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(expdir) eq 0 then begin
  print,'Syntax - make_measurement_catalog,expdir,version=version,redo=redo,stp=stp'
  return
endif

; Check if output file already exists
base = file_basename(expdir)
dir = file_dirname(expdir)

outfile = expdir+'/'+base+'_meas.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS already and /redo not set'
  return
endif

print,'Creating measurement catalog for exposure = ',base

;;  Load the exposure and metadata files
metafile = expdir+'/'+base+'_meta.fits'
catfile = expdir+'/'+base+'_cat.fits'
meta = MRDFITS(metafile,1,/silent)
meta.file = strtrim(meta.file,2)
meta.base = strtrim(meta.base,2)
meta.dateobs = strtrim(meta.dateobs,2)
nmeta = n_elements(meta)
cat = MRDFITS(catfile,1,/silent)

stop

;; Convert to final format


;; Get the OBJECTID from the combined healpix file IDSTR structure
;;  remove any sources that weren't used

idstr = MRDFITS(objfile,3,/silent)
idstr.sourceid = strtrim(idstr.sourceid,2)
idstr.exposure = strtrim(idstr.exposure,2)
idstr.expnum = strtrim(idstr.expnum,2)
idstr.objectid = strtrim(idstr.objectid,2)
nidstr = n_elements(idstr)







;; Output


stop

end
