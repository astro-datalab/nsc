pro create_v3_symlinks

;; Create symbolic links from v3 to v2 original source/measurement files

NSC_ROOTDIRS,dldir,mssdir,localdir
srcdir = '/dl1/users/dnidever/nsc/instcal/v2/'
destdir = '/dl1/users/dnidever/nsc/instcal/v3/'

;; Create the three instrument directories if they don't exist
if file_test(destdir+'c4d',/directory) eq 0 then file_mkdir,destdir+'c4d'
if file_test(destdir+'k4m',/directory) eq 0 then file_mkdir,destdir+'k4m'
if file_test(destdir+'ksb',/directory) eq 0 then file_mkdir,destdir+'ksb'

expstr = mrdfits('/dl1/users/dnidever/nsc/instcal/v2/lists/nsc_measure_summary.fits',1)
expstr.dir = strtrim(expstr.dir,2)
expstr.base = strtrim(expstr.base,2)
expstr.instrument = strtrim(expstr.instrument,2)
nexp = n_elements(expstr)

for i=0,nexp-1 do begin
  print,strtrim(i+1,2),' ',expstr[i].dir
  ;; Make the date directory if it doesn't exist
  arr = strsplit(expstr[i].dir,'/',/extract)
  narr = n_elements(arr)
  base = arr[narr-1]
  date = arr[narr-2]
  if file_test(destdir+expstr[i].instrument+'/'+date,/directory) eq 0 then file_mkdir,destdir+expstr[i].instrument+'/'+date
  ;; Make the exposure directory if it doesn't exist
  indir = expstr[i].dir+'/'
  outdir = destdir+expstr[i].instrument+'/'+date+'/'+base+'/'
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
  ;; Create symlinks for _cat.fits, _meta.fits, and .log files
  ;;   only create link if source file exists and destination does not
  catfile = base+'_cat.fits'
  metafile = base+'_meta.fits'
  logfile = base+'.log'
  if file_test(indir+catfile) eq 1 and file_test(outdir+catfile) eq 0 then file_link,indir+catfile,outdir+catfile
  if file_test(indir+metafile) eq 1 and file_test(outdir+metafile) eq 0then file_link,indir+metafile,outdir+metafile
  if file_test(indir+logfile) eq 1 and file_test(outdir+logfile) eq 0 then file_link,indir+logfile,outdir+logfile
endfor

stop

end
