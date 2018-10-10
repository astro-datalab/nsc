pro copytot3a

;; Copy exposure files for Healpix 148487 to t3a

meta = mrdfits('/dl1/users/dnidever/nsc/instcal/v2/combine/148/148487.fits.gz',1)
meta.base = strtrim(meta.base,2)
expstr = mrdfits('/dl1/users/dnidever/nsc/instcal/v2/lists/nscdr1_exposures.fits',1)
expstr.base = strtrim(expstr.base,2)
MATCH,expstr.base,meta.base,ind1,ind2,/sort
meta = meta[ind2]
expstr = expstr[ind1]
expstr.expdir = strtrim(expstr.expdir,2)

for i=0,n_elements(meta)-1 do begin

  print,strtrim(i+1,2),'/',strtrim(n_elements(meta),2),' ',expstr[i].expdir

  files = file_search(expstr[i].expdir+'/*',count=nfiles)
  if nfiles eq 0 then stop,'NO FILES FOUND FOR '+expstr[i].expdir

  ;; Remove calibration output files
  ;;-rw-r--r-- 1 dnidever users    16441920 Nov 24 02:06 c4d_120925_044741_ooi_r_a1_cat.fits
  ;;-rw-r--r-- 1 dnidever users       46080 Nov 24 02:06 c4d_120925_044741_ooi_r_a1_meta.fits
  ;;-rw-r--r-- 1 dnidever users       25337 Nov 24 02:06 c4d_120925_044741_ooi_r_a1_calib.log
  ;;-rw-r--r-- 1 dnidever users     9077760 Feb 13 15:32 c4d_120925_044741_ooi_r_a1_meas.fits
  bd = where(stregex(files,'_cat.fits',/boolean) eq 1 or stregex(files,'_meta.fits',/boolean) eq 1 or $
             stregex(files,'_calib.log',/boolean) eq 1 or stregex(files,'_meas.fits',/boolean) eq 1,nbd)
  if nbd gt 0 then remove,bd,files

  ;; Create output directory
  lo = strpos(expstr[i].expdir,'/v2/')
  dir1 = strmid(expstr[i].expdir,lo+4)
  newdir = '/dl1/users/dnidever/nsc/instcal/t3a/'+dir1
  if file_test(newdir,/directory) eq 0 then file_mkdir,newdir

  ;; Copy over the files
  file_copy,files,newdir

  ;stop

endfor


stop

end
