pro nsc_measure_summary,version,nosources=nosources

; Make a summary file of all of the exposures that were Source Extracted

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'

; Find all of the directories
print,'Getting the exposure directories'
c4d_expdirs = file_search(dir+'c4d/20??????/*',/test_directory,count=nc4d_expdirs)
if nc4d_expdirs gt 0 then push,expdirs,c4d_expdirs
k4m_expdirs = file_search(dir+'k4m/20??????/*',/test_directory,count=nk4m_expdirs)
if nk4m_expdirs gt 0 then push,expdirs,k4m_expdirs
ksb_expdirs = file_search(dir+'ksb/20??????/*',/test_directory,count=nksb_expdirs)
if nksb_expdirs gt 0 then push,expdirs,ksb_expdirs

nexpdirs = n_elements(expdirs)
print,strtrim(nexpdirs,2),' exposure directories'

; Create structure
expstr = replicate({dir:'',instrument:'',base:'',nchips:0L,nsources:0L,success:0,dt:0.0,chip1date:0LL,logdate:0LL},nexpdirs)
expstr.dir = expdirs
undefine,instrument
if nc4d_expdirs gt 0 then push,instrument,strarr(nc4d_expdirs)+'c4d'
if nk4m_expdirs gt 0 then push,instrument,strarr(nk4m_expdirs)+'k4m'
if nksb_expdirs gt 0 then push,instrument,strarr(nksb_expdirs)+'ksb'
expstr.instrument = instrument
expstr.base = file_basename(expdirs)

; Loop through the exposure directories
for i=0,nexpdirs-1 do begin
  if i mod 5000 eq 0 then print,i
  dir1 = expdirs[i]
  base1 = expstr[i].base
  ; Get chip files
  chipfiles = file_search(dir1+'/'+base1+'_*.fits',count=nchipfiles)
  expstr[i].nchips = nchipfiles
  ; First chip date
  if nchipfiles gt 0 then begin
    info1 = file_info(chipfiles[0])
    expstr[i].chip1date = info1.mtime
  endif
  ; It succeeded if the final log file exists
  logfile = dir1+'/'+base1+'.log'
  expstr[i].success = file_test(logfile)
  if expstr[i].success gt 0 then begin
    loginfo = file_info(logfile)
    expstr[i].logdate = loginfo.mtime
  endif
  ; dt, need chip files
  if nchipfiles gt 0 then begin
    info1 = file_info(chipfiles[0])
    info2 = file_info(chipfiles[nchipfiles-1])
    dt = info2.mtime-info1.mtime
    ; this is the time between the end of first chip and last chip
    ;  correct for that
    expstr[i].dt = dt * nchipfiles / (nchipfiles-1.0)
  endif
  ; Get sources
  ;for j=0,nchipfiles-1 do begin
  ;  hd = headfits(chipfiles[j],exten=2,errmsg=errmsg)
  ;  if errmsg eq '' then begin
  ;    nsources = sxpar(hd,'naxis2')
  ;    expstr[i].nsources += sxpar(hd,'naxis2')
  ;  endif
  ;endfor
  if expstr[i].success eq 1 and not keyword_set(nosources) then begin
    READLINE,logfile,lines
    g = where(stregex(lines,'sextracted',/boolean) eq 1,ng)
    for j=0,ng-1 do begin
      line1 = lines[g[j]]
      arr = strsplit(line1,'\',/extract)
      g2 = where(stregex(arr,'Objects: detected',/boolean) eq 1,ng2)
      if ng2 eq 0 then goto,BOMB1
      line2 = arr[g2[0]]
      pos = strpos(line2,'sextracted')
      if pos eq -1 then goto,BOMB1
      nsrc = long(strmid(line2,pos+10))
      expstr[i].nsources += nsrc
      BOMB1:
    endfor
  endif
endfor

; Write out the file
outfile = dir+'lists/nsc_measure_summary.fits'
print,'Writing summary file to ',outfile
MWRFITS,expstr,outfile,/create

stop

end
