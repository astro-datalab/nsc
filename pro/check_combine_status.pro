pro check_combine_status,pixinfo

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
if n_elements(version) eq 0 then version='v2'

files = file_search(localdir+'dnidever/nsc/instcal/'+version+'/tmp/nsccmb*batch',count=nfiles)
print,strtrim(nfiles,2),' command files'
cmd = strarr(nfiles)
for i=0,nfiles-1 do begin
  READLINE,files[i],lines
  cmd[i] = lines[0]
endfor
arr = strsplitter(cmd,',',/extract)
pix = reform(arr[1,*])

ui = uniq(pix,sort(pix))
upix = pix[ui]
nupix = n_elements(upix)
print,strtrim(nupix,2),' unique pixels'
outfile = dldir+'users/dnidever/nsc/instcal/'+version+'/combine/'+strtrim(long(upix)/1000,2)+'/'+strtrim(upix,2)+'.fits.gz'
test = file_test(outfile)
done = where(test eq 1,ndone,comp=left,ncomp=nleft)
if nleft eq 0 then begin
  print,'all done'
  stop
endif
leftpix = upix[left]
print,strtrim(nleft,2),' pixels with no output file'

; Check on the ones that are left
pixinfo = replicate({pix:0L,fitsexists:-1,nbatch:-1,startdate:0.0,enddate:0.0,dt:-1.0,nloglines:0L,nexp:-1L,nobj:-1L,$
                     wrotefile:-1,failure:-1,nosources:-1,missingfiles:-1,alreadyexists:-1,stillrunning:-1},nleft)
pixinfo.pix = leftpix
for i=0,nleft-1 do begin
  outfits = dldir+'users/dnidever/nsc/instcal/'+version+'/combine/'+strtrim(long(pixinfo[i].pix)/1000,2)+'/'+strtrim(pixinfo[i].pix,2)+'.fits'
  if file_test(outfits) eq 1 then pixinfo[i].fitsexists=1 else pixinfo[i].fitsexists=0
  MATCH,pix,pixinfo[i].pix,ind1,ind2,/sort,count=nmatch
  pixinfo[i].nbatch = nmatch
  if nmatch gt 0 then begin
    cmdfile = files[ind1]
    info = file_info(cmdfile)
    logfile = files[ind1]+'.log'
    loginfo = file_info(logfile)
    ; if multiple matches use the most recent one
    if nmatch gt 1 then begin
      si = reverse(sort(info.mtime))
      cmdfile = cmdfile[si[0]]
      info = info[si[0]]
      logfile = logfile[si[0]]
      loginfo = loginfo[si[0]]
    endif
    pixinfo[i].startdate = min([info.atime,info.mtime,info.ctime])
    pixinfo[i].enddate = loginfo.mtime
    pixinfo[i].dt = pixinfo[i].enddate - pixinfo[i].startdate
    READLINE,logfile[0],loglines,count=nloglines
    pixinfo[i].nloglines = nloglines
    ; get nexposures
    ind = where(stregex(loglines,'exposures that overlap this pixel and neighbors',/boolean) eq 1,nind)
    if nind gt 0 then begin
      line1 = loglines[ind[0]]
      nexp = long(first_el(strsplit(line1,' ',/extract)))
      pixinfo[i].nexp = nexp
    endif
    ; get nobjects
    ind = where(stregex(loglines,'final objects fall inside the pixel',/boolean) eq 1 and $
                stregex(loglines,'^None',/boolean) eq 0,nind)
    if nind gt 0 then begin
      line1 = loglines[ind[0]]
      nobj = long(first_el(strsplit(line1,' ',/extract)))
      pixinfo[i].nobj = nobj
    endif
    ; did it say it wrote the file
    ind = where(stregex(loglines,'Writing combined catalog to',/boolean) eq 1,nind)
    if nind gt 0 then pixinfo[i].wrotefile=1 else pixinfo[i].wrotefile=0
    ; check for "no sources" near the end
    ind = where(stregex(loglines,'^No sources in this pixel',/boolean) eq 1 or $
                stregex(loglines,'^None of the final objects fall inside the pixel',/boolean) eq 1,nind)
    if nind gt 0 then pixinfo[i].nosources=1 else pixinfo[i].nosources=0
    ; check for failure
    ind = where(stregex(loglines,'^% Execution halted at:',/boolean) eq 1,nind)
    if nind gt 0 then pixinfo[i].failure=1 else pixinfo[i].failure=0
    ; needed catalogs not found
    ind = where(stregex(loglines,'needed catalog files NOT FOUND',/boolean) eq 1,nind)
    if nind gt 0 then pixinfo[i].missingfiles=1 else pixinfo[i].missingfiles=0
    ; says output already exists
    ;  some say the catalog exists but it doesn't anymore  
    ind = where(stregex(loglines,'EXISTS already and /redo not set',/boolean) eq 1,nind)
    if nind gt 0 then pixinfo[i].alreadyexists=1 else pixinfo[i].alreadyexists=0
    ; still running
    if pixinfo[i].enddate gt systime(1)-3600. then pixinfo[i].stillrunning=1 else pixinfo[i].stillrunning=0

    ;if pixinfo[i].nosources eq 0 and pixinfo[i].failure eq 0 and pixinfo[i].missingfiles eq 0 then stop
  endif
endfor
; Other failures
; -some just end abruptly for no apparent reason, killed?

; Summary statistics
print,'Summary:'
ind = where(pixinfo.nosources eq 1,nind)
print,strtrim(nind,2),' pixels with no sources'
ind = where(pixinfo.failure eq 1,nind)
print,strtrim(nind,2),' pixels failed'
ind = where(pixinfo.missingfiles eq 1,nind)
print,strtrim(nind,2),' pixels had missing files'
ind = where(pixinfo.alreadyexists eq 1,nind)
print,strtrim(nind,2),' pixels used to have output file'
ind = where(pixinfo.stillrunning eq 1,nind)
print,strtrim(nind,2),' pixels are still running'

;stop

end
