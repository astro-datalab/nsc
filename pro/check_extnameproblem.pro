pro check_extnameproblem,dirs

;; Check if there is any problem with extension number issues
;; with the flux and mask files.

mssdir = '/mss1/'

;; Chip name and CCDNUM relations
;decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

undefine,sumstr

ndirs = n_elements(dirs)

;; Directory loop
For n=0,ndirs-1 do begin
  idir = dirs[n]
  inight = file_basename(idir)
  print,strtrim(n+1,2),' ',inight

  CD,idir

  outfile = '/net/dl1/users/dnidever/nsc/instcal/v3/extcheck/'+inight+'_summary.fits'
  if file_test(outfile) eq 1 and not keyword_set(redo) then begin
    print,outfile,' EXISTS and /red NOT set'
    goto,NIGHTBOMB
  endif

  ;; Get subdirectories
  expdirs = file_search('*',/test_directory,count=nexpdirs)

  undefine,sumstr

  ;; Loop over exposures
  For f=0,nexpdirs-1 do begin
    CD,expdirs[f]
    base = file_basename(expdirs[f])

    ;; Get mass more names from the log files
    logfile = base+'.log'
    if file_test(logfile) eq 0 then goto,EXPBOMB
    READLINE,logfile,loglines,nlineread=10
    ind = where(stregex(loglines,'Step #2: Copying InstCal images from mass store archive',/boolean) eq 1,nind)
    fline = loglines[ind[0]+1]
    lo = strpos(fline,'/archive')
    ;; make sure the mss1 directory is correct for this server
    fluxfile = mssdir+strtrim(strmid(fline,lo+1),2)
    wline = loglines[ind[0]+2]
    lo = strpos(wline,'/archive')
    wtfile = mssdir+strtrim(strmid(wline,lo+1),2)
    mline = loglines[ind[0]+3]
    lo = strpos(mline,'/archive')
    maskfile = mssdir+strtrim(strmid(mline,lo+1),2)

    ;; 2018-11-06 06:53:43,830 [INFO ]  Running SExtractor on c4d_161216_003410_ooi_r_v1 on host=gp06
    ;; 2018-11-06 06:53:43,830 [INFO ]    Temporary directory is: /data0/dnidever/nsc/instcal/v3/tmp/c4d_161216_003410_ooi_r_v1.1
    ;; 2018-11-06 06:53:43,831 [INFO ]  Step #2: Copying InstCal images from mass store archive
    ;; 2018-11-06 06:53:45,116 [INFO ]    /net/mss1/archive/pipe/20161215/ct4m/2012B-0001/c4d_161216_003410_ooi_r_v1.fits.fz
    ;; 2018-11-06 06:53:46,960 [INFO ]    /net/mss1/archive/pipe/20161215/ct4m/2012B-0001/c4d_161216_003410_oow_r_v1.fits.fz
    ;; 2018-11-06 06:53:47,089 [INFO ]    /net/mss1/archive/pipe/20161215/ct4m/2012B-0001/c4d_161216_003410_ood_r_v1.fits.fz
    ;; 2018-11-06 06:53:47,617 [INFO ]  Step #3: Running SExtractor on all subimages

    newstr = {night:'',exposure:'',fluxfile:'',maskfile:'',wtfile:'',okay:0}
    newstr.night = inight
    newstr.exposure = base
    newstr.okay = -1

    ;; Get fluxfile and maskfile header/extension information
    ;;   use symlink to make fits_open think it's a normal
    ;;   FITS file
    tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
    tmpfile += '.fits'
    FILE_LINK,fluxfile,tmpfile
    FITS_OPEN,tmpfile,fcb & FITS_CLOSE,fcb
    FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink 
    tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
    tmpfile += '.fits'
    FILE_LINK,maskfile,tmpfile
    FITS_OPEN,tmpfile,mfcb & FITS_CLOSE,mfcb
    FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink 

    newstr.fluxfile = fluxfile
    newstr.maskfile = maskfile
    newstr.wtfile = wtfile
    newstr.okay = 1  ;; good until found bad

    ;; Check if there are differences
    bd = where(fcb.extname ne mfcb.extname,nbd)
    if nbd gt 0 then begin
      print,inight,' ',base,'  PROBLEM'
      newstr.okay = 0
    endif else print,inight,' ',base,'  OK'

    PUSH,sumstr,newstr

    EXPBOMB:
    CD,idir
  Endfor ; exposure loop

  ;; Write out the summary information
  print,'Writing summary information to ',outfile
  MWRFITS,sumstr,outfile,/create

  NIGHTBOMB:
Endfor  ; night loop

;summary:
;
;for i=0,nnights-1 do begin
;  inight = nights[i]
;  outfile = delvedir+'exposures/rfile_check/'+inight+'_summary.fits'
;  if file_test(outfile) eq 1 then begin
;    sumstr = mrdfits(outfile,1,/silent)
;    bd = where(sumstr.okay eq 0,nbd)
;    if nbd eq 0 then print,strtrim(i+1,2),' ',inight else print,strtrim(i+1,2),' ',inight,' ',strtrim(nbd,2),' bad'
;    ;if nbd eq 0 then print,strtrim(i+1,2),' ',inight else print,strtrim(i+1,2),' ',inight,' ',strtrim(nbd,2),' bad: ',strjoin(sumstr[bd].expnum,' ')
;  endif else print,outfile,' NOT FOUND'
;endfor

;stop

end
