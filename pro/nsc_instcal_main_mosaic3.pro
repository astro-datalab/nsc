;+
;
; NSC_INSTCAL_MAIN_MOSAIC3
;
; This runs SExtractor on the Mosaic3 InstCal images.
; This is a wrapper around nsc_instcal.py which runs
; on one individual exposure.
;
; INPUTS:
;  =redo     Rerun on exposures that were previously processed.
;  =nmulti   The number of simultaneously jobs to run. Default is 30.
;  =maxjobs  The maximum number of exposures to attempt to process.
;              The default is 40,000.
;  /silent   Don't print much to the screen.
;  /unlock   Ignore the lock the files.
;
; OUTPUTS:
;  A log "journal" file is put in ROOTDIR+users/dnidever/nsc/instcal/logs/
;  as well as a structure with information on the jobs that were run.
;  The individual catalogs are put in ROOTDIR+users/dnidever/nsc/instcal/NIGHT/EXPOSURENAME/.
;
; USAGE:
;  IDL>nsc_instcal_main_mosaic3
;
; By D.Nidever  Feb 2017
;-

pro nsc_instcal_main_mosaic3,redo=redo,nmulti=nmulti,maxjobs=maxjobs,silent=silent,unlock=unlock

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
;dir = "/datalab/users/dnidever/decamcatalog/"
dir = dldir+'users/dnidever/nsc/instcal/'
if n_elements(maxjobs) eq 0 then maxjobs=4e4
if n_elements(nmulti) eq 0 then nmulti=30

t0 = systime(1)

; Log file
;------------------
; format is nsc_main_laf.DATETIME.log
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logtime = smonth+sday+syear+shour+sminute+ssecond
logfile = dir+'logs/nsc_instcal_main_mosaic3.'+logtime+'.log'
JOURNAL,logfile

print, "Running SExtractor on the Mosaic3 InstCal Images"

; Load all of the instcal exposures
;str = mrdfits(dir+"decam_instcal.fits.gz",1)
;n = n_elements(str)
;str.prodtype = strtrim(str.prodtype,2)
;; Only keep the newest versions for now
;gd = where(stregex(str.reference,'c4d',/boolean) eq 1,ngd)
;str = str[gd]
;add_tag,str,'base','',str
;str.base = strmid(strtrim(file_basename(str.reference),2),9)
;; Load the filenames
;mss = mrdfits(dir+"mss_filenames.fits.gz",1)
;mss.base = strtrim(mss.base,2)
;MATCH,mss.base,str.base,ind1,ind2,/sort,count=nmatch
;add_tag,str,'fluxfile','',str
;str[ind2].fluxfile = strtrim(mss[ind1].filename,2)
;; only keep the one with matches
;str = str[ind2]
;
;; Use dtacqnam to get EXPNUM
;
;; You can get EXPNUM from the pimary header
;
;; Use Date-OBS to find unique exposures
;; and use PLVER to get the most recent one
;dateobs = str.date_obs
;alldbl = doubles(dateobs,/all)
;dbl = doubles(dateobs)
;ndbl = n_elements(dbl)
;undefine,torem
;for i=0,ndbl-1 do begin
;  match,dateobs[alldbl],dateobs[dbl[i]],ind1,ind2,/sort,count=nmatch
;  dblind1 = alldbl[ind1]
;  plver = str[dblind1].plver
;  bestind = first_el(maxloc(plver))
;  left = dblind1
;  remove,bestind,left
;  push,torem,left  
;endfor
;; removing duplicate entires
;remove,torem,str
;
;; Get Mask and Weight filenames
;add_tag,str,'maskfile','',str
;add_tag,str,'wtfile','',str
;base = file_basename(str.fluxfile,'.fits.fz')
;mssbase = file_basename(mss.filename,'.fits.fz')
;mbase = repstr(base,'ooi','ood')
;MATCH,mssbase,mbase,ind1,ind2,/sort
;str[ind2].maskfile = mss[ind1].filename
;wbase = repstr(base,'ooi','oow')
;MATCH,mssbase,wbase,ind1,ind2,/sort
;str[ind2].wtfile = mss[ind1].filename
;; Only keeping ones with mask/weight files
;gd3 = where(str.fluxfile ne '' and str.maskfile ne '' and str.wtfile ne '',ngd3)
;str = str[gd3]
;MWRFITS,str,dir+'decam_instcal_list.fits',/create

str = MRDFITS(dir+'mosaic3_instcal_list.fits',1)
nstr = n_elements(str)
print,strtrim(nstr,2),' Mosaic3 InstCal images'

;stop

;sra = str[ui].ra
;ra = double(sexig2ten(sra)*15.0d0)
;sdec = str[ui].dec
;dec = double(sexig2ten(sdec))
glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
filt = strmid(str.filter,0,1)
exptime = str.exposure
;
;; The str.prodtype strings have extra ... at the end
;; which doesn't allow me to do the comparison properly
;prodtype = strarr(n)
;for i=0,n-1 do begin
;  if strcmp(str[i].prodtype,'expmap',6) eq 1 then prodtype[i]='expmap'
;  if strcmp(str[i].prodtype,'image',5) eq 1 then prodtype[i]='image'
;  if strcmp(str[i].prodtype,'wtmap',5) eq 1 then prodtype[i]='wtmap'
;  if strcmp(str[i].prodtype,'dqmask',6) eq 1 then prodtype[i]='dqmask'
;  if strcmp(str[i].prodtype,'image1',6) eq 1 then prodtype[i]='image1'
;endfor

gdexp = lindgen(nstr)
ngdexp = nstr

; Check the exposures
print,'Checking on the exposures'
expstr = replicate({fluxfile:'',wtfile:'',maskfile:'',allexist:0,outfile:'',done:0,locked:0,torun:0,cmd:'',cmddir:'',submitted:0},ngdexp)
for i=0,ngdexp-1 do begin
  if i mod 5000 eq 0 then print,i

  fluxfile = strtrim(str[gdexp[i]].fluxfile,2)
  wtfile = strtrim(str[gdexp[i]].wtfile,2)
  maskfile = strtrim(str[gdexp[i]].maskfile,2)
  ;wtfile = repstr(fluxfile,'ooi','oow')
  ;maskfile = repstr(fluxfile,'ooi','ood')
  base = file_basename(fluxfile)

  ; Change the root directory name
  ;  /net/mss1/blah/blah/
  fluxfile = mssdir+strmid(fluxfile,10)
  wtfile = mssdir+strmid(wtfile,10)
  maskfile = mssdir+strmid(maskfile,10)
  expstr[i].fluxfile = fluxfile
  expstr[i].wtfile = wtfile
  expstr[i].maskfile = maskfile

  ; Check if the output already exists.
  dateobs = str[gdexp[i]].date_obs
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  baseroot = file_basename(base,'.fits.fz')
  outfile = dldir+'users/dnidever/decamcatalog/instcal/'+night+'/'+baseroot+'/'+baseroot+'_'+strtrim(1,2)+'.fits'
  expstr[i].outfile = outfile

  ; Do all three files exist?
  if file_test(fluxfile) eq 1 and file_test(wtfile) eq 1 and file_test(maskfile) eq 1 then expstr[i].allexist=1
  ; Does the output file exist
  if file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1 then expstr[i].done = 1

  ; Not all three files exist
  if expstr[i].allexist eq 0 then begin
    if not keyword_set(silent) then print,'Not all three flux/wt/mask files found for ',fluxfile
    goto,BOMB
  endif

  ; Already done
  if (expstr[i].done eq 1) and not keyword_set(redo) then begin
    if not keyword_set(silent) then print,outfile,' EXISTS and /redo NOT set'
    goto,BOMB
  endif

  ;lock = djs_lockfile(outfile)
  lockfile = outfile+'.lock'
  testlock = file_test(lockfile)  

  ; No lock file
  ;if lock eq 1 or keyword_set(unlock) then begin
  if (testlock eq 0 or keyword_set(unlock)) then begin
    ;dum = djs_lockfile(outfile)  ; this is slow
    ;if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
    ;if testlock eq 0 then touchzero,outfile+'.lock'  ; this is fast
    expstr[i].cmd = '/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal.py '+fluxfile+' '+wtfile+' '+maskfile
    expstr[i].cmddir = localdir+'dnidever/nsc/instcal/tmp/'
    expstr[i].torun = 1
  ; Lock file exists
  endif else begin
    expstr[i].locked = 1
    expstr[i].torun = 0
    if not keyword_set(silent) then print,'Lock file exists ',outfile+'.lock'
  endelse
  BOMB:
endfor

torun = where(expstr.torun eq 1,ntorun)
stop
if ntorun eq 0 then begin
  print,'No exposures to process.'
  return
endif

; Pick the jobs to run
; MAXJOBS
if ntorun gt maxjobs then begin
  print,'More jobs than MAXJOBS.  Cutting down to ',strtrim(maxjobs,2),' jobs'
  expstr[torun[0:maxjobs-1]].submitted = 1
endif else expstr[torun].submitted = 1
tosubmit = where(expstr.submitted eq 1,ntosubmit)
cmd = expstr[tosubmit].cmd
cmddir = expstr[tosubmit].cmddir

; Lock the files that will be submitted
print,'Locking files to be submitted'
for i=0,ntosubmit-1 do begin
  outfile = expstr[tosubmit[i]].outfile
  if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
  lockfile = outfile+'.lock'
  testlock = file_test(lockfile)
  if testlock eq 0 then touchzero,outfile+'.lock'  ; this is fast
  expstr[tosubmit[i]].locked = 1
endfor

; Saving the structure of jobs to run
runfile = dir+'logs/nsc_instcal_main.'+logtime+'_run.fits'
print,'Writing running information to ',runfile
MWRFITS,expstr,runfile,/create

; Run PBS_DAEMON
stop
a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,cmddir,/hyperthread,prefix='nsc',wait=10,nmulti=nmulti

; Unlocking files
print,'Unlocking processed files'
file_delete,expstr[tosubmit].outfile+'.lock',/allow,/quiet
;for i=0,ntosubmit-1 do file_delete,expstr[tosubmit[i]].outfile+'.lock',/allow
;for i=0,ntosubmit-1 do djs_unlockfile,expstr[tosubmit[i]].outfile

print,'dt=',stringize(systime(1)-t0,ndec=2),' sec'

; End logfile
;------------
JOURNAL

stop

end
