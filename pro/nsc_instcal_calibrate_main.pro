pro nsc_instcal_calibrate_main,redo=redo

; Main NOAO DECam source catalog
dir = "/datalab/users/dnidever/decamcatalog/"

; Log file
;------------------
; format is nds_main_laf.DATETIME.log
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
logfile = dir+'nsc_instcal_calibrate_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

nmulti = 20

print, "Calibrating DECam InstCal SExtractor catalogs"

; Find all of the directories
expdirs = file_search(dir+'instcal/20??????/*',/test_directory,count=nexpdirs)
print,strtrim(nexpdirs,2),' exposure directories'

cmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(nexpdirs)+'/data0/dnidever/decamcatalog/instcal/tmp/'


; Run PBS_DAEMON
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='nsccalib',wait=10,nmulti=nmulti

; Load all the summary/metadata files
expstr = replicate({expdir:'',metafile:'',success:0,file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
                    airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,gaianmatch:0L,zpterm:0.0,zptermerr:0.0,zptermsig:0.0,$
                    nrefmatch:0L},nexpdirs)
expstr.expdir = expdirs
for i=0,nexpdirs-1 do begin
  base = file_basename(expdirs[i])
  metafile = expdirs[i]+'/'+base+'_meta.fits'
  expstr[i].metafile = metafile
  if file_test(metafile) eq 1 then begin
    expstr1 = MRDFITS(metafile,1,/silent)
    chstr1 = MRDFITS(metafile,2,/silent)
    temp = expstr[i]
    struct_assign,expstr1,temp,/nozero
    expstr[i] = temp
    ;expstr[i].rarms = median(chstr1.rarms)
    ;expstr[i].decrms = median(chstr1.decrms)
    ;expstr[i].gaianmatch = median(chstr1.gaianmatch)
    expstr[i].success = 1
  endif else expstr[i].success=0
endfor
gd = where(expstr.success eq 1,ngd)
print,strtrim(ngd,2),' exposure successfully calibrated'
print,'Writing summary file to ',dir+'instcal/nsc_instcal_calibrate.fits'
MWRFITS,expstr,dir+'instcal/nsc_instcal_calibrate.fits',/create

; End logfile
;------------
JOURNAL

;stop

end
