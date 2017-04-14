pro nsc_instcal_calibrate_main,nmulti=nmulti,redo=redo

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+"users/dnidever/nsc/"
;dir = dldir+"users/dnidever/decamcatalog/"
;dir = "/datalab/users/dnidever/decamcatalog/"

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

if n_elements(nmulti) eq 0 then nmulti = 20
wait = 1

print, "Calibrating DECam InstCal SExtractor catalogs"

; Find all of the directories
print,'Getting the exposure directories'
expdirs = file_search(dir+'instcal/20??????/*',/test_directory,count=nexpdirs)
; GDL file_search doesn't have /test_directory
;expdirs = file_search(dir+'instcal/20??????/*',count=nexpdirs)
;gddirs = where(file_test(expdirs,/directory) eq 1,ngddirs)
;expdirs = expdirs[gddirs]
;nexpdirs = n_elements(expdirs)
print,strtrim(nexpdirs,2),' exposure directories'

cmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(nexpdirs)+localdir+'/dnidever/nsc/instcal/tmp/'
;dirs = strarr(nexpdirs)+'/data0/dnidever/decamcatalog/instcal/tmp/'

; ----- Run the LAF and Stripe82 exposures ----
str = MRDFITS(dir+'decam_instcal_list.fits',1)
nstr = n_elements(str)
glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
filt = strmid(str.filter,0,1)
exptime = str.exposure

;; Stripe82, -60<RA<60 and -1.26 < DEC < 1.26
;gdexp82 = where((str.ra lt 61 or str.ra gt 259) and (str.dec ge -1.5 and str.dec le 1.5) and $
;                str.exposure ge 30,ngdexp82)
;print,strtrim(ngdexp82,2),' Stripe82 exposures';
;
; COSMOS field
;; RA +150.11916667 (10:00:28.600)
; DEC +2.20583333 (+02:12:21.00)
; covers a 2 square degree region around this center,
;gdexpcos = where((str.ra ge 148 and str.ra le 152) and (str.dec ge 0.0 and str.dec le 4.0) and $
;                 str.exposure ge 30,ngdexpcos)
;print,strtrim(ngdexpcos,2),' COSMOS exposures'
;
;; LAF region
;gdexplaf = where(mlon ge 40 and mlon le 80 and mlat ge -50 and mlat le 30 and $
;                 str.exposure gt 30,ngdexplaf)
;;                 (filt eq 'g' or filt eq 'i') and str.exposure gt 30,ngdexplaf)
;print,strtrim(ngdexplaf,2),' LAF exposures'
;
;; Combine them
;gdexp = [gdexp82, gdexpcos, gdexplaf]
;ngdexp = n_elements(gdexp)
;
;str.fluxfile = strtrim(str.fluxfile,2)
;base = file_basename(str[gdexp].fluxfile,'.fits.fz')
;expbase=file_basename(expdirs)
;match,expbase,base,ind1,ind2,/sort
;stop
;PBS_DAEMON,cmd[ind1],dirs[ind1],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; Now do the rest
;remove,ind1,cmd,dirs
;PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

stop

; Run PBS_DAEMON
;PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

;; Load all the summary/metadata files
;print,'Creating calibration summary file'
;expstr = replicate({expdir:'',metafile:'',success:0,file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
;                    airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,zpterm:0.0,zptermerr:0.0,zptermsig:0.0,$
;                    nrefmatch:0L},nexpdirs)
;expstr.expdir = expdirs
;for i=0,nexpdirs-1 do begin
;  if (i+1) mod 5000 eq 0 then print,i+1
;  base = file_basename(expdirs[i])
;  metafile = expdirs[i]+'/'+base+'_meta.fits'
;  expstr[i].metafile = metafile
;  if file_test(metafile) eq 1 then begin
;    expstr1 = MRDFITS(metafile,1,/silent)
;    ;chstr1 = MRDFITS(metafile,2,/silent)
;    temp = expstr[i]
;    struct_assign,expstr1,temp,/nozero
;    expstr[i] = temp
;    ;expstr[i].rarms = median(chstr1.rarms)
;    ;expstr[i].decrms = median(chstr1.decrms)
;    ;expstr[i].gaianmatch = median(chstr1.gaianmatch)
;    expstr[i].success = 1;
;
;    ; Fix missing DATE-OBS
;    if strtrim(expstr[i].dateobs,2) eq '' or strtrim(expstr[i].dateobs,2) eq '0' then begin
;      fluxfile = strtrim(expstr[i].file)
;      lo = strpos(fluxfile,'archive') 
;      fluxfile = mssdir+strmid(fluxfile,lo)
;      head = headfits(fluxfile,exten=0) 
;      expstr[i].dateobs = sxpar(head,'DATE-OBS') 
;    endif
;    ; Fix missing AIRMASS
;    if expstr[i].airmass lt 0.9 then begin
;      OBSERVATORY,'ctio',obs
;      lat = obs.latitude 
;      lon = obs.longitude
;      jd = date2jd(expstr[i].dateobs) 
;      ra = expstr[i].ra 
;      dec = expstr[i].dec
;      expstr[i].airmass = AIRMASS(jd,ra,dec,lat,lon)
;    endif
;
;  endif else expstr[i].success=0
;endfor
;gd = where(expstr.success eq 1,ngd)
;print,strtrim(ngd,2),' exposure successfully calibrated'
;print,'Writing summary file to ',dir+'instcal/nsc_instcal_calibrate.fits'
;MWRFITS,expstr,dir+'instcal/nsc_instcal_calibrate.fits',/create

; End logfile
;------------
JOURNAL

;stop

end
