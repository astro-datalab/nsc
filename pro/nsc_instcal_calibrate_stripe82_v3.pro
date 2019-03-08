pro nsc_instcal_calibrate_stripe82_v3

;; Run nsc_instcal_calibrate on Stripe82 exposures with the
;; initial model magnitude equations.

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='t3b'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(dir,/directory) eq 0 then file_mkdir,dir+'logs/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir

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
logfile = dir+'logs/nsc_instcal_calibrate_stripe82_v3.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile


listdir = '/dl1/users/dnidever/nsc/instcal/v3/lists/'

list = mrdfits(listdir+'decam_instcal_list.fits',1)
base = file_basename(strtrim(list.fluxfile,2),'.fits.fz')
sum = mrdfits(listdir+'nsc_measure_summary.fits',1)
sum.base = strtrim(sum.base,2)
match,sum.base,base,ind1,ind2,/sort,count=nmatch
sum2 = sum[ind1]
list2 = list[ind2]

;; Pick successful Stripe82 griz exposures
filter = strmid(list2.filter,0,1)
;gd = where((list2.ra lt 61 or list2.ra gt 259) and (list2.dec ge -1.5 and list2.dec le 1.5) and $
;           (filter eq 'g' or filter eq 'r' or filter eq 'i' or filter eq 'z') and $
;           list2.exposure ge 30 and sum2.success eq 1,ngd)

;; Add more low latitude fields, 349 exposures
glactc,ra1,dec1,2000.0,5,1,2,/deg
glactc,ra2,dec2,2000.0,5,5,2,/deg 
glactc,ra3,dec3,2000.0,5,10,2,/deg
glactc,ra4,dec4,2000.0,5,15,2,/deg
glactc,ra5,dec5,2000.0,10,20,2,/deg
glactc,ra6,dec6,2000.0,15,20,2,/deg
glactc,ra7,dec7,2000.0,30,25,2,/deg
ra = [ra1,ra2,ra3,ra4,ra5,ra6,ra7]
dec = [dec1,dec2,dec3,dec4,dec5,dec6,dec7]
for i=0,n_elements(ra)-1 do begin
  gd1 = where(list2.ra ge ra[i]-2 and list2.ra le ra[i]+2 and list2.dec ge dec[i]-2 and list2.dec le dec[i]+2 and $
              (filter eq 'g' or filter eq 'r' or filter eq 'i' or filter eq 'z') and $
              list2.exposure ge 30 and sum2.success eq 1,ngd1)
  if ngd1 eq 0 then stop,'no exposures'
  if ngd1 gt 0 then push,gd,gd1
endfor
ngd = n_elements(gd)

;;; Add the SMASH u-band exposures
;restore,'/dl1/users/dnidever/nsc/smash_matched_catalog_v3.dat'
;fieldstr.field = strtrim(fieldstr.field,2)
;smashinfo = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/check_calibrated_v6.fits',1)
;smashinfo.field = strtrim(smashinfo.field,2)
;MATCH,smashinfo.field,fieldstr.field,ind1,ind2,/sort
;smashinfo2 = smashinfo[ind1]
;undefine,sind
;for i=0,n_elements(smashinfo2)-1 do begin
;  dist = sphdist(smashinfo2[i].ra,smashinfo2[i].dec,list2.ra,list2.dec,/deg)
;  ind1 = where(dist lt 1.0 and filter eq 'u' and list2.exposure gt 10,nind1)
;  if nind1 gt 0 then push,sind,ind1
;endfor
;gd = [gd,sind]
;ngd = n_elements(gd)

;; 6746 exposures, about 1200-1900 in each band
;; 633 SMASH u-band exposures as well
expdirs = strtrim(sum2[gd].dir,2)
expdirs = repstr(expdirs,'/v3/','/t3b/')
bd = where(strmid(expdirs,0,4) eq '/net',nbd)
if nbd gt 0 then expdirs[bd] = strmid(expdirs[bd],4)

;; Make symlinks from t3b to v3
for i=0,ngd-1 do begin
  expdir1 = expdirs[i]
  print,strtrim(i+1,2),'/',strtrim(ngd,2),' ',expdir1
  if file_test(expdir1,/directory) eq 0 then begin
    file_mkdir,expdir1
    base = file_basename(expdir1)
    origdir1 = repstr(expdir1,'/t3b/','/v3/')
    fitsfiles = file_search(origdir1+'/'+base+'_*.fits',count=nfitsfiles)
    if nfitsfiles gt 0 then begin
      print,'Making fits symlinks for ',expdir1
      FILE_LINK,fitsfiles,expdir1
    endif else print,'No fits files for ',origdir1
    logfile = file_search(origdir1+'/'+base+'.log',count=nlogfile)
    if nlogfile gt 0 then begin
      print,'Making log symlinks for ',expdir1
      FILE_LINK,logfile,expdir1
    endif else print,'No log files for ',origdir1
  endif else print,expdir1,' already exists'
endfor

cmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(ngd)+tmpdir

stop

nmulti = 15
PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=1,nmulti=nmulti

; End logfile
;------------
JOURNAL

stop

end
