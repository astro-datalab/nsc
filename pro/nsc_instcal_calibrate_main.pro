pro nsc_instcal_calibrate_main,version,nmulti=nmulti,redo=redo

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
;dir = dldir+"users/dnidever/decamcatalog/"
;dir = "/datalab/users/dnidever/decamcatalog/"
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
logfile = dir+'logs/nsc_instcal_calibrate_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

if n_elements(nmulti) eq 0 then nmulti = 20
wait = 1

print, "Calibrating DECam/Mosiac3/Bok InstCal SExtractor catalogs"

; Find all of the directories
print,'Getting the exposure directories'
c4d_expdirs = file_search(dir+'c4d/20??????/*',/test_directory,count=nc4d_expdirs)
if nc4d_expdirs gt 0 then push,expdirs,c4d_expdirs
k4m_expdirs = file_search(dir+'k4m/20??????/*',/test_directory,count=nk4m_expdirs)
if nk4m_expdirs gt 0 then push,expdirs,k4m_expdirs
ksb_expdirs = file_search(dir+'ksb/20??????/*',/test_directory,count=nksb_expdirs)
if nksb_expdirs gt 0 then push,expdirs,ksb_expdirs
; GDL file_search doesn't have /test_directory
;expdirs = file_search(dir+'instcal/20??????/*',count=nexpdirs)
;gddirs = where(file_test(expdirs,/directory) eq 1,ngddirs)
;expdirs = expdirs[gddirs]
nexpdirs = n_elements(expdirs)
print,strtrim(nexpdirs,2),' exposure directories'

cmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(nexpdirs)+tmpdir
;dirs = strarr(nexpdirs)+'/data0/dnidever/decamcatalog/instcal/tmp/'

; ----- Run the LAF and Stripe82 exposures ----
;str = MRDFITS(dir+'decam_instcal_list.fits',1)
;nstr = n_elements(str)
;glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
;gal2mag,glon,glat,mlon,mlat
;filt = strmid(str.filter,0,1)
;exptime = str.exposure

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

;sum=mrdfits(dldir+'users/dnidever/nsc/instcal/v2/lists/nsc_instcal_calibrate.fits',1)
;sum.filter = strtrim(sum.filter,2)
;sum.expdir = strtrim(sum.expdir,2)
;g=where(sum.instrument eq 'c4d' and (sum.filter eq 'g' or sum.filter eq 'r') and sum.dec lt -29 and sum.nrefmatch lt 10,ng)
;cmd = 'nsc_instcal_calibrate,"'+sum[g].expdir+'",/redo'
;dirs = strarr(ng)+tmpdir
;
;; hulk
;PBS_DAEMON,cmd[0:6000],dirs[0:6000],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; Rerun all u-band exposures
;sum=mrdfits(dldir+'users/dnidever/nsc/instcal/v2/lists/nsc_instcal_calibrate.fits',1)
;sum.filter = strtrim(sum.filter,2)
;sum.expdir = strtrim(sum.expdir,2)
;g=where(sum.filter eq 'u',ng)
;expdir = sum[g].expdir
;bd = where(strmid(expdir,0,4) eq '/net',nbd)
;expdir[bd] = strmid(expdir[bd],4)
;cmd = 'nsc_instcal_calibrate,"'+expdir+'",/redo'
;dirs = strarr(ng)+tmpdir

; Run self-cal on failed u-band exposures
sum=mrdfits(dldir+'users/dnidever/nsc/instcal/v2/lists/nsc_instcal_calibrate.fits',1)
sum.filter = strtrim(sum.filter,2)
sum.instrument = strtrim(sum.instrument,2)
sum.expdir = strtrim(sum.expdir,2)
instfilt = sum.instrument+'-'+sum.filter
ifilt = 'c4d-u'
bd = where(instfilt eq ifilt and sum.nrefmatch eq 0,nbd)
gd = where(instfilt eq ifilt and sum.nrefmatch gt 10,ngd)
gdsum = sum[gd]

refmatch = lonarr(nbd)
for i=0,nbd-1 do begin
  minradius = 0.43
  if sum[bd[i]].instrument eq 'ksb' then minradius=minradius > 0.75
  if sum[bd[i]].instrument eq 'c4d' then minradius=minradius > 1.1
  MATCH,gdsum.pix,sum[bd[i]].pix,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    gddec = gdsum[ind1].dec
    dist = sphdist(gdsum[ind1].ra,gdsum[ind1].dec,sum[bd[i]].ra,sum[bd[i]].dec,/deg)
    ;gddist = where(dist lt minradius or gddec gt -30,ngddist)
    gddist = where(dist lt minradius,ngddist)
    refmatch[i] = ngddist
 endif
endfor

g = where(refmatch gt 0,ng)
expdir = sum[bd[g]].expdir
bb = where(strmid(expdir,0,4) eq '/net',nbb)
expdir[bb] = strmid(expdir[bb],4)
cmd = 'nsc_instcal_calibrate,"'+expdir+'",/redo,/selfcal'
dirs = strarr(ng)+tmpdir
;PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti 

stop

; hulk
;PBS_DAEMON,cmd[0:3383],dirs[0:3383],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; thing
;PBS_DAEMON,cmd[3384:6767],dirs[3384:6767],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; gp09
;PBS_DAEMON,cmd[6768:*],dirs[6768:*],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; Rerun failures with selfcal
;sum=mrdfits(dldir+'users/dnidever/nsc/instcal/v2/lists/nsc_instcal_calibrate.fits',1)
;sum.filter = strtrim(sum.filter,2)
;sum.expdir = strtrim(sum.expdir,2)
;g=where(sum.success eq 1 and sum.nrefmatch le 5,ng)
;expdir = sum[g].expdir
;bd = where(strmid(expdir,0,4) eq '/net',nbd)
;expdir[bd] = strmid(expdir[bd],4)
;cmd = 'nsc_instcal_calibrate,"'+expdir+'",/redo,/selfcal'
;dirs = strarr(ng)+tmpdir
;;lockfiles = expdir+file_basename(expdir)+'.lock'
;;donefiles = expdir+file_basename(expdir)+'.caldone'
;;done = where(file_test(donefiles) eq 1,ndone)
;;if ndone gt 0 then remove,done,cmd,dirs,expdir

; hulk
;PBS_DAEMON,cmd[0:3586],dirs[0:3586],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; thing
;PBS_DAEMON,cmd[3587:7173],dirs[3587:7173],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; gp09
;PBS_DAEMON,cmd[7174:*],dirs[7174:*],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti


; which ones still need to be done
;outfile = expdir+'/'+file_basename(expdir)+'_meta.v1.fits'
;otest = file_test(outfile)
;done2 = where(otest eq 1,ndone2)
;if ndone2 gt 0 then remove,done2,cmd,dirs

; hulk
;PBS_DAEMON,cmd[0:1143],dirs[0:1143],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; thing
;PBS_DAEMON,cmd[1144:2287],dirs[1144:2287],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

; gp09
;PBS_DAEMON,cmd[2288:*],dirs[2288:*],jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

stop

; Run PBS_DAEMON
;PBS_DAEMON,cmd,dirs,lockfiles=lockfiles,donefiles=donefiles,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti


;stop

; USE NSC_CALIBRATE_SUMMARY.PRO TO CREATE SUMMARY FILE NOW!!!

; Load all the summary/metadata files
;print,'Creating calibration summary file'
;expstr = replicate({expdir:'',instrument:'',metafile:'',success:0,file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
;                    airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,zpterm:0.0,zptermerr:0.0,zptermsig:0.0,$
;                    zpspatialvar_rms:999999.0,zpspatialvar_range:999999.0,zpspatialvar_nccd:0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nexpdirs)
;expstr.expdir = expdirs
;chstr = replicate({expdir:'',instrument:'',success:0,filename:'',ccdnum:0L,nsources:0L,cenra:999999.0d0,cendec:999999.0d0,$
;                   gaianmatch:0L,rarms:999999.0,racoef:dblarr(4),decrms:999999.0,$
;                   deccoef:dblarr(4),vra:dblarr(4),vdec:dblarr(4),zpterm:999999.0,zptermerr:999999.0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nexpdirs*61)
;chcnt = 0LL
;for i=0,nexpdirs-1 do begin
;  if (i+1) mod 5000 eq 0 then print,i+1
;  base = file_basename(expdirs[i])
;  type = ['c4d','k4m','ksb']
;  obs = ['ctio','kpno','kpno']
;  obsname = 'ctio'    ; by default
;  instrument = 'c4d'  ; by default
;  for j=0,n_elements(type)-1 do begin
;    if stregex(expdirs[i],'/'+type[j]+'/',/boolean) eq 1 then begin
;      instrument = type[j]
;      obsname = obs[j]
;    endif
;  endfor
;
;  metafile = expdirs[i]+'/'+base+'_meta.fits'
;  expstr[i].metafile = metafile
;  if file_test(metafile) eq 1 then begin
;    ; Exposure structure
;    expstr1 = MRDFITS(metafile,1,/silent)
;    temp = expstr[i]
;    struct_assign,expstr1,temp,/nozero
;    expstr[i] = temp
;    ;expstr[i].rarms = median(chstr1.rarms)
;    ;expstr[i].decrms = median(chstr1.decrms)
;    ;expstr[i].gaianmatch = median(chstr1.gaianmatch)
;    expstr[i].instrument = instrument
;    expstr[i].success = 1
;    ; Chip structure
;    chstr1 = MRDFITS(metafile,2,/silent)
;    nchstr1 = n_elements(chstr1)
;    temp = chstr[chcnt:chcnt+nchstr1-1]
;    struct_assign,chstr1,temp,/nozero
;    chstr[chcnt:chcnt+nchstr1-1] = temp
;    chstr[chcnt:chcnt+nchstr1-1].expdir = expdirs[i]
;    chstr[chcnt:chcnt+nchstr1-1].instrument = instrument
;    chstr[chcnt:chcnt+nchstr1-1].success = 1
;    chcnt += nchstr1
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
;      OBSERVATORY,obsname,obs
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
;; Trim CHSTR structure
;chstr = chstr[0:chcnt-1]
;gd = where(expstr.success eq 1,ngd)
;print,strtrim(ngd,2),' exposure successfully calibrated'
;print,'Writing summary file to ',dir+'lists/nsc_instcal_calibrate.fits'
;MWRFITS,expstr,dir+'lists/nsc_instcal_calibrate.fits',/create
;MWRFITS,chstr,dir+'lists/nsc_instcal_calibrate.fits',/silent

; End logfile
;------------
JOURNAL

stop

end
