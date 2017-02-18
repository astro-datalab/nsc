pro nsc_instcal_main,redo=redo

; Main NOAO DECam source catalog
dir = "/datalab/users/dnidever/decamcatalog/"

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
logfile = dir+'nsc_instcal_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

print, "Running SExtractor on the DECam InstCal Images"

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
str = MRDFITS(dir+'decam_instcal_list.fits',1)

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

; Select the good groups
; NO release year!!??
;ryear = long(strmid(str.release_date,0,4))
;rmonth long(strmid(str.release_date,5,2))
;public = (ryear lt 2017 or (ryear eq 2017 and rmonth le 2))

; Stripe82, -60<RA<60 and -1.26 < DEC < 1.26
gdexp82 = where((str.ra lt 61 or str.ra gt 259) and (str.dec ge -1.5 and str.dec le 1.5) and $
                exptime ge 30,ngdexp82)
print,strtrim(ngdexp82,2),' Stripe82 exposures'

; COSMOS field
; RA +150.11916667 (10:00:28.600)
; DEC +2.20583333 (+02:12:21.00)
; covers a 2 square degree region around this center,
gdexpcos = where((str.ra ge 148 and str.ra le 152) and (str.dec ge 0.0 and str.dec le 4.0) and $
                 exptime ge 30,ngdexpcos)
print,strtrim(ngdexpcos,2),' COSMOS exposures'

; LAF region
gdexplaf = where(mlon ge 40 and mlon le 80 and mlat ge -50 and mlat le 30 and $
                 (filt eq 'g' or filt eq 'i') and exptime gt 30,ngdexplaf)
;              (filt eq 'g' or filt eq 'i') and (public eq 1) and str.exposure gt 50,ngdgrp)
print,strtrim(ngdexplaf,2),' LAF exposures'

; Combine them
gdexp = [gdexp82, gdexpcos, gdexplaf]
ngdexp = n_elements(gdexp)

;stop

; Use jobs_daemon to run on all the unique stacks.
;for i=0,ngroups-1 do begin
for i=0,ngdexp-1 do begin

  fluxfile = strtrim(str[gdexp[i]].fluxfile,2)
  wtfile = strtrim(str[gdexp[i]].wtfile,2)
  maskfile = strtrim(str[gdexp[i]].maskfile,2)
  ;wtfile = repstr(fluxfile,'ooi','oow')
  ;maskfile = repstr(fluxfile,'ooi','ood')
  base = file_basename(fluxfile)

  ;; No wtfile in the same directory
  ;if file_test(wtfile) eq 0  then begin
  ;  MATCH,mss.base,file_basename(wtfile),ind1,ind2,/sort,count=nmatch
  ;  if nmatch eq 0 then begin
  ;    print,'Weight file missing for ',fluxfile
  ;    goto,BOMB
  ;  endif
  ;  wtfile = mss[ind1[0]].filename
  ;endif
  ;
  ;; No maskfile in the same directory
  ;if file_test(maskfile) eq 0  then begin
  ;  MATCH,mss.base,file_basename(maskfile),ind1,ind2,/sort,count=nmatch
  ;  if nmatch eq 0 then begin
  ;    print,'Mask file missing for ',fluxfile
  ;    goto,BOMB
  ;  endif
  ;  maskfile = mss[ind1[0]].filename
  ;endif

  ; Check if the output already exists.
  ;fhead = headfits(fluxfile,exten=0)
  ;dateobs = sxpar(fhead,'DATE-OBS')
  dateobs = str[gdexp[i]].date_obs
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  baseroot = file_basename(base,'.fits.fz')
  outfile = '/datalab/users/dnidever/decamcatalog/instcal/'+night+'/'+baseroot+'/'+baseroot+'_'+strtrim(1,2)+'.fits'
  if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
    print,outfile,' EXISTS and /redo NOT set'
    goto,BOMB
  endif

  if file_test(fluxfile) eq 1 and file_test(wtfile) eq 1 and file_test(maskfile) eq 1 then begin
    push,cmd,'/home/dnidever/projects/noaosourcecatalog/python/nsc_instcal.py '+fluxfile+' '+wtfile+' '+maskfile
    push,dirs,'/data0/dnidever/decamcatalog/instcal/tmp/'
    ; Run nsc_fullstack.py
   ;   retcode = subprocess.call(["nsc_fullstack.py",fluxfile,wtfile,maskfile])
    ;retcode = subprocess.call(["sex","flux.fits","-c","default.config"],stdout=sf,stderr=subprocess.STDOUT)
  endif else begin
    print,'Not all three flux/wt/mask files found for ',fluxfile
  endelse
  BOMB:
endfor

;stop

; Run PBS_DAEMON
; this works
PBS_DAEMON,cmd,dirs,/hyperthread,prefix='nsc',wait=10,nmulti=30
;PBS_DAEMON,cmd,dirs,/hyperthread,prefix='nsc',wait=10,nmulti=20

; End logfile
;------------
JOURNAL

stop

end
