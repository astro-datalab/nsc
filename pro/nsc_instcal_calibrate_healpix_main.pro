pro nsc_instcal_calibrate_healpix_main

; Drive for NSC_INSTCAL_CALIBRATE_HEALPIX

if n_elements(nside) eq 0 then nside=64
radeg = 180.0d0 / !dpi

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(dir,/directory) eq 0 then file_mkdir,dir+'logs/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir

t0 = systime(1)

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
logfile = dir+'logs/nsc_instcal_calibrate_healpix_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

if n_elements(nmulti) eq 0 then nmulti = 20
wait = 1

print, "Calibrating DECam/Mosiac3/Bok InstCal SExtractor catalogs"

; ---- Create the list ----

; Load the list of exposures
list1 = MRDFITS(dir+'/lists/decam_instcal_list.fits',1)
list2 = MRDFITS(dir+'/lists/mosaic3_instcal_list.fits',1)
list3 = MRDFITS(dir+'/lists/bok90prime_instcal_list.fits',1)
str = [list1,list2,list3]
undefine,list1,list2,list3
nstr = n_elements(str)
str.fluxfile = strtrim(str.fluxfile,2)
str.maskfile = strtrim(str.maskfile,2)
str.wtfile = strtrim(str.wtfile,2)
print,strtrim(nstr,2),' InstCal images'

; Start output file
list = replicate({expdir:'',pix:-1L,instrument:'',filter:''},nstr)
list.instrument = strtrim(str.instrument,2)
list.filter = strtrim(strmid(str.filter,0,2),2)
g = where(strmid(str.filter,0,4) eq 'bokr',ng)
if ng gt 0 then list[g].filter = 'r'
g = where(str.instrument eq 'k4m' and strmid(str.filter,0,2) eq 'zd',ng)
if ng gt 0 then list[g].filter = 'z'

; Calculate the output directory
for i=0,nstr-1 do begin
  dateobs = strtrim(str[i].date_obs,2)
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  baseroot = file_basename(strtrim(str[i].fluxfile,2),'.fits.fz')
  expdir = dir+str[i].instrument+'/'+night+'/'+baseroot+'/'
  list[i].expdir = expdir
endfor

; 318 exposures have RA=DEC=NAN
;  get their coordinates from the fluxfile
bdra = where(finite(str.ra) eq 0 or finite(str.dec) eq 0,nbdra)
for i=0,nbdra-1 do begin
  fluxfile = str[bdra[i]].fluxfile
  lo = strpos(fluxfile,'/mss1')
  fluxfile = mssdir+strmid(fluxfile,lo+6)
  head = headfits(fluxfile,exten=0)
  str[bdra[i]].ra = sxpar(head,'ra')
  str[bdra[i]].dec = sxpar(head,'dec')
endfor

; Calculate the healpix
theta = (90-str.dec)/radeg
phi = str.ra/radeg
ANG2PIX_RING,nside,theta,phi,ipring
list.pix = ipring

; Save the list
print,'Writing list to ',dir+'/lists/nsc_calibrate_healpix_list.fits'
MWRFITS,list,dir+'/lists/nsc_calibrate_healpix_list.fits',/create


;================
; RUN THE JOBS
;================
; Get the unique healpix pixels
uipix = uniq(list.pix,sort(list.pix))
upix = list[uipix].pix
npix = n_elements(upix)
print,strtrim(npix,2),' healpix to run'
; Create the commands
cmd = 'nsc_instcal_calibrate_healpix,'+strtrim(upix,2)
;if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(npix)+tmpdir

; RANDOMIZE
rnd = sort(randomu(0,npix))
cmd = cmd[rnd]
dirs = dirs[rnd]

stop

a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalibhpix',wait=wait,nmulti=nmulti 


; Load all the summary/metadata files
print,'Creating calibration summary file'
expstr = replicate({expdir:'',instrument:'',metafile:'',success:0,file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
                    airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,zpterm:0.0,zptermerr:0.0,zptermsig:0.0,$
                    zpspatialvar_rms:999999.0,zpspatialvar_range:999999.0,zpspatialvar_nccd:0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nstr)
expstr.expdir = list.expdir
chstr = replicate({expdir:'',instrument:'',success:0,filename:'',ccdnum:0L,nsources:0L,cenra:999999.0d0,cendec:999999.0d0,$
                   gaianmatch:0L,rarms:999999.0,racoef:dblarr(4),decrms:999999.0,$
                   deccoef:dblarr(4),vra:dblarr(4),vdec:dblarr(4),zpterm:999999.0,zptermerr:999999.0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nexpdirs*61)
chcnt = 0LL
for i=0,nstr-1 do begin
  if (i+1) mod 5000 eq 0 then print,i+1
  base = file_basename(expstr[i].expdir)
  type = ['c4d','k4m','ksb']
  obs = ['ctio','kpno','kpno']
  obsname = 'ctio'    ; by default
  instrument = 'c4d'  ; by default

  for j=0,n_elements(type)-1 do begin
    if stregex(expdirs[i],'/'+type[j]+'/',/boolean) eq 1 then begin
      instrument = type[j]
      obsname = obs[j]
    endif
  endfor

  metafile = expstr[i].expdir+'/'+base+'_meta.fits'
  expstr[i].metafile = metafile
  if file_test(metafile) eq 1 then begin
    ; Exposure structure
    expstr1 = MRDFITS(metafile,1,/silent)
    temp = expstr[i]
    struct_assign,expstr1,temp,/nozero
    expstr[i] = temp
    expstr[i].instrument = instrument
    expstr[i].success = 1
    ; Chip structure
    chstr1 = MRDFITS(metafile,2,/silent)
    nchstr1 = n_elements(chstr1)
    temp = chstr[chcnt:chcnt+nchstr1-1]
    struct_assign,chstr1,temp,/nozero
    chstr[chcnt:chcnt+nchstr1-1] = temp
    chstr[chcnt:chcnt+nchstr1-1].expdir = expdirs[i]
    chstr[chcnt:chcnt+nchstr1-1].instrument = instrument
    chstr[chcnt:chcnt+nchstr1-1].success = 1
    chcnt += nchstr1
    
    ; Fix missing DATE-OBS
    if strtrim(expstr[i].dateobs,2) eq '' or strtrim(expstr[i].dateobs,2) eq '0' then begin
      fluxfile = strtrim(expstr[i].file)
      lo = strpos(fluxfile,'archive')
      fluxfile = mssdir+strmid(fluxfile,lo)
      head = headfits(fluxfile,exten=0)
      expstr[i].dateobs = sxpar(head,'DATE-OBS')
    endif
    ; Fix missing AIRMASS
    if expstr[i].airmass lt 0.9 then begin
      OBSERVATORY,obsname,obs
      lat = obs.latitude
      lon = obs.longitude
      jd = date2jd(expstr[i].dateobs)
      ra = expstr[i].ra
      dec = expstr[i].dec
      expstr[i].airmass = AIRMASS(jd,ra,dec,lat,lon)
   endif

 endif else expstr[i].success=0
endfor
; Trim CHSTR structure
chstr = chstr[0:chcnt-1]
gd = where(expstr.success eq 1,ngd)
print,strtrim(ngd,2),' exposure successfully calibrated'
print,'Writing summary file to ',dir+'lists/nsc_instcal_calibrate.fits'
MWRFITS,expstr,dir+'lists/nsc_instcal_calibrate.fits',/create
MWRFITS,chstr,dir+'lists/nsc_instcal_calibrate.fits',/silent


print,'dt=',stringize(systime(1)-t0,ndec=2),' sec'

; End logfile
;------------
JOURNAL

stop

end
