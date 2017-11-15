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

; Get good RA/DEC
;  this was obtained with grabcoords_all.pro
coords = mrdfits(dir+'lists/allcoords.fits',1)
coords.file = strtrim(coords.file,2)
MATCH,str.fluxfile,'/net'+coords.file,ind1,ind2,/sort
dist=sphdist(str[ind1].ra,str[ind1].dec,coords[ind2].ra,coords[ind2].dec,/deg)
str[ind1].ra = coords[ind2].ra
str[ind1].dec = coords[ind2].dec
; 7 didn't match b/c the fluxfiles aren't found, trim them out
str = str[ind1]
nstr = n_elements(str)

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
  if i mod 1000 eq 0 then print,i
  dateobs = strtrim(str[i].date_obs,2)
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  baseroot = file_basename(strtrim(str[i].fluxfile,2),'.fits.fz')
  expdir = dir+str[i].instrument+'/'+night+'/'+baseroot+'/'
  list[i].expdir = expdir

  ;; Get RA/DEC values from the FLUXFILE header
  ;; Many of the RA/DEC values from the metadata database are incorrect
  ;fluxfile = str[i].fluxfile
  ;lo = strpos(fluxfile,'/mss1')
  ;fluxfile = mssdir+strmid(fluxfile,lo+6)
  ;head = headfits(fluxfile,exten=1)
  ;str[i].ra = sxpar(head,'crval1')
  ;str[i].dec = sxpar(head,'crval2')
endfor

;; 318 exposures have RA=DEC=NAN
;;  get their coordinates from the fluxfile
;bdra = where(finite(str.ra) eq 0 or finite(str.dec) eq 0,nbdra)
;for i=0,nbdra-1 do begin
;  fluxfile = str[bdra[i]].fluxfile
;  lo = strpos(fluxfile,'/mss1')
;  fluxfile = mssdir+strmid(fluxfile,lo+6)
;  head = headfits(fluxfile,exten=0)
;  str[bdra[i]].ra = sxpar(head,'ra')
;  str[bdra[i]].dec = sxpar(head,'dec')
;endfor

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

; gp09, run 1st batch, 8755
;cmd = cmd[0:8754]
;dirs = dirs[0:8754]

; Hulk, run 2nd batch, 8755
;cmd = cmd[8755:17509]
;dirs = dirs[8755:17509]

; Thing, run 3rd batch, 8755
;cmd = cmd[17510:*]
;dirs = dirs[17510:*]

stop


a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalibhpix',wait=wait,nmulti=nmulti 

stop

; Make the summary file with nsc_calibrate_summary.pro

print,'dt=',stringize(systime(1)-t0,ndec=2),' sec'

; End logfile
;------------
JOURNAL

stop

end
