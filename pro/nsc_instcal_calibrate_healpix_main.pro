pro nsc_instcal_calibrate_healpix_main,version=version,nmulti=nmulti,redo=redo

; Drive for NSC_INSTCAL_CALIBRATE_HEALPIX

if n_elements(nside) eq 0 then nside=64
radeg = 180.0d0 / !dpi

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir,host
hostname = first_el(strsplit(host,'.',/extract))
if n_elements(version) eq 0 then version='v3'
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
logtime = smonth+sday+syear+shour+sminute+ssecond
logfile = dir+'logs/nsc_instcal_calibrate_healpix_main.'+logtime+'.log'
JOURNAL,logfile

if n_elements(nmulti) eq 0 then nmulti = 10
wait = 1

print, "Calibrating DECam/Mosiac3/Bok InstCal SExtractor catalogs"

; ---- Create the list ----

; Load the list of exposures
list1 = MRDFITS(dir+'/lists/decam_instcal_list.fits.gz',1)
list2 = MRDFITS(dir+'/lists/mosaic3_instcal_list.fits.gz',1)
list3 = MRDFITS(dir+'/lists/bok90prime_instcal_list.fits.gz',1)
str = [list1,list2,list3]
undefine,list1,list2,list3
nstr = n_elements(str)
str.filter = strtrim(str.filter,2)
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
;; 13 didn't match b/c the fluxfiles aren't found, trim them out
;str = str[ind1]
;nstr = n_elements(str)

; Get good RA/DEC
bd = where(finite(str.ra) eq 0 or finite(str.dec) eq 0,nbd)
print,'Getting RA/DEC for ',strtrim(nbd,2),' exposures'
for i=0,nbd-1 do begin
  if i mod 50 eq 0 then print,i
  file = strtrim(str[bd[i]].fluxfile)
  if strmid(file,0,4) eq '/net' then file=strmid(file,4)
  head = headfits(file,exten=1,errmsg=errmsg)
  if errmsg eq '' then begin
    ra = sxpar(head,'crval1',count=nra)
    if nra gt 0 then str[bd[i]].ra = ra
    dec = sxpar(head,'crval2',count=ndec)
    if ndec gt 0 then str[bd[i]].dec = dec
  endif
endfor

; Reruning exposures that had coordinates problems but finished
; successfully originally
;;test = file_test(outfile)
;oldsum = mrdfits(dir+'lists/nsc_instcal_calibrate.fits.bak111217',1)
;bad = where(dist gt 0.1 and oldsum.success eq 1,nbad)

;cmd = 'nsc_instcal_calibrate,"'+strtrim(oldsum[bad].expdir,2)+'",/redo'
;dirs = strarr(nbad)+tmpdir
;stop
;PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti

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

;; Rerun u-band exposures with dec<0 with Skymapper u-band for calibration
bd = where(strmid(str.filter,0,1) eq 'u' and str.dec lt 0,nbd)
list = list[bd]
str = str[bd]

;; Only rerunning on failed exposures
;sumstr = mrdfits(dir+'lists/nsc_instcal_calibrate_failures.fits',1)
;bd = where(sumstr.nsources gt 100 and sumstr.fwhm le 2 and sumstr.exptime ge 30 and sumstr.meta_exists eq 0,nbd)
;failed_expdirs = strtrim(sumstr[bd].expdir,2)
;failed_expdirs = trailingslash(repstr(failed_expdirs,'/net/dl1/','/dl1/'))
;list.expdir = repstr(list.expdir,'/net/dl1/','/dl1/')
;MATCH,list.expdir,failed_expdirs,ind1,ind2,/sort,count=nmatch
;print,'Only keeping ',strtrim(nmatch,2),' failed exposures'
;list = list[ind1]
;str = str[ind1]

;failed = mrdfits(dir+'lists/nsc_instcal_calibrate_failures.fits',1)
;failed.expdir = strtrim(failed.expdir,2)
;list.expdir = repstr(list.expdir,'/net/dl1/','/dl1/')
;MATCH,list.expdir,failed.expdir,ind1,ind2,/sort,count=nmatch
;print,'Only keeping ',strtrim(nmatch,2),' failed exposures'
;list = list[ind1]
;str = str[ind1]

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
allcmd = 'nsc_instcal_calibrate_healpix,'+strtrim(upix,2)
if keyword_set(redo) then allcmd+=',/redo'
alldirs = strarr(npix)+tmpdir
allupix = upix

;; Only run healpix with DEC>-29
;g = where(str.dec gt -29,ng)
;uipix = uniq(list[g].pix,sort(list[g].pix))
;upix = list[g[uipix]].pix
;npix = n_elements(upix)
;print,strtrim(npix,2),' healpix to run'
;; Create the commands
;cmd = 'nsc_instcal_calibrate_healpix,'+strtrim(upix,2)
;if keyword_set(redo) then cmd+=',/redo'
;dirs = strarr(npix)+tmpdir

; RANDOMIZE
rnd = sort(randomu(0,npix))
allcmd = allcmd[rnd]
alldirs = alldirs[rnd]
allupix = allupix[rnd]

;; Hulk finished the first 3622 jobs (7/24/19)
;print,'Removing the first 3622 jobs that hulk already finished'
;remove,lindgen(3622),allcmd,alldirs,allupix
;npix = n_elements(allcmd)

;; Parcel out the jobs
hosts = ['gp09','hulk','thing']
nhosts = n_elements(hosts)
torun = lindgen(npix)
nperhost = npix/nhosts
for i=0,nhosts-1 do $
  if stregex(host,hosts[i],/boolean) eq 1 then torun=torun[i*nperhost:(i+1)*nperhost-1]
ntorun = n_elements(torun)
cmd = allcmd[torun]
dirs = alldirs[torun]
upix = allupix[torun]
print,'Running ',strtrim(n_elements(torun),2),' on ',hostname

; Saving the structure of jobs to run
runfile = dir+'lists/nsc_instcal_calibrate_healpix_main.'+hostname+'.'+logtime+'_run.fits'
print,'Writing running information to ',runfile
runstr = replicate({pix:0L,cmd:'',dir:''},n_elements(cmd))
runstr.pix = upix
runstr.cmd = cmd
runstr.dir = dirs
MWRFITS,runstr,runfile,/create


; 29,495 healpix, 14748 each

; ALL HEALPIX
; gp09, run 1st batch, 14748
;cmd = cmd[0:14747]
;dirs = dirs[0:14747]

; Thing, run 2nd batch, 14747
;cmd = cmd[14748:*]
;dirs = dirs[14748:*]

;; Hulk, run the gp09/thing failed jobs
;readline,dir+'/lists/nsccalibhpix_failed_gp09_hpix.txt',failed_hpix1
;readline,dir+'/lists/nsccalibhpix_failed_thing_hpix.txt',failed_hpix2
;failed_hpix = [failed_hpix1,failed_hpix2]
;cmd = 'nsc_instcal_calibrate_healpix,'+strtrim(failed_hpix,2)
;dirs = strarr(n_elements(cmd))+tmpdir

; rerun exposures that had coordinate problems and finished
; successfully the first time

; Run healpix with failures
;sum = mrdfits(dir+'lists/nsc_instcal_calibrate.fits',1)
;bad = where(sum.success eq 0,nbad)
;ui = uniq(sum[bad].pix,sort(sum[bad].pix))
;upix = sum[bad[ui]].pix
;npix = n_elements(upix)
;; Create the commands
;cmd = 'nsc_instcal_calibrate_healpix,'+strtrim(upix,2)
;;if keyword_set(redo) then cmd+=',/redo'
;dirs = strarr(npix)+tmpdir

; Thing, run first half
;cmd = cmd[0:317]
;dirs = dirs[0:317]

; Hulk, run second half
;cmd = cmd[318:*]
;dirs = dirs[318:*]


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
