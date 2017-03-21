pro nsc_instcal_combine_main,redo=redo,nmulti=nmulti

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
;nside = 256
nside = 128
nmulti = 30
radeg = 180.0d0 / !dpi

; Log file
;------------------
; format is nsc_combine_main.DATETIME.log
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
logfile = dir+'combine/nsc_instcal_combine_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

print, "Combining DECam InstCal catalogs"

; Fix MJD values from night numbers to full MJD
;  and bad AIRMASS values
;str = mrdfits(dir+'nsc_instcal_calibrate.fits',1,/silent)
;nstr = n_elements(str)
;OBSERVATORY,'ctio',obs
;for i=0,nstr-1 do begin
;  if (i+1) mod 10000 eq 0 then print,i+1
;  if str[i].success eq 1 then begin
;    dateobs = strtrim(str[i].dateobs,2)
;    if dateobs eq '' or dateobs eq '0' then begin
;      fluxfile = strtrim(str[i].file,2)
;      lo = strpos(fluxfile,'archive')
;      fluxfile = mssdir+strmid(fluxfile,lo)
;      head = headfits(fluxfile,exten=0)
;      dateobs = sxpar(head,'DATE-OBS',count=ndateobs)
;    if ndateobs eq 0 then stop,'no dateobs'
;      str[i].dateobs = dateobs
;    endif
;    str[i].mjd=date2jd(str[i].dateobs,/mjd)
;    ; Fix airmass values
;    if str[i].airmass lt 0.9 then begin
;      jd = date2jd(str[i].dateobs)
;      ra = str[i].ra
;      dec = str[i].dec
;      str[i].airmass = AIRMASS(jd,ra,dec,obs.latitude,obs.longitude)
;    endif
;  endif
;endfor
;;mwrfits,str,dir+'nsc_instcal_calibrate.fits',/create
;stop


; Restore the calibration summary file
str0 = mrdfits(dir+'nsc_instcal_calibrate.fits',1,/silent)
gd = where(str0.success eq 1,ngd)
str = str0[gd]
nstr = n_elements(str)
str.expdir = strtrim(str.expdir,2)
str.metafile = strtrim(str.metafile,2)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)
str.filter = strtrim(str.filter,2)
file = str.file
g1 = where(stregex(file,'/net/mss1/',/boolean) eq 1,ng1)
file[g1] = strmid(file[g1],10)
g2 = where(stregex(file,'/mss1/',/boolean) eq 1,ng2)
file[g2] = strmid(file[g2],6)

; APPLY RELEASE-DATE CUTS
list = MRDFITS(dir+'decam_instcal_list.fits',1)
list.fluxfile = strtrim(list.fluxfile,2)
fluxfile = strmid(list.fluxfile,10)
MATCH,fluxfile,file,ind1,ind2,/sort,count=nmatch
; some don't match because they were from a previous version
;  of the input list
release_date = strarr(n_elements(str))+'2020-01-01 00:00:00'
release_date[ind2] = list[ind1].release_date
release_year = long(strmid(release_date,0,4))
release_month = long(strmid(release_date,5,2))
release_day = long(strmid(release_date,8,2))
release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
release_cutoff = [2017,4,1]  ; April 1, 2017 for now
release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
gdrelease = where(release_mjd le release_cutoff_mjd,ngdrelease)
print,strtrim(ngdrelease,2),' exposures are PUBLIC'
str = str[gdrelease]  ; impose the public data cut

;; Fix bad MJD values
;;  my MJD values are off from the archive ones.
;bd = where(str.mjd lt 1000,nbd)
;print,'Fixing ',strtrim(nbd,2),' bad MJDs'
;for i=0,nbd-1 do begin
;  if (i+1) mod 100 eq 0 then print,i+1
;  fluxfile = str[bd[i]].file
;  lo = strpos(fluxfile,'archive')
;  fluxfile = mssdir+strmid(fluxfile,lo)
;  head = headfits(fluxfile,exten=0)
;  str[bd[i]].dateobs = sxpar(head,'DATE-OBS')
;  str[bd[i]].mjd = date2jd(str[bd[i]].dateobs,/mjd)
;endfor
;; Fix bad AIRMASS
;bd = where(str.airmass lt 0.9,nbd)
;for i=0,nbd-1 do begin
;  OBSERVATORY,'ctio',obs
;  lat = obs.latitude
;  lon = obs.longitude
; ;jd = str[bd].mjd + 2400000.5d0
;  jd = date2jd(str[bd[i]].dateobs)
;  ra = str[bd[i]].ra
;  dec = str[bd[i]].dec
;  str[bd[i]].airmass = AIRMASS(jd,ra,dec,lat,lon)
;endfor

;stop

; Airmass extinction term structure
amstr = replicate({filter:'',amcoef:fltarr(2)},7)
amstr.filter = ['u','g','r','i','z','Y','VR']
;ZP(u) = -0.38340071*X -1.6202186
amstr[0].amcoef = [-1.6202186,-0.38340071]
;ZP(g) = -0.19290602*X + 0.25785509
amstr[1].amcoef = [0.25785509,-0.19290602]
;ZP(r) = -0.17150109*X + 0.58451093
amstr[2].amcoef = [0.58451093,-0.17150109]
;ZP(i) = -0.10241476*X + 0.39506315
amstr[3].amcoef = [0.39506315,-0.10241476]
;ZP(z) = -0.07762681*X + 0.097165843
amstr[4].amcoef = [0.097165843,-0.07762681]
;ZP(Y) = -0.03525268*X -1.0885863
amstr[5].amcoef = [-1.0885863,-0.03525268]
;ZP(VR) = -0.091242237*X + 1.0202652
amstr[6].amcoef = [1.0202652,-0.091242237]
namstr = n_elements(amstr)

; APPLY QA CUTS IN ZEROPOINT AND SEEING
print,'APPLY QA CUTS IN ZEROPOINT AND SEEING'
fwhmthresh = 3.0  ; arcsec
;filters = ['u','g','r','i','z','Y','VR']
;nfilters = n_elements(filters)
;zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
badmask = lonarr(n_elements(str))

for i=0,namstr-1 do begin
  ind = where(str.filter eq amstr[i].filter and str.success eq 1,nind)
  print,amstr[i].filter,' ',strtrim(nind,2),' exposures'
  if nind gt 0 then begin
    zpterm = str[ind].zpterm
    am = str[ind].airmass
    mjd = str[ind].mjd
    ; Correct for airmass extinction effect
    zpterm -= poly(am,amstr[i].amcoef)
    ; Remove temporal variations
    ;gd1 = where(abs(zpterm) lt 2,ngd1)
    ;si = sort(mjd[gd1])
    ;BINDATA,mjd[gd1[si]],zpterm[gd1[si]],xbin,ybin,binsize=60,/med,gdind=gdind
    ;xbin = xbin[gdind] & ybin=ybin[gdind]
    ; this doesn't work so well

    ;sigzp = mad(str[ind].zpterm)
    ; We are using ADDITIVE zpterm 
    ;  calmag = instmag + zpterm
    ; if there are clouds then instmag is larger/fainter
    ;  and zpterm is smaller (more negative)
    ;bdind = where(str[ind].zpterm-medzp lt -zpthresh[i],nbdind)
    bdind = where(zpterm lt -zpthresh[i],nbdind)
    print,'  ',strtrim(nbdind,2),' exposures with ZPTERM below the threshold'    
    if nbdind gt 0 then badmask[ind[bdind]] = 1
  endif
endfor
; Many of the short u-band exposures have weird ZPTERMs, not sure why
bdexp = where(str.fwhm gt fwhmthresh or badmask eq 1,nbdexp)
print,'QA cuts remove ',strtrim(nbdexp,2),' exposures'
REMOVE,bdexp,str
nstr = n_elements(str)

; Which healpix pixels have data
print,'Finding the Healpix pixels with data'
radius = 1.1
healstr = replicate({file:'',base:'',pix:0L},1e5)
nhealstr = n_elements(healstr)
cnt = 0LL
for i=0,nstr-1 do begin
  if i mod 1e3 eq 0 then print,i
  ;head = headfits(str[i].file,exten=0)
  ;sra = sxpar(head,'ra')
  ;sdec = sxpar(head,'dec')
  ;ra = sexig2ten(sra)*15.0d0
  ;dec = sexig2ten(sdec)
  theta = (90-str[i].dec)/radeg
  phi = str[i].ra/radeg
  ANG2VEC,theta,phi,vec
  QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

  ; Add new elements to array
  if cnt+nlistpix gt nhealstr then begin
    old = healstr
    healstr = replicate({file:'',base:'',pix:0L},nhealstr+1e4)
    healstr[0:nhealstr-1] = old
    nhealstr += 1e4
    undefine,old
  endif

  healstr[cnt:cnt+nlistpix-1].file = str[i].expdir+'/'+str[i].base+'_cat.fits'
  healstr[cnt:cnt+nlistpix-1].base = str[i].base
  healstr[cnt:cnt+nlistpix-1].pix = listpix
  cnt += nlistpix

endfor
; Trim extra elements
healstr = healstr[0:cnt-1]
nhealstr = n_elements(healstr)

; Get uniq pixels
ui = uniq(healstr.pix,sort(healstr.pix))
upix = healstr[ui].pix
nupix = n_elements(upix)
print,strtrim(nupix,2),' Healpix pixels have overlapping data'

; Get start/stop indices for each pixel
idx = sort(healstr.pix)
healstr = healstr[idx]
q = healstr.pix
lo = where(q ne shift(q,1),nlo)
;hi = where(q ne shift(q,-1))
hi = [lo[1:nlo-1]-1,nhealstr-1]
nexp = hi-lo+1
index = replicate({pix:0L,lo:0L,hi:0L,nexp:0L},nupix)
index.pix = upix
index.lo = lo
index.hi = hi
index.nexp = nexp

; Write the full list plus an index
print,'Writing list to ',dir+'combine/nsc_healpix_list.fits'
MWRFITS,healstr,dir+'combine/nsc_healpix_list.fits',/create
MWRFITS,index,dir+'combine/nsc_healpix_list.fits',/silent
; Copy to local directory for faster reading speed
file_copy,dir+'combine/nsc_healpix_list.fits',localdir+'dnidever/nsc/instcal/',/over
; PUT NSIDE IN HEADER!!

;; Loop over each pixel and get list of overlapping exposures
;for i=0,nupix-1 do begin
;  healstr1 = healstr[idx[lo[i]:hi[i]]]
;  ; Create a file with this list
;  ;lines = healstr1.file+'  '+healstr1.base+'  '+healstr1.pix
;  WRITELINE,dir+'combine/lists/'+strtrim(upix[i],2)+'.lst',healstr1.file
;  ;MWRFITS,healstr1,dir+'combine/lists/'+strtrim(upix[i],2)+'.fits',/create
;endfor

; Make the commands
cmd = "nsc_instcal_combine,"+strtrim(upix,2)+",nside="+strtrim(nside,2)
if keyword_set(redo) then cmd+=',/redo'
cmddir = strarr(nupix)+localdir+'dnidever/nsc/instcal/tmp/'
;dirs = strarr(nupix)+'/data0/dnidever/decamcatalog/tmp/'

; Check if the output file exists
if not keyword_set(redo) then begin
  outfiles = dir+'combine/'+strtrim(upix,2)+'.fits'
  test = file_test(outfiles)
  gd = where(test eq 0,ngd,comp=bd,ncomp=nbd)
  if nbd gt 0 then begin
    print,strtrim(nbd,2),' files already exist and /redo not set.'
  endif 
  if ngd eq 0 then begin
    print,'No files to process'
    return
  endif
  print,strtrim(ngd,2),' files left to process'
  cmd = cmd[gd]
  cmddir = cmddir[gd]
endif

stop

; Now run the combination program on each healpix pixel
PBS_DAEMON,cmd,cmddir,/hyperthread,/idle,prefix='nsccmb',jobs=jobs,nmulti=nmulti,wait=1


; Load all the summary/metadata files
print,'Creating Healpix summary file'
sumstr = replicate({pix:0L,nexposures:0L,nobjects:0L,success:0},nupix)
sumstr.pix = upix
for i=0,nupix-1 do begin
  if (i+1) mod 5000 eq 0 then print,i+1
  file = dir+'combine/'+strtrim(upix[i],2)+'.fits'
  if file_test(file) eq 1 then begin
    meta = MRDFITS(file,1,/silent)
    sumstr[i].nexposures = n_elements(meta)
    hd = headfits(file,exten=2)
    sumstr[i].nobjects = sxpar(hd,'naxis2')
    sumstr[i].success = 1
  endif else begin
    sumstr[i].success = 0
  endelse
endfor
gd = where(sumstr.success eq 1,ngd)
print,strtrim(ngd,2),' Healpix successfully processed'
print,'Writing summary file to ',dir+'combine/nsc_instcal_combine.fits'
MWRFITS,expstr,dir+'combine/nsc_instcal_combine.fits',/create

; End logfile
;------------
JOURNAL



stop

end
