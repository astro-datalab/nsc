pro nsc_calibrate_summary

; Create calibration summary file

radeg = 180.0d0 / !dpi

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'

t0 = systime(1)

print, "Creating calibration summary file"

; Load the list
list = MRDFITS(dir+'/lists/nsc_calibrate_healpix_list.fits',1,/silent)
nlist = n_elements(list)
list.expdir = strtrim(list.expdir,2)
list.instrument = strtrim(list.instrument,2)
list.filter = strtrim(list.filter,2)

; Load all the summary/metadata files
expstr = replicate({expdir:'',instrument:'',pix:0L,metafile:'',metadate:0LL,success:0,file:'',wtfile:'',maskfile:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d0,dateobs:'',$
                    mjd:0.0d,filter:'',exptime:0.0,airmass:0.0,nsources:0L,ngoodsources:-1L,fwhm:0.0,nchips:0L,rarms:999999.0,rastderr:999999.0,decrms:999999.0,decstderr:999999.0,$
                    ebv:0.0,ngaiamatch:0L,ngoodgaiamatch:-1L,zptype:1,zpterm:999999.0,zptermerr:99999.0,zptermsig:999999.0,zpspatialvar_rms:999999.0,$
                    zpspatialvar_range:999999.0,zpspatialvar_nccd:0,nrefmatch:0L,ngoodrefmatch:0L,depth95:99.99,depth10sig:99.99},nlist)
expstr.expdir = list.expdir
expstr.instrument = list.instrument
expstr.filter = list.filter
expstr.pix = list.pix
chstr = replicate({expdir:'',instrument:'',success:0,filename:'',ccdnum:0L,nsources:0L,cenra:999999.0d0,cendec:999999.0d0,$
                   gaianmatch:0L,gaiagoodnmatch:0L,rarms:999999.0,rastderr:999999.0,racoef:dblarr(4),decrms:999999.0,$
                   decstderr:999999.0,deccoef:dblarr(4),vra:dblarr(4),vdec:dblarr(4),zpterm:999999.0,$
                   zptermerr:999999.0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nlist*61)
chcnt = 0LL
for i=0,nlist-1 do begin
  if (i+1) mod 5000 eq 0 then print,i+1
  base = file_basename(expstr[i].expdir)

  ; Getting observatory name
  case expstr[i].instrument of
  'c4d': obsname='ctio'
  'k4m': obsname='kpno'
  'ksb': obsname='kpno'
  else: obsname='ctio'
  endcase

  metafile = expstr[i].expdir+'/'+base+'_meta.fits'
  lo = strpos(metafile,'/dl1')
  metafile = dldir + strmid(metafile,lo+5)
  expstr[i].metafile = metafile
  if file_test(metafile) eq 1 then begin
    metainfo = file_info(metafile)
    expstr[i].metadate = metainfo.ctime
    ; Exposure structure
    expstr1 = MRDFITS(metafile,1,/silent)
    temp = expstr[i]
    STRUCT_ASSIGN,expstr1,temp,/nozero
    ; old versions
    if tag_exist(expstr1,'gaianmatch') then temp.ngaiamatch=expstr1.gaianmatch
    if tag_exist(expstr1,'gaiagoodnmatch') then temp.ngoodgaiamatch=expstr1.gaiagoodnmatch
    if tag_exist(expstr1,'nrefgdmatch') then temp.ngoodrefmatch=expstr1.nrefgdmatch
    expstr[i] = temp
    expstr[i].success = 1
    ; Chip structure
    chstr1 = MRDFITS(metafile,2,/silent)
    nchstr1 = n_elements(chstr1)
    temp = chstr[chcnt:chcnt+nchstr1-1]
    STRUCT_ASSIGN,chstr1,temp,/nozero
    chstr[chcnt:chcnt+nchstr1-1] = temp
    chstr[chcnt:chcnt+nchstr1-1].expdir = expstr[i].expdir
    chstr[chcnt:chcnt+nchstr1-1].instrument = expstr[i].instrument
    chstr[chcnt:chcnt+nchstr1-1].success = 1
    chcnt += nchstr1
    ; Add RASTDERR/DECSTDERR to expstr
    if tag_exist(chstr1,'RASTDERR') then begin
      expstr[i].rastderr = median(chstr1.rastderr)
      expstr[i].decstderr = median(chstr1.decstderr)
    endif    

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

stop

end
