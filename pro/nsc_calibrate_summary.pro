pro nsc_calibrate_summary,version

; Create calibration summary file

radeg = 180.0d0 / !dpi

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v3'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'

t0 = systime(1)

print, "Creating calibration summary file"

; Load the list
list = MRDFITS(dir+'/lists/nsc_calibrate_healpix_list.fits',1,/silent)
nlist = n_elements(list)
print,strtrim(nlist,2),' exposures to check'
list.expdir = strtrim(list.expdir,2)
list.instrument = strtrim(list.instrument,2)
list.filter = strtrim(list.filter,2)

;; Chip offsets
chipoff = mrdfits(dldir+'users/dnidever/nsc/instcal/decam_chip_xyoff.fits',1)

; Load all the summary/metadata files
expstr = replicate({expdir:'',instrument:'',pix:0L,metafile:'',metadate:0LL,success:0,file:'',wtfile:'',maskfile:'',base:'',expnum:'',ra:0.0d0,dec:0.0d0,$
                    dateobs:'',mjd:0.0d,filter:'',exptime:0.0,plver:'',airmass:0.0,meas_logfile_success:0,meas_success:0,meas_dt:0.0,meas_chip1date:0LL,$
                    meas_logdate:0LL,nsources:0L,ngoodsources:-1L,fwhm:999.0,nchips:0L,rarms:999999.0,$
                    rastderr:999999.0,decrms:999999.0,decstderr:999999.0,ebv:0.0,ngaiamatch:0L,ngoodgaiamatch:-1L,zptype:1,zpterm:999999.0,$
                    zptermerr:99999.0,zptermsig:999999.0,zpspatialvar_rms:999999.0,zpspatialvar_range:999999.0,zpspatialvar_nccd:0,nrefmatch:0L,$
                    ngoodrefmatch:0L,depth95:99.99,depth10sig:99.99},nlist)
expstr.expdir = list.expdir
expstr.base = file_basename(list.expdir)
expstr.instrument = list.instrument
expstr.filter = list.filter
expstr.pix = list.pix

; Load the list of exposures
explist1 = MRDFITS(dir+'/lists/decam_instcal_list.fits.gz',1)
explist2 = MRDFITS(dir+'/lists/mosaic3_instcal_list.fits.gz',1)
explist3 = MRDFITS(dir+'/lists/bok90prime_instcal_list.fits.gz',1)
explist = [explist1,explist2,explist3]
undefine,explist1,explist2,explist3
explist.fluxfile = strtrim(explist.fluxfile,2)
MATCH,expstr.base,file_basename(explist.fluxfile,'.fits.fz'),ind1,ind2,/sort,count=nmatch
expstr[ind1].file = explist[ind2].fluxfile
expstr[ind1].wtfile = explist[ind2].wtfile
expstr[ind1].maskfile = explist[ind2].maskfile
expstr[ind1].expnum = strtrim(explist[ind2].expnum,2)
expstr[ind1].ra = explist[ind2].ra
expstr[ind1].dec = explist[ind2].dec
expstr[ind1].exptime = explist[ind2].exposure
expstr[ind1].dateobs = explist[ind2].date_obs
expstr[ind1].mjd = explist[ind2].mjd_obs
expstr[ind1].plver = explist[ind2].plver

;; Match to MEAS SUMMARY FILE
meas = mrdfits(dir+'lists/nsc_measure_summary.fits',1)
meas.base = strtrim(meas.base,2)
MATCH,meas.base,expstr.base,ind1,ind2,/sort
expstr[ind2].nchips = meas[ind1].nchips
expstr[ind2].nsources = meas[ind1].nsources
expstr[ind2].meas_logfile_success = meas[ind1].logfile_success
expstr[ind2].meas_success = meas[ind1].success
expstr[ind2].meas_dt = meas[ind1].dt
expstr[ind2].meas_chip1date = meas[ind1].chip1date
expstr[ind2].meas_logdate = meas[ind1].logdate

;; Chip structure
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

  ;; No metafile
  endif else begin
    expstr[i].success = 0

    ;; Calculate FWHM using one chip catalog
    ;;  28 and 35 are at the center, use them to central position
    chips = [28,35,1,indgen(58)+3,62]
    medfwhm = 999.0
    flag = 0
    ch = 0
    undefine,cat
    ncat = 0
    while (flag eq 0) do begin
      chipfile = expstr[i].expdir+expstr[i].base+'_'+strtrim(chips[ch],2)+'.fits'
      if file_test(chipfile) eq 1 then begin
        hd = headfits(chipfile,exten=2)
        ncat = sxpar(hd,'naxis2')
        if ncat gt 5 then cat=mrdfits(chipfile,2,/silent)
      endif
      ch++
      if n_elements(cat) gt 0 or ch gt 59 then flag=1
    endwhile
    if n_elements(cat) gt 0 then begin
      gdcat = where(cat.mag_auto lt 50 and cat.magerr_auto lt 0.05 and cat.class_star gt 0.8,ngdcat)
      if ngdcat eq 0 then gdcat = where(cat.mag_auto lt 50 and cat.magerr_auto lt 0.10,ngdcat)
      if ngdcat gt 5 then begin
        medfwhm = median(cat[gdcat].fwhm_world*3600.)
        expstr[i].fwhm = medfwhm
        raoff = 0.0
        decoff = 0.0
        MATCH,chipoff.chip,chips[ch],ind1,ind2,/sort,count=nmatch
        if expstr[i].instrument eq 'c4d' and nmatch gt 0 then begin
          raoff = chipoff[ind1].xoff
          decoff = chipoff[ind1].yoff
        endif 
        cenra = median(cat[gdcat].alpha_j2000)+raoff
        cendec = median(cat[gdcat].delta_j2000)+decoff
        ;; only use if archive value seems off
        if sphdist(expstr[i].ra,expstr[i].dec,cenra,cendec,/deg) gt 0.5 then begin
          expstr[i].ra = cenra
          expstr[i].dec = cendec
        endif
      endif
    endif
  endelse  ; no metafile
endfor
; Trim CHSTR structure
chstr = chstr[0:chcnt-1]
gd = where(expstr.success eq 1,ngd)
print,strtrim(ngd,2),' exposure successfully calibrated'
print,'Writing summary file to ',dir+'lists/nsc_instcal_calibrate.fits'
MWRFITS,expstr,dir+'lists/nsc_instcal_calibrate.fits',/create
MWRFITS,chstr,dir+'lists/nsc_instcal_calibrate.fits',/silent
if file_test(dir+'lists/nsc_instcal_calibrate.fits.gz') eq 1 then file_delete,dir+'lists/nsc_instcal_calibrate.fits.gz',/allow
spawn,['gzip',dir+'lists/nsc_instcal_calibrate.fits'],/noshell

print,'dt=',stringize(systime(1)-t0,ndec=2),' sec'

stop

end
