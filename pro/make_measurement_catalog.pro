pro make_measurement_catalog,expdir,version=version,redo=redo,stp=stp

;; Make the measurement catalog for a single exposure

t0 = systime(1)

nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(expdir) eq 0 then begin
  print,'Syntax - make_measurement_catalog,expdir,version=version,redo=redo,stp=stp'
  return
endif

; Check if output file already exists
base = file_basename(expdir)
outfile = expdir+'/'+base+'_meas.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS already and /redo not set'
  return
endif

print,'Creating measurement catalog for exposure = ',base

;;  Load the exposure and metadata files
metafile = expdir+'/'+base+'_meta.fits'
catfile = expdir+'/'+base+'_cat.fits'
meta = MRDFITS(metafile,1,/silent)
meta.file = strtrim(meta.file,2)
meta.base = strtrim(meta.base,2)
meta.dateobs = strtrim(meta.dateobs,2)
nmeta = n_elements(meta)
cat = MRDFITS(catfile,1,/silent)
ncat = n_elements(cat)

;; Convert to final format
schema = {id:'',objectid:'',instrument:'',exposure:'',ccdnum:0L,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,raerr:0.0,dec:0.0d0,decerr:0.0,$
          mag_auto:0.0,magerr_auto:0.0,mag_aper1:0.0,magerr_aper1:0.0,mag_aper2:0.0,magerr_aper2:0.0,$
          mag_aper4:0.0,magerr_aper4:0.0,mag_aper6:0.0,magerr_aper6:0.0,mag_aper8:0.0,magerr_aper8:0.0,kron_radius:0.0,$
          background:0.0,threshold:0.0,isoarea:0.0,x2:0.0,x2err:0.0,y2:0.0,y2err:0.0,xy:0.0,xyerr:0.0,$
          asemi:0.0,asemierr:0.0,bsemi:0.0,bsemierr:0.0,theta:0.0,thetaerr:0.0,fwhm:0.0,flags:0,imaflags_iso:0L,nimaflags_iso:0L,class_star:0.0}
meas = replicate(schema,ncat)
STRUCT_ASSIGN,cat,meas,/nozero
meas.id = strtrim(cat.sourceid,2)
meas.instrument = strtrim(meta.instrument,2)
meas.exposure = base
meas.x = cat.x_image
meas.y = cat.y_image
meas.mag_auto = cat.cmag
meas.magerr_auto = cat.cerr
meas.mag_aper1 = cat.mag_aper[0]
meas.magerr_aper1 = cat.magerr_aper[0]
meas.mag_aper2 = cat.mag_aper[1]
meas.magerr_aper2 = cat.magerr_aper[1]
meas.mag_aper4 = cat.mag_aper[2]
meas.magerr_aper4 = cat.magerr_aper[2]
meas.mag_aper6 = cat.mag_aper[3]
meas.magerr_aper6 = cat.magerr_aper[3]
meas.mag_aper8 = cat.mag_aper[4]
meas.magerr_aper8 = cat.magerr_aper[4]
meas.isoarea = cat.isoarea_world * 3600.^2  ; convert to arcsec^2
meas.x2 = cat.x2_world * 3600.^2            ; convert to arcsec^2
meas.x2err = cat.errx2_world * 3600.^2      ; convert to arcsec^2
meas.y2 = cat.y2_world * 3600.^2            ; convert to arcsec^2
meas.y2err = cat.erry2_world * 3600.^2      ; convert to arcsec^2
meas.xy = cat.xy_world * 3600.^2            ; convert to arcsec^2
meas.xyerr = cat.errxy_world * 3600.^2      ; convert to arcsec^2
meas.asemi = cat.a_world * 3600.            ; convert to arcsec
meas.asemierr = cat.erra_world * 3600.      ; convert to arcsec
meas.bsemi = cat.b_world * 3600.            ; convert to arcsec
meas.bsemierr = cat.errb_world * 3600.      ; convert to arcsec
meas.theta = 90-cat.theta_world             ; make CCW E of N
meas.thetaerr = cat.errtheta_world
meas.fwhm = cat.fwhm_world * 3600.          ; convert to arcsec

; possible columns to remove: background, threshold, flags,
; imaflags_iso, nimaflags_iso

;; Calibrate the aperture photometry
tags = tag_names(schema)
apertags = where(strmid(tags,0,8) eq 'MAG_APER',napertags)
for i=0,napertags-1 do begin
  ;; only calibrate "good" measurements, no 99.99
  gd = where(meas.(apertags[i]) lt 50,ngd)
  if ngd gt 0 then meas[gd].(apertags[i]) += meta.zpterm
endfor

;; Get the OBJECTID from the combined healpix file IDSTR structure
;;  remove any sources that weren't used

;; Figure out which healpix this figure overlaps

theta = (90-meas.dec)/radeg
phi = meas.ra/radeg
ANG2PIX_RING,nside,theta,phi,pix
ui = uniq(pix,sort(pix))
upix = pix[ui]
npix = n_elements(upix)

;; Load the healpix list
;listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
;if file_test(listfile) eq 0 then begin
;  print,listfile,' NOT FOUND'
;  return
;endif
;healstr = MRDFITS(listfile,1,/silent)  ; takes ~20s to load
;healstr.file = strtrim(healstr.file,2)
;healstr.base = strtrim(healstr.base,2)
;ind = where(healstr.base eq base,nind)  ; takes ~2s
;upix = healstr[ind].pix

;; Loop over the pixels
for i=0,npix-1 do begin
  objfile = dir+'combine/'+strtrim(long(upix[i])/1000,2)+'/'+strtrim(upix[i],2)+'.fits.gz'
  if file_test(objfile) eq 1 then begin
    idstr = MRDFITS(objfile,3,/silent)
    idstr.sourceid = strtrim(idstr.sourceid,2)
    nidstr = n_elements(idstr)
    MATCH,idstr.sourceid,meas.id,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then meas[ind2].objectid=strtrim(idstr[ind1].objectid,2)
    print,i+1,upix[i],nmatch,format='(I5,I10,I7)'
  endif else print,objfile,' NOT FOUND'
endfor

;; Only keep sources with an objectid
ind = where(meas.objectid ne '',nind)
if nind gt 0 then begin
  print,'Keeping ',strtrim(nind,2),' of ',strtrim(ncat,2),' sources'
  meas = meas[ind]
endif else begin
  print,'No sources to keep'
  return
endelse

;; Output
print,'Writing measurement catalog to ',outfile
MWRFITS,meas,outfile,/create

if keyword_set(stp) then stop

end
