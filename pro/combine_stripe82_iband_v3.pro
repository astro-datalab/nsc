pro combine_stripe82_iband_v3

; Combine together the exposure catalogs in Stripe82 for a single band

filter = 'i'

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+"users/dnidever/nsc/"
str = mrdfits(dir+'instcal/t3b/lists/nsc_instcal_calibrate.fits',1)
;str.expdir = file_dirname(strtrim(str.expdir,2))
str.expdir = strtrim(str.expdir,2)
str.filter = strtrim(str.filter,2)
str.expnum = strtrim(str.expnum,2)

;outfile = str.expdir+'/'+file_basename(str.expdir)+'_cat.fits'
;ind0 = where((str.ra lt 61 or str.ra gt 299) and abs(str.dec) lt 3.0 and str.zpterm ne 0 and str.filter eq filter,nind0)
;medzpterm = median(str[ind0].zpterm)
;sigzpterm = mad(str[ind0].zpterm)
;ind = where((str.ra lt 61 or str.ra gt 299) and abs(str.dec) lt 3.0 and str.zpterm ne 0 and str.filter eq filter and $
;            abs(str.zpterm-medzpterm) lt 3*sigzpterm,nind)
ind = where(str.filter eq filter,nind)

;test = file_test(outfile[gd])
print,strtrim(nind,2),' exposures for BAND=',filter

;; Load the reference data
ref1 = mrdfits('/dl1/users/dnidever/nsc/Stripe82_v3_ejk.fits',1)
ref2 = mrdfits('/dl1/users/dnidever/nsc/Stripe82_v3_ejk_midplane.fits',1)
;ref = [ref,ref2]


for i=0,nind-1 do begin
  expdir = str[ind[i]].expdir
  base = file_basename(str[ind[i]].expdir)
  ;file = repstr(strtrim(str[ind[i]].metafile,2),'meta','cat')
  metafile = str[ind[i]].expdir+'/'+base+'_meta.fits'
  catfile = str[ind[i]].expdir+'/'+base+'_cat.fits'

  ;; Individual meas files
  if file_test(metafile) eq 1 and file_test(catfile) eq 0 then begin
    measfiles = file_search(expdir+'/'+base+'_*_meas.fits',count=nmeasfiles)
    undefine,cat
    ;; Get the number of sources first
    ncat = 0L
    for j=0,nmeasfiles-1 do begin
      head = headfits(measfiles[j],exten=1)
      ncat += sxpar(head,'NAXIS2')
    endfor
    cat0 = MRDFITS(measfiles[0],1,/silent)
    schema = cat0[0]
    struct_assign,{dum:''},schema
    cat = REPLICATE(schema,ncat)
    catcnt = 0LL
    ;; Now load the data
    for j=0,nmeasfiles-1 do begin
      cat1 = MRDFITS(measfiles[j],1,/silent)
      ncat1 = n_elements(cat1)
      cat[catcnt:catcnt+ncat1-1] = cat1
      catcnt += ncat1
    endfor

  ;; Single old calibrated photometry file
  endif else begin
    cat1 = mrdfits(catfile,1,/silent)
    ncat1 = n_elements(cat1)
    expstr = mrdfits(metafile,1,/silent)
    
    ;; Make a cut on quality mask flag (IMAFLAGS_ISO)
    bdcat = where(cat1.imaflags_iso gt 0,nbdcat)
    if nbdcat gt 0 then begin
      if nbdcat eq ncat1 then goto,BOMB
      REMOVE,bdcat,cat1
      ncat1 = n_elements(cat1)
    endif  
    ;; Make cuts on SE FLAGS
    ;;   this removes problematic truncatd sources near chip edges
    bdseflags = where( ((cat1.flags and 8) eq 8) or $             ; object truncated
                     ((cat1.flags and 16) eq 16),nbdseflags)    ; aperture truncated
    if nbdseflags gt 0 then begin
      if nbdseflags eq ncat1 then goto,BOMB
      REMOVE,bdseflags,cat1
      ncat1 = n_elements(cat1)
    endif
    ;; Removing low-S/N sources
    ;;  snr = 1.087/err
    snrcut = 5.0
    bdsnr = where(1.087/cat1.magerr_auto lt snrcut,nbdsnr)
    if nbdsnr gt 0 then begin
      if nbdsnr eq ncat1 then goto,BOMB
      REMOVE,bdsnr,cat1
      ncat1 = n_elements(cat1)
    endif
    ;; Convert to measurement format
    schema = {measid:'',objectid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,raerr:0.0,dec:0.0d0,decerr:0.0,$
              mag_auto:0.0,magerr_auto:0.0,mag_aper1:0.0,magerr_aper1:0.0,mag_aper2:0.0,magerr_aper2:0.0,$
              mag_aper4:0.0,magerr_aper4:0.0,mag_aper8:0.0,magerr_aper8:0.0,kron_radius:0.0,$
              asemi:0.0,asemierr:0.0,bsemi:0.0,bsemierr:0.0,theta:0.0,thetaerr:0.0,fwhm:0.0,flags:0}
    cat = replicate(schema,ncat1)
    STRUCT_ASSIGN,cat1,cat,/nozero
    cat.measid = strtrim(cat1.sourceid,2)
    cat.exposure = base
    cat.x = cat1.x_image
    cat.y = cat1.y_image
    cat.mag_auto = cat1.cmag
    cat.magerr_auto = cat1.cerr
    cat.mag_aper1 = cat1.mag_aper[0] + 2.5*alog10(expstr.exptime) + expstr.zpterm
    cat.magerr_aper1 = cat1.magerr_aper[0]
    cat.mag_aper2 = cat1.mag_aper[1] + 2.5*alog10(expstr.exptime) + expstr.zpterm
    cat.magerr_aper2 = cat1.magerr_aper[1]
    cat.mag_aper4 = cat1.mag_aper[2] + 2.5*alog10(expstr.exptime) + expstr.zpterm
    cat.magerr_aper4 = cat1.magerr_aper[2]
    cat.mag_aper8 = cat1.mag_aper[4] + 2.5*alog10(expstr.exptime) + expstr.zpterm
    cat.magerr_aper8 = cat1.magerr_aper[4]
    cat.asemi = cat1.a_world * 3600.            ; convert to arcsec
    cat.asemierr = cat1.erra_world * 3600.      ; convert to arcsec
    cat.bsemi = cat1.b_world * 3600.            ; convert to arcsec
    cat.bsemierr = cat1.errb_world * 3600.      ; convert to arcsec
    cat.theta = 90-cat1.theta_world             ; make CCW E of N
    cat.thetaerr = cat1.errtheta_world
    cat.fwhm = cat1.fwhm_world * 3600.          ; convert to arcsec
    ;add_tag,cat,'expnum','',cat
    ;cat.exposure = str[ind[i]].expnum
  endelse
  ncat = n_elements(cat)
  print,strtrim(i+1,2),' ',base,' ',str[ind[i]].base,' ',strtrim(ncat,2)

  cendec = mean(minmax(cat.dec))
  cenra = mean(minmax(cat.ra))
  if range(cat.ra) gt 100 then begin
    ra = cat.ra
    bd = where(ra gt 180,nbd)
    if nbd gt 0 then ra[bd]-=360
    cenra = mean(minmax(ra))
    if cenra lt 0 then cenra+=360
  endif
  glactc,cenra,cendec,2000.0,cengl,cengb,1,/deg


  ; Matching
  dcr = 1.0
  if abs(cendec) lt 5 then begin
    SRCMATCH,ref1.ra,ref1.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
    print,'  ',strtrim(nmatch,2),' matches to reference data'
    if nmatch eq 0 then goto,BOMB
    newref = ref1[ind1]
    newcat = cat[ind2]
  endif else begin
    SRCMATCH,ref2.ra,ref2.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
    print,'  ',strtrim(nmatch,2),' matches to reference data'
    if nmatch eq 0 then goto,BOMB
    newref = ref2[ind1]
    newcat = cat[ind2]
  endelse

  if n_elements(allcat) eq 0 then begin
    cat0 = cat[0]
    struct_assign,{dum:''},cat0
    allcat = replicate(cat0,2e7)
    ref0 = ref1[0]
    struct_assign,{dum:''},ref0
    allref = replicate(ref0,2e7)
    cnt = 0LL
  endif
  tempcat = allcat[cnt:cnt+nmatch-1]
  struct_assign,newcat,tempcat
  allcat[cnt:cnt+nmatch-1] = tempcat
  tempref = allref[cnt:cnt+nmatch-1]
  struct_assign,newref,tempref
  allref[cnt:cnt+nmatch-1] = tempref
  cnt += nmatch

  BOMB:
endfor
; Trim extra elements
allcat = allcat[0:cnt-1]
allref = allref[0:cnt-1]

; Save the matched catalogs
;save,allcat,allref,file='/dl1/users/dnidever/nsc/instcal/t3b/combine/combine_stripe82_iband_v3.dat'

stop


; Make the plot
!p.font = 0
setdisp
plotdir = '/dl1/users/dnidever/nsc/instcal/t3b/plots/'

file = plotdir+'stripe82_iband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = allref.tmass_jmag-allref.tmass_kmag-0.17*allcat.ebv
; SM_IMAG+0.041*COLOR+0.010*EBV-0.003
model_mag = allref.sm_imag + 0.041*jk0 + 0.010*allref.ebv - 0.003
gd = where(allcat.class_star gt 0.8 and allref.tmass_phqual eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='i-band'
bindata,jk0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.3 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
; -0.0970436    0.0542384
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.3,0.3],[-2,2],linestyle=1,co=255
oplot,[0.7,0.7],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; versus EBV
file = plotdir+'stripe82_iband_magdiff_ebv'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
hess,allcat[gd].ebv,model_mag[gd]-allcat[gd].cmag,dx=0.01,dy=0.02,xr=[0,0.8],yr=[-1,1],/log,xtit='E(B-V)',ytit='Model-Mag',tit='i-band'
oplot,[-1,3],[0,0],linestyle=2,co=255
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

;; Now use the new equation
file = plotdir+'stripe82_iband_magdiff_color_adjusted'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = allref.tmass_jmag-allref.tmass_kmag-0.17*allcat.ebv
; ORIGINAL: SM_IMAG+0.041*COLOR+0.010*EBV-0.003
; ADJUSTED: SM_IMAG-0.0132*COLOR+0.010*EBV+0.0940
model_mag = allref.sm_imag - 0.0132*jk0 + 0.010*allref.ebv + 0.0940
gd = where(allcat.class_star gt 0.8 and allref.tmass_phqual eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='ADJUSTED i-band'
bindata,jk0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.3 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.3,0.3],[-2,2],linestyle=1,co=255
oplot,[0.7,0.7],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

pdfcombine,plots+'.pdf',plotdir+'stripe82_iband_combine.pdf',/clobber

stop

end
