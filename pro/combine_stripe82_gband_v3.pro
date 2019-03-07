pro combine_stripe82_gband_v3

; Combine together the exposure catalogs in Stripe82 for a single band

filter = 'g'

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+"users/dnidever/nsc/"
str = mrdfits(dir+'instcal/t3b/lists/nsc_instcal_calibrate.fits',1)
str.expdir = file_dirname(strtrim(str.expdir,2))
;str.expdir = strtrim(str.expdir,2)
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
ref = mrdfits('/dl1/users/dnidever/nsc/Stripe82_v3.fits',1)

for i=0,nind-1 do begin
  base = file_basename(str[ind[i]].expdir)
  ;file = repstr(strtrim(str[ind[i]].metafile,2),'meta','cat')
  file = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir,'_meta.fits')+'_cat.fits'
  if file_test(file) eq 0 then begin
    print,file,' NOT FOUND'
    goto,BOMB
  endif
  cat = mrdfits(file,1,/silent)
  ncat = n_elements(cat)
  add_tag,cat,'expnum','',cat
  cat.expnum = str[ind[i]].expnum
  print,strtrim(i+1,2),' ',base,' ',str[ind[i]].expnum,' ',strtrim(ncat,2)

  cendec = mean(minmax(cat.dec))
  cenra = mean(minmax(cat.ra))
  if range(cat.ra) gt 100 then begin
    ra = cat.ra
    bd = where(ra gt 180,nbd)
    if nbd gt 0 then ra[bd]-=360
    cenra = mean(minmax(ra))
    if cenra lt 0 then cenra+=360
  endif

  ;; Load the Gaia file
  ;gaiafile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_GAIA.fits'
  ;if file_test(gaiafile) eq 0 then goto,BOMB
  ;gaia = MRDFITS(gaiafile,1,/silent)
  ;
  ;; Load the 2MASS file
  ;tmassfile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_TMASS.fits'
  ;if file_test(tmassfile) eq 0 then goto,BOMB
  ;tmass = MRDFITS(tmassfile,1,/silent)
  ;
  ;; Load the APASS file
  ;apassfile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_APASS.fits'
  ;if file_test(apassfile) eq 0 then goto,BOMB
  ;apass = MRDFITS(apassfile,1,/silent)

  ; Matching
  index = lonarr(ncat,3)-1
  dcr = 1.0
  ;SRCMATCH,gaia.ra_icrs,gaia.de_icrs,cat.ra,cat.dec,dcr,gind1,gind2,/sph,count=ngmatch
  ;if ngmatch gt 0 then index[gind2,0] = gind1
  ;SRCMATCH,tmass.raj2000,tmass.dej2000,cat.ra,cat.dec,dcr,tind1,tind2,/sph,count=ntmatch
  ;if ntmatch gt 0 then index[tind2,1] = tind1
  ;SRCMATCH,apass.raj2000,apass.dej2000,cat.ra,cat.dec,dcr,aind1,aind2,/sph,count=namatch
  ;if namatch gt 0 then index[aind2,2] = aind1
  ;gd = where(total(index gt -1,2) eq 3,ngd)
  ;print,'  ',strtrim(ngd,2),' matches to GAIA, 2MASS and APASS'
  ;if ngd eq 0 then begin
  ;  print,'No matches to GAIA, 2MASS and APASS'
  ;  goto,BOMB
  ;endif
  ;cat1 = cat[gd]
  ;gaia1 = gaia[index[gd,0]]
  ;tmass1 = tmass[index[gd,1]]
  ;apass1 = apass[index[gd,2]]
  SRCMATCH,ref.ra,ref.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
  print,'  ',strtrim(nmatch,2),' matches to reference data'

  ;; Get reference data
  if nmatch eq 0 then begin
    filters='c4d-'+['u','g','r','i','z','Y','VR']
    ref0 = GETREFDATA(filters,cenra,cendec,1.8)
    ;; gets ps, gaia, tmass, allwise, galex, ejk
    ;; put ref data in correct format
    refschema = ref[0]
    struct_assign,{dum:''},refschema
    ;; make all 99.99 by default
    for k=0,n_tags(refschema)-1 do begin
      type = size(refschema.(k),/type)
      if type eq 4 or type eq 5 then refschema.(k)=99.99
    endfor
    ref1 = replicate(refschema,n_elements(ref0))
    ref1.ra = ref0.ra
    ref1.dec = ref0.dec
    ref1.ps1_gmag = ref0.ps_gmag
    ref1.ps1_rmag = ref0.ps_rmag
    ref1.ps1_imag = ref0.ps_imag
    ref1.ps1_zmag = ref0.ps_zmag
    ref1.ps1_ymag = ref0.ps_ymag
    ref1.gaia_gmag = ref0.gmag
    ref1.gaia_gerr = ref0.e_gmag
    ref1.gaia_bp = ref0.bp
    ref1.gaia_bperr = ref0.e_bp
    ref1.gaia_rp = ref0.rp
    ref1.gaia_rperr = ref0.e_rp
    ref1.tmass_jmag = ref0.jmag
    ref1.tmass_jerr = ref0.e_jmag
    ref1.tmass_hmag = ref0.hmag
    ref1.tmass_herr = ref0.e_hmag
    ref1.tmass_kmag = ref0.kmag
    ref1.tmass_kerr = ref0.e_kmag
    ref1.tmass_phqual = ref0.qflg
    ref1.galex_nuv = ref0.nuv
    ref1.galex_nuverr = ref0.e_nuv
    ;; apass, sm, allwise
    ref1.allwise_w1mag = ref0.w1mag
    ref1.allwise_w1err = ref0.e_w1mag
    ref1.allwise_w2mag = ref0.w2mag
    ref1.allwise_w2err = ref0.e_w2mag
    ref1.ebv = ref0.ebv_sfd
    ref1.ejk = ref0.ejk
    ref1.e_ejk = ref0.e_jk
    ref1.ext_type = ref0.ext_type
    ;; Add Skymapper data
    ;; need to get apass and skymapper
    ;apass0 = GETREFCAT(cenra,cendec,'APASS')
    sm = GETREFCAT(cenra,cendec,1.8,'SKYMAPPER')
    SRCMATCH,ref1.ra,ref1.dec,sm.raj2000,sm.dej2000,1.0,ind1,ind2,/sph,coun=nmatch
    if nmatch gt 0 then begin
      ref1[ind1].sm_gmag = sm[ind2].sm_gmag
      ref1[ind1].sm_gerr = sm[ind2].e_sm_gmag
      ref1[ind1].sm_rmag = sm[ind2].sm_rmag
      ref1[ind1].sm_rerr = sm[ind2].e_sm_rmag
      ref1[ind1].sm_imag = sm[ind2].sm_imag
      ref1[ind1].sm_ierr = sm[ind2].e_sm_imag
      ref1[ind1].sm_zmag = sm[ind2].sm_zmag
      ref1[ind1].sm_zerr = sm[ind2].e_sm_zmag
    endif

    ;; Now match to the DECam exposure
    SRCMATCH,ref1.ra,ref1.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
    if nmatch eq 0 then goto,BOMB

stop
  endif
  ref1 = ref[ind1]
  cat1 = cat[ind2]


  if n_elements(allcat) eq 0 then begin
    cat0 = cat[0]
    struct_assign,{dum:''},cat0
    allcat = replicate(cat0,2e7)
    ref0 = ref[0]
    struct_assign,{dum:''},ref0
    allref = replicate(ref0,2e7)
    cnt = 0LL
  endif
  tempcat = allcat[cnt:cnt+nmatch-1]
  struct_assign,cat1,tempcat
  allcat[cnt:cnt+nmatch-1] = tempcat
  tempref = allref[cnt:cnt+nmatch-1]
  struct_assign,ref1,tempref
  allref[cnt:cnt+nmatch-1] = tempref
  cnt += nmatch

  ;stop
  BOMB:
endfor
; Trim extra elements
allcat = allcat[0:cnt-1]
allref = allref[0:cnt-1]

; Save the matched catalogs
;save,allcat,allref,file='/dl1/users/dnidever/nsc/instcal/t3b/combine/combine_stripe82_gband_v3.dat'

stop


; Make the plot
!p.font = 0
setdisp
plotdir = '/dl1/users/dnidever/nsc/instcal/t3b/plots/'

file = plotdir+'stripe82_gband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = allref.tmass_jmag-allref.tmass_kmag-0.17*allcat.ebv
; old version
;model_mag = allapass.g_mag - 0.1433*jk0 - 0.05*allcat.ebv - 0.0138
; new version
;model_mag = allapass.g_mag - 0.0421*jk0 - 0.05*allcat.ebv - 0.0620
; SM_GMAG+0.229*COLOR+0.150*EBV-0.013
model_mag = allref.sm_gmag + 0.229*jk0 + 0.150*allref.ebv - 0.013
gd = where(allcat.class_star gt 0.8 and allref.tmass_phqual eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='g-band'
bindata,jk0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.3 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
; 0.0598187   -0.0947644
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
file = plotdir+'stripe82_gband_magdiff_ebv'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
hess,allcat[gd].ebv,model_mag[gd]-allcat[gd].cmag,dx=0.01,dy=0.02,xr=[0,0.8],yr=[-1,1],/log,xtit='E(B-V)',ytit='Model-Mag',tit='g-band'
oplot,[-1,3],[0,0],linestyle=2,co=255
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

;; Now use the new equation
file = plotdir+'stripe82_gband_magdiff_color_adjusted'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = allref.tmass_jmag-allref.tmass_kmag-0.17*allcat.ebv
; ORIGINAL: SM_GMAG+0.229*COLOR+0.150*EBV-0.013
;model_mag = allref.sm_gmag + 0.229*jk0 + 0.150*allref.ebv - 0.013
; ADJUSTED: SM_GMAG+0.324*COLOR+0.150*EBV-0.073
model_mag = allref.sm_gmag + 0.324*jk0 + 0.150*allref.ebv - 0.073
gd = where(allcat.class_star gt 0.8 and allref.tmass_phqual eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='ADJUSTED g-band'
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

pdfcombine,plots+'.pdf',plotdir+'stripe82_gband_combine.pdf',/clobber

stop

end
