pro combine_stripe82_uband_v3

; Combine together the exposure catalogs in Stripe82 for a single band

filter = 'u'

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
restore,'/dl1/users/dnidever/nsc/smash_matched_catalog_v3.dat'
ref = mobj
undefine,mobj

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
  if nmatch eq 0 then goto,BOMB
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
;save,allcat,allref,file='/dl1/users/dnidever/nsc/instcal/t3b/combine/combine_stripe82_uband_v3.dat'

stop


; Make the plot
!p.font = 0
setdisp
plotdir = '/dl1/users/dnidever/nsc/instcal/t3b/plots/'

file = plotdir+'stripe82_uband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
gj0 = allref.gaia_gmag - allref.tmass_jmag - 1.12*allref.ebv
model_mag = 0.2452*allref.galex_nuv + 0.7486*allref.gaia_gmag + 0.3717*gj0 + 0.350*allref.ebv + 0.1352
gd = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,gj0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.5,2.0],yr=[-1,1],/log,xtit='(G-J)o',ytit='Model-Mag',tit='u-band'
bindata,gj0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0.5,max=1.5
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.8 and xbin le 1.1,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
; -0.0262557    0.0196494
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.8,0.8],[-2,2],linestyle=1,co=255
oplot,[1.1,1.1],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(G-J)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Get extinction part
gd1 = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0 and gj0 ge 0.8 and gj0 lt 1.2 and allcat.ebv gt 0.15,ngd1)
a = dblarr(4,ngd1)
a[0,*] = allref[gd1].galex_nuv 
a[1,*] = allref[gd1].gaia_gmag
a[2,*] = allcat[gd1].ebv
a[3,*] = 1
b=allcat[gd1].cmag
SVDC, A, W, U, V
factor = SVSOL(U, W, V, B)
print,factor
; 0.27682490      0.73827996      0.39430181     0.052752100

; Now nail down extinction and get the other terms
gd1 = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0 and gj0 ge 0.8 and gj0 lt 1.2 and allcat.ebv lt 0.1,ngd1)
a = dblarr(4,ngd1)
a[0,*] = allref[gd1].galex_nuv
a[1,*] = allref[gd1].gaia_gmag
a[2,*] = gj0[gd1]
a[3,*] = 1
b=allcat[gd1].cmag-0.3943*allcat[gd1].ebv
SVDC, A, W, U, V
factor = SVSOL(U, W, V, B)
print,factor
; 0.23011191      0.76160673      0.49370716      0.13443835

; versus EBV
file = plotdir+'stripe82_uband_magdiff_ebv'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
gd2 = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0 and gj0 lt 1.2,ngd2)
hess,allcat[gd2].ebv,model_mag[gd2]-allcat[gd2].cmag,dx=0.01,dy=0.02,xr=[0,0.8],yr=[-1,1],/log,xtit='E(B-V)',ytit='Model-Mag',tit='u-band'
oplot,[-1,3],[0,0],linestyle=2,co=255
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

;; Now use the new equation
file = plotdir+'stripe82_uband_magdiff_color_adjusted'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
gj0 = allref.gaia_gmag - allref.tmass_jmag - 1.12*allref.ebv
; ORIGINAL: 0.2452*NUV+0.7486*GMAG+0.3717*COLOR+0.350*EBV+0.1352
; ADJUSTED: 0.2301*NUV+0.7616*GMAG+0.4937*COLOR+0.3943*EBV+0.1344
model_mag = 0.2246*allref.galex_nuv + 0.7703*allref.gaia_gmag + 0.4420*gj0 + 0.3799*allref.ebv + 0.1624
gd = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,gj0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.5,2.0],yr=[-1,1],/log,xtit='(G-J)o',ytit='Model-Mag',tit='u-band'
bindata,gj0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0.5,max=1.5
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.8 and xbin le 1.1,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;  0.16465859     -0.21008923
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.8,0.8],[-2,2],linestyle=1,co=255
oplot,[1.1,1.1],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(G-J)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

pdfcombine,plots+'.pdf',plotdir+'stripe82_uband_combine.pdf',/clobber

stop

end
