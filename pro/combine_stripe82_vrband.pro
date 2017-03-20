pro combine_stripe82_vrband

; Combine together the exposure catalogs in Stripe82 for a single band

filter = 'VR'

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+"users/dnidever/nsc/"
str = mrdfits(dir+'instcal/nsc_instcal_calibrate.fits',1)
str.expdir = strtrim(str.expdir,2)
str.filter = strtrim(str.filter,2)
str.expnum = strtrim(str.expnum,2)

outfile = str.expdir+'/'+file_basename(str.expdir)+'_cat.fits'
ind0 = where(str.zpterm ne 0 and str.filter eq filter,nind0)
medzpterm = median(str[ind0].zpterm)
sigzpterm = mad(str[ind0].zpterm)
ind = where(abs(str.dec) lt 23 and str.exptime gt 60 and str.filter eq filter and abs(str.zpterm-medzpterm) lt 3*sigzpterm,nind)
; only 19 exposures in stripe82, so expanding the area

;test = file_test(outfile[gd])
print,strtrim(nind,2),' exposures for BAND=',filter

for i=0,nind-1 do begin
  base = file_basename(str[ind[i]].expdir)
  ;file = repstr(strtrim(str[ind[i]].metafile,2),'meta','cat')
  file = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_cat.fits'
  if file_test(file) eq 0 then begin
    print,file,' NOT FOUND'
    goto,BOMB
  endif
  cat = mrdfits(file,1,/silent)
  ncat = n_elements(cat)
  add_tag,cat,'expnum','',cat
  cat.expnum = str[ind[i]].expnum
  print,strtrim(i+1,2),' ',base,' ',str[ind[i]].expnum,' ',strtrim(ncat,2)

  ; Load the Gaia file
  gaiafile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_GAIA.fits'
  if file_test(gaiafile) eq 0 then goto,BOMB
  gaia = MRDFITS(gaiafile,1,/silent)

  ; Load the 2MASS file
  tmassfile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_TMASS.fits'
  if file_test(tmassfile) eq 0 then goto,BOMB
  tmass = MRDFITS(tmassfile,1,/silent)

  ; Matching
  index = lonarr(ncat,2)-1
  dcr = 1.0
  SRCMATCH,gaia.ra_icrs,gaia.de_icrs,cat.ra,cat.dec,dcr,gind1,gind2,/sph,count=ngmatch
  if ngmatch gt 0 then index[gind2,0] = gind1
  SRCMATCH,tmass.raj2000,tmass.dej2000,cat.ra,cat.dec,dcr,tind1,tind2,/sph,count=ntmatch
  if ntmatch gt 0 then index[tind2,1] = tind1
  gd = where(total(index gt -1,2) eq 2,ngd)
  print,'  ',strtrim(ngd,2),' matches to GAIA and 2MASS'
  if ngd eq 0 then begin
    print,'No matches to GAIA and 2MASS'
    goto,BOMB
  endif
  cat1 = cat[gd]
  gaia1 = gaia[index[gd,0]]
  tmass1 = tmass[index[gd,1]]

  if n_elements(allcat) eq 0 then begin
    cat0 = cat[0]
    struct_assign,{dum:''},cat0
    allcat = replicate(cat0,1e7)
    gaia0 = gaia[0]
    struct_assign,{dum:''},gaia0
    allgaia = replicate(gaia0,1e7)
    tmass0 = tmass[0]
    struct_assign,{dum:''},tmass0
    alltmass = replicate(tmass0,1e7)
    cnt = 0LL
  endif
  tempcat = allcat[cnt:cnt+ngd-1]
  struct_assign,cat1,tempcat
  allcat[cnt:cnt+ngd-1] = tempcat
  tempgaia = allgaia[cnt:cnt+ngd-1]
  struct_assign,gaia1,tempgaia
  allgaia[cnt:cnt+ngd-1] = tempgaia
  temptmass = alltmass[cnt:cnt+ngd-1]
  struct_assign,tmass1,temptmass
  alltmass[cnt:cnt+ngd-1] = temptmass
  cnt += ngd

  ;stop
  BOMB:
endfor
; Trim extra elements
allcat = allcat[0:cnt-1]
allgaia = allgaia[0:cnt-1]
alltmass = alltmass[0:cnt-1]
; Maybe match to PS1 as well

; Save the matched catalogs
;save,allcat,allgaia,alltmass,file='combine_stripe82_vrband.dat'  

; Make the plot
!p.font = 0
setdisp
file = 'stripe82_vrband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
model_mag = allgaia.gmag
gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-allcat[gd].cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='VR-band'
bindata,jk0[gd],model_mag[gd]-allcat[gd].cmag,xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.2 and xbin le 0.6,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.2,0.2],[-2,2],linestyle=1,co=255
oplot,[0.6,0.6],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; LEAVE THE EQUATION AS IS!!  ADJUST COLOR RANGE TO 0.2<JK0<0.60

; This is the "corrected" relation!!!
;model_mag2 = allgaia.gmag

; versus EBV
file = 'stripe82_vrband_magdiff_ebv'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
gd2 = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0 and jk0 lt 0.6,ngd2)
hess,allcat[gd2].ebv,model_mag[gd2]-allcat[gd2].cmag,dx=0.01,dy=0.02,xr=[0,0.8],yr=[-1,1],/log,xtit='E(B-V)',ytit='Model-Mag',tit='VR-band'
oplot,[-1,3],[0,0],linestyle=2,co=255
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; don't see much EBV dependence

stop

end
