pro combine_stripe82_zband

; Combine together the exposure catalogs in Stripe82 for a single band

filter = 'z'

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+"users/dnidever/nsc/"
str = mrdfits(dir+'instcal/nsc_instcal_calibrate.fits',1)
str.expdir = strtrim(str.expdir,2)
str.filter = strtrim(str.filter,2)
str.expnum = strtrim(str.expnum,2)

outfile = str.expdir+'/'+file_basename(str.expdir)+'_cat.fits'
ind0 = where((str.ra lt 61 or str.ra gt 299) and abs(str.dec) lt 3.0 and str.zpterm ne 0 and str.filter eq filter,nind0)
medzpterm = median(str[ind0].zpterm)
sigzpterm = mad(str[ind0].zpterm)
ind = where((str.ra lt 61 or str.ra gt 299) and abs(str.dec) lt 3.0 and str.zpterm ne 0 and str.filter eq filter and $
            abs(str.zpterm-medzpterm) lt 3*sigzpterm,nind)

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

  ; Load the 2MASS file
  tmassfile = str[ind[i]].expdir+'/'+file_basename(str[ind[i]].expdir)+'_TMASS.fits'
  if file_test(tmassfile) eq 0 then goto,BOMB
  tmass = MRDFITS(tmassfile,1,/silent)

  ; Matching
  index = lonarr(ncat)-1
  dcr = 1.0
  SRCMATCH,tmass.raj2000,tmass.dej2000,cat.ra,cat.dec,dcr,tind1,tind2,/sph,count=ntmatch
  if ntmatch gt 0 then index[tind2] = tind1
  gd = where(index gt -1,ngd)
  print,'  ',strtrim(ngd,2),' matches to 2MASS'
  if ngd eq 0 then begin
    print,'No matches to 2MASS'
    goto,BOMB
  endif
  cat1 = cat[gd]
  tmass1 = tmass[index[gd]]

  if n_elements(allcat) eq 0 then begin
    cat0 = cat[0]
    struct_assign,{dum:''},cat0
    allcat = replicate(cat0,5e6)
    tmass0 = tmass[0]
    struct_assign,{dum:''},tmass0
    alltmass = replicate(tmass0,5e6)
    cnt = 0LL
  endif
  tempcat = allcat[cnt:cnt+ngd-1]
  struct_assign,cat1,tempcat
  allcat[cnt:cnt+ngd-1] = tempcat
  temptmass = alltmass[cnt:cnt+ngd-1]
  struct_assign,tmass1,temptmass
  alltmass[cnt:cnt+ngd-1] = temptmass
  cnt += ngd

  ;stop
  BOMB:
endfor
; Trim extra elements
allcat = allcat[0:cnt-1]
alltmass = alltmass[0:cnt-1]
; Maybe match to PS1 as well

; Make the plot
!p.font = 0
setdisp
file = 'stripe82_zband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
model_mag = alltmass.jmag + 0.765720*jk0 + 0.40*allcat.ebv +  0.605658
hess,jk0,model_mag-allcat.cmag,dx=0.02,dy=0.02,xr=[-0.1,1.3],yr=[-1,1],/log,xtit='(J-Ks)o',ytit='Model-Mag',tit='z-band'
bindata,jk0,model_mag-allcat.cmag,xbin,ybin,binsize=0.05,/med,gdind=gdind,min=0,max=1.2
oplot,xbin[gdind],ybin[gdind],ps=-1,co=255
gd = where(xbin ge 0.2 and xbin le 0.8,ngd)
coef = robust_poly_fitq(xbin[gd],ybin[gd],1)
xx = scale_vector(findgen(100),-1,3)
oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.3,0.3],[-2,2],linestyle=1,co=255
oplot,[0.7,0.7],[-2,2],linestyle=1,co=255
al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!dn0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; This is the "corrected" relation!!!
;model_mag = alltmass.jmag + 0.765720*jk0 + 0.40*allcat.ebv +  0.605658

; versus EBV
file = 'stripe82_zband_magdiff_ebv'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
hess,allcat.ebv,model_mag-allcat.cmag,dx=0.02,dy=0.02,xr=[0,0.8],yr=[-1,1],/log,xtit='E(B-V)',ytit='Model-Mag',tit='z-band'
;bindata,jk0,model_mag-allcat.cmag,xbin,ybin,binsize=0.05,/med,gdind=gdind,min=0,max=1.2
;oplot,xbin[gdind],ybin[gdind],ps=-1,co=255
;gd = where(xbin ge 0.2 and xbin le 0.8,ngd)
;coef = robust_poly_fitq(xbin[gd],ybin[gd],1)
;xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!dn0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

stop

end
