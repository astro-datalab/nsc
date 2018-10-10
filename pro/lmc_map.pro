pro lmc_map

;;  Make density map of LMC periphery with NSC data

dir = '/dl1/users/dnidever/nsc/instcal/v2/'
outdir = '/dl1/users/dnidever/nsc/lmc/'

str = mrdfits(dir+'lists/nsc_instcal_combine.fits',1)
str.file = strtrim(str.file,2)
cel2lmc,str.ra,str.dec,lmcpa,lmcrad
cel2smc,str.ra,str.dec,smcpa,smcrad
glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
gd = where(lmcrad gt 6 and lmcrad lt 30 and smcrad gt 2.0 and abs(glat) gt 15 and str.exptimefilt[1] gt 0 and $
           str.exptimefilt[2] gt 0,ngd)
;gd = where(rad gt 6 and rad lt 20 and abs(glat) gt 15 and str.exptimefilt[1] gt 0 and $
;           str.exptimefilt[2] gt 0,ngd)
;           (str.exptimefilt[2] gt 0 or str.exptimefilt[3] gt 0.0),ngd)
print,strtrim(ngd,2),' HEALPix pixels around LMC'

pix = str[gd].pix
outfiles = outdir+strtrim(pix,2)+'_lmc.fits'
test = file_test(outfiles)
bd = where(test eq 0,nbd)

cmd = 'get_lmc_data,'+strtrim(pix[bd],2)+',/redo'
cmddir = '/data0/dnidever/nsc/lmc/'+strarr(nbd)

stop

;pbs_daemon,cmd,cmddir,jobs=jobs,nmulti=20,/hyper,/idle,prefix='nsclmc',wait=0.1

;print,'Loading the data'
;all = mrdfits(outdir+'lmc_map.fits',1)
;goto,plotting

; Load the files
cnt = 0LL
for i=0,ngd-1 do begin
  file = outdir+strtrim(str[gd[i]].pix,2)+'_lmc.fits'
  if file_test(file) eq 1 then begin
    cat = mrdfits(file,1,/silent)
    ncat = n_elements(cat)
     print,strtrim(i+1,2),' ',file,' ',strtrim(ncat,2)

    if n_elements(all) eq 0 then begin
      schema = cat[0]
      struct_assign,{dum:''},schema
      all = replicate(schema,30e6)
      nall = n_elements(all)
    endif

    if cnt+ncat gt nall then begin
      print,'Adding new elements'
      old = all
      all = replicate(schema,nall+5e6)
      all[0:nall-1] = old
      nall = n_elements(all)
      undefine,old 
    endif

    ; add new catalog
    all[cnt:cnt+ncat-1] = cat
    cnt += ncat

  endif else print,file,' NOT FOUND'

endfor
; trim extra elements
all = all[0:cnt-1]

;mwrfits,all,outdir+'lmc_map.fits',/create
;mwrfits,all,outdir+'lmc_map.v1.fits',/create

stop

plotting:
setdisp

glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
add_tag,str,'mlon',0.0,str
add_tag,str,'mlat',0.0,str
str.mlon = mlon
str.mlat = mlat
add_tag,str,'ebv',0.0,str
str.ebv = dust_getval(glon,glat,/noloop,/interp)

nside = 128
radeg=180.0d0 / !dpi
theta=(90-all.dec)/radeg    
phi=all.ra/radeg
ang2pix_ring,nside,theta,phi,ipring


;gd = where(rad gt 6 and rad lt 20 and abs(glat) gt 15 and str.exptimefilt[1] gt 0 and $
;           str.exptimefilt[2] gt 0 and str.depth95filt[1] gt 22.8 and str.ebv lt 0.5,ngd)
gd = where(lmcrad gt 6 and lmcrad lt 30 and smcrad gt 3.0 and abs(glat) gt 15 and str.exptimefilt[1] gt 0 and $
           str.exptimefilt[2] gt 0 and str.depth95filt[1] gt 22.8 and str.ebv lt 0.5,ngd)

; deredden
mag = all.gmag
col = all.gmag-all.rmag
ag = all.ebv*3.303
ar = all.ebv*2.285
mag0 = all.gmag-ag
col0 = (all.gmag-ag)-(all.rmag-ar)

dens = fltarr(196608)
pix = lindgen(196608)
xcut = [ 0.0038121236,    0.0038121236,      0.41930888,      0.30972732]
;xcut = [ -0.0561879,   -0.0561879,     0.359309,     0.249727]
ycut = [ 21.8, 22.8, 22.8, 21.8]
;roi_cut,xcut,ycut,col0,mag0,ind,cutind,fac=100
roi_cut,xcut,ycut,all.gmag-all.rmag,all.gmag,ind,cutind,fac=100
dens[ipring[cutind]]++
add_tag,str,'dens',0.0,str
match,str.pix,pix,ind1,ind2,/sort
str[ind1].dens = dens[ind2]

; get depth for each star
pixdepth = fltarr(196608)
pixdepth[str[gd].pix] = str[gd].depth95filt[1]
depth = pixdepth[ipring]
; mask
pixmask = intarr(196608)
pixmask[str[gd].pix] = 1
mask = pixmask[ipring]

; LMC MSTO luminosity function
deepind=where(all.mlon ge 8.714 and all.mlon lt 9.9 and all.mlat gt 1.73 and all.mlat lt 2.83)
all2=all[deepind]           
gg=where(all2.gmag-all2.rmag gt 0.11 and all2.gmag-all2.rmag lt 0.433)
hist=histogram(all2[gg].gmag,bin=0.2,locations=xhist)
;plot,xhist,hist
;coef=robust_poly_fit(xhist[38:48],hist[38:48],3)
coef = [ -1.90342e+06,      246319.,     -10635.1,      153.300]

; correct to depth of 23.0
lum23 = poly(23.0,coef)
lumfrac = lum23/poly(str.depth95filt[1],coef)

;before scaling for the variable depth I need to SUBTRACT the background for that depth
;(mainly MW halo stars).  The MW halo is pretty flat with magnitude, so the background
;should just scale with magnitude range used.

plotc,str[gd].mlon,str[gd].mlat,str[gd].dens*lumfrac[gd],ps=8,/xflip,/log,sym=2,min=15,/trim
ml0 = 0.21550614
mb0 = 2.3321053
pa = scale_vector(findgen(1000),0,2*!dpi)
oplot,16*sin(pa)+ml0,16*cos(pa)+mb0
oplot,14*sin(pa)+ml0,14*cos(pa)+mb0
oplot,12*sin(pa)+ml0,12*cos(pa)+mb0
oplot,10*sin(pa)+ml0,10*cos(pa)+mb0


; 0.1x0.1 deg maps
hess,all[cutind].mlon,all[cutind].mlat,dum,im1,dx=0.1,dy=0.1,xr=[-25,20],yr=[-20,23],xarr=xarr1,yarr=yarr1,/noplot
hess,all.mlon,all.mlat,all.ebv,ebv1,dx=0.1,dy=0.1,/mean,xr=[-25,20],yr=[-20,23],/noplot
hess,all.mlon,all.mlat,depth,depth1,dx=0.1,dy=0.1,/mean,xr=[-25,20],yr=[-20,23],/noplot
hess,all.mlon,all.mlat,mask,mask1,dx=0.1,dy=0.1,/tot,xr=[-25,20],yr=[-20,23],/noplot
; 0.2x0.2 deg maps
hess,all[cutind].mlon,all[cutind].mlat,dum,im2,dx=0.2,dy=0.2,xr=[-25,20],yr=[-20,23],xarr=xarr2,yarr=yarr2,/noplot
hess,all.mlon,all.mlat,all.ebv,ebv2,dx=0.2,dy=0.2,/mean,xr=[-25,20],yr=[-20,23],/noplot
hess,all.mlon,all.mlat,depth,depth2,dx=0.2,dy=0.2,/mean,xr=[-25,20],yr=[-20,23],/noplot
hess,all.mlon,all.mlat,mask,mask2,dx=0.2,dy=0.2,/tot,xr=[-25,20],yr=[-20,23],/noplot
;save,im1,xarr1,yarr1,ebv1,depth1,mask1,im2,xarr2,yarr2,ebv2,depth2,mask2,file=outdir+'lmc_map.dat' 


stop

end
