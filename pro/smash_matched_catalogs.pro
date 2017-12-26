pro smash_matched_catalogs

; Get calibrated fields with full ugriz
; and not in the LMC/SMC main bodies
rootdir = smashred_rootdir()
info = mrdfits(rootdir+'cp/red/photred/catalogs/final/v5/check_calibrated_v5.fits',1)
;restore,'/data/smash/cp/red/photred/catalogs/pro/check_calibrated.dat'
gd = where(info.calflag eq 2 and info.nuchips gt 0 and long(strmid(info.field,5)) gt 60,ngd)
fields = info[gd].field
nfields = n_elements(fields)

catdir = rootdir+'cp/red/photred/catalogs/final/v5/'
gaiadir = rootdir+'cp/red/photred/gaia/'
galexdir = rootdir+'cp/red/photred/galex/'
;fields = ['Field101','Field110','Field134']
;nfields = n_elements(fields)

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

;; Set up the fitting parameters
;;gaiastr = replicate({filter:'',bestcoef:dblarr(4),bestgirange:fltarr(2),decentcoef:dblarr(5),decentgirange:fltarr(2)},5)
;gaiastr = replicate({filter:'',colnames:strarr(2),norder:0,bestcoef:dblarr(4),bestcolrange:fltarr(2),decentcoef:dblarr(5),decentcolrange:fltarr(2),yr:fltarr(2)},5)
;gaiastr.filter = filters
;for i=0,nfilters-1 do begin
;  filt = filters[i]
;  ; Insert the color names and ranges
;  case filt of
;  'u': begin
;    ; Decent fit for 0.30 < g-i < 2.7
;    gaiastr[i].decentcolrange = [0.30,2.7]
;    ; Best fit for 0.60< g-i < 2.2
;    gaiastr[i].bestcolrange = [0.60,2.2]
;    yr = [-1.0,6.0]
;    colnames = ['g','i']
;    norder = 4
;  end
;  'g': begin
;    ; Becent fit for 0.0 < g-i < 2.8
;    ;gaiastr[i].decentcolrange = [0.0,2.8]
;    gaiastr[i].decentcolrange = [0.15,1.0]
;    ; Best fit for 0.50 < g-i < 1.7
;    ;gaiastr[i].bestcolrange = [0.5,1.7]
;    gaiastr[i].bestcolrange = [0.2,0.6]
;    yr = [-0.5,3.0]
;    ;colnames = ['g','i']
;    colnames = ['r','z']
;    norder = 3
;  end
;  'r': begin
;    ; Decent fit for 0.25 < g-i < 2.60
;    gaiastr[i].decentcolrange = [0.25,2.60]
;    ; Best fit for 0.50 < g-i < 2.0
;    gaiastr[i].bestcolrange = [0.5,2.0]
;    yr = [-1.0,1.5]
;    colnames = ['g','i']
;    norder = 3
;  end
;  'i': begin
;    ; Decent fit for 0.0 < g-i < 2.8
;    gaiastr[i].decentcolrange = [0.0,2.8]
;    ; Best fit for  0.40 < g-i < 1.9 
;    gaiastr[i].bestcolrange = [0.4,1.9]
;    yr = [-1.5,1.0]
;    colnames = ['g','i']
;    norder = 3
;  end
;  'z': begin
;    ; Decent fit for 0.15 < g-i < 2.9
;    gaiastr[i].decentcolrange = [0.15,2.9]
;    ; Best fit for 0.50 < g-i < 2.3
;    gaiastr[i].bestcolrange = [0.5,2.3]
;    yr = [-1.5,0.5]
;    colnames = ['g','i']
;    norder = 3
;  end
;  else: stop,' not an option'
;  endcase
;  gaiastr[i].colnames = colnames
;  gaiastr[i].norder = norder
;  gaiastr[i].yr = yr
;endfor


;restore,'/data/smash/cp/red/photred/gaia/gaiasmash_matchedcats.dat'
;goto,fitcoef


; Load the photometric fields and gaia data
undefine,mobj,mgaia
fieldstr = replicate({field:'',fieldid:0,nmatch:0L,ucoef:dblarr(5),urms:99.99,ubin:fltarr(31),gcoef:dblarr(4),grms:99.99,$
                      gbin:fltarr(31),rcoef:dblarr(4),rrms:99.99,rbin:fltarr(31),icoef:dblarr(4),irms:99.99,ibin:fltarr(31),$
                      zcoef:dblarr(5),zrms:99.99,zbin:fltarr(31)},nfields)
ftags = tag_names(fieldstr)
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfields,2),' ',fields[i]

  ;; Load the GALEX catalog
  galexfile = galexdir+fields[i]+'_galex.fits'
  if file_test(galexfile) eq 0 then begin
    print,galexfile,' NOT FOUND'
    goto,BOMB
  endif
  galex = mrdfits(galexfile,1)

  obj = mrdfits(catdir+fields[i]+'_combined_allobj.fits.gz',1)
  obj.id = strtrim(obj.id,2)
  ;gaia = mrdfits(gaiadir+fields[i]+'_gaia.fits',1)
  otags = tag_names(obj)

  ;; Load the XMATCH catalog
  xmatch = mrdfits(catdir+fields[i]+'_combined_allobj_xmatch.fits.gz',1)
  xmatch.id = strtrim(xmatch.id,2)

  ;; Only keep objects with UBAND data and GALEX, 2MASS, and GAIA matches!!!
  MATCH,obj.id,xmatch.id,ind1,ind2,/sort,count=nmatch
  if nmatch eq 0 then begin
    print,'No XMATCH matches'
    goto,BOMB
  endif
  obj = obj[ind1]
  xmatch = xmatch[ind2]

  ; Now crossmatch with GALEX
  srcmatch,obj.ra,obj.dec,galex.raj2000,galex.dej2000,0.5,ind1,ind2,/sph,count=nmatch
  print,strtrim(nmatch,2),' matches'
  if nmatch eq 0 then begin
    print,'No Galex matches'
    goto,BOMB
  endif
  obj = obj[ind1]
  xmatch = xmatch[ind1]
  galex = galex[ind2]

  ; ONLY KEEP OJECTS WITH GALEX, 2MASS AND GAIA MATCHES
  gdkeep = where(xmatch.tmass_match eq 1 and xmatch.gaia_match eq 1,ngdkeep)
  print,strtrim(ngdkeep,2),' final sources with SMASH, Galex, 2MASS and Gaia matches'
  if ngdkeep eq 0 then goto,BOMB
  obj = obj[gdkeep]
  xmatch = xmatch[gdkeep]
  galex = galex[gdkeep]

  ; Convert to common obj schema without indices
  ;  gaia too
  schema_obj = {id:'',fieldid:0,ra:0.d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,$
                i:0.0,ierr:0.0,z:0.0,zerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ebv:0.0,$
                gaia_source:0L,gaia_gmag:0.0,gaia_gerr:0.0,tmass_id:'',tmass_jmag:0.0,tmass_jerr:0.0,$
                tmass_hmag:0.0,tmass_herr:0.0,tmass_kmag:0.0,tmass_kerr:0.0,tmass_qflg:'',$
                galex_id:0L,galex_fuv:0.0,galex_fuverr:0.0,galex_nuv:0.0,galex_nuverr:0.0}
  newobj = replicate(schema_obj,ngdkeep)
  struct_assign,obj,newobj
  newobj.fieldid = long(strmid(fields[i],5))
  struct_assign,xmatch,newobj,/nozero
  newobj.galex_id = galex.objid
  newobj.galex_fuv = galex.fuv
  newobj.galex_fuverr = galex.e_fuv
  newobj.galex_nuv = galex.nuv
  newobj.galex_nuverr = galex.e_nuv

  fieldstr[i].field = fields[i]
  fieldstr[i].fieldid = long(strmid(fields[i],5))
  fieldstr[i].nmatch = ngdkeep
  ;for j=0,nfilters-1 do begin
  ;  magind = where(otags eq strupcase(gaiastr[j].filter),nmagind)
  ;  colmagind1 = where(otags eq strupcase(gaiastr[j].colnames[0]),ncolmagind1)
  ;  colmagind2 = where(otags eq strupcase(gaiastr[j].colnames[1]),ncolmagind2)
  ;  coef = robust_poly_fitq(obj1.(colmagind1)-obj1.(colmagind2),obj1.(magind)-gaia1._gmag_,gaiastr[j].norder)
  ;  rms = mad(obj1.(magind)-gaia1._gmag_-poly(obj1.(colmagind1)-obj1.(colmagind2),coef))
  ;  bindata,obj1.(colmagind1)-obj1.(colmagind2),obj1.(magind)-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  ;  coefind = where(ftags eq strupcase(gaiastr[j].filter)+'COEF',ncoefind)
  ;  rmsind = where(ftags eq strupcase(gaiastr[j].filter)+'RMS',nrmsind)
  ;  binind = where(ftags eq strupcase(gaiastr[j].filter)+'BIN',nbinind)
  ;  fieldstr[i].(coefind) = coef
  ;  fieldstr[i].(rmsind) = rms
  ;  fieldstr[i].(binind) = ybin
  ;endfor

  push,mobj,newobj

  BOMB:

  ;stop
endfor

stop

;; RMS in the binned values
;urms = mad(fieldstr.ubin,dim=2)  ; ~5%
;grms = mad(fieldstr.gbin,dim=2)  ; ~2-5%
;rrms = mad(fieldstr.rbin,dim=2)  ; ~1%
;irms = mad(fieldstr.ibin,dim=2)  ; ~1%
;zrms = mad(fieldstr.zbin,dim=2)  ; ~1.5%

;save,fieldstr,mobj,file='/datalab/users/dnidever/nsc/smash_matched_catalog.dat'
stop

; Do the fitting of the color-color relations
FITCOEF:

; Field132 shows higher scatter and multiple "sequences"
;  gaia or smash issue?
; Field171 as well
setdisp
!p.font = 0
nobj = n_elements(mobj)
otags = tag_names(mobj)
ftags = tag_names(fieldstr)
;colr = mobj.g-mobj.i
x = scale_vector(findgen(100),0.0,3.0)
nfields = n_elements(fieldstr)


; Filter loop
for i=0,4 do begin
  filt = filters[i]

  if keyword_set(save) then begin
    file = '/data/smash/cp/red/photred/gaia/plots/gaiasmash_colorterms_'+filt
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=9.5
  endif

  colnames = gaiastr[i].colnames
  yr = gaiastr[i].yr
  colmagind1 = where(otags eq strupcase(colnames[0]),ncolmagind1)
  colmagind2 = where(otags eq strupcase(colnames[1]),ncolmagind2)
  colr = mobj.(colmagind1) - mobj.(colmagind2)

  ; Get good fields and their stars
  rmsind = where(ftags eq strupcase(filt)+'RMS',nrmsind)
  gdfields = where(fieldstr.(rmsind) lt 0.5,ngdfields)
  goodmask = intarr(nobj)
  for j=0,ngdfields-1 do begin
    MATCH,mobj.fieldid,fieldstr[gdfields[j]].fieldid,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then goodmask[ind1]=1
  endfor
  magind = where(otags eq strupcase(filt),nmagind)
  ; Get good stars to use
  ;gdobj = where(goodmask eq 1 and mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21 and mobj.(magind) lt 21,ngdobj)
  gdobj = where(goodmask eq 1 and mobj.(colmagind1) lt 21 and mobj.(colmagind2) lt 21 and mgaia._gmag_ lt 21 and mobj.(magind) lt 21 and $
                abs(mobj.sharp) lt 1 and mobj.chi lt 4,ngdobj)

  ; Hess diagram
  ;hess,mobj[gdobj].g-mobj[gdobj].i,mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,dx=0.01,dy=0.01,xr=[0,3.2],yr=yr,$
  ;     xtit='g-i',ytit=filt+'-G',tit=filt+'-G',/log
  hess,colr[gdobj],mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,dx=0.01,dy=0.01,xr=[0,3.2],yr=yr,$
       xtit=colnames[0]+'-'+colnames[1],ytit=filt+'-G',tit=filt+'-G',/log

  ; Overplot the individual fits
  coefind = where(ftags eq strupcase(filt)+'COEF',ncoefind)
  ;for j=0,ngdfields-1 do oplot,x,poly(x,fieldstr[gdfields[j]].(coefind)),co=255
  ; Get values for individual fits and calculate RMS
  fit = fltarr(100,ngdfields)
  for j=0,ngdfields-1 do fit[*,j]=poly(x,fieldstr[gdfields[j]].(coefind))
  ;plot,x,mad(fit,dim=2),ps=-1

  ; Poly fit to all stars
  coef_all = robust_poly_fitq(colr[gdobj],mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,3)
  ;oplot,x,poly(x,coef_all),co=250

  ; Bin the data and poly fitting
  bindata,colr[gdobj],mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.0,max=3.0,gdind=gdind
  ;oplot,xbin,ybin,ps=-1,sym=3,co=80
  coef_bin_all = robust_poly_fitq(xbin[gdind],ybin[gdind],4)
  ;oplot,x,poly(x,coef_bin_all),co=150

  ; Get decent fit
  decentind = where(xbin ge gaiastr[i].decentcolrange[0] and xbin le gaiastr[i].decentcolrange[1] and $
                    finite(ybin) eq 1,ndecentind)
  coef_bin_decent = robust_poly_fitq(xbin[decentind],ybin[decentind],4)
  gaiastr[i].decentcoef = coef_bin_decent
  oplot,xbin[decentind],poly(xbin[decentind],coef_bin_decent),co=255

  ; Get best fit
  bestind = where(xbin ge gaiastr[i].bestcolrange[0] and xbin le gaiastr[i].bestcolrange[1] and $
                  finite(ybin) eq 1,nbestind)
  coef_bin_best = robust_poly_fitq(xbin[bestind],ybin[bestind],3)
  gaiastr[i].bestcoef[0:3] = coef_bin_best
  oplot,xbin[bestind],poly(xbin[bestind],coef_bin_best),co=250,linestyle=2

  ; Legend
  ;al_legend,['Binned','Decent','Best'],textcolor=[80,200,250],psym=[1,0,0],charsize=1.5,/top,/left
  al_legend,['Decent','Best'],textcolor=[255,250],charsize=1.5,/top,/left

  if keyword_set(save) then begin
    ps_close
    ps2png,file+'.eps',/eps
  endif

  ; Now test it on the binned values per field
  xbin = findgen(31)*0.1+0.05  ; midpoint
  binind = where(tag_names(fieldstr) eq strupcase(filt)+'BIN',nbinind)
  medarr = fltarr(nfields)
  binarr = fltarr(31,nfields)
  for j=0,nfields-1 do begin
    ymodel = poly(xbin,gaiastr[i].bestcoef)   ; model for all points
    ybin = fieldstr[j].(binind)
    gdind = where(xbin ge gaiastr[i].bestcolrange[0] and xbin le gaiastr[i].bestcolrange[1] and $
                  finite(ybin) eq 1,ngdind)

    med = median(ybin[gdind]-ymodel[gdind])
    binarr[*,j] = ybin-med
    medarr[j] = med
  endfor
  rms = mad(binarr,dim=2)
  orms = mad(fieldstr.(binind),dim=2)
  plot,xbin,rms,ps=-1,yr=[0,0.1]
  oplot,xbin,orms,ps=-1,co=250
  ; u-band, best RMS point is now 0.013 when before it was 0.033
  ;   but overall the RMS is not great
  ; g-band, best RMS is now 0.02 while before it was 0.02-0.04
  ; r-band, best RMS is now 0.002 while before it was 0.01
  ; i-band, best RMS is now 0.002 while before it was 0.01
  ; z-band, best RMS is now 0.004 while before it was 0.015

  ;stop

endfor

; MWRFITS,gaiastr,gaiadir+'gaiasmash_colorterms.fits',/create

stop

end
