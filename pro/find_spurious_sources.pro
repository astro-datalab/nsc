pro find_spurious_sources,expname

; Find spurious sources in a single exposure catalog

; Criteria for spurious near-bright-star sources
; -within a certain distance of bright, saturated star
; -internal SE flags set (1 or 2)  (only good in non-crowded regions, in crowded regions many have flags=3)
; -delta mag > some number
; -low S/N???
; -many elliptical (but not all)

; NEEDS TO WORK IN BOTH CROWDED AND UNCROWDED REGIONS!!!

; a large majority of sources in crowded fields have flags=3, while
; most flags=3 sources in uncrowded fields are spurious.
; HOWEVER, essentially all spurious sources around bright stars
; (whether in crowded or uncrowded fields) have flags=3
; So it's still a useful criteria (mainly for uncrowded fields)


; To find the bright, saturated stars:
; -find all Gaia/2MASS sources brighter than max(cat.cmag)+0.5 or so
; -match all of them with our catalog
; -the ones that DO NOT match should be saturated

rootdir = '/Users/nidever/datalab/nsc/qa/'
expname = 'c4d_140108_022046_ooi_z_v1'

; Load the catalog
cat = mrdfits(rootdir+expname+'_cat.fits.gz',1)

; Get the reference catalogs
tmass = mrdfits(rootdir+expname+'_TMASS.fits.gz',1)
gaia = mrdfits(rootdir+expname+'_GAIA.fits.gz',1)

; Only keep sources in the observed region
cenra = mean([min(cat.ra),max(cat.ra)])
if range(cat.ra) gt 100 then begin
  ra = cat.ra
  bd = where(ra gt 180,nbd)
  if nbd gt 0 then ra[bd]-=360
  cenra = mean([min(ra),max(ra)])
  if cenra lt 0 then cenra+=360
endif
cendec = mean([min(cat.dec),max(cat.dec)])
rotsphcen,cat.ra,cat.dec,cenra,cendec,lon,lat,/gnomic
rotsphcen,tmass.raj2000,tmass.dej2000,cenra,cendec,tlon,tlat,/gnomic
rotsphcen,gaia.ra_icrs,gaia.de_icrs,cenra,cendec,glon,glat,/gnomic
rlon = [min(lon)-0.05,max(lon)+0.05]
rlat = [min(lat)-0.05,max(lat)+0.05]
gdtmass = where(tlon ge rlon[0] and tlon le rlon[1] and $
                tlat ge rlat[0] and tlat le rlat[1],ngdtmass)
tmass = tmass[gdtmass]
gdgaia = where(glon ge rlon[0] and glon le rlon[1] and $
               glat ge rlat[0] and glat le rlat[1],ngdgaia)
gaia = gaia[gdgaia]

; Find bright, saturated sources
;-------------------------------
; Use brightest observed magnitudes to figure out saturation point
minmag = min(cat.cmag)
; Get all reference sources brighter than this
brtgaiaind = where(gaia.gmag le (minmag+1.0),nbrtgaia)
brtgaia = gaia[brtgaiaind]
print,strtrim(nbrtgaia,2),' bright Gaia sources'

; Crossmatch and remove any that were detected (couldn't be saturated)
srcmatch,brtgaia.ra_icrs,brtgaia.de_icrs,cat.ra,cat.dec,1.0,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' bright sources observed. ',strtrim(nbrtgaia-nmatch,2),' left'
if nmatch eq nbrtgaia then begin
  print,'No bright Gaia sources left'
  return 
endif else begin
  remove,ind1,brtgaia
  nbrtgaia = n_elements(brtgaia) 
endelse

; Use MATCHALL_SPH to find the closest matches
dcr = 5.0       ; arcsec
res = matchall_sph(brtgaia.ra_icrs,brtgaia.de_icrs,cat.ra,cat.dec,dcr/3600,nmatch,distance=distance)
; res gives reverse indices

;plot,cat.ra,cat.dec,ps=3
;oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=8,co=250

; Loop through bright sources and find close neighbors
For i=0,nbrtgaia-1 do begin
  brt1 = brtgaia[i]
   
  ;plot,cat.ra,cat.dec,ps=3
  ;oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=8,co=250
  ;oplot,[brtgaia[i].ra_icrs],[brtgaia[i].de_icrs],ps=4,co=150,sym=4
  
  if res[i] ne res[i+1] then begin
    ind = res[res[i]:res[i+1]-1]
    nind = n_elements(ind)
    print,'Gmag=',stringize(brt1.gmag,ndec=2)
    print,'nmatches = ',strtrim(nind,2)
    writecol,-1,cat[ind].cmag,cat[ind].flags,1/cat[ind].cerr,cat[ind].ellipticity

    plot,cat.ra,cat.dec,ps=3,xr=[-0.005,0.005]/cos(cendec/!radeg)+brtgaia[i].ra_icrs,yr=[-0.005,0.005]+brtgaia[i].de_icrs,xs=1,ys=1
    oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=1,co=250,sym=2
    oplot,[cat[ind].ra],[cat[ind].dec],ps=4,co=150
    gd = where(cat.ra ge brtgaia[i].ra_icrs-0.005/cos(cendec/!radeg) and cat.ra le brtgaia[i].ra_icrs+0.005/cos(cendec/!radeg) and $
               cat.dec ge brtgaia[i].de_icrs-0.005 and cat.dec le brtgaia[i].de_icrs+0.005,ngd)
    for j=0,ngd-1 do begin
      cat1 = cat[gd[j]]
      coords = ellipsecoords(cat1.a_world,cat1.b_world,cat1.ra,cat1.dec,cat1.theta_world)
      oplot,coords[0,*],coords[1,*]
    endfor
    
    stop
  endif else print,'no matches'

Endfor


stop

end
