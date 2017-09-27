pro find_spurious_sources,expdir

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

;expname = 'c4d_140108_022046_ooi_z_v1'
expdir = '/dl1/users/dnidever/nsc/instcal/c4d/20140108/c4d_140108_022046_ooi_z_v1/'

NSC_ROOTDIRS,dldir,mssdir,localdir

; Not enough inputs
if n_elements(expdir) eq 0 then begin
  print,'Syntax - find_spurious_sources,expdir'
  return
endif

; Make sure the directory exists
if file_test(expdir,/directory) eq 0 then begin
  print,expdir,' NOT FOUND'
  return
endif

t00 = systime(1)

; Setting pool thread values
if n_elements(ncpu) eq 0 then ncpu=1
CPU, TPOOL_NTHREADS = ncpu

base = file_basename(expdir)
logf = -1
;outfile = expdir+'/'+base+'_cat.fits'

printlog,logf,'Find spurious sources in exposure ',base,' in ',expdir

; What instrument is this?
instrument = 'c4d'  ; by default
if stregex(expdir,'/k4m/',/boolean) eq 1 then instrument='k4m'
if stregex(expdir,'/ksb/',/boolean) eq 1 then instrument='ksb'
printlog,logf,'This is a '+instrument+' exposure'


;rootdir = '/Users/nidever/datalab/nsc/qa/'
;expname = 'c4d_140108_022046_ooi_z_v1'

; Load the catalog
cat = mrdfits(expdir+'/'+base+'_cat.fits',1)
ncat = n_elements(cat)

;; SAME EXACT CUTS AS IN NSC_INSTCAL_COMBINE.PRO!!!!
;; Make a cut on quality mask flag (IMAFLAGS_ISO)
;;   Don't use "difference image masking" for pre-V3.5 or so version
;;   because it had problems.  We don't have PLVER but just use
;;   the maximum IMAFLAGS_ISO value
;if max(cat.imaflags_iso) gt 10 then begin
;  bdcat = where(cat.imaflags_iso gt 0 and cat.imaflags_iso lt 120,nbdcat)
;;THIS IS A BITMASK, make sure none of the lower values are set
;; bdcat = where(cat.imaflags_iso gt 0 and cat.imagflags_iso ne 128,nbdcat)
;endif else begin
;  bdcat = where(cat.imaflags_iso gt 0,nbdcat)
;endelse
;if nbdcat gt 0 then begin
;  print,'  Removing ',strtrim(nbdcat,2),' sources contaminated by bad pixels.'
;  if nbdcat eq ncat then return
;  REMOVE,bdcat,cat
;  ncat = n_elements(cat)
;endif


; --CP bit masks, Pre-V3.5.0 (PLVER)
; Bit   DQ Type  PROCTYPE
; 1  detector bad pixel          InstCal
; 1  detector bad pixel/no data  Resampled
; 1  No data                     Stacked
; 2  saturated                   InstCal/Resampled
; 4  interpolated                InstCal/Resampled
; 16  single exposure cosmic ray InstCal/Resampled
; 64  bleed trail                InstCal/Resampled
; 128  multi-exposure transient  InstCal/Resampled  I turned this OFF
; --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks
;  1 = bad (in static bad pixel mask)
;  2 = no value (for stacks)
;  3 = saturated
;  4 = bleed mask
;  5 = cosmic ray
;  6 = low weight
;  7 = diff detect  I turned this OFF

; I THINK WE CAN GET MOST OF THE SATURATED SOURCES FROM THE CATALOG
; MAYBE COMBINE WITH BRIGHT GAIA SOURCES

b = where(((cat.imaflags_iso and 3) eq 3) or ((cat.imaflags_iso and 4) eq 4))

cendec = mean(minmax(cat.dec))
decrange = range(cat.dec)
cenra = mean(minmax(cat.ra))
rarange = range(cat.ra)*cos(cendec/!radeg)
; Wrapping around RA=0
if range(minmax(cat.ra)) gt 100 then begin
 ra = cat.ra
 bdra = where(ra gt 180,nbdra)
 if nbdra gt 0 then ra[bdra]-=360
 cenra = mean(minmax(ra))
 if cenra lt 0 then cenra+=360
 rarange = range(ra)*cos(cendec/!radeg)
 rawrap = 1
endif else rawrap=0                                
; Search radius                                                                                                                                                                    
radius = 1.1 * sqrt( (0.5*rarange)^2 + (0.5*decrange)^2 )

; Get the reference catalogs
tmass = getrefcat(cenra,cendec,radius,'2MASS')
gaia = getrefcat(cenra,cendec,radius,'Gaia')
;tmass = mrdfits(rootdir+expname+'_TMASS.fits.gz',1)
;gaia = mrdfits(rootdir+expname+'_GAIA.fits.gz',1)

; Match the two catalogs
srcmatch,gaia.ra_icrs,gaia.de_icrs,tmass.raj2000,tmass.dej2000,1.0,ind1,ind2,/sph,count=nmatch
add_tag,gaia,'jmag',99.99,gaia
add_tag,gaia,'kmag',99.99,gaia
if nmatch gt 0 then begin
  gaia[ind1].jmag = tmass[ind2].jmag
  gaia[ind1].kmag = tmass[ind2].kmag
endif

; Only keep sources in the observed region
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
gstar = where(cat.class_star ge 0.7,ngstar)
if ngstar eq 0 then gstar=where(cat.class_star ge 0.5,ngstar)
minmag = min(cat[gstar].cmag)
; Get all reference sources brighter than this
brtgaiaind = where(gaia.gmag le (minmag+1.0),nbrtgaia)
brtgaia = gaia[brtgaiaind]
print,strtrim(nbrtgaia,2),' bright Gaia sources'
; THIS IS A VERY **BAD** CRITERIA FOR SATURATED STARS
; THERE CAN BE SATURATED STARS ~4 MAGS FAINTER THAN THIS!!!

; Add reddening
add_tag,gaia,'ebv',0.0,gaia
glactc,gaia.ra_icrs,gaia.de_icrs,2000.0,glon,glat,1,/deg
gaia.ebv = dust_getval(glon,glat,/noloop,/interp)

; SOMETIMES THE BRIGHT SATURATED STARS **ARE** DETECTED!!!
; Some sources have bad CP flags.  Why are those kept??
; Not all spurious sources have FLAGS=2 or 3

;; Crossmatch and remove any that were detected (couldn't be saturated)
;srcmatch,brtgaia.ra_icrs,brtgaia.de_icrs,cat.ra,cat.dec,1.0,ind1,ind2,/sph,count=nmatch
;print,strtrim(nmatch,2),' bright sources observed. ',strtrim(nbrtgaia-nmatch,2),' left'
;if nmatch gt 0 then begin
;  if nmatch eq nbrtgaia then begin
;    print,'No bright Gaia sources left'
;    return 
;  endif else begin
;    remove,ind1,brtgaia
;    nbrtgaia = n_elements(brtgaia) 
;  endelse
;endif

; Copy mask file to local directory
maskfile = '/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
tempmaskfile = localdir+'dnidever/nsc/instcal/'+file_basename(maskfile)
file_copy,maskfile,tempmaskfile,/over
spawn,['funpack','-D',tempmaskfile],out,errout,/noshell
tempmaskfile = strmid(tempmaskfile,0,strlen(tempmaskfile)-3)

head0 = headfits(tempmaskfile,exten=0)
filter = sxpar(head0,'filter')
filt = strmid(filter,0,1)

; Loop over all chips
for i=0,62 do begin
  fits_read,tempmaskfile,mask,mhead,exten=i+1,/no_abort,message=message
  if message ne '' then goto,BOMB
  nx = sxpar(mhead,'NAXIS1')
  ny = sxpar(mhead,'NAXIS2')
  satmask = fix( ((mask and 3) eq 3) )
  satmask3 = convol(float(satmask),[[1,1,1],[1,1,1],[1,1,1]])

  imfile = '/dl1/users/dnidever/nsc/instcal/c4d/20140108/c4d_140108_022046_ooi_z_v1/c4d_140108_022046_ooi_z_v1.fits'
  fits_read,imfile,im,head,exten=i+1

  ; Get Gaia sources in this this
  head_xyad,mhead,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/deg
  rra = minmax(vra)
  rdec = minmax(vdec)
  gd1 = where(gaia.ra_icrs ge rra[0]-0.001 and gaia.ra_icrs le rra[1]+0.001 and $
             gaia.de_icrs ge rdec[0]-0.001 and gaia.de_icrs le rdec[1]+0.001,ngd1)
  gaia1 = gaia[gd1]
  head_adxy,mhead,gaia1.ra_icrs,gaia1.de_icrs,gx1,gy1,/deg
  gd2 = where(gx1 ge 0 and gx1 le (nx-1) and gy1 ge 0 and gy1 le (ny-1),ngd2)
  gaia2 = gaia1[gd2]
  gx2 = 0 > round(gx1[gd2]) < (nx-1)
  gy2 = 0 > round(gy1[gd2]) < (ny-1)
  satval = satmask[gx2,gy2]
  satval3 = satmask3[gx2,gy2]

  ; require that to be 3 neighboring pixels need to saturated as well
  gbad = where(satval eq 1 and satval3 ge 4,ngbad)
  if ngbad gt 0 then begin
    badgaia1 = gaia2[gbad]
    push,badgaia,badgaia1
  endif  

  ; THIS IS PICKING UP SOME FAINT STARS THAT HAPPEN TO LAND ON A CR AND
  ; AND OTHER SATURATED PIXELS, bleed trails/columns

  ;stop
  BOMB:

endfor

; Get model magnitudes
model_mag = getmodelmag(gaia,filt)

srcmatch,cat.ra,cat.dec,gaia.ra_icrs,gaia.de_icrs,1.0,ind1,ind2,/sph
cat2=cat[ind1]
gaia2=gaia[ind2]
model_mag2=model_mag[ind2]
plotc,cat2.cmag,cat2.cmag-model_mag2,cat2.ellipticity,ps=3,xr=[10,20],yr=[-5,5]
; a bunch of bright source have residuals that are way off,
; are these bright galaxies?
; Almost all of those are detected and in the catalog.  maybe
; that's why the residuals are off, the measured values are BAD
; they have BLENDED flags, but not many have saturated ones.
; These are real sources but close to the broad wings of very bright
; stars (high background) and I think that flux is getting
; incorporated and giving MUCH brighter mags.
; good example is source near X=1074, Y=2910 in ccdnum=1
; c4d_140108_022046_ooi_z_v1
; BACKGROUND isn't much help here.
; MAG_APER-MAG_AUTO gives LARGE differences, ~3 mags, same for
; MAG_APER-MAG_ISO, up to ~4 mags difference
; MAG_ISO and MAG_AUTO are pretty close though

; MAG_ISO  flux inside detection isophote
; MAG_APER flux inside given aperture (10 pixels)
; MAG_AUTO uses moments to get elliptical aperture and calculate total flux

; Maybe set a flag if MAG_APER and MAG_AUTO are VERY different.
; but what to do with extended sources where you expect there to be
; differences?

plotc,cat2.cmag,cat2.cmag-model_mag2,cat2.fwhm_world*3600*4,ps=3,xr=[10,20],yr=[-5,5],max=40
; FWHM is also really HUGE!

; ELLIPTICITY/ELONGATION is high as well!
; It must think the background flux is part of the object or
; something.

; NO, NO, NO!!  THESE ARE JUST THE SATURATED STARS!!!
; The measured magnitudes are much fainter because the saturated
; pixels are interpolated over with lower values. IMAFLAGS_ISO=7
; I think the 7 is a bitwise OR of 3 and 4
plotc,gaia2.gmag,cat2.cmag-model_mag2,cat2.imaflags_iso,ps=3,xr=[8,21],yr=[-5,5]
; Yep, these are the bad ones

;I should have used
;BACKPHOTO_TYPE LOCAL instead of GLOGAL and that would have
;calculated a more local backgrounds.

; DO WE REALLY NEED TO FIND THE SATURATED *SOURCES*??  Why not just
; look for detected sources near saturated pixels and grow the bad
; pixels and figure out which detected sources it affects.

stop

; OLD STUFF

; Use MATCHALL_SPH to find the closest matches
dcr = 10.0       ; arcsec
res = matchall_sph(brtgaia.ra_icrs,brtgaia.de_icrs,cat.ra,cat.dec,dcr/3600,nmatch,distance=distance)
; res gives reverse indices

;plot,cat.ra,cat.dec,ps=3
;oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=8,co=250

setdisp

ccdstr = importascii('~/projects/noaosourcecatalog/pro/decam_chips.txt',/header)

; Loop through bright sources and find close neighbors
For i=0,nbrtgaia-1 do begin
  brt1 = brtgaia[i]
   
  ;plot,cat.ra,cat.dec,ps=3
  ;oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=8,co=250
  ;oplot,[brtgaia[i].ra_icrs],[brtgaia[i].de_icrs],ps=4,co=150,sym=4
  
  if res[i] ne res[i+1] then begin
    ind = res[res[i]:res[i+1]-1]
    nind = n_elements(ind)
    dist = distance[res[i]-res[0]:res[i+1]-1-res[0]]*3600
    print,'Gmag=',stringize(brt1.gmag,ndec=2)
    print,'nmatches = ',strtrim(nind,2)
    ;print,'     NUMBER     DISTANCE      CMAG      FLAGS      SNR      ELLIPTICITY'
    print,'  NUM   DIST  CMAG  CPFLAGS SEFLAGS   SNR   ELLIPTICITY'
    writecol,-1,indgen(nind)+1,dist,cat[ind].cmag,cat[ind].imaflags_iso,$
             cat[ind].flags,1/cat[ind].cerr,cat[ind].ellipticity,fmt='(I5,F7.2,F7.2,2I7,F9.1,F9.2)'

    ; Load the image
    ;imfile = '/mss1/archive/pipeline/Q20141118/DEC13B/20140105/c4d_140108_022046_ooi_z_v1.fits.fz'
    imfile = expdir+'/'+base+'.fits'
    ccdnum = cat[ind[0]].ccdnum
    ccdind = where(ccdstr.ccdnum eq ccdnum,nccdind)
    fits_read,imfile,im,head,extname=ccdstr[ccdind].detpos

    maskfile='/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
    tempfile = mktemp('mask',outdir=expdir)
    file_delete,tempfile
    spawn,['funpack','-E',ccdstr[ccdind].detpos,'-O',tempfile,maskfile],out,errout,/noshell
    fits_read,tempfile,mask,mhead
    file_delete,tempfile

    head_adxy,head,brtgaia[i].ra_icrs,brtgaia[i].de_icrs,xbrt,ybrt,/deg
    head_adxy,head,cat.ra,cat.dec,x,y,/deg
    buff = 100
    xr = [xbrt-buff,xbrt+buff]
    yr = [ybrt-buff,ybrt+buff]
    loadcol,3
    displayc,im,xr=xr,yr=yr,/z
    ;displayc,im,xr=xr,yr=yr,/log,min=median(im)*0.9,max=median(im)*1.2
    ;plot,cat.ra,cat.dec,ps=3,xr=[-0.005,0.005]/cos(cendec/!radeg)+brtgaia[i].ra_icrs,$
    ;  yr=[-0.005,0.005]+brtgaia[i].de_icrs,xs=1,ys=1
    ;oplot,brtgaia.ra_icrs,brtgaia.de_icrs,ps=1,co=250,sym=2
    ;oplot,[cat[ind].ra],[cat[ind].dec],ps=4,co=150
    loadct,39,/silent
    oplot,[xbrt],[ybrt],ps=1,co=250,sym=2
    oplot,[x[ind]],[y[ind]],ps=4,co=150
    oplot,[x],[y],ps=3

    ;gd = where(cat.ra ge brtgaia[i].ra_icrs-0.005/cos(cendec/!radeg) and cat.ra le brtgaia[i].ra_icrs+0.005/cos(cendec/!radeg) and $
    ;           cat.dec ge brtgaia[i].de_icrs-0.005 and cat.dec le brtgaia[i].de_icrs+0.005,ngd)
    gd = where(x ge xbrt-buff and x le xbrt+buff and $
               y ge ybrt-buff and y le ybrt+buff,ngd)
    for j=0,ngd-1 do begin
      cat1 = cat[gd[j]]
      coords = ellipsecoords(cat1.a_world,cat1.b_world,cat1.ra,cat1.dec,cat1.theta_world,head=head,/world2pix)
      ;coords = ellipsecoords(cat1.a_world,cat1.b_world,cat1.ra,cat1.dec,cat1.theta_world)
      oplot,coords[0,*],coords[1,*]
    endfor
    
    stop
  endif else print,'no matches'

Endfor


stop

end
