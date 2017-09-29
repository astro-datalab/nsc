pro find_spurious_sources,expdir,spurind,usemask=usemask

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

t00 = systime(1)

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
; Load the metadata information
meta = mrdfits(expdir+'/'+base+'_meta.fits',1)
filter = meta.filter
fwhm = meta.fwhm  ; arcsec
; Load one chip-level catalog
cat1 = mrdfits(expdir+'/'+base+'_1.fits',1)
head1 = cat1.field_header_card
pixscale = sxpar(head1,'SEXPXSCL')
fwhmpix = fwhm/pixscale

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

; Get central coordinates and radius
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
ref = GETREFDATA(filter,cenra,cendec,radius)

; Only keep sources in the observed region
rotsphcen,cat.ra,cat.dec,cenra,cendec,lon,lat,/gnomic
rotsphcen,ref.ra,ref.dec,cenra,cendec,glon,glat,/gnomic
rlon = [min(lon)-0.05,max(lon)+0.05]
rlat = [min(lat)-0.05,max(lat)+0.05]
gdref = where(glon ge rlon[0] and glon le rlon[1] and $
              glat ge rlat[0] and glat le rlat[1],ngdref)
ref = ref[gdref]


; CALCULATE SATURATION LEVEL
;----------------------------
zp = median(cat.cmag-cat.mag_auto)
flux = 10^((25.0-ref.model_mag+zp)/2.5)  ; in counts
add_tag,ref,'flux',0.0,ref
ref.flux = flux
; 2D Gaussian profile
gmaxflux = flux/(1.138*fwhmpix^2)        ; in counts
add_tag,ref,'gmaxflux',0.0,ref
ref.gmaxflux = gmaxflux
; 2D Moffat profile
; beta = 2.5 or 3.0
; profile = (beta-1)/(pi*alpha^2)*(1+(r/alpha)^2)^(-beta) 
; max flux = flux * (beta-1)/(pi*alpha^2)
beta = 3.0d0
alpha = double( fwhmpix/(2*sqrt(2^(1.0/beta)-1)) )
mmaxflux = flux * (beta-1)/(!dpi*alpha^2)
add_tag,ref,'mmaxflux',0.0,ref
ref.mmaxflux = mmaxflux

; Find the saturated stars
satlevel = 50000.
badsat = where(ref.model_mag lt 50 and gmaxflux gt satlevel,nbadsat)
satmag = max(ref[badsat].model_mag)
print,'Saturation magnitude = ',stringize(satmag,ndec=2),' mag'
satref0 = ref[badsat]
print,strtrim(nbadsat,2),' saturated Ref stars'
; saturation level for stars, galaxies can be brighter because
;  they aren't as peaky, so check measured FWHM if they are detected
; If they are in the 2MASS-PSC then they should be STARS!!!

; Get unique chip names
uiccd = uniq(cat.ccdnum,sort(cat.ccdnum))
uccdnum = cat[uiccd].ccdnum
nccd = n_elements(uccdnum)
print,strtrim(nccd,2),' unique ccds'

; Use the mask image to make sure they are saturated
if keyword_set(usemask) then begin

  ; Copy mask file to local directory
  maskfile = '/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
  tempmaskfile = localdir+'dnidever/nsc/instcal/'+file_basename(maskfile)
  file_copy,maskfile,tempmaskfile,/over
  spawn,['funpack','-D',tempmaskfile],out,errout,/noshell
  tempmaskfile = strmid(tempmaskfile,0,strlen(tempmaskfile)-3)

  ; Loop over all chips
  for i=0,nccd-1 do begin
    fits_read,tempmaskfile,mask,mhead,exten=i+1,/no_abort,message=message
    if message ne '' then goto,BOMB
    nx = sxpar(mhead,'NAXIS1')
    ny = sxpar(mhead,'NAXIS2')
    satmask = fix( ((mask and 3) eq 3) )
    ;satmask3 = convol(float(satmask),[[1,1,1],[1,1,1],[1,1,1]])

    ;imfile = '/dl1/users/dnidever/nsc/instcal/c4d/20140108/c4d_140108_022046_ooi_z_v1/c4d_140108_022046_ooi_z_v1.fits'
    ;fits_read,imfile,im,head,exten=i+1

    ; Get ref sources in this chip
    head_xyad,mhead,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/deg
    rra = minmax(vra)
    rdec = minmax(vdec)
    gd1 = where(satref0.ra ge rra[0]-0.001 and satref0.ra le rra[1]+0.001 and $
               satref0.dec ge rdec[0]-0.001 and satref0.dec le rdec[1]+0.001,ngd1)
    satref1 = satref0[gd1]
    head_adxy,mhead,satref1.ra,satref1.dec,gx1,gy1,/deg
    gd2 = where(gx1 ge 0 and gx1 le (nx-1) and gy1 ge 0 and gy1 le (ny-1),ngd2)
    satref2 = satref1[gd2]
    gx2 = 0 > round(gx1[gd2]) < (nx-1)
    gy2 = 0 > round(gy1[gd2]) < (ny-1)
    satval = satmask[gx2,gy2]
    ;satval3 = satmask3[gx2,gy2]

    ; require that to be 3 neighboring pixels need to saturated as well
    ;gbad = where(satval eq 1 and satval3 ge 4,ngbad)
    gbad = where(satval eq 1,ngbad)
    if ngbad gt 0 then begin
      ;badref1 = ref2[gbad]
      ;push,badref,badref1
      push,badrefsource,satref2.source
    endif  
    print,strtrim(ngbad,2),'/',strtrim(ngd2,2)

    ; THIS IS PICKING UP SOME FAINT STARS THAT HAPPEN TO LAND ON A CR AND
    ; AND OTHER SATURATED PIXELS, bleed trails/columns

    ;stop
    BOMB:

  endfor
  ; this loop takes about 10s and 7s is reading the fits file
  ; this is without the mask convolution

  MATCH,ref.source,badrefsource,ind1,ind2,/sort,count=nsatref
  satref = ref[ind1]
  print,'Final saturated sources = ',strtrim(nsatref,2)

  ; Almost all sources that FALL in the chip end up being kept as saturated

; Don't use mask, only keep source inside chips
endif else begin

  ; Remove "saturated" sources that were observed and have OKAY CP flags
  srcmatch,satref0.ra,satref0.dec,cat.ra,cat.dec,0.5,ind1,ind2,/sph,count=nsatmatch
  if nsatmatch gt 0 then begin
    okayind = where(cat[ind2].imaflags_iso eq 0,nokayind)
    ; Some okay ones to remove from SATREF0
    if nokayind gt 0 then begin
      if nokayind eq n_elements(satref0) then begin
        print,'No saturated stars in image'
        return
      endif
      remove,ind1[okayind],satref0
    endif
  endif

  ; Loop over all chips
  undefine,satind
  for i=0,nccd-1 do begin
    cathd = mrdfits(expdir+'/'+base+'_'+strtrim(uccdnum[i],2)+'.fits',1,/silent)
    head = cathd.field_header_card    
    nx = sxpar(head,'NAXIS1')
    ny = sxpar(head,'NAXIS2')
    ; Get ref sources in this chip
    head_xyad,head,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/deg
    rra = minmax(vra)
    rdec = minmax(vdec)
    gd1 = where(satref0.ra ge rra[0]-0.001 and satref0.ra le rra[1]+0.001 and $
                satref0.dec ge rdec[0]-0.001 and satref0.dec le rdec[1]+0.001,ngd1)
    satref1 = satref0[gd1]
    head_adxy,head,satref1.ra,satref1.dec,gx1,gy1,/deg
    gd2 = where(gx1 ge 0 and gx1 le (nx-1) and gy1 ge 0 and gy1 le (ny-1),ngd2)
    push,satind,gd1[gd2]
  endfor
  if n_elements(satind) gt 0 then begin
    satref = satref0[satind]
  endif else begin
    print,'No saturated stars in image'
    return
  endelse
  nsatref = n_elements(satref)
endelse

; Find candidate spurious nearby sources
;  use MATCHALL_SPH to find the closest matches
dcr = 10.0       ; arcsec
; THIS MIGHT BE TOO SMALL FOR SOME VERY SATURATED STARS
res = matchall_sph(satref.ra,satref.dec,cat.ra,cat.dec,dcr/3600,nmatch,distance=distance)
; res gives reverse indices

ccdstr = importascii('~/projects/noaosourcecatalog/pro/decam_chips.txt',/header)

; Loop through bright sources and find close neighbors
undefine,allspur
For i=0,nsatref-1 do begin
  sat1 = satref[i]

  ; Some matches  
  if res[i] ne res[i+1] then begin
    ind = res[res[i]:res[i+1]-1]
    nind = n_elements(ind)
    dist = distance[res[i]-res[0]:res[i+1]-1-res[0]]*3600
    ccdnum = cat[ind[0]].ccdnum
    ;print,'Gmag=',stringize(sat1.gmag,ndec=2)
    print,strtrim(i+1,2),'/',strtrim(nsatref,2),' gmag=',stringize(sat1.gmag,ndec=2),$
          '  nmatches=',strtrim(nind,2)
    ;print,'     NUMBER     DISTANCE      CMAG      FLAGS      SNR      ELLIPTICITY'
    verbose = 0
    if keyword_set(verbose) then begin
      print,'  NUM   DIST  CMAG  CPFLAGS SEFLAGS   SNR   FWHM  ELLIPTICITY'
      writecol,-1,indgen(nind)+1,dist,cat[ind].cmag,cat[ind].imaflags_iso,$
               cat[ind].flags,1/cat[ind].cerr,cat[ind].fwhm_world*3600,cat[ind].ellipticity,$
               fmt='(I5,F7.2,F7.2,2I7,F9.1,F6.1,F9.2)'
    endif

    ; Load the image
    ;imfile = '/mss1/archive/pipeline/Q20141118/DEC13B/20140105/c4d_140108_022046_ooi_z_v1.fits.fz'
    ;imfile = expdir+'/'+base+'.fits'
    ;ccdind = where(ccdstr.ccdnum eq ccdnum,nccdind)
    ;fits_read,imfile,im,head,extname=ccdstr[ccdind].detpos

    ;maskfile='/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
    ;tempfile = mktemp('mask',outdir=expdir)
    ;file_delete,tempfile
    ;spawn,['funpack','-E',ccdstr[ccdind].detpos,'-O',tempfile,maskfile],out,errout,/noshell
    ;fits_read,tempfile,mask,mhead
    ;file_delete,tempfile

    ; I need to figure out what the "effective" radius of the star is in
    ; the image.
    ; Maybe where the model 2D Gaussian flux is some sigma above the background?
    ; you get Lorentzian wings
    ; Moffat function
    ; https://www.gnu.org/software/gnuastro/manual/html_node/PSF.html
    ; beta = 2.5 or 3.0
    ; alpha = fwhm/(2*sqrt(2^(1.0/beta)-1))
    ; profile = h*(1+((x-x0)/alpha)^2)^(-beta) + background
    ; normalized 2D Moffat function is:
    ; profile = (beta-1)/(pi*alpha^2)*(1+(r/alpha)^2)^(-beta) 
    ; smaller beta is broader wings

    ; At what point does the profile drop to a certain absolute height
    ; above the background, maybe do it in terms of sigma of the image
    ; sigma = background noise and readnoise
    ; sigma = sqrt(background*gain + rdnoise^2)/gain
    ;background = median(im)
    ;rdnoise = sxpar(head,'rdnoisea')
    ;gain = sxpar(head,'gaina')
    ;noise = sqrt(background*gain + rdnoise^2)/gain
    ; can we get these values from the SE catalog meta-data?
    ; HDU1 in single-chip files, got gain, rdnoise, fwhm
    ; can get background from cat.background
    cathd = mrdfits(expdir+base+'_'+strtrim(ccdnum,2)+'.fits',1,/silent)
    head = cathd.field_header_card
    ;background = cat[ind].background
    background = sxpar(head,'SEXBKGND')
    rdnoise = sxpar(head,'rdnoisea')
    gain = sxpar(head,'SEXGAIN')
    noise = sqrt(background*gain + rdnoise^2)/gain

    ; The radius of height h is
    ; r_h = alpha*sqrt( (Io/h)^(1/beta) - 1 )
    ; where Io is the constant in front
    hlim = noise  ; 3*noise
    I0 = sat1.flux * (beta-1)/(!dpi*alpha^2)
    rlim = alpha*sqrt( (I0/hlim)^(1/beta) - 1)
    pixscale = 0.25  ; for decam
    print,'Limiting radius = ',stringize(rlim,ndec=1),' pixels = ',stringize(rlim*pixscale,ndec=2),' arcsec'

    ; Select spurious sources
    spurind = where(dist le rlim*pixscale and (cat[ind].cmag-sat1.model_mag) > 1 and $
                    1.0/cat[ind].cerr lt 20,nspur)
    if nspur gt 0 then begin
      spur = ind[spurind]
      push,allspur,spur
    endif
    ; elliptical, fwhm, SE flags
    print,strtrim(nspur,2),' spurious source(s) found'


    ; Plotting
    pl = 0  ;1
    if keyword_set(pl) then begin
      ; Load the flux file
      ccdind = where(ccdstr.ccdnum eq ccdnum,nccdind)
      imfile = '/mss1/archive/pipeline/Q20141118/DEC13B/20140105/c4d_140108_022046_ooi_z_v1.fits.fz'
      fits_read,imfile,im,head,extname=ccdstr[ccdind].detpos

      ; Load the mask image
      maskfile='/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
      tempfile = mktemp('mask',outdir=expdir)
      file_delete,tempfile
      spawn,['funpack','-E',ccdstr[ccdind].detpos,'-O',tempfile,maskfile],out,errout,/noshell
      fits_read,tempfile,mask,mhead
      file_delete,tempfile

      head_adxy,head,satref[i].ra,satref[i].dec,xsat,ysat,/deg
      head_adxy,head,cat.ra,cat.dec,x,y,/deg
      buff = 100
      xr = [xsat-buff,xsat+buff]
      yr = [ysat-buff,ysat+buff]
      setdisp,/silent
      loadcol,3
      satmask = fix( ((mask and 3) eq 3) or ((mask and 4) eq 4) )
      displayc,im*(1-satmask),xr=xr,yr=yr,/z
      ;displayc,im,xr=xr,yr=yr,/log,min=median(im)*0.9,max=median(im)*1.2
      ;plot,cat.ra,cat.dec,ps=3,xr=[-0.005,0.005]/cos(cendec/!radeg)+satref[i].ra,$
      ;  yr=[-0.005,0.005]+satref[i].dec,xs=1,ys=1
      ;oplot,satref.ra,satref.dec,ps=1,co=250,sym=2
      ;oplot,[cat[ind].ra],[cat[ind].dec],ps=4,co=150
      loadct,39,/silent
      oplot,[xsat],[ysat],ps=1,co=250,sym=2
      oplot,[x[ind]],[y[ind]],ps=4,co=150
      if nspur gt 0 then oplot,[x[spur]],[y[spur]],ps=6,co=250
      oplot,[x],[y],ps=1,sym=0.5,co=80
      phi = scale_vector(findgen(50),0,2*!dpi)
      oplot,rlim*sin(phi)+xsat,rlim*cos(phi)+ysat

      ;gd = where(cat.ra ge satref[i].ra-0.005/cos(cendec/!radeg) and cat.ra le satref[i].ra+0.005/cos(cendec/!radeg) and $
      ;           cat.dec ge satref[i].dec-0.005 and cat.dec le satref[i].dec+0.005,ngd)
      gd = where(x ge xsat-buff and x le xsat+buff and $
                 y ge ysat-buff and y le ysat+buff,ngd)
      for j=0,ngd-1 do begin
        cat1 = cat[gd[j]]
        coords = ellipsecoords(cat1.a_world,cat1.b_world,cat1.ra,cat1.dec,cat1.theta_world,head=head,/world2pix)
        ;coords = ellipsecoords(cat1.a_world,cat1.b_world,cat1.ra,cat1.dec,cat1.theta_world)
        oplot,coords[0,*],coords[1,*],co=80
      endfor
    endif ; plotting    

    ;stop
  endif ;else print,'no matches'

Endfor

; Get unique elements
ui = uniq(allspur,sort(allspur))
spurind = allspur[ui]
nspur = n_elements(spurind)
print,strtrim(nspur,2),' final spurious sources'
spurious = cat[spurind]

; Final good indices
gdind = lindgen(n_elements(cat))
remove,spurind,gdind

print,'dt = ',systime(1)-t00,' sec'

; 20s on a "normal" uncrowded field
; 10s without using the mask

; Should we do a final cross-match of spurious sources with Gaia to
; see if we have matches.  If there are matches, then only toss them
; if they have bad CPFLAGS.  My worry is accidentally throwing out
; good bright but not saturated stars.
; YEP, 142 matches in 0.5".
;srcmatch,ref.ra,ref.dec,cat[allspur].ra,cat[allspur].dec,0.5,ind1,ind2,/sph,count=ngmatch
; most (113/142) have bad flags, but 29 do not.  they might be "good".
; double-check

stop

end
