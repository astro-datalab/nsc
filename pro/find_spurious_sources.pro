;+
;
; FIND_SPURIOUS_SOURCES
;
; Find spurious sources in a single exposure catalog 
;
; INPUTS:
;  expdir     The exposure directory.
;  /usemask   Use the quality mask image information to make sure
;               the bright stars are saturated.
;  /verbose   Print a lot of information to the screen.  verbose=2 is
;               even more verbose.
;  =logfile   Name of log file.
;  /diagplot  Show diagnostic plots.
;  /stp       Stop at end of program.
;
; OUTPUTS:
;  spurind    Index of spurious sources in catalog.
;  =count     Number of spurious sources.
;  =spurious  Structure of spurious catalog sources.
;  =cat       Observed catalog.
;  =ref       Reference catalog
;
; USAGE:
;  IDL>find_spurious_sources,expdir,spurind,usemask=usemask,spurious=spurious
; 
; By D. Nidever  Oct 2017
;-

pro find_spurious_sources,expdir,spurind,count=count,usemask=usemask,spurious=spurious,cat=cat,ref=ref,$
                          verbose=verbose,logfile=logfile,diagplot=diagplot,stp=stp

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

; Initialize output values
undefine,spurind,spurious
cound = 0

t00 = systime(1)

;expdir = '/dl1/users/dnidever/nsc/instcal/c4d/20140108/c4d_140108_022046_ooi_z_v1/'

; Not enough inputs
if n_elements(expdir) eq 0 then begin
  print,'Syntax - find_spurious_sources,expdir,spurind,usemask=usemask,spurious=spurious,cat=cat,ref=ref'
  return
endif

; Default settings
if n_elements(verbose) eq 0 then verbose=0
if n_elements(logfile) eq 0 then logf=-1 else logf=logfile
if n_elements(diagplot) eq 0 then diagplot=0

; Getting directory names
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

; Setting pool thread values
if n_elements(ncpu) eq 0 then ncpu=1
CPU, TPOOL_NTHREADS = ncpu

base = file_basename(expdir)


printlog,logf,'Find spurious sources in exposure ',base,' in ',expdir

; What instrument is this?
instrument = 'c4d'  ; by default
if stregex(expdir,'/k4m/',/boolean) eq 1 then instrument='k4m'
if stregex(expdir,'/ksb/',/boolean) eq 1 then instrument='ksb'
printlog,logf,'This is a '+instrument+' exposure'

; Get extension names for CCDs
case instrument of
'c4d': chipstr = importascii('~/projects/noaosourcecatalog/params/decam_chips.txt',/header,/silent)
'k4m': chipstr = importascii('~/projects/noaosourcecatalog/params/k4m_chips.txt',/header,/silent)
'ksb': chipstr = importascii('~/projects/noaosourcecatalog/params/ksb_chips.txt',/header,/silent)
else: stop,instrument+' not supported'
endcase


; --- Load the catalogs ---
cat = mrdfits(expdir+'/'+base+'_cat.fits',1,/silent)
ncat = n_elements(cat)
; Load the metadata information
meta = mrdfits(expdir+'/'+base+'_meta.fits',1,/silent)
filter = meta.filter
fwhm = meta.fwhm  ; arcsec
; Load one chip-level catalog
cat1 = mrdfits(expdir+'/'+base+'_1.fits',1,/silent)
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

; Add some columns to REF structure
; flux, gmaxflux, mmaxflux, ccdnum
schema = ref[0]
struct_assign,{dum:''},schema
schema = create_struct(schema,'flux',0.0,'gmaxflux',0.0,'mmaxflux',0.0,'ccdnum',0L,'rlim',0.0)
old = ref
ref = replicate(schema,n_elements(ref))
struct_assign,old,ref,/nozero
undefine,old

; CALCULATE SATURATION LEVEL
;----------------------------
zp = median(cat.cmag-cat.mag_auto)
printlog,logf,'Zero-point offset = ',stringize(zp,ndec=2),' mag'
flux = 10^((25.0-ref.model_mag+zp)/2.5)  ; in counts
ref.flux = flux
; 2D Gaussian profile
gmaxflux = flux/(1.138*fwhmpix^2)        ; in counts
ref.gmaxflux = gmaxflux
; 2D Moffat profile
; beta = 2.5 or 3.0
; profile = (beta-1)/(pi*alpha^2)*(1+(r/alpha)^2)^(-beta) 
; max flux = flux * (beta-1)/(pi*alpha^2)
beta = 3.0d0
alpha = double( fwhmpix/(2*sqrt(2^(1.0/beta)-1)) )
mmaxflux = flux * (beta-1)/(!dpi*alpha^2)
ref.mmaxflux = mmaxflux

; Find the saturated stars
satlevel = 50000.
badsat = where(ref.model_mag lt 50 and gmaxflux gt satlevel,nbadsat)
satmag = max(ref[badsat].model_mag)
printlog,logf,'Saturation magnitude = ',stringize(satmag,ndec=2),' mag'
satref0 = ref[badsat]
printlog,logf,strtrim(nbadsat,2),' saturated reference stars'
; saturation level for stars, galaxies can be brighter because
;  they aren't as peaky, so check measured FWHM if they are detected
; If they are in the 2MASS-PSC then they should be STARS!!!

; Get unique chip names
uiccd = uniq(cat.ccdnum,sort(cat.ccdnum))
uccdnum = cat[uiccd].ccdnum
nccd = n_elements(uccdnum)
printlog,logf,strtrim(nccd,2),' unique ccds'

; Get chip information and SExtractor info
ccdstr = replicate({ccdnum:0L,nx:0L,ny:0L,cenra:0.0d0,cendec:0.0d0,$
                    pixscale:0.0,vra:dblarr(4),vdec:dblarr(4),background:0.0,$
                    rdnoise:0.0,gain:0.0,noise:0.0,alpha:float(alpha),$
                    beta:float(beta),head:ptr_new()},nccd)
for i=0,nccd-1 do begin
  cathd = mrdfits(expdir+'/'+base+'_'+strtrim(uccdnum[i],2)+'.fits',1,/silent)
  head = cathd.field_header_card    
  nx = sxpar(head,'NAXIS1')
  ny = sxpar(head,'NAXIS2')
  ; Get central coordinates
  head_xyad,head,nx/2,ny/2,cenra1,cendec1,/deg
  ; Get vertices
  head_xyad,head,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/deg
  ; Get pixscale
  ;head_xyad,head,0.0,0.0,ra1,dec1,/degree
  ;head_xyad,head,1.0,0.0,ra2,dec2,/degree
  ;pixscale = sphdist(ra1,dec1,ra2,dec2,/deg)*3600.

  ; Stuff the information in the structure
  ccdstr[i].ccdnum = sxpar(head,'ccdnum')
  ccdstr[i].nx = nx
  ccdstr[i].ny = ny
  ccdstr[i].cenra = cenra1
  ccdstr[i].cendec = cendec1
  ;ccdstr[i].pixscale = pixscale
  ccdstr[i].pixscale = sxpar(head,'SEXPXSCL')
  ccdstr[i].vra = vra
  ccdstr[i].vdec = vdec
  ccdstr[i].background = sxpar(head,'SEXBKGND')
  rdnoise = sxpar(head,'rdnoise',count=nrdnoise)
  if rdnoise eq 0 then rdnoise = sxpar(head,'rdnoisea',count=nrdnoise)
  if nrdnoise eq 0 then stop,'NO READNOISE'
  ccdstr[i].rdnoise = rdnoise
  ccdstr[i].gain = sxpar(head,'SEXGAIN')
  ccdstr[i].noise = sqrt(ccdstr[i].background*ccdstr[i].gain + ccdstr[i].rdnoise^2)/ccdstr[i].gain
  ccdstr[i].head = ptr_new(head)
endfor
; this takes ~0.4s to run

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
;cathd = mrdfits(expdir+base+'_'+strtrim(ccdnum,2)+'.fits',1,/silent)
;head = cathd.field_header_card
;background = sxpar(head,'SEXBKGND')
;rdnoise = sxpar(head,'rdnoisea')
;gain = sxpar(head,'SEXGAIN')
;noise = sqrt(background*gain + rdnoise^2)/gain
;; The radius of height h is
;; r_h = alpha*sqrt( (Io/h)^(1/beta) - 1 )
;; where Io is the constant in front
;hlim = noise  ; 3*noise
;I0 = sat1.flux * (beta-1)/(!dpi*alpha^2)
;rlim = alpha*sqrt( (I0/hlim)^(1/beta) - 1)
;pixscale = 0.25  ; for decam
;printlog,logf,'Limiting radius = ',stringize(rlim,ndec=1),' pixels = ',stringize(rlim*pixscale,ndec=2),' arcsec'


; Use the mask image to make sure they are saturated
if keyword_set(usemask) then begin

  ; Filenames
  fluxfile = meta.file
  if strmid(mssdir,0,4) eq '/net' and strmid(fluxfile,0,4) ne '/net' then fluxfile='/net'+fluxfile
  ; assume mask file is in same location for now
  maskfile = repstr(fluxfile,'_ooi_','_ood_')
  if file_test(maskfile) eq 0 then stop,'MASK FILE NOT FOUND'

  ; Copy mask file to local directory
  ;maskfile = '/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
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
    ccdnum = long(sxpar(mhead,'ccdnum'))
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
    satref0[gd1].ccdnum = long(sxpar(head,'ccdnum'))

    ; Calculate source size on image in pixels
    MATCH,ccdstr.ccdnum,sat1.ccdnum,ccdind,ind2,/sort,count=nmatch
    ccdstr1 = ccdstr[ccdind[0]]
    head = *ccdstr1.head
    ; The radius of height h is
    ; r_h = alpha*sqrt( (Io/h)^(1/beta) - 1 )
    ; where Io is the constant in front
    hlim = ccdstr1.noise  ; 3*noise
    I0 = satref0[gd1].flux * (ccdstr1.beta-1)/(!dpi*ccdstr1.alpha^2)
    satref0[gd1].rlim = ccdstr1.alpha*sqrt( (I0/hlim)^(1/ccdstr1.beta) - 1)

    ; More selection based on mask
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
    printlog,logf,strtrim(ngbad,2),'/',strtrim(ngd2,2)

    ; THIS IS PICKING UP SOME FAINT STARS THAT HAPPEN TO LAND ON A CR AND
    ; AND OTHER SATURATED PIXELS, bleed trails/columns

    ;stop
    BOMB:

  endfor
  ; this loop takes about 10s and 7s is reading the fits file
  ; this is without the mask convolution

  MATCH,ref.source,badrefsource,ind1,ind2,/sort,count=nsatref
  satref = ref[ind1]
  printlog,logf,'Final saturated sources = ',strtrim(nsatref,2)

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
        printlog,logf,'No saturated stars in image'
        count = 0
        return
      endif
      remove,ind1[okayind],satref0
    endif
  endif

  ; Loop over all chips
  undefine,satind
  for i=0,nccd-1 do begin
    MATCH,uccdnum[i],ccdstr.ccdnum,ind1,ccdind,/sort,count=nmatch
    if nmatch eq 0 then stop,'NO CCD MATCH'
    ccdstr1 = ccdstr[ccdind[0]]
    head = *ccdstr1.head
    vra = ccdstr1.vra
    vdec = ccdstr1.vdec
    ; Get ref sources in this chip
    rra = minmax(vra)
    rdec = minmax(vdec)
    gd1 = where(satref0.ra ge rra[0]-0.001 and satref0.ra le rra[1]+0.001 and $
                satref0.dec ge rdec[0]-0.001 and satref0.dec le rdec[1]+0.001,ngd1)
    satref1 = satref0[gd1]
    head_adxy,head,satref1.ra,satref1.dec,gx1,gy1,/deg
    gd2 = where(gx1 ge 0 and gx1 le (nx-1) and gy1 ge 0 and gy1 le (ny-1),ngd2)
    gdind = gd1[gd2]
    satref0[gdind].ccdnum = uccdnum[i]

    ; Calculate source size on image in pixels
    ; The radius of height h is
    ; r_h = alpha*sqrt( (Io/h)^(1/beta) - 1 )
    ; where Io is the constant in front
    hlim = ccdstr1.noise  ; 1xnoise
    I0 = satref0[gdind].flux * (ccdstr1.beta-1)/(!dpi*ccdstr1.alpha^2)
    satref0[gdind].rlim = ccdstr1.alpha*sqrt( (I0/hlim)^(1/ccdstr1.beta) - 1)

    ; Save the indices
    push,satind,gd1[gd2]
  endfor
  if n_elements(satind) gt 0 then begin
    satref = satref0[satind]
  endif else begin
    printlog,logf,'No saturated stars in image'
    count = 0
    return
  endelse
  nsatref = n_elements(satref)
endelse

; Find candidate spurious nearby sources
;  use MATCHALL_SPH to find the closest matches
dcr = max(satref.rlim)
printlog,logf,'Maximum radius is ',stringize(dcr,ndec=2),' arcsec'
res = matchall_sph(satref.ra,satref.dec,cat.ra,cat.dec,dcr/3600,nmatch,distance=distance)
; res gives reverse indices


; Loop through bright sources and find close neighbors
;-----------------------------------------------------
undefine,allspur
For i=0,nsatref-1 do begin
  sat1 = satref[i]

  MATCH,sat1.ccdnum,ccdstr.ccdnum,ind1,ccdind,/sort,count=nmatch
  ccdstr1 = ccdstr[ccdind[0]]

  ; Some matches  
  if res[i] ne res[i+1] then begin
    ind = res[res[i]:res[i+1]-1]
    nind = n_elements(ind)
    dist = distance[res[i]-res[0]:res[i+1]-1-res[0]]*3600
    ccdnum = cat[ind[0]].ccdnum
    ;printlog,logf,'Gmag=',stringize(sat1.gmag,ndec=2)
    ;printlog,logf,'     NUMBER     DISTANCE      CMAG      FLAGS      SNR      ELLIPTICITY'
    if keyword_set(verbose) then begin
      printlog,logf,strtrim(i+1,2),'/',strtrim(nsatref,2),' gmag=',stringize(sat1.gmag,ndec=2),$
            '  nmatches=',strtrim(nind,2)
      printlog,logf,'Limiting radius = ',stringize(sat1.rlim,ndec=1),' pixels = ',stringize(sat1.rlim*ccdstr1.pixscale,ndec=2),' arcsec'
      if verbose eq 2 then begin
        printlog,logf,'  NUM   DIST  CMAG  CPFLAGS SEFLAGS   SNR   FWHM  ELLIPTICITY'
        writecol,-1,indgen(nind)+1,dist,cat[ind].cmag,cat[ind].imaflags_iso,$
                 cat[ind].flags,1/cat[ind].cerr,cat[ind].fwhm_world*3600,cat[ind].ellipticity,$
                 fmt='(I5,F7.2,F7.2,2I7,F9.1,F6.1,F9.2)'
      endif
    endif

    ; Select spurious sources
    ;spurind = where(dist le sat1.rlim*ccdstr1.pixscale and (cat[ind].cmag-sat1.model_mag) > 1 and $
    ;                1.0/cat[ind].cerr lt 20,nspur)
    spurind = where(dist le sat1.rlim*ccdstr1.pixscale,nspur)
    if nspur gt 0 then begin
      spur = ind[spurind]
      push,allspur,spur
    endif
    ; elliptical, fwhm, SE flags
    if keyword_set(verbose) then $
      printlog,logf,strtrim(nspur,2),' spurious source(s) found'

    ; Plotting
    if keyword_set(diagplot) then begin
      ; Filenames
      fluxfile = meta.file
      if strmid(mssdir,0,4) eq '/net' and strmid(fluxfile,0,4) ne '/net' then fluxfile='/net'+fluxfile
      ; assume mask file is in same location for now
      maskfile = repstr(fluxfile,'_ooi_','_ood_')
      if file_test(maskfile) eq 0 then stop,'MASK FILE NOT FOUND'

      ; Load the flux file
      chipind = where(chipstr.ccdnum eq ccdnum,nchipind)
      ;imfile = '/mss1/archive/pipeline/Q20141118/DEC13B/20140105/c4d_140108_022046_ooi_z_v1.fits.fz'
      FITS_READ,fluxfile,im,head,extname=chipstr[chipind].extname

      ; Load the mask image
      ;maskfile='/net/mss1/archive/pipeline/Q20141119/DEC13B/20140105/c4d_140108_022046_ood_z_v1.fits.fz'
      tempfile = mktemp('mask',outdir=expdir)
      file_delete,tempfile
      spawn,['funpack','-E',chipstr[chipind].extname,'-O',tempfile,maskfile],out,errout,/noshell
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
      oplot,sat1.rlim*sin(phi)+xsat,sat1.rlim*cos(phi)+ysat

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
      dum = ''
      read,'Hit any key to continue',dum
    endif ; plotting    
  endif ;else printlog,logf,'no matches'

Endfor  ; bright, saturated star loop

; Some spurious sources found
if n_elements(allspur) gt 0 then begin
  ; Get unique elements
  ui = uniq(allspur,sort(allspur))
  spurind = allspur[ui]
  nspur = n_elements(spurind)
  printlog,logf,strtrim(nspur,2),' final spurious sources'
  spurious = cat[spurind]

  ; Final good indices
  gdind = lindgen(n_elements(cat))
  remove,spurind,gdind

; No spurious sources
endif else begin
  printlog,logf,'No spurious sources found'
endelse
count = n_elements(spurind)

printlog,logf,'dt = ',stringize(systime(1)-t00,ndec=4),' sec'

; 20s on a "normal" uncrowded field
; 10s without using the mask
; Sped up to ~6s now

if keyword_set(stp) then stop

end
