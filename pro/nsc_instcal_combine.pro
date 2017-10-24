pro nsc_instcal_combine,pix,version=version,nside=nside,redo=redo,stp=stp,outdir=outdir

t0 = systime(1)

; Combine all the exposures that fall in this healpix pixel
;nside = 256
if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(pix) eq 0 then begin
  print,'Syntax - nsc_instcal_combine,pix,version=version,nside=nside,redo=redo,stp=stp'
  return
endif

; Check if output file already exists
if n_elements(outdir) eq 0 then outdir=dir+'combine/'
subdir = strtrim(long(pix)/1000,2)    ; use the thousands to create subdirectory grouping
outfile = outdir+'/'+subdir+'/'+strtrim(pix,2)+'.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,outfile,' EXISTS already and /redo not set'
  return
endif

print,'Combining InstCal SExtractor catalogs for Healpix pixel = ',strtrim(pix,2)

; Load the list
listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
if file_test(listfile) eq 0 then begin
  print,listfile,' NOT FOUND'
  return
endif
healstr = MRDFITS(listfile,1,/silent)
healstr.file = strtrim(healstr.file,2)
healstr.base = strtrim(healstr.base,2)
index = MRDFITS(listfile,2,/silent)
; Find our pixel
ind = where(index.pix eq pix,nind)
if nind eq 0 then begin
  print,'No entries for Healpix pixel "',strtrim(pix,2),'" in the list'
  return
endif
ind = ind[0]
list = healstr[index[ind].lo:index[ind].hi]
nlist = n_elements(list)

; GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
;  so we can deal with the edge cases
NEIGHBOURS_RING,nside,pix,neipix,nneipix
for i=0,nneipix-1 do begin
  ind = where(index.pix eq neipix[i],nind)
  if nind gt 0 then begin
    ind = ind[0]
    list1 = healstr[index[ind].lo:index[ind].hi]
    push,list,list1
  endif
endfor
; Get unique values
ui = uniq(list.file,sort(list.file))
list = list[ui]
nlist = n_elements(list)
print,strtrim(nlist,2),' exposures that overlap this pixel and neighbors'

; Get the boundary coordinates
;   healpy.boundaries but not sure how to do it in IDL
;   pix2vec_ring/nest can optionally return vertices but only 4
;     maybe subsample myself between the vectors
; Expand the boundary to include a "buffer" zone
;  to deal with edge cases
;PIX2VEC_RING,nside,pix,vec,vertex

; Use python code to get the boundary
;  this takes ~2s mostly from import statements
tempfile = MKTEMP('bnd')
file_delete,tempfile+'.fits',/allow
step = 100
pylines = 'python -c "from healpy import boundaries; from astropy.io import fits;'+$
          ' v=boundaries('+strtrim(nside,2)+','+strtrim(pix,2)+',step='+strtrim(step,2)+');'+$
          " fits.writeto('"+tempfile+".fits'"+',v)"'
spawn,pylines,out,errout
vecbound = MRDFITS(tempfile+'.fits',0,/silent)
file_delete,[tempfile,tempfile+'.fits'],/allow
VEC2ANG,vecbound,theta,phi
rabound = phi*radeg
decbound = 90-theta*radeg

; Expand the boundary by the buffer size
PIX2ANG_RING,nside,pix,centheta,cenphi
cenra = cenphi*radeg
cendec = 90-centheta*radeg
; reproject onto tangent plane
ROTSPHCEN,rabound,decbound,cenra,cendec,lonbound,latbound,/gnomic
; expand by a fraction, it's not an extact boundary but good enough
buffsize = 10.0/3600. ; in deg
radbound = sqrt(lonbound^2+latbound^2)
frac = 1.0 + 1.5*max(buffsize/radbound)
lonbuff = lonbound*frac
latbuff = latbound*frac
buffer = {cenra:cenra,cendec:cendec,lon:lonbuff,lat:latbuff}

; Initialize the ID structure
;  this will contain the SourceID, Exposure name, ObjectID
schema_idstr = {sourceid:'',exposure:'',expnum:'',objectid:'',objectindex:0LL}
idstr = replicate(schema_idstr,1e6)
nidstr = n_elements(idstr)
idcnt = 0LL

; Initialize the object structure
schema_obj = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,raerr:0.0d0,decerr:0.0d0,pmra:0.0d0,$
              pmdec:0.0d0,pmraerr:0.0d0,pmdecerr:0.0d0,mjd:0.0d0,deltamjd:0.0,ndet:0L,nphot:0L,$
              ndetu:0,nphotu:0,umag:0.0,urms:0.0,uerr:0.0,uasemi:0.0,ubsemi:0.0,utheta:0.0,$
              ndetg:0,nphotg:0,gmag:0.0,grms:0.0,gerr:0.0,gasemi:0.0,gbsemi:0.0,gtheta:0.0,$
              ndetr:0,nphotr:0,rmag:0.0,rrms:0.0,rerr:0.0,rasemi:0.0,rbsemi:0.0,rtheta:0.0,$
              ndeti:0,nphoti:0,imag:99.9,irms:0.0,ierr:0.0,iasemi:0.0,ibsemi:0.0,itheta:0.0,$
              ndetz:0,nphotz:0,zmag:0.0,zrms:0.0,zerr:0.0,zasemi:0.0,zbsemi:0.0,ztheta:0.0,$
              ndety:0,nphoty:0,ymag:0.0,yrms:0.0,yerr:0.0,yasemi:0.0,ybsemi:0.0,ytheta:0.0,$
              ndetvr:0,nphotvr:0,vrmag:0.0,vrrms:0.0,vrerr:0.0,vrasemi:0.0,vrbsemi:0.0,vrtheta:0.0,$
              x2:0.0,x2err:0.0,y2:0.0,y2err:0.0,xy:0.0,xyerr:0.0,asemi:0.0,asemierr:0.0,bsemi:0.0,$
              bsemierr:0.0,theta:0.0,thetaerr:0.0,ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
;              x2:0.0,x2err:0.0,y2:0.0,y2err:0.0,xy:0.0,xyerr:0.0,cxx:0.0,cxxerr:0.0,$
;              cxy:0.0,cxyerr:0.0,cyy:0.0,cyyerr:0.0,asemi:0.0,asemierr:0.0,bsemi:0.0,$
;              bsemierr:0.0,theta:0.0,thetaerr:0.0,elongation:0.0,$
;              ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
tags = tag_names(schema_obj)
obj = replicate(schema_obj,5e5)
nobj = n_elements(obj)
schema_totobj = {ra:0.0d0,dec:0.0d0,ramjd:0.0d0,decmjd:0.0d0,ramjd2:0.0d0,decmjd2:0.0d0,minmjd:999999.0d0,maxmjd:-999999.0d0,$
                 umag2:0.0d0,gmag2:0.0d0,rmag2:0.0d0,imag2:0.0d0,zmag2:0.0d0,ymag2:0.0d0,vrmag2:0.0d0,$
                 utot:0.0d0,gtot:0.0d0,rtot:0.0d0,itot:0.0d0,ztot:0.0d0,ytot:0.0d0,vrtot:0.0d0}
tottags = tag_names(schema_totobj)
totobj = replicate(schema_totobj,nobj)
cnt = 0LL

; Loop over the exposures
undefine,allmeta
FOR i=0,nlist-1 do begin
  print,strtrim(i+1,2),' Loading ',list[i].file

  ; Load the exposure catalog
  cat1 = MRDFITS(list[i].file,1,/silent)
  ncat1 = n_elements(cat1)
  print,'  ',strtrim(ncat1,2),' sources'

  ; Make sure it's in the right format
  if n_tags(cat1) ne 45 then begin   ; 49 for v1
    print,'  This catalog does not have the right format. Skipping'
    goto,BOMB
  endif

  metafile = repstr(list[i].file,'_cat','_meta')
  meta = MRDFITS(metafile,1,/silent)
  meta.base = strtrim(meta.base)
  meta.expnum = strtrim(meta.expnum)
  meta.mjd = date2jd(meta.dateobs,/mjd)  ; recompute because some MJD are bad
  ;head = headfits(list[i].file,exten=0)
  ;filtername = sxpar(head,'filter')
  ;if strmid(filtername,0,2) eq 'VR' then filter='VR' else filter=strmid(filtername,0,1)
  ;exptime = sxpar(head,'exptime')
  print,'  FILTER=',meta.filter,'  EXPTIME=',stringize(meta.exptime,ndec=1),' sec'

  ; Remove bad chip data
  ; Half of chip 31 for MJD>56660
  ;  c4d_131123_025436_ooi_r_v2 with MJD=56619 has problems too
  ;  if the background b/w the left and right half is large then BAd
  lft31 = where(cat1.x_image lt 1024 and cat1.ccdnum eq 31,nlft31)
  rt31 = where(cat1.x_image ge 1024 and cat1.ccdnum eq 31,nrt31)
  if nlft31 gt 10 and nrt31 gt 10 then begin
    lftback = median([cat1[lft31].background])
    rtback = median([cat1[rt31].background])
    mnback = 0.5*(lftback+rtback)
    sigback = mad(cat1.background)
    if abs(lftback-rtback) gt (sqrt(mnback)>sigback) then jump31=1 else jump31=0
    print,'  Big jump in CCDNUM 31 background levels'
  endif else jump31=0
  if meta.mjd gt 56600 or jump31 eq 1 then begin  
    badchip31 = 1
    ; Remove bad measurements
    ; X: 1-1024 okay
    ; X: 1025-2049 bad
    ; use 1000 as the boundary since sometimes there's a sharp drop
    ; at the boundary that causes problem sources with SExtractor
    bdind = where(cat1.x_image gt 1000 and cat1.ccdnum eq 31,nbdind,comp=gdind,ncomp=ngdind)
    if nbdind gt 0 then begin   ; some bad ones found
      if ngdind eq 0 then begin   ; all bad
        print,'NO useful measurements in ',fitsfile
        undefine,cat1
        ncat1 = 0
        goto,BOMB
      endif else begin
        print,'  Removing '+strtrim(nbdind,2)+' bad chip 31 measurements, '+strtrim(ngdind,2)+' left.'
        REMOVE,bdind,cat1
        ncat1 = n_elements(cat1)
      endelse
    endif  ; some bad ones to remove
  endif else badchip31=0     ; chip 31

  ; Removing BAD sources
  ; Source Extractor FLAGS
  ; 1   The object has neighbours, bright and close enough to significantly bias the MAG AUTO photometry,
  ;       or bad pixels (more than 10% of the integrated area affected),
  ; 2   The object was originally blended with another one,
  ; 4   At least one pixel of the object is saturated (or very close to),
  ; 8   The object is truncated (too close to an image boundary),
  ; 16  Object’s aperture data are incomplete or corrupted,
  ; 32  Object’s isophotal data are incomplete or corrupted,
  ; 64  A memory overflow occurred during deblending,
  ; 128 A memory overflow occurred during extraction.
  ; Don't use any of these.  That saturation one is already
  ;  covered in the CP bitmask flag

  ; --CP bit masks, Pre-V3.5.0 (PLVER)  
  ; Bit     DQ Type
  ; 1      detector bad pixel
  ; 2      saturated
  ; 4      interpolated
  ; 16     single exposure cosmic ray
  ; 64     bleed trail
  ; 128    multi-exposure transient
  ; --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks  
  ;  1 = bad (in static bad pixel mask)
  ;  2 = no value (for stacks)
  ;  3 = saturated
  ;  4 = bleed mask
  ;  5 = cosmic ray
  ;  6 = low weight
  ;  7 = diff detect   
  ; OR: the result is an arithmetic (bit-to-bit) OR of flag-map pixels.
  ; The NIMAFLAGS ISO catalog parameter contains a number of relevant
  ; flag-map pixels: the number of non-zero flag-map pixels in the
  ; case of an OR or AND FLAG TYPE.

  ; Make a cut on quality mask flag (IMAFLAGS_ISO)
  bdcat = where(cat1.imaflags_iso gt 0,nbdcat)
  if nbdcat gt 0 then begin
    print,'  Removing ',strtrim(nbdcat,2),' sources with bad CP flags.'
    if nbdcat eq ncat1 then goto,BOMB
    REMOVE,bdcat,cat1
    ncat1 = n_elements(cat1)
  endif

  ; Make cuts on SE FLAGS
  ;   this removes problematic truncatd sources near chip edges
  bdseflags = where( ((cat1.flags and 8) eq 8) or $             ; object truncated
                     ((cat1.flags and 16) eq 16),nbdseflags)    ; aperture truncate
  if nbdseflags gt 0 then begin
    print,'  Removing ',strtrim(nbdseflags,2),' truncated sources'
    if nbdseflags eq ncat1 then goto,BOMB
    REMOVE,bdseflags,cat1
    ncat1 = n_elements(cat1)
  endif

  ; Removing low-S/N sources
  ;  snr = 1.087/err
  snrcut = 5.0
  bdcatsnr = where(1.087/cat1.magerr_auto lt snrcut,nbdcatsnr)
  if nbdcatsnr gt 0 then begin
    print,'  Removing ',strtrim(nbdcatsnr,2),' sources with S/N<',strtrim(snrcut,2)
    if nbdcatsnr eq ncat1 then goto,BOMB
    REMOVE,bdcatsnr,cat1
    ncat1 = n_elements(cat1)
  endif

  ; Only include sources inside Boundary+Buffer zone
  ;  -use ROI_CUT
  ;  -reproject to tangent plane first so we don't have to deal
  ;     with RA=0 wrapping or pol issues
  ROTSPHCEN,cat1.ra,cat1.dec,buffer.cenra,buffer.cendec,lon,lat,/gnomic
  ROI_CUT,buffer.lon,buffer.lat,lon,lat,ind0,ind1,fac=100,/silent
  nmatch = n_elements(ind1)

  ; Only want source inside this pixel
  ;theta = (90-cat1.dec)/radeg
  ;phi = cat1.ra/radeg
  ;ANG2PIX_RING,nside,theta,phi,ipring
  ;MATCH,ipring,pix,ind1,ind2,/sort,count=nmatch
  if nmatch eq 0 then begin
    print,'  No sources inside this pixel'
    goto,BOMB
  endif
  print,'  ',strtrim(nmatch,2),' sources are inside this pixel'
  cat = cat1[ind1]
  ncat = nmatch

  ; Add metadata to ALLMETA
  ;  Make sure it's in the right format
  newmeta = {file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d0,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
             airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,badchip31:0B,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,$
             zpterm:0.0,zptermerr:0.0,zptermsig:0.0,nrefmatch:0L}
  STRUCT_ASSIGN,meta,newmeta
  newmeta.badchip31 = badchip31
  PUSH,allmeta,newmeta

  ; NDETX is good "detection" and morphology for this filter
  ; NPHOTX is good "photometry" for this filter
  detind = where(tags eq 'NDET'+strupcase(meta.filter),ndetind)
  magind = where(tags eq strupcase(meta.filter)+'MAG',nmagind)
  errind = where(tags eq strupcase(meta.filter)+'ERR',nerrind)
  totind = where(tottags eq strupcase(meta.filter)+'TOT',ntotind)
  mag2ind = where(tottags eq strupcase(meta.filter)+'MAG2',nmag2ind)
  detphind = where(tags eq 'NPHOT'+strupcase(meta.filter),nphotind)
  asemiind = where(tags eq strupcase(meta.filter)+'ASEMI',nasemiind)
  bsemiind = where(tags eq strupcase(meta.filter)+'BSEMI',nbsemiind)
  thetaind = where(tags eq strupcase(meta.filter)+'THETA',nthetaind)

  ; Combine the data
  ;-----------------
  ; First catalog
  If cnt eq 0 then begin

    ; Copy to final structure
    newexp = replicate(schema_obj,ncat)
    newexp.id = strtrim(pix,2)+'.'+strtrim(lindgen(ncat)+1,2)
    newexp.pix = pix
    newexp.ra = cat.ra
    newexp.dec = cat.dec
    newexp.raerr = 1.0/cat.raerr^2                           ; total(ra_wt)
    newexp.decerr = 1.0/cat.decerr^2                         ; total(dec_wt)
    newexp.pmra = (1.0/cat.raerr^2) * newmeta.mjd*cat.ra     ; total(wt*mjd*ra)
    newexp.pmdec = (1.0/cat.decerr^2) * newmeta.mjd*cat.dec  ; total(wt*mjd*dec)
    newexp.mjd = newmeta.mjd                                 ; total(mjd)
    newexp.ndet = 1
    ; Detection and morphology parameters for this FILTER
    newexp.(detind) = 1
    newexp.(asemiind) = cat.a_world
    newexp.(bsemiind) = cat.b_world
    newexp.(thetaind) = cat.theta_world
    ; Good photometry for this FILTER
    gdmag = where(cat.cmag lt 50,ngdmag)
    if ngdmag gt 0 then begin
      newexp[gdmag].(magind) = 2.5118864d^cat[gdmag].cmag * (1.0d0/cat[gdmag].cerr^2)
      newexp[gdmag].(errind) = 1.0d0/cat[gdmag].cerr^2
      newexp[gdmag].(detphind) = 1
    endif
    newexp.x2 = cat.x2_world
    newexp.x2err = cat.errx2_world^2
    newexp.y2 = cat.y2_world
    newexp.y2err = cat.erry2_world^2
    newexp.xy = cat.xy_world
    newexp.xyerr = cat.errxy_world^2
    ;newexp.cxx = cat.cxx_world
    ;newexp.cxxerr = cat.errcxx_world^2
    ;newexp.cyy = cat.cyy_world
    ;newexp.cyyerr = cat.errcyy_world^2
    ;newexp.cxy = cat.cxy_world
    ;newexp.cxyerr = cat.errcxy_world^2
    newexp.asemi = cat.a_world
    newexp.asemierr = cat.erra_world^2
    newexp.bsemi = cat.b_world
    newexp.bsemierr = cat.errb_world^2
    newexp.theta = cat.theta_world
    newexp.thetaerr = cat.errtheta_world^2
    ;newexp.elongation = cat.elongation
    newexp.ellipticity = cat.ellipticity
    newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
    newexp.flags = cat.flags
    newexp.class_star = cat.class_star
    obj[0:ncat-1] = newexp
    totobj[0:ncat-1].ra = cat.ra * (1.0/cat.raerr^2)               ; total(ra*wt)
    totobj[0:ncat-1].dec = cat.dec * (1.0/cat.decerr^2)            ; total(dec*wt)
    totobj[0:ncat-1].ramjd = (1.0/cat.raerr^2) * newmeta.mjd       ; total(wt_ra*mjd)
    totobj[0:ncat-1].decmjd = (1.0/cat.decerr^2) * newmeta.mjd     ; total(wt_dec*mjd)
    totobj[0:ncat-1].ramjd2 = (1.0/cat.raerr^2) * newmeta.mjd^2    ; total(wt_ra*mjd^2)
    totobj[0:ncat-1].decmjd2 = (1.0/cat.decerr^2) * newmeta.mjd^2  ; total(wt_dec*mjd^2)
    totobj[0:ncat-1].minmjd = newmeta.mjd
    totobj[0:ncat-1].maxmjd = newmeta.mjd
    totobj[gdmag].(totind) = cat[gdmag].cmag                       ; sum(mag)
    totobj[gdmag].(mag2ind) = double(cat[gdmag].cmag)^2            ; sum(mag^2), need dbl to precent underflow
    cnt += ncat

    ; Add new elements to IDSTR
    if idcnt+ncat gt nidstr then begin
      old = idstr
      nnew = 3e5 > ncat
      idstr = replicate(schema_idstr,nidstr+nnew)
      idstr[0:nidstr-1] = old
      nidstr = n_elements(idstr)
      undefine,old
    endif
    ; Add to IDSTR
    sourceid = strtrim(meta.expnum,2)+'.'+strtrim(cat.ccdnum,2)+'.'+strtrim(cat.number,2)
    idstr[idcnt:idcnt+ncat-1].sourceid = sourceid
    idstr[idcnt:idcnt+ncat-1].exposure = meta.base
    idstr[idcnt:idcnt+ncat-1].expnum = meta.expnum
    idstr[idcnt:idcnt+ncat-1].objectid = newexp.id
    idstr[idcnt:idcnt+ncat-1].objectindex = lindgen(ncat)
    idcnt += ncat

  ; Second and up
  Endif else begin

    ; Match new sources to the objects
    SRCMATCH,obj[0:cnt-1].ra,obj[0:cnt-1].dec,cat.ra,cat.dec,0.5,ind1,ind2,count=nmatch,/sph,/usehist  ; use faster histogram_nd method
    if nmatch gt 0 then if (max(ind2) gt ncat-1) or (max(ind1) gt cnt-1) then begin
      ; sometimes the histogram_nd method gives indices that are too large
      ; not sure why, use SRCOR method.
      print,'  Problem with HISTOGRAM_ND SRCMATCH indices, using SRCOR method instead.'
      SRCMATCH,obj[0:cnt-1].ra,obj[0:cnt-1].dec,cat.ra,cat.dec,0.5,ind1,ind2,count=nmatch,/sph,usehist=0
    endif
    print,'  ',strtrim(nmatch,2),' matched sources'
    ; Some matches, add data to existing record for these sources
    if nmatch gt 0 then begin

      ; When CMAG=99.99 the morphology parameters are still okay

      ; Combine the information
      cmb = obj[ind1]
      totcmb = totobj[ind1]
      cmb.ndet++
      newcat = cat[ind2]
      ; Coordinate errors
      cmb.raerr += (1.0/newcat.raerr^2)    ; total(ra_wt)
      cmb.decerr += (1.0/newcat.decerr^2)  ; total(dec_wt)
      ; Proper motion stuff
      cmb.pmra += (1.0/newcat.raerr^2) * newmeta.mjd*newcat.ra     ; total(wt*mjd*ra)
      cmb.pmdec += (1.0/newcat.decerr^2) * newmeta.mjd*newcat.dec  ; total(wt*mjd*dec)
      cmb.mjd += newmeta.mjd                                       ; total(mjd^2)
      ; Detection and morphology parameters for this FILTER
      cmb.(detind)++
      cmb.(asemiind) += newcat.a_world
      cmb.(bsemiind) += newcat.b_world
      cmb.(thetaind) += newcat.theta_world
      ; Good photometry for this FILTER
      gdmag = where(newcat.cmag lt 50,ngdmag)
      if ngdmag gt 0 then begin
        cmb[gdmag].(magind) += 2.5118864d^newcat[gdmag].cmag * (1.0d0/newcat[gdmag].cerr^2)
        cmb[gdmag].(errind) += 1.0d0/newcat[gdmag].cerr^2
        cmb[gdmag].(detphind) += 1
        ; NPHOTX means good PHOT detection
      endif
      cmb.x2 += newcat.x2_world
      cmb.x2err += newcat.errx2_world^2
      cmb.y2 += newcat.y2_world
      cmb.y2err += newcat.erry2_world^2
      cmb.xy += newcat.xy_world
      cmb.xyerr += newcat.errxy_world^2
      ;cmb.cxx += newcat.cxx_world
      ;cmb.cxxerr += newcat.errcxx_world^2
      ;cmb.cyy += newcat.cyy_world
      ;cmb.cyyerr += newcat.errcyy_world^2
      ;cmb.cxy += newcat.cxy_world
      ;cmb.cxyerr += newcat.errcxy_world^2
      cmb.asemi += newcat.a_world
      cmb.asemierr += newcat.erra_world^2
      cmb.bsemi += newcat.b_world
      cmb.bsemierr += newcat.errb_world^2
      cmb.theta += newcat.theta_world
      cmb.thetaerr += newcat.errtheta_world^2
      ;cmb.elongation += newcat.elongation
      cmb.ellipticity += newcat.ellipticity
      cmb.fwhm += newcat.fwhm_world*3600  ; in arcsec
      cmb.flags OR= newcat.flags
      cmb.class_star += newcat.class_star
      totcmb.ra += newcat.ra * (1.0/newcat.raerr^2)             ; total(wt*ra)
      totcmb.dec += newcat.dec * (1.0/newcat.decerr^2)          ; total(wt*dec)
      totcmb.ramjd +=  (1.0/newcat.raerr^2) * newmeta.mjd       ; total(wt_ra*mjd)
      totcmb.decmjd +=  (1.0/newcat.decerr^2) * newmeta.mjd     ; total(wt_dec*mjd)
      totcmb.ramjd2 +=  (1.0/newcat.raerr^2) * newmeta.mjd^2    ; total(wt_ra*mjd^2)
      totcmb.decmjd2 +=  (1.0/newcat.decerr^2) * newmeta.mjd^2  ; total(wt_dec*mjd^2)
      totcmb.minmjd <= newmeta.mjd
      totcmb.maxmjd >= newmeta.mjd
      totcmb[gdmag].(totind) += newcat[gdmag].cmag              ; sum(mag)
      totcmb[gdmag].(mag2ind) += double(newcat[gdmag].cmag)^2   ; sum(mag^2), need dbl to prevent underflow
      obj[ind1] = cmb  ; stuff it back in
      totobj[ind1] = totcmb

      ; Add new elements to IDSTR
      if idcnt+nmatch gt nidstr then begin
        old = idstr
        nnew = 3e5 > nmatch
        idstr = replicate(schema_idstr,nidstr+nnew)
        idstr[0:nidstr-1] = old
        nidstr = n_elements(idstr)
        undefine,old
      endif
      ; Add to IDSTR
      sourceid = strtrim(meta.expnum,2)+'.'+strtrim(newcat.ccdnum,2)+'.'+strtrim(newcat.number,2)
      idstr[idcnt:idcnt+nmatch-1].sourceid = sourceid
      idstr[idcnt:idcnt+nmatch-1].exposure = meta.base
      idstr[idcnt:idcnt+nmatch-1].expnum = meta.expnum
      idstr[idcnt:idcnt+nmatch-1].objectid = cmb.id
      idstr[idcnt:idcnt+nmatch-1].objectindex = ind1
      idcnt += nmatch

      ; Remove stars
      if nmatch lt n_elements(cat) then remove,ind2,cat else undefine,cat
      ncat = n_elements(cat)
    endif

    ; Some left, add records for these sources
    if n_elements(cat) gt 0 then begin
      print,'  ',strtrim(ncat,2),' sources left to add'

      ; Add new elements
      if cnt+ncat gt nobj then begin
        old = obj
        nnew = 3e5 > ncat
        obj = replicate(schema_obj,nobj+nnew)
        obj[0:nobj-1] = old
        oldtot = totobj
        totobj = replicate(schema_totobj,nobj+nnew)
        totobj[0:nobj-1] = oldtot
        nobj = n_elements(obj)
        undefine,old
      endif

      ; Copy to final structure
      newexp = replicate(schema_obj,ncat)
      newexp.id = strtrim(pix,2)+'.'+strtrim(cnt+lindgen(ncat)+1,2)
      newexp.pix = pix
      newexp.ra = cat.ra
      newexp.dec = cat.dec
      newexp.raerr = 1.0/cat.raerr^2                          ; total(ra_wt)
      newexp.decerr = 1.0/cat.decerr^2                        ; total(dec_wt)
      newexp.pmra = (1.0/cat.raerr^2) * newmeta.mjd*cat.ra    ; total(wt*mjd*ra)
      newexp.pmdec = (1.0/cat.decerr^2) * newmeta.mjd*cat.dec ; total(wt*mjd*dec)
      newexp.mjd = newmeta.mjd                                ; total(mjd) 
      newexp.ndet = 1
      ; Detection and morphology parameters for this FILTER
      newexp.(detind) = 1
      newexp.(asemiind) = cat.a_world
      newexp.(bsemiind) = cat.b_world
      newexp.(thetaind) = cat.theta_world
      gdmag = where(cat.cmag lt 50,ngdmag)
      if ngdmag gt 0 then begin
        newexp[gdmag].(magind) = 2.5118864d^cat[gdmag].cmag * (1.0d0/cat[gdmag].cerr^2)
        newexp[gdmag].(errind) = 1.0d0/cat[gdmag].cerr^2
        newexp[gdmag].(detphind) = 1
      endif
      newexp.x2 = cat.x2_world
      newexp.x2err = cat.errx2_world^2
      newexp.y2 = cat.y2_world
      newexp.y2err = cat.erry2_world^2
      newexp.xy = cat.xy_world
      newexp.xyerr = cat.errxy_world^2
      ;newexp.cxx = cat.cxx_world
      ;newexp.cxxerr = cat.errcxx_world^2
      ;newexp.cyy = cat.cyy_world
      ;newexp.cyyerr = cat.errcyy_world^2
      ;newexp.cxy = cat.cxy_world
      ;newexp.cxyerr = cat.errcxy_world^2
      newexp.asemi = cat.a_world
      newexp.asemierr = cat.erra_world^2
      newexp.bsemi = cat.b_world
      newexp.bsemierr = cat.errb_world^2
      newexp.theta = cat.theta_world
      newexp.thetaerr = cat.errtheta_world^2
      ;newexp.elongation = cat.elongation
      newexp.ellipticity = cat.ellipticity
      newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
      newexp.flags = cat.flags
      newexp.class_star = cat.class_star
      obj[cnt:cnt+ncat-1] = newexp   ; stuff it in
      totobj[cnt:cnt+ncat-1].ra = cat.ra * (1.0/cat.raerr^2)               ; total(ra*wt)
      totobj[cnt:cnt+ncat-1].dec = cat.dec * (1.0/cat.decerr^2)            ; total(dec*wt)
      totobj[cnt:cnt+ncat-1].ramjd = (1.0/cat.raerr^2) * newmeta.mjd       ; total(wt_ra*mjd)
      totobj[cnt:cnt+ncat-1].decmjd = (1.0/cat.decerr^2) * newmeta.mjd     ; total(wt_dec*mjd)
      totobj[cnt:cnt+ncat-1].ramjd2 = (1.0/cat.raerr^2) * newmeta.mjd^2    ; total(wt_ra*mjd^2)
      totobj[cnt:cnt+ncat-1].decmjd2 = (1.0/cat.decerr^2) * newmeta.mjd^2  ; total(wt_dec*mjd^2)
      totobj[cnt:cnt+ncat-1].minmjd = newmeta.mjd
      totobj[cnt:cnt+ncat-1].maxmjd = newmeta.mjd
      totobj[cnt+gdmag].(totind) = cat[gdmag].cmag                         ; sum(mag)
      totobj[cnt+gdmag].(mag2ind) = double(cat[gdmag].cmag)^2              ; sum(mag^2), need dbl to prevent underflow
      objectindex = lindgen(ncat)+cnt
      cnt += ncat

      ; Add new elements to IDSTR
      if idcnt+ncat gt nidstr then begin
        old = idstr
        nnew = 3e5 > ncat
        idstr = replicate(schema_idstr,nidstr+nnew)
        idstr[0:nidstr-1] = old
        nidstr = n_elements(idstr)
        undefine,old
      endif
      ; Add to IDSTR
      sourceid = strtrim(meta.expnum,2)+'.'+strtrim(cat.ccdnum,2)+'.'+strtrim(cat.number,2)
      idstr[idcnt:idcnt+ncat-1].sourceid = sourceid
      idstr[idcnt:idcnt+ncat-1].exposure = meta.base
      idstr[idcnt:idcnt+ncat-1].expnum = strtrim(meta.expnum,2)
      idstr[idcnt:idcnt+ncat-1].objectid = newexp.id
      idstr[idcnt:idcnt+ncat-1].objectindex = objectindex
      idcnt += ncat
    endif

  Endelse  ; second exposure and up
  BOMB:
ENDFOR  ; exposure loop
; No sources
if cnt eq 0 then begin
  print,'No sources in this pixel'
  return
endif
; Trim off the excess elements
obj = obj[0:cnt-1]
totobj = totobj[0:cnt-1]
nobj = n_elements(obj)
print,strtrim(nobj,2),' final objects'
idstr = idstr[0:idcnt-1]

; Make NPHOT from NPHOTX
obj.nphot = obj.nphotu+obj.nphotg+obj.nphotr+obj.nphoti+obj.nphotz+obj.nphoty+obj.nphotvr

; Convert total(mjd*ra) to true proper motion values
;  the slope of RA vs. MJD is
;  pmra=(total(wt*mjd*ra)/total(wt)-<mjd>*<ra>)/(total(wt*mjd^2)/total(wt)-<mjd>^2)
;  we are creating the totals cumulatively as we go
totobj.ra /= obj.raerr        ; wt mean RA (totalrawt/totalwt)
totobj.dec /= obj.decerr      ; wt mean DEC (totaldecwt/totalwt)
obj.mjd /= obj.ndet           ; mean MJD
totobj.ramjd /= obj.raerr     ; wt_ra mean MJD
totobj.decmjd /= obj.decerr   ; wt_dec mean MJD
pmra = (obj.pmra/obj.raerr-totobj.ramjd*totobj.ra)/(totobj.ramjd2/obj.raerr-totobj.ramjd^2)   ; deg[ra]/day
pmra *= (3600*1e3)*365.2425     ; mas/year
pmra *= cos(obj.dec/radeg)      ; mas/year, true angle
pmdec = (obj.pmdec/obj.decerr-totobj.decmjd*totobj.dec)/(totobj.decmjd2/obj.decerr-totobj.decmjd^2)  ; deg/day
pmdec *= (3600*1e3)*365.2425    ; mas/year
; Proper motion errors
; pmerr = 1/sqrt( sum(wt*mjd^2) - <mjd>^2 * sum(wt) )
;   if wt=1/err^2 with err in degrees, but we are using arcsec
;   Need to divide by 3600 for PMDECERR and 3600*cos(dec) for PMRAERR
pmraerr = 1.0/sqrt( totobj.ramjd2 - totobj.ramjd^2 * obj.raerr )
pmraerr /= (3600*cos(totobj.dec/!radeg))    ; correction for raerr in arcsec
pmraerr *= (3600*1e3)*365.2425     ; mas/year
pmraerr *= cos(obj.dec/radeg)      ; mas/year, true angle
pmdecerr = 1.0/sqrt( totobj.decmjd2 - totobj.decmjd^2 * obj.decerr )
pmdecerr /= 3600                   ; correction for decerr in arcsec
pmdecerr *= (3600*1e3)*365.2425    ; mas/year
gdet = where(obj.ndet gt 1,ngdet)
if ngdet gt 0 then begin
  obj[gdet].pmra = pmra[gdet]
  obj[gdet].pmdec = pmdec[gdet]
  obj[gdet].pmraerr = pmraerr[gdet]
  obj[gdet].pmdecerr = pmdecerr[gdet]
endif
bdet = where(obj.ndet lt 2 or finite(pmra) eq 0,nbdet)
; sometimes it happens that the denominator is 0.0 
;  when there are few closely spaced points
;  nothing we can do, just mark as bad
if nbdet gt 0 then begin
  obj[bdet].pmra = 999999.0
  obj[bdet].pmdec = 999999.0
  obj[bdet].pmraerr = 999999.0
  obj[bdet].pmdecerr = 999999.0
endif
obj.deltamjd = totobj.maxmjd-totobj.minmjd
; Average coordinates
obj.ra = totobj.ra   ; now stuff in the average coordinates
obj.dec = totobj.dec
obj.raerr = sqrt(1.0/obj.raerr)    ; err in wt mean RA, arcsec
obj.decerr = sqrt(1.0/obj.decerr)  ; err in wt mean DEC, arcsec

; Convert totalwt and totalfluxwt to MAG and ERR
;  and average the morphology parameters PER FILTER
filters = ['u','g','r','i','z','y','vr']
nfilters = n_elements(filters)
for i=0,nfilters-1 do begin
  ; NDETX is good "detection" and morphology for this filter
  ; NPHOTX is good "photometry" for this filter
  detind = where(tags eq 'NDET'+strupcase(filters[i]),ndetind)
  magind = where(tags eq strupcase(filters[i])+'MAG',nmagind)
  errind = where(tags eq strupcase(filters[i])+'ERR',nerrind)
  rmsind = where(tags eq strupcase(filters[i])+'RMS',nrmsind)
  totind = where(tottags eq strupcase(filters[i])+'TOT',ntotind)
  mag2ind = where(tottags eq strupcase(filters[i])+'MAG2',nmag2ind)
  detphind = where(tags eq 'NPHOT'+strupcase(filters[i]),nphotind)
  asemiind = where(tags eq strupcase(filters[i])+'ASEMI',nasemiind)
  bsemiind = where(tags eq strupcase(filters[i])+'BSEMI',nbsemiind)
  thetaind = where(tags eq strupcase(filters[i])+'THETA',nthetaind)
  
  newflux = obj.(magind) / obj.(errind)
  newmag = 2.50*alog10(newflux)
  newerr = sqrt(1.0/obj.(errind))
  obj.(magind) = newmag
  obj.(errind) = newerr
  bdmag = where(finite(newmag) eq 0,nbdmag)
  if nbdmag gt 0 then begin
    obj[bdmag].(magind) = 99.99
    obj[bdmag].(errind) = 9.99
  endif

  ; Calculate RMS scatter
  ;  RMS^2 * N = sum(mag^2) - 2*<mag>*sum(mag) + N*<mag>^2
  ;   where <mag> is a weighted average
  ;  RMS = sqrt( sum(mag^2)/N - 2*<mag>*sum(mag)/N + <mag>^2 )
  ;  sum(mag^2) is in the MAG2 column and sum(mag) is in TOT
  rms = fltarr(nobj)
  gdrms = where(obj.(detphind) gt 1,ngdrms,comp=bdrms,ncomp=nbdrms)
  if ngdrms gt 0 then $
    rms[gdrms] = sqrt( totobj[gdrms].(mag2ind)/obj[gdrms].(detphind) - $
                       2*obj[gdrms].(magind)*totobj[gdrms].(totind)/obj[gdrms].(detphind) + double(obj[gdrms].(magind))^2 )
  if nbdrms gt 0 then rms[bdrms] = 999999.
  obj.(rmsind) = rms

  ; Average the morphology parameters PER FILTER
  gdet = where(obj.(detind) gt 0,ngdet,comp=bdet,ncomp=nbdet)
  if ngdet gt 0 then begin
    obj[gdet].(asemiind) /= obj[gdet].(detind)
    obj[gdet].(bsemiind) /= obj[gdet].(detind)
    obj[gdet].(thetaind) /= obj[gdet].(detind)
  endif
  if nbdet gt 0 then begin
    obj[bdet].(asemiind) = 99.99
    obj[bdet].(bsemiind) = 99.99
    obj[bdet].(thetaind) = 99.99
  endif
endfor

; Average the morphology parameters, Need a separate counter for that maybe?
;mtags = ['x2','y2','xy','cxx','cyy','cxy','asemi','bsemi','theta','elongation','ellipticity','fwhm','class_star']
mtags = ['x2','y2','xy','asemi','bsemi','theta','ellipticity','fwhm','class_star']
nmtags = n_elements(mtags)
gdet = where(obj.ndet gt 0,ngdet,comp=bdet,ncomp=nbdet)
for i=0,nmtags-1 do begin
  ind = where(tags eq strupcase(mtags[i]),nind)
  ; Divide by the number of detections
  if ngdet gt 0 then obj[gdet].(ind) /= obj[gdet].ndet
  if nbdet gt 0 then obj[bdet].(ind) = 99.99   ; no good detections
endfor

; Get the average error
;metags = ['x2err','y2err','xyerr','cxxerr','cyyerr','cxyerr','asemierr','bsemierr','thetaerr']
metags = ['x2err','y2err','xyerr','asemierr','bsemierr','thetaerr']
nmetags = n_elements(metags)
for i=0,nmetags-1 do begin
  ind = where(tags eq strupcase(metags[i]),nind)
  ; Just take the sqrt to complete the addition in quadrature
  if ngdet gt 0 then obj[gdet].(ind) = sqrt(obj[gdet].(ind))
  if nbdet gt 0 then obj[bdet].(ind) = 99.99
endfor

; Add E(B-V)
print,'Getting E(B-V)'
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
obj.ebv = DUST_GETVAL(glon,glat,/noloop,/interp)

; ONLY INCLUDE OBJECTS WITH AVERAGE RA/DEC
; WITHIN THE BOUNDARY OF THE HEALPIX PIXEL!!!
theta = (90-obj.dec)/radeg
phi = obj.ra/radeg
ANG2PIX_RING,nside,theta,phi,ipring
MATCH,ipring,pix,ind1,ind2,/sort,count=nmatch
if nmatch eq 0 then begin
  print,'None of the final objects fall inside the pixel'
  return
endif
; Get trimmed objects and indices
objtokeep = lonarr(nobj)         ; boolean to keep or trim objects
objtokeep[ind1] = 1
if nmatch lt nobj then begin  ; some to trim
  trimind = lindgen(nobj)
  REMOVE,ind1,trimind
  trimobj = obj[trimind]          ; trimmed objects
endif
newobjindex = lonarr(nobj)-1    ; new indices
newobjindex[ind1] = lindgen(nmatch)
; Keep the objects inside the Healpix
obj = obj[ind1]
print,strtrim(nmatch,2),' final objects fall inside the pixel'

; Remove trimmed objects from IDSTR
totrim = where(objtokeep[idstr.objectindex] eq 0,ntotrim)  ;using old index
if ntotrim gt 0 then begin
  ; Trim objects
  REMOVE,totrim,idstr
  ; Update IDSTR.objectindex
  old_idstr_objectindex = idstr.objectindex
  idstr.objectindex = newobjindex[old_idstr_objectindex]
endif

; Create final summary structure from ALLMETA
;  get exposures that are in IDSTR
;  sometimes EXPNUM numbers have the leading 0s removed
;  and sometimes not, so turn to LONG to match
uiexpnum = uniq(long(idstr.expnum),sort(long(idstr.expnum)))
uexpnum = long(idstr[uiexpnum].expnum)
nuexpnum = n_elements(uexpnum)
MATCH,long(allmeta.expnum),uexpnum,ind1,ind2,/sort,count=nmatch
sumstr = allmeta[ind1]
add_tag,sumstr,'nobjects',0L,sumstr
add_tag,sumstr,'healpix',long(pix),sumstr
; get number of objects per exposure
expnum = long(idstr.expnum)
siexp = sort(expnum)
expnum = expnum[siexp]
if nuexpnum gt 1 then begin
  brklo = where(expnum ne shift(expnum,1),nbrk)
  brkhi = [brklo[1:nbrk-1]-1,n_elements(expnum)-1]
  numobjexp = brkhi-brklo+1
endif else numobjexp=n_elements(expnum)
MATCH,long(sumstr.expnum),uexpnum,ind1,ind2,/sort,count=nmatch
sumstr[ind1].nobjects = numobjexp

; Write the output file
print,'Writing combined catalog to ',outfile
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
if file_test(outdir+'/'+subdir,/directory) eq 0 then file_mkdir,outdir+'/'+subdir
MWRFITS,sumstr,outfile,/create      ; first, summary table
MWRFITS,obj,outfile,/silent         ; second, catalog
MWRFITS,idstr,outfile,/silent       ; third, ID table
if file_test(outfile+'.gz') eq 1 then file_delete,outfile+'.gz',/allow
spawn,['gzip',outfile],/noshell     ; compress final catalog

dt = systime(1)-t0
print,'dt = ',stringize(dt,ndec=2),' sec.'

if keyword_set(stp) then stop

end
