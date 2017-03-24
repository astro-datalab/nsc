pro nsc_instcal_combine,pix,nside=nside,redo=redo,stp=stp

t0 = systime(1)

; Combine all the exposures that fall in this healpix pixel
;nside = 256
if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(pix) eq 0 then begin
  print,'Syntax - nsc_instcal_combine,pix,nside=nside,redo=redo,stp=stp'
  return
endif

; Check if output file already exists
outfile = dir+'combine/'+strtrim(pix,2)+'.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS already and /redo not set'
  return
endif

print,'Combining InstCal SExtractor catalogs or Healpix pixel = ',strtrim(pix,2)

; Load the list
listfile = localdir+'dnidever/nsc/instcal/nsc_healpix_list.fits'
;listfile = dir+'combine/nsc_healpix_list.fits'
;listfile = dir+'combine/lists/'+strtrim(pix,2)+'.fits'
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
schema_obj = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,ndet:0L,$
              ndetu:0,nphotu:0,umag:0.0,uerr:0.0,uasemi:0.0,ubsemi:0.0,utheta:0.0,$
              ndetg:0,nphotg:0,gmag:0.0,gerr:0.0,gasemi:0.0,gbsemi:0.0,gtheta:0.0,$
              ndetr:0,nphotr:0,rmag:0.0,rerr:0.0,rasemi:0.0,rbsemi:0.0,rtheta:0.0,$
              ndeti:0,nphoti:0,imag:99.9,ierr:0.0,iasemi:0.0,ibsemi:0.0,itheta:0.0,$
              ndetz:0,nphotz:0,zmag:0.0,zerr:0.0,zasemi:0.0,zbsemi:0.0,ztheta:0.0,$
              ndety:0,nphoty:0,ymag:0.0,yerr:0.0,yasemi:0.0,ybsemi:0.0,ytheta:0.0,$
              ndetvr:0,nphotvr:0,vrmag:0.0,vrerr:0.0,vrasemi:0.0,vrbsemi:0.0,vrtheta:0.0,$
              x2:0.0,x2err:0.0,y2:0.0,y2err:0.0,xy:0.0,xyerr:0.0,cxx:0.0,cxxerr:0.0,$
              cxy:0.0,cxyerr:0.0,cyy:0.0,cyyerr:0.0,asemi:0.0,asemierr:0.0,bsemi:0.0,$
              bsemierr:0.0,theta:0.0,thetaerr:0.0,elongation:0.0,$
              ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
tags = tag_names(schema_obj)
obj = replicate(schema_obj,1e5)
nobj = n_elements(obj)
cnt = 0LL

; Loop over the exposures
undefine,allmeta
FOR i=0,nlist-1 do begin
  print,strtrim(i+1,2),' Loading ',list[i].file

  ; Load the exposure catalog
  cat1 = MRDFITS(list[i].file,1,/silent)
  ncat1 = n_elements(cat1)
  print,'  ',strtrim(ncat1,2),' sources'

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
  ;   Don't use "difference image masking" for pre-V3.5 or so version
  ;   because it had problems.  We don't have PLVER but just use
  ;   the maximum IMAFLAGS_ISO value
  if max(cat1.imaflags_iso) gt 10 then begin
    bdcat = where(cat1.imaflags_iso gt 0 and cat1.imaflags_iso lt 120,nbdcat)
  endif else begin
    bdcat = where(cat1.imaflags_iso gt 0,nbdcat)
  endelse
  if nbdcat gt 0 then begin
    print,'  Removing ',strtrim(nbdcat,2),' sources contaminated by bad pixels.'
    if nbdcat eq ncat1 then goto,BOMB
    REMOVE,bdcat,cat1
    ncat1 = n_elements(cat1)
  endif

  metafile = repstr(list[i].file,'_cat','_meta')
  meta = MRDFITS(metafile,1,/silent)
  meta.base = strtrim(meta.base)
  meta.expnum = strtrim(meta.expnum)
  ;head = headfits(list[i].file,exten=0)
  ;filtername = sxpar(head,'filter')
  ;if strmid(filtername,0,2) eq 'VR' then filter='VR' else filter=strmid(filtername,0,1)
  ;exptime = sxpar(head,'exptime')
  print,'  FILTER=',meta.filter,'  EXPTIME=',stringize(meta.exptime,ndec=1),' sec'

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
  push,allmeta,meta

  ; NDETX is good "detection" and morphology for this filter
  ; NPHOTX is good "photometry" for this filter
  detind = where(tags eq 'NDET'+strupcase(meta.filter),ndetind)
  magind = where(tags eq strupcase(meta.filter)+'MAG',nmagind)
  errind = where(tags eq strupcase(meta.filter)+'ERR',nerrind)
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
    newexp.cxx = cat.cxx_world
    newexp.cxxerr = cat.errcxx_world^2
    newexp.cyy = cat.cyy_world
    newexp.cyyerr = cat.errcyy_world^2
    newexp.cxy = cat.cxy_world
    newexp.cxyerr = cat.errcxy_world^2
    newexp.asemi = cat.a_world
    newexp.asemierr = cat.erra_world^2
    newexp.bsemi = cat.b_world
    newexp.bsemierr = cat.errb_world^2
    newexp.theta = cat.theta_world
    newexp.thetaerr = cat.errtheta_world^2
    newexp.elongation = cat.elongation
    newexp.ellipticity = cat.ellipticity
    newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
    newexp.flags = cat.flags
    newexp.class_star = cat.class_star
    obj[0:ncat-1] = newexp
    cnt += ncat

    ; Add new elements to IDSTR
    if idcnt+ncat gt nidstr then begin
      old = idstr
      idstr = replicate(schema_idstr,nidstr+1e5)
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
      cmb.ndet++
      newcat = cat[ind2]
      ; Detection and morphology parameters for this FILTER
      cmb.(detind)++
      cmb.(asemiind) += newcat.a_world
      cmb.(bsemiind) += newcat.b_world
      cmb.(thetaind) += newcat.theta_world
      ; Good photometry for this FILTER
      gdmag = where(newcat.cmag lt 50,ngdmag)
      if ngdmag gt 0 then begin
        cmb[gdmag].(magind) = 2.5118864d^newcat[gdmag].cmag * (1.0d0/newcat[gdmag].cerr^2)
        cmb[gdmag].(errind) = 1.0d0/newcat[gdmag].cerr^2
        cmb[gdmag].(detphind) += 1
        ; NPHOTX means good PHOT detection
      endif
      cmb.x2 += newcat.x2_world
      cmb.x2err += newcat.errx2_world^2
      cmb.y2 += newcat.y2_world
      cmb.y2err += newcat.erry2_world^2
      cmb.xy += newcat.xy_world
      cmb.xyerr += newcat.errxy_world^2
      cmb.cxx += newcat.cxx_world
      cmb.cxxerr += newcat.errcxx_world^2
      cmb.cyy += newcat.cyy_world
      cmb.cyyerr += newcat.errcyy_world^2
      cmb.cxy += newcat.cxy_world
      cmb.cxyerr += newcat.errcxy_world^2
      cmb.asemi += newcat.a_world
      cmb.asemierr += newcat.erra_world^2
      cmb.bsemi += newcat.b_world
      cmb.bsemierr += newcat.errb_world^2
      cmb.theta += newcat.theta_world
      cmb.thetaerr += newcat.errtheta_world^2
      cmb.elongation += newcat.elongation
      cmb.ellipticity += newcat.ellipticity
      cmb.fwhm += newcat.fwhm_world*3600  ; in arcsec
      cmb.flags OR= newcat.flags
      cmb.class_star += newcat.class_star
      obj[ind1] = cmb  ; stuff it back in

      ; Add new elements to IDSTR
      if idcnt+nmatch gt nidstr then begin
        old = idstr
        idstr = replicate(schema_idstr,nidstr+1e5)
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
        obj = replicate(schema_obj,nobj+1e5)
        obj[0:nobj-1] = old
        nobj = n_elements(obj)
        undefine,old
      endif

      ; Copy to final structure
      newexp = replicate(schema_obj,ncat)
      newexp.id = strtrim(pix,2)+'.'+strtrim(cnt+lindgen(ncat)+1,2)
      newexp.pix = pix
      newexp.ra = cat.ra
      newexp.dec = cat.dec
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
      newexp.cxx = cat.cxx_world
      newexp.cxxerr = cat.errcxx_world^2
      newexp.cyy = cat.cyy_world
      newexp.cyyerr = cat.errcyy_world^2
      newexp.cxy = cat.cxy_world
      newexp.cxyerr = cat.errcxy_world^2
      newexp.asemi = cat.a_world
      newexp.asemierr = cat.erra_world^2
      newexp.bsemi = cat.b_world
      newexp.bsemierr = cat.errb_world^2
      newexp.theta = cat.theta_world
      newexp.thetaerr = cat.errtheta_world^2
      newexp.elongation = cat.elongation
      newexp.ellipticity = cat.ellipticity
      newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
      newexp.flags = cat.flags
      newexp.class_star = cat.class_star
      obj[cnt:cnt+ncat-1] = newexp   ; stuff it in
      objectindex = lindgen(ncat)+cnt
      cnt += ncat

      ; Add new elements to IDSTR
      if idcnt+ncat gt nidstr then begin
        old = idstr
        idstr = replicate(schema_idstr,nidstr+1e5)
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
nobj = n_elements(obj)
print,strtrim(nobj,2),' final objects'
idstr = idstr[0:idcnt-1]

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
mtags = ['x2','y2','xy','cxx','cyy','cxy','asemi','bsemi','theta','elongation','ellipticity','fwhm','class_star']
nmtags = n_elements(mtags)
gdet = where(obj.ndet gt 0,ngdet,comp=bdet,ncomp=nbdet)
for i=0,nmtags-1 do begin
  ind = where(tags eq strupcase(mtags[i]),nind)
  ; Divide by the number of detections
  if ngdet gt 0 then obj[gdet].(ind) /= obj[gdet].ndet
  if nbdet gt 0 then obj[bdet].(ind) = 99.99   ; no good detections
endfor

; Get the average error
metags = ['x2err','y2err','xyerr','cxxerr','cyyerr','cxyerr','asemierr','bsemierr','thetaerr']
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
trimind = lindgen(nobj)
REMOVE,ind1,trimind
trimobj = obj[trimind]          ; trimmed objects
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
MWRFITS,sumstr,outfile,/create      ; first, summary table
MWRFITS,obj,outfile,/silent         ; second, catalog
MWRFITS,idstr,outfile,/silent       ; third, ID table

dt = systime(1)-t0
print,'dt = ',stringize(dt,ndec=2),' sec.'

if keyword_set(stp) then stop

end
