pro nsc_instcal_combine,pix,nside=nside,redo=redo,stp=stp

; Combine all the exposures that fall in this healpix pixel
;nside = 256
if n_elements(nside) eq 0 then nside = 128
dir = '/datalab/users/dnidever/decamcatalog/instcal/'
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

; Load the list
listfile = dir+'combine/healpix_list.fits'
;listfile = dir+'combine/lists/'+strtrim(pix,2)+'.fits'
if file_test(listfile) eq 0 then begin
  print,listfile,' NOT FOUND'
  return
endif
healstr = MRDFITS(listfile,1,/silent)
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

print,'Combining InstCal SExtractor catalogs or Healpix pixel = ',strtrim(pix,2)
print,strtrim(nlist,2),' overlapping exposure(s)'

; Initialize the object structure
schema_obj = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,ndet:0L,umag:0.0,uerr:0.0,ndetu:0,gmag:0.0,$
              gerr:0.0,ndetg:0,rmag:0.0,rerr:0.0,ndetr:0,imag:99.9,ierr:0.0,ndeti:0,$
              zmag:0.0,zerr:0.0,ndetz:0,ymag:0.0,yerr:0.0,ndety:0,vrmag:0.0,vrerr:0.0,$
              ndetvr:0,x2:0.0,y2:0.0,xy:0.0,cxx:0.0,cxy:0.0,cyy:0.0,asemi:0.0,bsemi:0.0,theta:0.0,$
              elongation:0.0,ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
tags = tag_names(schema_obj)
obj = replicate(schema_obj,1e5)
nobj = n_elements(obj)
cnt = 0LL

; Loop over the exposures
for i=0,nlist-1 do begin
  print,strtrim(i+1,2),' Loading ',list[i].file

  ; Load the exposure catalog
  cat1 = MRDFITS(list[i].file,1,/silent)
  ncat1 = n_elements(cat1)
  print,'  ',strtrim(ncat1,2),' sources'

  ; Remove sources near bad pixels
  bdcat = where(cat1.imaflags_iso gt 10,nbdcat)
  if nbadcat gt 0 then begin
    print,'  Removing ',strtrim(nbdcat,2),' sources contaminated by bad pixels.'
    if nbdcat eq ncat1 then goto,BOMB
    REMOVE,bdcat,cat1
    ncat1 = n_elements(cat1)
  endif

  metafile = repstr(list[i].file,'_cat','_meta')
  meta = MRDFITS(metafile,1,/silent)
  ;head = headfits(list[i].file,exten=0)
  ;filtername = sxpar(head,'filter')
  ;if strmid(filtername,0,2) eq 'VR' then filter='VR' else filter=strmid(filtername,0,1)
  ;exptime = sxpar(head,'exptime')
  print,'  FILTER=',meta.filter,'  EXPTIME=',stringize(meta.exptime,ndec=1),' sec'

  ; Only want source inside this pixel
  theta = (90-cat1.dec)/radeg
  phi = cat1.ra/radeg
  ANG2PIX_RING,nside,theta,phi,ipring
  MATCH,ipring,pix,ind1,ind2,/sort,count=nmatch
  if nmatch eq 0 then begin
    print,'  No sources inside this pixel'
    goto,BOMB
  endif
  print,'  ',strtrim(nmatch,2),' sources are inside this pixel'
  cat = cat1[ind1]
  ncat = nmatch

  magind = where(tags eq strupcase(meta.filter)+'MAG',nmagind)
  errind = where(tags eq strupcase(meta.filter)+'ERR',nerrind)
  detind = where(tags eq 'NDET'+strupcase(meta.filter),ndetind)

  ; Combine the data
  ;-----------------
  ; First catalog
  If cnt eq 0 then begin

    ; Copy to final structure
    newexp = replicate(schema_obj,ncat)
    newexp.id = lindgen(ncat)+1
    newexp.pix = pix
    newexp.ra = cat.ra
    newexp.dec = cat.dec
    newexp.ndet = 1
    gdmag = where(cat.cmag lt 50,ngdmag)
    ;print,ngdmag,' good phot'
    if ngdmag gt 0 then begin
      newexp[gdmag].(magind) = 2.5118864d^cat[gdmag].cmag * (1.0d0/cat[gdmag].cerr^2)
      newexp[gdmag].(errind) = 1.0d0/cat[gdmag].cerr^2
      newexp[gdmag].(detind) = 1
    endif
    newexp.x2 = cat.x2_world
    newexp.y2 = cat.y2_world
    newexp.xy = cat.xy_world
    newexp.cxx = cat.cxx_world
    newexp.cyy = cat.cyy_world
    newexp.cxy = cat.cxy_world
    newexp.asemi = cat.a_world
    newexp.bsemi = cat.b_world
    newexp.theta = cat.theta_world
    newexp.elongation = cat.elongation
    newexp.ellipticity = cat.ellipticity
    newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
    newexp.flags = cat.flags
    newexp.class_star = cat.class_star
    obj[0:ncat-1] = newexp
    cnt += ncat

  ; Second and up
  Endif else begin

    ; Match new sources to the objects
    SRCMATCH,obj[0:cnt-1].ra,obj[0:cnt-1].dec,cat.ra,cat.dec,0.5,ind1,ind2,count=nmatch,/sph,/usehist  ; use faster histogram_nd method
    print,'  ',strtrim(nmatch,2),' matched sources'
    ; Some matches, add data to existing record for these sources
    if nmatch gt 0 then begin

      ; When CMAG=99.99 the morphology parameters are still okay

      ; Combine the information
      cmb = obj[ind1]
      cmb.ndet++
      newcat = cat[ind2]
      gdmag = where(newcat.cmag lt 50,ngdmag)
      ;print,ngdmag,' good phot'
      if ngdmag gt 0 then begin
        cmb[gdmag].(magind) = 2.5118864d^newcat[gdmag].cmag * (1.0d0/newcat[gdmag].cerr^2)
        cmb[gdmag].(errind) = 1.0d0/newcat[gdmag].cerr^2
        cmb[gdmag].(detind) += 1
        ; NDETX means good PHOT detection
      endif
      cmb.x2 += newcat.x2_world
      cmb.y2 += newcat.y2_world
      cmb.xy += newcat.xy_world
      cmb.cxx += newcat.cxx_world
      cmb.cyy += newcat.cyy_world
      cmb.cxy += newcat.cxy_world
      cmb.asemi += newcat.a_world
      cmb.bsemi += newcat.b_world
      cmb.theta += newcat.theta_world
      cmb.elongation += newcat.elongation
      cmb.ellipticity += newcat.ellipticity
      cmb.fwhm += newcat.fwhm_world*3600  ; in arcsec
      cmb.flags OR= newcat.flags
      cmb.class_star += newcat.class_star
      obj[ind1] = cmb  ; stuff it back in

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
      newexp.id = cnt+lindgen(ncat)+1
      newexp.pix = pix
      newexp.ra = cat.ra
      newexp.dec = cat.dec
      newexp.ndet = 1
      gdmag = where(cat.cmag lt 50,ngdmag)
      ;print,ngdmag,' good phot'
      if ngdmag gt 0 then begin
        newexp[gdmag].(magind) = 2.5118864d^cat[gdmag].cmag * (1.0d0/cat[gdmag].cerr^2)
        newexp[gdmag].(errind) = 1.0d0/cat[gdmag].cerr^2
        newexp[gdmag].(detind) = 1
      endif
      newexp.x2 = cat.x2_world
      newexp.y2 = cat.y2_world
      newexp.xy = cat.xy_world
      newexp.cxx = cat.cxx_world
      newexp.cyy = cat.cyy_world
      newexp.cxy = cat.cxy_world
      newexp.asemi = cat.a_world
      newexp.bsemi = cat.b_world
      newexp.theta = cat.theta_world
      newexp.elongation = cat.elongation
      newexp.ellipticity = cat.ellipticity
      newexp.fwhm = cat.fwhm_world*3600  ; in arcsec
      newexp.flags = cat.flags
      newexp.class_star = cat.class_star
      obj[cnt:cnt+ncat-1] = newexp   ; stuff it in
      cnt += ncat
    endif

  Endelse
  BOMB:
endfor
; No sources
if cnt eq 0 then begin
  print,'No sources in this pixel'
  return
endif
; Trim off the excess elements
obj = obj[0:cnt-1]
nobj = n_elements(obj)
print,strtrim(nobj,2),' final objects'

; Convert totalwt and totalfluxwt to MAG and ERR
filters = ['u','g','r','i','z','y','vr']
nfilters = n_elements(filters)
for i=0,nfilters-1 do begin
  magind = where(tags eq strupcase(filters[i])+'MAG',nmagind)
  errind = where(tags eq strupcase(filters[i])+'ERR',nerrind)
  detind = where(tags eq 'NDET'+strupcase(filters[i]),ndetind)
  
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
endfor

; do we use NDET for the morphology paramters???

; Average the morphology paramters, Need a separate counter for that maybe?
mtags = ['x2','y2','xy','cxx','cyy','cxy','asemi','bsemi','theta','elongation','ellipticity','fwhm','class_star']
nmtags = n_elements(mtags)
for i=0,nmtags-1 do begin
  ind = where(tags eq strupcase(mtags[i]),nind)
  ; Divide by the number of detections
  gdet = where(obj.ndet gt 0,ngdet)
  if ngdet gt 0 then obj[gdet].(ind) /= obj[gdet].ndet
endfor

; Add E(B-V)
print,'Getting E(B-V)'
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
obj.ebv = DUST_GETVAL(glon,glat,/noloop,/interp)

; Write the output file
print,'Writing combined catalog to ',outfile
MWRFITS,obj,outfile,/create

if keyword_set(stp) then stop

end
