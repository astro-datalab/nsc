pro nsc_combine_getgals

; Make galaxies catalogs for each healpix

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
nside = 128
radeg = 180.0d0 / !dpi

index = mrdfits(dir+'combine/nsc_healpix_list.fits',2,/silent)
npix = n_elements(index)

sumstr = mrdfits(dir+'combine/gal/nsc_combine_gal.fits',1,/silent)
sumstr.galfile = strtrim(sumstr.galfile,2)
goto,combine

; Load all the summary/metadata files
print,'Creating Healpix summary file'
sumstr = replicate({file:'',galfile:'',exists:0,pix:0L,ra:0.0d0,dec:0.0d0,nexposures:0L,nfilters:0L,filters:'',nobjects:0L,ngals:0L,success:0,size:0LL,mtime:0LL},npix)
sumstr.pix = index.pix
PIX2ANG_RING,nside,sumstr.pix,theta,phi
sumstr.ra = phi*radeg
sumstr.dec = 90-theta*radeg
files = dir+'combine/'+strtrim(sumstr.pix,2)+'.fits'
info = file_info(files)
sumstr.file = info.name
sumstr.exists = info.exists
sumstr.size = info.size
sumstr.mtime = info.mtime
gdexist = where(sumstr.exists eq 1,ngdexist)
print,strtrim(ngdexist,2),'/',strtrim(npix,2),' exist'
for i=0,npix-1 do begin
  if (i+1) mod 100 eq 0 then print,i+1
  if sumstr[i].exists eq 1 then begin
    hd = headfits(sumstr[i].file,exten=2,errmsg=errmsg)
    if errmsg eq '' then begin  ; check that we can read it
      meta = MRDFITS(sumstr[i].file,1,/silent,status=status)
      if status lt 0 then goto,BOMB 
     if tag_exist(meta,'filter') eq 0 then begin  ; old file
        sumstr[i].success = 0
        goto,BOMB
      endif
      sumstr[i].nexposures = n_elements(meta)
      filter = meta.filter
      ui = uniq(filter,sort(filter))
      sumstr[i].nfilters = n_elements(ui)
      sumstr[i].filters = strjoin(filter[ui],',')
      sumstr[i].nobjects = sxpar(hd,'naxis2')
      sumstr[i].success = 1
      ; Get the catalog
      cat = MRDFITS(sumstr[i].file,2,/silent,status=status2)
      gdgal = where(cat.fwhm gt 3.0 and cat.class_star lt 0.1,ngdgal)
      ;print,'  ',strtrim(ngdgal,2),' galaxies'
      if ngdgal gt 0 then begin
        gal = cat[gdgal]
        outfile = file_dirname(sumstr[i].file)+'/gal/'+strtrim(sumstr[i].pix,2)+'_gal.fits'
        MWRFITS,gal,outfile,/create
        sumstr[i].ngals = ngdgal
      endif
    endif
  endif
  BOMB:
endfor
gd = where(sumstr.success eq 1,ngd)
print,strtrim(ngd,2),' Healpix successfully processed'
print,'Writing summary file to ',dir+'combine/gal/nsc_combine_gal.fits'
MWRFITS,sumstr,dir+'combine/gal/nsc_combine_gal.fits',/create

stop

; Now sum them all up, keep on what's necessary
combine:

schema_gal = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,filter:'',mag:99.0,fwhm:0.0,asemi:0.0,bsemi:0.0,class_star:0.0,ellipticity:0.0}
allgal = replicate(schema_gal,1e7)
nallgal = n_elements(allgal)
cnt = 0LL

; id, ra, dec, fiducial_mag, fwhm, class_star, ellipticity
filters = ['r','i','g','z','Y','VR','u']
nfilters = n_elements(filters)
magind = lonarr(nfilters)
for i=0,npix-1 do begin
  if (i+1) mod 100 eq 0 then print,i+1,cnt
  if sumstr[i].ngals gt 0 then begin
    outfile = file_dirname(sumstr[i].file)+'/gal/'+strtrim(sumstr[i].pix,2)+'_gal.fits'
    sumstr[i].galfile = outfile
    gal = mrdfits(sumstr[i].galfile,1,/silent)
    ngal = n_elements(gal)
    gal.id = strtrim(gal.id,2)

    ; Get tags and magnitude indices
    tags = tag_names(gal)
    for j=0,nfilters-1 do magind[j]=where(tags eq strupcase(filters[j])+'MAG') 

    ; Put in new structure
    newgal = replicate(schema_gal,ngal)
    struct_assign,gal,newgal,/nozero

    ; Get fiducial magnitude
    nhavemag = 0
    for j=0,nfilters-1 do begin
      if nhavemag lt ngal then begin
        ind = where(gal.(magind[j]) lt 50 and newgal.mag gt 50,nind)
        if nind gt 0 then begin
          newgal[ind].mag = gal[ind].(magind[j])
          newgal[ind].filter = filters[j]
          nhavemag += nind
        endif
      endif
    endfor

    ; Add new elements
    if cnt+ngal gt nallgal then begin
      orig = allgal
      allgal = replicate(schema_gal,nallgal+5e6)
      allgal[0:nallgal-1] = orig
      nallgal = n_elements(allgal)
    endif

    ; Add to the final structure
    allgal[cnt:cnt+ngal-1] = newgal
    cnt += ngal
  endif
endfor
; Trim allgal
allgal = allgal[0:cnt-1]

; save the output
print,'Writing combined catalog to ',dir+'combine/gal/nsc_combine_gal_allcat.fits'
MWRFITS,allgal,dir+'combine/gal/nsc_combine_gal_allcat.fits',/create

stop

end
