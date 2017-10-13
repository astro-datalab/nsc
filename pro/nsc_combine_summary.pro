pro nsc_combine_summary,version=version

; Create a summary file of all the Healpix combined files

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
nside = 128
radeg = 180.0d0 / !dpi

index = mrdfits(dir+'lists/nsc_healpix_list.fits',2,/silent)
npix = n_elements(index)

; Load all the summary/metadata files
print,'Creating Healpix summary file'
; u, g, r, i, z, Y, VR
sumstr = replicate({file:'',exists:0,pix:0L,ra:0.0d0,dec:0.0d0,nexposures:0L,nfilters:0L,filters:'',nobjects:0L,nsources:-1L,nobjfilt:lonarr(7),$
                    nexpfilt:lonarr(7),exptimefilt:fltarr(7),depth95filt:fltarr(7)+99.99,depth10sigfilt:fltarr(7)+99.99,success:0,size:0LL,mtime:0LL},npix)
sumstr.pix = index.pix
PIX2ANG_RING,nside,sumstr.pix,theta,phi
sumstr.ra = phi*radeg
sumstr.dec = 90-theta*radeg
files = dir+'combine/'+strtrim(sumstr.pix/1000,2)+'/'+strtrim(sumstr.pix,2)+'.fits.gz'
info = file_info(files)
sumstr.file = info.name
sumstr.exists = info.exists
sumstr.size = info.size
sumstr.mtime = info.mtime
gdexist = where(sumstr.exists eq 1,ngdexist)
print,strtrim(ngdexist,2),'/',strtrim(npix,2),' exist'
filters = ['u','g','r','i','z','Y','VR']
nfilters = n_elements(filters)
for i=0,npix-1 do begin
  if (i+1) mod 100 eq 0 then print,i+1
  if sumstr[i].exists eq 1 then begin
    hd2 = headfits(sumstr[i].file,exten=2,errmsg=errmsg2)
    if errmsg2 eq '' then begin  ; check that we can read it
      ; Load the summary catalog
      meta = MRDFITS(sumstr[i].file,1,/silent,status=status1)
      if status1 lt 0 then goto,BOMB 
      if tag_exist(meta,'filter') eq 0 then begin  ; old file
        sumstr[i].success = 0
        goto,BOMB
      endif
      sumstr[i].nexposures = n_elements(meta)
      filter = meta.filter
      ui = uniq(filter,sort(filter))
      sumstr[i].nfilters = n_elements(ui)
      sumstr[i].filters = strjoin(filter[ui],',')
      sumstr[i].nobjects = sxpar(hd2,'naxis2')
      ; Load the catalog
      cat = MRDFITS(sumstr[i].file,2,/silent,status=status2)
      tags = tag_names(cat)
      for j=0,nfilters-1 do begin
        ; Total exposure time in this filter
        filtind = where(meta.filter eq filters[j],nfiltind)
        sumstr[i].nexpfilt[j] = nfiltind
        if nfiltind gt 0 then sumstr[i].exptimefilt[j]=total(meta[filtind].exptime)
        ; Check photometry
        magind = where(tags eq strupcase(filters[j])+'MAG',nmagind)
        errind = where(tags eq strupcase(filters[j])+'ERR',nerrind)
        ; Nobjects per filter
        ind = where(cat.(magind) lt 50,nind)
        sumstr[i].nobjfilt[j] = nind
        ; Depth per filter
        ; 95th percentile
        cmag = cat[ind].(magind)
        si = sort(cmag)
        cmag = cmag[si]
        depth95 = cmag[round(0.95*nind)-1]
        sumstr[i].depth95filt[j] = depth95
        ; 10 sigma
        ;  S/N = 1.087/err
        ;  so S/N=5 is for err=1.087/5=0.2174
        ;  S/N=10 is for err=1.087/10=0.1087
        depth10sig = 99.99
        depind = where(cat.(magind) lt 50 and cat.(magind) gt depth95-3.0 and cat.(errind) ge 0.0987 and cat.(errind) le 0.1187,ndepind)
        if ndepind lt 5 then depind = where(cat.(magind) lt 50 and cat.(magind) gt depth95-3.0 and cat.(errind) ge 0.0787 and cat.(errind) le 0.1387,ndepind)
        if ndepind gt 5 then begin
          depth10sig = median([cat[depind].(magind)])
        endif else begin
          depind = where(cat.(magind) lt 50,ndepind)
          if ndepind gt 0 then depth10sig=max([cat[depind].(magind)])
        endelse
        sumstr[i].depth10sigfilt[j] = depth10sig
      endfor
      ; Load the source table header
      hd3 = headfits(sumstr[i].file,exten=3,errmsg=errmsg3)
      if errmsg3 eq '' then sumstr[i].nsources = sxpar(hd3,'naxis2')
      sumstr[i].success = 1
    endif
  endif
  BOMB:
endfor
gd = where(sumstr.success eq 1,ngd)
print,strtrim(ngd,2),' Healpix successfully processed'
print,'Writing summary file to ',dir+'lists/nsc_instcal_combine.fits'
MWRFITS,sumstr,dir+'lists/nsc_instcal_combine.fits',/create

stop

end
