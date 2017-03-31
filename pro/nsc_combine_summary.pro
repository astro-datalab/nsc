pro nsc_combine_summary

; Create a summary file of all the Healpix combined files

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
nside = 128
radeg = 180.0d0 / !dpi

index = mrdfits(dir+'combine/nsc_healpix_list.fits',2,/silent)
npix = n_elements(index)

; Load all the summary/metadata files
print,'Creating Healpix summary file'
; u, g, r, i, z, Y, VR
sumstr = replicate({file:'',exists:0,pix:0L,ra:0.0d0,dec:0.0d0,nexposures:0L,nfilters:0L,filters:'',nobjects:0L,nsources:-1L,nobjfilt:lonarr(7),$
                    exptimefilt:fltarr(7),depthfilt:fltarr(7)+99.99,success:0,size:0LL,mtime:0LL},npix)
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
        if nfiltind gt 0 then sumstr[i].exptimefilt[j]=total(meta[filtind].exptime)
        ; Check photometry
        magind = where(tags eq strupcase(filters[j])+'MAG',nmagind)
        errind = where(tags eq strupcase(filters[j])+'ERR',nerrind)
        ; Nobjects per filter
        ind = where(cat.(magind) lt 50,nind)
        sumstr[i].nobjfilt[j] = nind
        ; Depth per filter
        ;  S/N = 1.087/err
        ;  so S/N=5 is for err=1.087/5=0.2174
        ;  S/N=10 is for err=1.087/10=0.1087
        ind2 = where(cat.(magind) lt 50 and cat.(errind) ge 0.0987 and cat.(errind) le 0.1187,nind2)
        if nind2 lt 5 then ind2 = where(cat.(magind) lt 50 and cat.(errind) ge 0.0787 and cat.(errind) le 0.1387,nind2)
        if nind2 gt 5 then begin
          depth = median([cat[ind2].(magind)])
          sumstr[i].depthfilt[j] = depth
        endif
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
print,'Writing summary file to ',dir+'combine/nsc_instcal_combine.fits'
MWRFITS,sumstr,dir+'combine/nsc_instcal_combine.fits',/create

stop

end
