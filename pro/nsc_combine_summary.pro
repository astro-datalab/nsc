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
sumstr = replicate({file:'',exists:0,pix:0L,ra:0.0d0,dec:0.0d0,nexposures:0L,nfilters:0L,filters:'',nobjects:0L,success:0,size:0LL,mtime:0LL},npix)
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
    meta = MRDFITS(sumstr[i].file,1,/silent)
    if tag_exist(meta,'filter') eq 0 then begin  ; old file
      sumstr[i].success = 0
      goto,BOMB
    endif
    sumstr[i].nexposures = n_elements(meta)
    filter = meta.filter
    ui = uniq(filter,sort(filter))
    sumstr[i].nfilters = n_elements(ui)
    sumstr[i].filters = strjoin(filter[ui],',')
    hd = headfits(sumstr[i].file,exten=2)
    sumstr[i].nobjects = sxpar(hd,'naxis2')
    sumstr[i].success = 1
  endif else begin
    sumstr[i].success = 0
  endelse
  BOMB:
endfor
gd = where(sumstr.success eq 1,ngd)
print,strtrim(ngd,2),' Healpix successfully processed'
print,'Writing summary file to ',dir+'combine/nsc_instcal_combine.fits'
MWRFITS,expstr,dir+'combine/nsc_instcal_combine.fits',/create

stop

end
