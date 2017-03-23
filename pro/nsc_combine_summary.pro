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
sumstr = replicate({pix:0L,ra:0.0d0,dec:0.0d0,nexposures:0L,nfilters:0L,filters:'',nobjects:0L,success:0},npix)
sumstr.pix = index.pix
pix2ang_ring,nside,sumstr.pix,theta,phi
sumstr.ra = phi*radeg
sumstr.dec = 90-theta*radeg
for i=0,npix-1 do begin
  if (i+1) mod 5000 eq 0 then print,i+1
  file = dir+'combine/'+strtrim(sumstr[i].pix,2)+'.fits'
  if file_test(file) eq 1 then begin
    meta = MRDFITS(file,1,/silent)
    sumstr[i].nexposures = n_elements(meta)
    filter = meta.filter
    ui = uniq(filter,sort(filter))
    sumstr[i].nfilters = n_elements(ui)
    sumstr[i].filters = strjoin(filter[ui],',')
    hd = headfits(file,exten=2)
    sumstr[i].nobjects = sxpar(hd,'naxis2')
    sumstr[i].success = 1
  endif else begin
    sumstr[i].success = 0
  endelse
endfor
gd = where(sumstr.success eq 1,ngd)
print,strtrim(ngd,2),' Healpix successfully processed'
print,'Writing summary file to ',dir+'combine/nsc_instcal_combine.fits'
MWRFITS,expstr,dir+'combine/nsc_instcal_combine.fits',/create

stop

end
