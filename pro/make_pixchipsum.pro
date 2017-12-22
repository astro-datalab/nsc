pro make_pixchipsum,pix,redo=redo

; Make small summary table of the number of chip sources in each pixel
  
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
cmbdir = dir +'combine/'

if n_elements(pix) eq 0 then begin
  print,'syntax - make_pixchipsum,pix'
  return
endif

pixfile = cmbdir+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits.gz'
if file_test(pixfile) eq 0 then begin
  print,pixfile,' NOT FOUND'
  return
endif

outfile = cmbdir+'chipsum/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'_chipsum.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS and /redo NOT SET'
  return
endif

; Load the IDSTR structure
idstr = MRDFITS(pixfile,3,/silent)
idstr.exposure = strtrim(idstr.exposure,2)
idstr.sourceid = strtrim(idstr.sourceid,2)
nidstr = n_elements(idstr)
uiidstr = uniq(idstr.sourceid,sort(idstr.sourceid))
; Get unique exposure+ccdnum IDs and add up the sources
dum = strsplitter(idstr.sourceid,'.',/extract)
ccdnums = reform(dum[2,*])
chids = idstr.exposure+'-'+ccdnums
uichids = uniq(chids,sort(chids))
uchids = chids[uichids]
nuchids = n_elements(uchids)
si = sort(chids)
; find the breaks
brklo = where(chids[si] ne shift(chids[si],1),nbrk)
if nbrk gt 0 then begin
  brkhi = [brklo[1:nbrk-1]-1,n_elements(chids)-1]
  nchsrc = brkhi-brklo+1
endif else begin
  brklo = 0
  brkhi = n_elements(chids)-1
  nchsrc = n_elements(chids)
endelse
uidstrindex = {id:uchids,exposure:strarr(nuchids),ccdnum:lonarr(nuchids),lo:brklo,hi:brkhi,num:nchsrc,index:si}
dum2 = strsplitter(uchids,'-',/extract)
uchids_exposure_all = reform(dum2[0,*])  ; all
uidstrindex.exposure = reform(dum2[0,*])
uidstrindex.ccdnum = long(reform(dum2[1,*]))
; Make the summary structure
pixchipsum = replicate({exposure:'',ccdnum:0L,nsources:0L},nuchids)
pixchipsum.exposure = uidstrindex.exposure
pixchipsum.ccdnum = uidstrindex.ccdnum
pixchipsum.nsources = uidstrindex.num
; Save the file
print,'Writing to ',outfile
if file_test(cmbdir+'chipsum/'+strtrim(long(pix)/1000,2)+'/',/directory) eq 0 then file_mkdir,cmbdir+'chipsum/'+strtrim(long(pix)/1000,2)+'/'
MWRFITS,pixchipsum,outfile,/create

;stop

end
