pro nsc_instcal_combine_exposures_needed,allpix,expdir,list=list

; The exposures needed to run certain healpix pixels

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
radeg = 180.0d0 / !dpi

undefine,expdir,list

nallpix = n_elements(allpix)
if nallpix eq 0 then begin
  print,'Syntax - nsc_instcal_combine_exposues_needed,allpix,expdir'
  return
endif

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

; LOOP OVER THE PIXELS
undefine,list
FOR p=0L,nallpix-1 do begin
  pix = allpix[p]

  ; Find our pixel
  ind = where(index.pix eq pix,nind)
  if nind eq 0 then begin
    print,'No entries for Healpix pixel "',strtrim(pix,2),'" in the list'
    goto,bomb
  endif
  ind = ind[0]
  list1 = healstr[index[ind].lo:index[ind].hi]
  nlist1 = n_elements(list1)

  ; GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
  ;  so we can deal with the edge cases
  NEIGHBOURS_RING,nside,pix,neipix,nneipix
  for i=0,nneipix-1 do begin
    ind = where(index.pix eq neipix[i],nind)
    if nind gt 0 then begin
      ind = ind[0]
      neilist = healstr[index[ind].lo:index[ind].hi]
      push,list1,neilist
    endif
  endfor
  ; Get unique values
  ui = uniq(list1.file,sort(list1.file))
  list1 = list1[ui]
  nlist1 = n_elements(list1)
  print,strtrim(p+1,2),'/',strtrim(nallpix,2),'  ',strtrim(pix,2),' ',strtrim(nlist1,2),' exposures that overlap this pixel and neighbors'

  ; Add to final list
  push,list,list1

  BOMB:
ENDFOR  ; pixel loop

; Get unique values
ui = uniq(list.file,sort(list.file))
list = list[ui]
nlist = n_elements(list)
print,strtrim(nlist,2),' exposures needed to process these healpix'
expdir = file_dirname(list.file)

; Fix the directories if needed
;  /net/dl1
if strmid(dldir,0,4) eq '/net' then begin
  bd = where(strmid(expdir,0,4) ne '/net',nbd)
  if nbd gt 0 then expdir[bd]='/net'+expdir[bd]
;  /dl1
endif else begin
  bd = where(strmid(expdir,0,4) eq '/net',nbd)
  if nbd gt 0 then expdir[bd]=strmid(expdir[bd],4)
endelse

;stop

end
