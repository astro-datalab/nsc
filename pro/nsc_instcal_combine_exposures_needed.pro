pro nsc_instcal_combine_exposures_needed,pixlist,expdir,list=list

; The exposures needed to run certain healpix pixels

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
radeg = 180.0d0 / !dpi

undefine,expdir,list

npixlist = n_elements(pixlist)
if npixlist eq 0 then begin
  print,'Syntax - nsc_instcal_combine_exposues_needed,pixlist,expdir'
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

; Get all of the neighbors
undefine,neipix
For p=0L,npixlist-1 do begin
  pix = pixlist[p]
  NEIGHBOURS_RING,nside,pix,neipix1,nneipix1
  push,neipix,neipix1
Endfor

; Combine pixels list and get unique ones
allpix = [pixlist,neipix]
ui = uniq(allpix,sort(allpix))
allpix = allpix[ui]
nallpix = n_elements(allpix)

; MATCH to index
MATCH,index.pix,allpix,ind1,ind2,/sort,count=nmatch
if nmatch lt nallpix then print,strtrim(nallpix-nmatch,2),' pixels or neighbors not covered'
allpix = allpix[ind2]
allpix_lo = index[ind1].lo
allpix_hi = index[ind1].hi
nallpix = n_elements(allpix)

; LOOP OVER THE PIXELS
schema = healstr[0]
struct_assign,{dum:''},schema
list = replicate(schema,10*npixlist)
nlist = n_elements(list)
cnt = 0LL
FOR p=0L,nallpix-1 do begin
  pix = allpix[p]

  ;; Find our pixel
  ;ind = where(index.pix eq pix,nind)
  ;if nind eq 0 then begin
  ;  print,'No entries for Healpix pixel "',strtrim(pix,2),'" in the list'
  ;  goto,bomb
  ;endif
  ;ind = ind[0]
  ;list1 = healstr[index[ind].lo:index[ind].hi]
  list1 = healstr[allpix_lo[p]:allpix_hi[p]]  
  nlist1 = n_elements(list1)
  print,strtrim(p+1,2),'/',strtrim(nallpix,2),'  ',strtrim(pix,2),' ',strtrim(nlist1,2),' exposures that overlap this pixel'

  ; Add more elements
  if cnt+nlist1 gt nlist then begin
    orig = list
    list = replicate(schema,nlist+((2*npixlist>1000)>nlist1))
    list[0:nlist-1] = orig
    nlist = n_elements(list)
    undefine,orig
  endif

  ; Add to final list
  ;push,list,list1
  list[cnt:cnt+nlist1-1] = list1
  cnt += nlist1

  BOMB:
ENDFOR  ; pixel loop
; Trim the structure
list = list[0:cnt-1]

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
