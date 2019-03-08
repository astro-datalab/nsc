;+
;
; GETEXTTYPE
;
;  Get the extinction type for the NSC.
;
; INPUTS:
;  cenra    Right Ascension at the center of the image.
;  cendec   Declination at the center of the image.
;  radius   Radius of the image.
;
; OUTPUTS:
;  exttype  Extinction type:
;            1 - SFD
;            2 - RJCE ALLWISE
;            3 - RJCE GlIMPSE
;            4 - RJCE SAGE
;
; USAGE:
;  IDL>ext_type = getexttype(cenra,cendec,radius)
;
; By D. Nidever  Feb 2019
;-

function getexttype,cenra,cendec,radius

;; Not enough inputs
if n_elements(cenra) eq 0 or n_elements(cendec) eq 0 or n_elements(radius) eq 0 then begin
  print,'Syntax - ext_type = getexttype(cenra,cendec,radius)'
  return,-1
endif

  
;; Figure out what reddening method we are using
;;----------------------------------------------
;; Extinction types:
;; 1 - SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2
;; 2 - RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE
;;     data available
;; 3 - RJCE GLIMPSE, GLIMPSE data available
;; 4 - RJCE SAGE, SAGE data available
ext_type = 0
GLACTC,cenra,cendec,2000.0,cengl,cengb,1,/deg

;; Get grid of SFD EBV values across the field to see
;;  what regime we are working in
x = scale_vector(findgen(100),-radius,radius)#replicate(1,100)
y = replicate(1,100)#scale_vector(findgen(100),-radius,radius)
ROTSPHCEN,x,y,cenra,cendec,rr,dd,/gnomic,/reverse
GLACTC,rr,dd,2000.0,gl,gb,1,/deg
ebv_grid = dust_getval(gl,gb,/noloop,/interp)
maxebv = max(ebv_grid)

;; Check if there is any GLIMPSE data available
;  do 3x3 grid with small matching radius
if abs(cengb) lt 5 and (cengl lt 65 or cengl gt 290) then begin
  x = scale_vector(findgen(3),-radius,radius)#replicate(1,3)
  y = replicate(1,3)#scale_vector(findgen(3),-radius,radius)
  ROTSPHCEN,x,y,cenra,cendec,rr,dd,/gnomic,/reverse
  rr = reform(rr)
  dd = reform(dd)
  ncat = 0
  cnt = 0L
  while (ncat eq 0) and (cnt lt 9) do begin
    cat = QUERYVIZIER('II/293/glimpse',[rr[cnt],dd[cnt]],1.0,count=ncat,/silent)
    cnt++
  endwhile
  if ncat gt 0 then ext_type = 3
endif

;; Check if there is any SAGE data available
;  do 3x3 grid with small matching radius
lmcrad = sphdist(81.9,-69.867,cenra,cendec,/deg)
smcrad = sphdist(13.183,-72.8283,cenra,cendec,/deg)
if lmcrad lt 5.0 or smcrad lt 4.0 then begin
  x = scale_vector(findgen(3),-radius,radius)#replicate(1,3)
  y = replicate(1,3)#scale_vector(findgen(3),-radius,radius)
  ROTSPHCEN,x,y,cenra,cendec,rr,dd,/gnomic,/reverse
  rr = reform(rr)
  dd = reform(dd)
  ncat = 0
  cnt = 0L
  while (ncat eq 0) and (cnt lt 9) do begin
    cat = QUERYVIZIER('II/305/archive',[rr[cnt],dd[cnt]],1.0,count=ncat,/silent)
    cnt++
  endwhile
  if ncat gt 0 then ext_type = 4
endif

;; Use RJCE ALLWISE, |b|<16 or max(EBV)>0.2 and not GLIMPSE or SAGE
;;     data available
if ext_type eq 0 and (abs(cengb) lt 16 or maxebv gt 0.2) then ext_type=2

;; SFD, |b|>16 and RLMC>5.0 and RSMC>4.0 and max(EBV)<0.2
if ext_type eq 0 then ext_type=1

return,ext_type

end
