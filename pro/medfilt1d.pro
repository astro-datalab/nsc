function medfilt1d,array,width,edge_copy=edge_copy,even=even,stp=stp

;+
;
; MEDFILT1D
;
; This is a wrapper for MEDIAN in the 1D case.  It additionally
; can fix the edges by repeating the last good median filtered
; value.
;
; INPUTS:
;  array    A 1D array to be median filtered
;  width    The sie of the median smoothing filter
;  /edge_copy  For the first and last width/2 pixel, where
;               there aren't enough pixels to do the median,
;               the last "good" median value will be used.
;               Otherwise these pixels are not filtered and
;               have their original values..
;  /even    Return the average of the two middle number if there
;              are an even number of elements in array.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;   result  The 1D array after median filtering.
;
; USAGE:
;  IDL>result = medfilt1d(array,10)
;
; By D. Nidever  Jan 2010
;-

; Not enough inputs
if n_elements(array) eq 0 or n_elements(width) eq 0 then begin
  print,'Syntax - result = medfilt1d(array,width,edge_copy=edge_copy,even=even,stp=stp)'
  return,-1
endif

; Check width
if n_elements(width) gt 1 then begin
  print,'WIDTH must be a scalar'
  return,-1
endif

sz = size(reform(array))

; Check that the array is 1D
if sz[0] ne 1 then begin
  print,'ARRAY must be 1D'
  return,-1
endif

if width gt sz[1] then begin
  print,'WIDTH mut be <= the number of elements in the first dimension'
  return,-1
endif

; Regular median
if not keyword_set(edge_copy) then return,MEDIAN(reform(array),width,even=even)

; Perform the median filtering
result = MEDIAN(reform(array),width,even=even)

; Deal with the edges
wid2b = width/2      ; edge pixels at beginning
wid2e = width-wid2b  ; edge pixels at end, odd if width odd

if edge_copy eq 1 then begin
  ; replicate the last median out to the edges
  ; if the number of pixels in the array (along this dimension)
  ; are even then wid2e=wid2b, not sure why
  if sz[1] mod 2 eq 0 then wid2e=wid2b
  
  ;   Copy the last "good" median value
  result[0:wid2b-1] = REPLICATE(result[wid2b],wid2b)
  result[sz[1]-1-wid2e+1:sz[1]-1] = REPLICATE(result[sz[1]-1-wid2e],wid2e)
endif else begin
  ; reflect the data for median of the edge
  tmp=0.*result
  for i=0,wid2b-1 do begin
    left=0
    right=i+width/2
    ngood=right-left+1
    nleft = width - ngood
    tmp[i] = MEDIAN([array[0:ngood],array[0:nleft-1]])
  endfor
  result[0:wid2b-1]=tmp[0:wid2b-1]
  for i=sz[1]-wid2e,sz[1]-1 do begin
    left=i-width/2
    right=sz[1]-1
    ngood=right-left+1
    nleft = width - ngood
    tmp[i] = MEDIAN([array[left:right],array[sz[1]-1-nleft:sz[1]-1]])
  endfor
  result[sz[1]-wid2e:sz[1]-1]=tmp[sz[1]-wid2e:sz[1]-1]
endelse

if keyword_set(stp) then stop

return,result

end
