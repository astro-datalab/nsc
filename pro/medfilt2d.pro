function medfilt2d,array,width,dim=dim,edge_copy=edge_copy,even=even,stp=stp

;+
;
; MEDFILT2D
;
; This performs 1D median filtering of a 2D array.
; Unfortunately the built-in MEDIAN function doesn't
; do this already.
;
; The basic idea is to collapse the 2D array into
; a 1D array and then use MEDIAN to do the 1D median
; filtering.  Then the result is unwrapped and the
; edges are fixed.  The first and last width/2 pixels
; are not filtered because they aren't enough pixels
; to do the filtering.
; Setting /edge_copy will fill these width/2 pixels
; with the last "good" median value.
; The number of "bad" pixels at the end depends on
; whether width is even/odd and the number of elements
; in array along the filtering dimension is even/odd
; At the beginning it is always width/2 and at the
; end it is either width/2 or width/2+1.
;
; INPUTS:
;  array    A 2D array to be 1D median smoothed
;  width    The sie of the median smoothing filter
;  =dim     The dimension over which to do the median
;             smoothing starting with 1 (i.e. dim=1 or 2)
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
;   result  The 2D array after 1D median filtering.
;
; USAGE:
;  IDL>result = medfilt2d(array,10,dim=1)
;
; By D. Nidever  Jan 2010
;-

; Not enough inputs
if n_elements(array) eq 0 or n_elements(width) eq 0 or n_elements(dim) eq 0 then begin
  print,'Syntax - result = medfilt2d(array,width,dim=dim,edge_copy=edge_copy,even=even,stp=stp)'
  return,-1
endif

; Check that array is 2D
sz = size(array)
if sz[0] ne 2 then begin
  print,'ARRAY must be 2D'
  return,-1
endif

; Check dim
if n_elements(dim) gt 1 then begin
  print,'DIM must be a scalar'
  return,-1
endif
if dim ne 1 and dim ne 2 then begin
  print,'DIM must be 1 or 2'
  return,-1
endif

; Check width
if n_elements(width) gt 1 then begin
  print,'WIDTH must be a scalar'
  return,-1
endif

; Median over FIRST dimension
;-----------------------------
if (dim eq 1) then begin

  if width gt sz[1] then begin
    print,'WIDTH mut be <= the number of elements in the first dimension'
    return,-1
  endif

  ; Make array 1D
  array1d = (array)(*)

  ; Perform the median filtering
  med = MEDIAN(array1d,width,even=even)

  ; Unwrap the result to 2D
  result = array*0    ; initialize result
  result[*,*] = med

  ; Deal with the edges
  wid2b = width/2      ; edge pixels at beginning
  wid2e = width-wid2b  ; edge pixels at end, odd if width odd
  ; if the number of pixels in the array (along this dimension)
  ; are even then wid2e=wid2b, not sure why
  if sz[1] mod 2 eq 0 then wid2e=wid2b


  ;   Leave the pixels unfiltered, their original values
  if not keyword_set(edge_copy) then begin
    result[0:wid2b-1,*] = array[0:wid2b-1,*]
    result[sz[1]-1-wid2e+1:sz[1]-1,*] = array[sz[1]-1-wid2e+1:sz[1]-1,*]
  ;   Copy the last "good" median value
  endif else begin
    result[0:wid2b-1,*] = (fltarr(wid2b)+1.0)#reform(result[wid2b,*])
    result[sz[1]-1-wid2e+1:sz[1]-1,*] = (fltarr(wid2e)+1.0)#reform(result[sz[1]-1-wid2e,*])
  endelse


; Median over SECOND dimension
;-----------------------------
endif else begin

  if width gt sz[2] then begin
    print,'WIDTH mut be <= the number of elements in the second dimension'
    return,-1
  endif

  ; Make array 1D, after transposing
  array1d = (transpose(array))(*)

  ; Perform the median filtering
  med = MEDIAN(array1d,width,even=even)

  ; Unwrap the result to 2D
  result = transpose(array*0)  ; initialize result
  result[*,*] = med

  ; Transpose back
  result = transpose(result)

  ; Deal with the edges
  wid2b = width/2      ; edge pixels at beginning
  wid2e = width-wid2b  ; edge pixels at end, odd if width odd

  ;   Leave the pixels unfiltered, their original values
  if not keyword_set(edge_copy) then begin
    result[*,0:wid2b-1] = array[*,0:wid2b-1]
    result[*,sz[2]-1-wid2e+1:sz[2]-1] = array[*,sz[2]-1-wid2e+1:sz[2]-1]
  ;   Copy the last "good" median value
  endif else begin
    result[*,0:wid2b-1] = result[*,wid2b]#(fltarr(wid2b)+1.0)
    result[*,sz[2]-1-wid2e+1:sz[2]-1] = result[*,sz[2]-1-wid2e]#(fltarr(wid2e)+1.0)
  endelse


endelse ; median filtering over second dimension

if keyword_set(stp) then stop

return,result

end
