pro grabcoords,files

; Grab the RA/DEC values out of the header and stuff them in a small
; file 

for i=0,n_elements(files)-1 do begin
  file = files[i]
  if file_test(file) eq 0 then goto,bomb
  head = headfits(file,exten=1,errmsg=errmsg)
  if errmsg eq '' then begin
    ra = sxpar(head,'crval1',count=nra)
    if nra eq 0 then ra=999999.
    dec = sxpar(head,'crval2',count=ndec)
    if ndec eq 0 then dec=999999.
  endif else begin
    ra = 999999.
    dec = 999999.
  endelse
  line = file+'  '+strtrim(ra,2)+'  '+strtrim(dec,2)
  base = file_basename(file,'.fits.fz')
  outfile = '/d0/dnidever/nsc/instcal/v2/tmp/coords/'+base+'_coords.txt'
  writeline,outfile,line
  BOMB:
endfor

end
