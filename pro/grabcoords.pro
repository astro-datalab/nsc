pro grabcoords,files

; Grab the RA/DEC values out of the header and stuff them in a small file 
; upgraded to also grab the WCSCAL and TELSTAT from header
; to tell if it was tracking or not

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

  head0 = headfits(file,exten=0,errmsg=errmsg0)
  if errmsg0 eq '' then begin
    wcscal = sxpar(head0,'wcscal',count=nwcscal)
    if nwcscal eq 0 then wcscal='NAN'
    telstat = sxpar(head0,'telstat',count=ntelstat)
    if ntelstat eq 0 then telstat='NAN'
  endif else begin
    wcscal = 'NAN'
    telstat = 'NAN'
  endelse

  line = file+'  '+strtrim(ra,2)+'  '+strtrim(dec,2)+'  '+strtrim(wcscal,2)+'  '+strtrim(telstat,2)
  base = file_basename(file,'.fits.fz')
  outfile = '/d0/dnidever/nsc/instcal/v2/tmp/coords/'+base+'_coords.txt'
  writeline,outfile,line
  BOMB:
endfor

end
