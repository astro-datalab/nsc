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
    ; add calibration file stuff
    ; biasfil, bpmfil, bpm, pupfil, flatfil, illumfil, illumcor, starfil
    biasfil = sxpar(head,'biasfil',count=nbiasfil)
    if nbiasfil eq 0 then biasfil=''
    bpmfil = sxpar(head,'bpmfil',count=nbpmfil)
    if nbpmfil eq 0 then bpmfil=''
    bpm = sxpar(head,'bpm',count=nbpm)
    if nbpm eq 0 then bpm=''
    pupfil = sxpar(head,'pupfil',count=npupfil)
    if npupfil eq 0 then pupfil=''
    flatfil = sxpar(head,'flatfil',count=nflatfil)
    if nflatfil eq 0 then flatfil=''
    illumfil = sxpar(head,'illumfil',count=nillumfil)
    if nillumfil eq 0 then illumfil=''
    illumcor = sxpar(head,'illumcor',count=nillumcor)
    if nillumcor eq 0 then illumcor=''
    starfil = sxpar(head,'starfil',count=nstarfil)
    if nstarfil eq 0 then starfil=''
  endif else begin
    ra = 999999.
    dec = 999999.
    biasfil = ''
    bpmfil = ''
    bpm = ''
    pupfil = ''
    flatfil = ''
    illumfil = ''
    illumcor = ''
    starfil = ''
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

  sep = ' | '
  line = file+sep+strtrim(ra,2)+sep+strtrim(dec,2)+sep+strtrim(wcscal,2)+sep+strtrim(telstat,2)+sep+$
         biasfil+sep+bpmfil+sep+bpm+sep+pupfil+sep+flatfil+sep+illumfil+sep+illumcor+sep+starfil
  base = file_basename(file,'.fits.fz')
  ;outfile = '/d0/dnidever/nsc/instcal/v2/tmp/coords/'+base+'_coords.txt'
  outfile = '/data0/dnidever/nsc/instcal/v3/tmp/coords/'+base+'_coords.txt'
  writeline,outfile,line
  BOMB:
endfor

end
