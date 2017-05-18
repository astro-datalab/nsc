pro merge_ps1_cats

; Merge Arjun/Frank's PS1 catalogs, one file per filter right now

;dir = '/net/dl1/data/cp_calib_cats/ps1/'
dir = '/data0/dnidever/nsc/ps1/'
outdir = dir+'merged/'

;filters = ['g','r','i','z','Y','V']
filters = ['g','r','i','z','Y']
nfilters = n_elements(filters)

; Each directory is organized as files such as ps1_g-20_000.fits for
; PS1, g-band, -20.5 to -19.5 dec, 359.5 to 0.5 ra.  Unfortunately,
; there is a bug in the naming of the ps1 files for negative decs in g,
; r, and z that is offset by a degree.
;
; Each file is a FITS binary table of ra, dec, and mag.  I can't recall
; for sure right now but I believe the PS1 files have ra and dec
; substituted from Gaia.

; 000 001 002 003 004 
nra = 360
ra0 = 0
ra = string(lindgen(nra)+ra0,format='(i03)')
ndec = 122
dec0 = -31
;  +00, +01, -21, -22
dec = string(lindgen(ndec)+dec0,format='(i03)')
b = where(long(dec) ge 0,nb)
dec[b] = '+'+string(long(dec[b]),format='(i02)')

schema = {ra:0.0d0,dec:0.0d0,g:99.99,r:99.99,i:99.99,z:99.99,y:99.99}
tags = tag_names(schema)

; Loop over RA
for i=0,nra-1 do begin
;for i=10,nra-1 do begin
  ; Loop over DEC
  for j=0,ndec-1 do begin
;  for j=10,ndec-1 do begin
    undefine,final
    ; Loop through filters
    for k=0,nfilters-1 do begin
      ; Each directory is organized as files such as ps1_g-20_000.fits for
      ; PS1, g-band, -20.5 to -19.5 dec, 359.5 to 0.5 ra
      file = dir+filters[k]+'/ps1_'+filters[k]+dec[j]+'_'+ra[i]+'.fits'
      ; fix for g/r/z
      if long(dec[j]) lt 0 and (filters[k] eq 'g' or filters[k] eq 'r' or filters[k] eq 'z') then $
        file = dir+filters[k]+'/ps1_'+filters[k]+dec[j+1]+'_'+ra[i]+'.fits'
      if file_test(file) eq 1 then begin
        str = MRDFITS(file,1,/silent)
        nstr = n_elements(str)
        print,ra[i],' ',dec[j],' ',filters[k],' ',file,' ',nstr

        ; The number of stars each filter file are DIFFERENT
        ; Need to crossmatch them.
        magind = where(tags eq strupcase(filters[k]),nmagind)
        if n_elements(final) eq 0 then begin
          final = replicate(schema,nstr)
          final.ra = str.ra
          final.dec = str.dec
          final.(magind) = str.mag
        endif else begin
          SRCMATCH,final.ra,final.dec,str.ra,str.dec,0.5,ind1,ind2,/sph,count=nmatch
          print,' ',strtrim(nmatch,2),' matches'
          if nmatch gt 0 then final[ind1].(magind)=str[ind2].mag
          ;  Some "new" ones
          if nmatch lt nstr then begin
            left = str
            remove,ind2,left
            nleft = n_elements(left)
            print,' Adding ',strtrim(nleft,2),' new stars'
            new = replicate(schema,nleft)
            new.ra = left.ra
            new.dec = left.dec
            new.(magind) = left.mag
            final = [final,new]
          endif
        endelse
      endif else print,file,' NOT FOUND'
    endfor ; filter loop
    si = sort(final.ra)
    final = final[si]  ; sort by RA

    outfile = outdir+'ps1_'+dec[j]+'_'+ra[i]+'.fits'
    print,'Writing to ',outfile
    MWRFITS,final,outfile,/create

    ;stop
  endfor  ; dec loop
endfor  ; ra loop

stop

end
