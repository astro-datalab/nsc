pro checkcmblogs

; Check the combine logs
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = localdir+'dnidever/nsc/instcal/tmp/'
cd,dir

files = file_search('nsccmb*.batch.log',count=nfiles)
info = file_info(files)
si = sort(info.mtime)
info = info[si]
print,strtrim(nfiles,2),' nsccmb logs found'

; cutoff
cutoff = 1.48924e+09  ; March 11, 2017
gd = where(info.mtime gt cutoff,ngd)
info = info[gd]
ninfo = n_elements(info)
print,strtrim(ninfo,2),' log files after the cutoff'

str = replicate({file:'',size:0LL,mtime:0LL,nexp:-1L,pix:-1L,noverlap:-1L,dt:-1.0,nobjects:-1L,empty:-1,bugerror:0,badformaterror:0,recent:0},ninfo)
str.file = info.name
str.size = info.size
str.mtime = info.mtime
; Check if this log was modified recently, might still be in progress
bd = where(str.mtime ge (systime(1)-3600.),nbd)
if nbd gt 0 then str[bd].recent=1
; Loop through the files
for i=0,ninfo-1 do begin                                             
  ; Get head
  spawn,['head','-100',str[i].file],hout,herrout,/noshell
  if n_elements(hout) gt 0 then begin
    ind1 = where(stregex(hout,'Healpix pixel = ',/boolean) eq 1,nind1)
    if nind1 gt 0 then begin
      dum = strsplit(hout[ind1[0]],' ',/extract)
      ; Combining InstCal SExtractor catalogs or Healpix pixel = 47400
      pix1 = long(dum[8])
      str[i].pix = pix1
    endif
    ind2 = where(stregex(hout,'exposures that overlap this pixel and neighbors',/boolean) eq 1,nind2)
    if nind2 gt 0 then begin
      dum = strsplit(hout[ind2[0]],' ',/extract)
      ; 1 exposures that overlap this pixel and neighbors
      noverlap = long(dum[0])
      str[i].noverlap = noverlap
    endif
  endif
  ; Get tail
  spawn,['tail','-100',str[i].file],tout,terrout,/noshell
  if n_elements(tout) gt 0 then begin
    ind3 = where(stregex(tout,'dt = ',/boolean) eq 1,nind3)
    if nind3 gt 0 then begin
      dum2 = strsplit(tout[ind3[0]],' ',/extract)
      sdt = dum2[2]
      sdt = repstr(sdt,',','')
      dt = float(sdt)
      str[i].dt = dt
    endif
    ind4 = where(stregex(tout,'final objects fall inside the pixel',/boolean) eq 1,nind4)
    if nind4 gt 0 then begin
      ; 699 final objects fall inside the pixel
      ; None of the final objects fall inside the pixel
      dum3 = strsplit(tout[ind4[0]],' ',/extract)
      if dum3[0] eq 'None' then begin
        str[i].empty = 1
      endif else begin
        nobjects = long(dum3[0])
        str[i].nobjects = nobjects
        str[i].empty = 0
      endelse
    endif
    ; If there are no sources in this pixel then it
    ; ends with "No sources in this pixel" and that's it
    ind5 = where(stregex(tout,'No sources in this pixel',/boolean) eq 1,nind5)
    if nind5 gt 0 then str[i].empty = 1
    ; Check for error due to bug for single exposure at end
    ind6 = where(stregex(tout,'Illegal subscript range: BRKLO.',/boolean) eq 1,nind6)
    if nind6 gt 0 then str[i].bugerror = 1
    ; Check for error due to bad catalog format
    ind7 = where(stregex(tout,'Illegal subscript range: BRKLO.',/boolean) eq 1,nind7)
    if nind6 gt 0 then str[i].badformaterror = 1
  endif  
  print,strtrim(i+1,2),' ',str[i].file,' ',strtrim(str[i].pix,2),' ',strtrim(str[i].nobjects,2),' ',strtrim(str[i].dt,2),' ',strtrim(str[i].empty,2),$
        ' ',strtrim(str[i].bugerror,2),' ',strtrim(str[i].badformaterror,2)
endfor

index=mrdfits(localdir+'dnidever/nsc/instcal/nsc_healpix_list.fits',2)
match,index.pix,str.pix,ind1,ind2,/sort,count=nmatch
str[ind2].nexp = index[ind1].nexp

outfile = localdir+'dnidever/nsc/instcal/nsccmb_summary.fits'
print,'Writing results to ',outfile
mwrfits,str,outfile,/create

stop

end
