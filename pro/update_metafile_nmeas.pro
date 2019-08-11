pro update_metafile_nmeas,expdirs

;; Add NMEAS to the exposure structure in the metadatafile

if n_elements(expdirs) eq 0 then begin
  print,'Syntax - update_metafile_nmeas,expdirs'
  return
endif

;; Loop over exposure directories
nexpdir = n_elements(expdirs)
for i=0,nexpdir-1 do begin
  expdir = expdirs[i]
  print,strtrim(i+1,2),'/',strtrim(nexpdir,2),' ',expdir

  base = file_basename(expdir)
  metafile = expdir+'/'+base+'_meta.fits'
  if file_test(metafile) eq 0 then begin
    print,metafile,' NOT FOUND'
    goto,BOMB
  endif

  ;; Check if NMEAS was already added
  hd = headfits(metafile,exten=1)
  bd = where(stregex(hd,'TTYPE',/boolean) eq 1 and stregex(hd,'NMEAS',/boolean) eq 1,nbd)
  if nbd gt 0 then begin
    print,'NMEAS already added'
    goto,BOMB
  endif

  expstr = mrdfits(metafile,1)
  chstr = mrdfits(metafile,2)
  add_tag,expstr,'nmeas',0L,expstr
  expstr.nmeas = long(total(chstr.nmeas))
  print,'NMEAS = ',strtrim(expstr.nmeas,2)

  file_move,metafile,metafile+'.orig',/over,/allow
  MWRFITS,expstr,metafile,/create
  MWRFITS,chstr,metafile,/silent
  if file_test(metafile) eq 1 then file_delete,metafile+'.orig'

  BOMB:
endfor

end
