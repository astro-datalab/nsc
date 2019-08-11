pro update_metafile_nmeas,expdir

;; Add NMEAS to the exposure structure in the metadatafile

if n_elements(expdir) eq 0 then begin
  print,'Syntax - update_metafile_nmeas,expdir'
  return
endif

base = file_basename(expdir)
metafile = expdir+'/'+base+'_meta.fits'
if file_test(metafile) eq 0 then begin
  print,metafile,' NOT FOUND'
  return
endif

;; Check if NMEAS was already added
hd = headfits(metafile,exten=1)
bd = where(stregex(hd,'NMEAS',/boolean,/fold_case) eq 1,nbd)
if nbd gt 0 then begin
  print,'NMEAS already added'
  return
endif

expstr = mrdfits(metafile,1)
chstrstr = mrdfits(metafile,2)
add_tag,expstr,'nmeas',0L,expstr
expstr.nmeas = long(total(chstr.nmeas))

stop

file_move,metafile,metafile+'.orig',/over,/allow
MWRFITS,expstr,metafile,/create
MWRFITS,chstr,metafile,/silent
if file_test(metafile) eq 1 then file_delete,metafile+'.orig'

stop

end
