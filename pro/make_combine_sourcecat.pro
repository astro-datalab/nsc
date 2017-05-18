pro make_combine_sourcecat,pix,source,redo=redo,stp=stp,nooutput=nooutput

; Make the source catalog for a single combined healpix

t0 = systime(1)

nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(pix) eq 0 then begin
  print,'Syntax - make_combine_sourcecat,pix,source,redo=redo,stp=stp,nooutput=nooutput'
  return
endif

; Check if output file already exists
outfile = dir+'combine/source/'+strtrim(pix,2)+'_source.fits'
if not keyword_set(nooutput) and (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,outfile,' EXISTS already and /redo not set'
  return
endif

print,'Creating source catalog for Healpix pixel = ',strtrim(pix,2)

; Does the object file exist
objfile = dir+'combine/'+strtrim(pix,2)+'.fits.gz'
if file_test(objfile) eq 0 then begin
  print,objfile,' NOT FOUND'
  return
endif

; Load the data
meta = MRDFITS(objfile,1,/silent)
meta.file = strtrim(meta.file,2)
meta.base = strtrim(meta.base,2)
meta.dateobs = strtrim(meta.dateobs,2)
nmeta = n_elements(meta)
print,strtrim(nmeta,2),' exposure(s)'
idstr = MRDFITS(objfile,3,/silent)
idstr.sourceid = strtrim(idstr.sourceid,2)
idstr.exposure = strtrim(idstr.exposure,2)
idstr.expnum = strtrim(idstr.expnum,2)
idstr.objectid = strtrim(idstr.objectid,2)
nidstr = n_elements(idstr)

; Initializing the source table


; Load the catalogs
for i=0,nmeta-1 do begin
  dateobs = meta[i].dateobs
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  catfile = dir+night+'/'+meta[i].base+'/'+meta[i].base+'_cat.fits'
  cat = MRDFITS(catfile,1,/silent)
  ncat = n_elements(cat)
  ; Initializing the source schema and catalog
  if n_elements(schema) eq 0 then begin
    schema = cat[0]
    struct_assign,{dum:''},schema
    add_tag,schema,'objectid','',schema
    nsource = 100000L
    source = replicate(schema,nsource)
    cnt = 0LL
  endif
  catid = strtrim(cat.ccdnum,2)+'.'+strtrim(cat.number,2)
  ind = where(idstr.exposure eq meta[i].base,nind)
  idstr1 = idstr[ind]
  srcid0 = idstr1.sourceid
  ; missing the instrument tag and expnum not correct for bok
  ;  just use ccdnum.number
  dum = strsplitter(srcid0,'.',/extract)
  srcid = reform(dum[1,*])+'.'+reform(dum[2,*])
  MATCH,catid,srcid,ind1,ind2,/sort,count=nmatch
  if nmatch ne nind then stop,'not all matched'
  print,strtrim(i+1,2),' ',meta[i].base,' ',strtrim(ncat,2),' ',strtrim(nmatch,2)

  ; Add new elements
  if cnt+nmatch gt nsource then begin
    old = source
    nnew = 100000L
    source = replicate(schema,nsource+nnew)
    source[0:nsource-1] = old
    nsource += nnew
    undefine,old
  endif

  ; Add to the source catalog
  temp = source[cnt:cnt+nmatch-1]
  struct_assign,cat[ind1],temp,/nozero
  temp.objectid = idstr1[ind2].objectid
  source[cnt:cnt+nmatch-1] = temp
  cnt += nmatch
endfor
; Trim any extra elements
source = source[0:cnt-1]
nsource = n_elements(source)
print,strtrim(nsource,2),' final sources'

; Write out
if not keyword_set(nooutput) then begin
  print,'Writing to ',outfile
  MWRFITS,source,outfile,/create
endif

if keyword_set(stp) then stop

end
