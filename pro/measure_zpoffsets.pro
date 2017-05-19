pro measure_zpoffsets,meta,obj,src

; Measure zero-point offsets for fields with many exposures.

dir = '/dl1/users/dnidever/nsc/instcal/combine/'
print,'Loading files'
if n_elements(meta) eq 0 then begin
  meta = mrdfits(dir+'stripe82_ra220dec0_expmeta.fits',1)
  meta.filter = strtrim(meta.filter,2)
  meta.base = strtrim(meta.base,2)
endif
nmeta = n_elements(meta)
if n_elements(obj) eq 0 then begin
  obj = mrdfits(dir+'stripe82_ra220dec0_object.fits',1)
  obj.id = strtrim(obj.id,2)
endif
nobj = n_elements(obj)
if n_elements(src) eq 0 then begin
  src = mrdfits(dir+'stripe82_ra220dec0_source.fits',1)
  src.sourceid = strtrim(src.sourceid,2)
  src.objectid = strtrim(src.objectid,2)
endif
nsrc = n_elements(src)
otag = tag_names(obj)

dum = strsplitter(src.sourceid,'.',/extract)
src_expname = reform(dum[0,*])+'.'+reform(dum[1,*])
indx = CREATE_INDEX(src_expname)

offsets = fltarr(nmeta)
rms = fltarr(nmeta)
print,strtrim(nmeta,2),' exposures'
for i=0,nmeta-1 do begin
  base = meta[i].base
  expnum = meta[i].expnum
  omagind = where(otag eq strupcase(meta[i].filter)+'MAG',nomagind)
  onumind = where(otag eq 'NDET'+strupcase(meta[i].filter),nonumind)
  ; get instrument
  instrument = ''
  if strmid(base,0,2) eq 'tu' or strmid(base,0,3) eq 'c4d' then instrument='c4d'
  if strmid(base,0,3) eq 'ksb' then instrument='ksb'
  if strmid(base,0,3) eq 'k4m' then instrument='k4m'
  if instrument eq '' then stop,'instrument not known'
  
  expname = instrument+'.'+strtrim(expnum,2)
  ;gd = where(src_expname eq expname,ngd)
  ;src1 = src[gd]
  expind = where(indstr.name eq expname,nexpind)
  src1 = src[indx.index[indx.lo[expind[0]]:indx.hi[expind[0]]]]
  ; only keep ones with good photometry
  gd1 = where(src1.cmag lt 50,ngd1)
  src2 = src1[gd1]
  ; now match objects
  MATCH,obj.id,src2.objectid,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    omag = obj[ind1].(omagind)
    onum = obj[ind1].(onumind)
    smag = src2[ind2].cmag
    ind = where(onum gt 1,nind)  ; only want ones with other detections
    diff = smag[ind]-omag[ind]
    offsets[i] = median([diff])
    rms[i] = mad([diff])
    print,strtrim(i+1,2),' ',expname,' ',strtrim(nind,2),' ',stringize(offsets[i],ndec=4),' ',stringize(rms[i],ndec=4),' mag'
  endif else print,'No matches'

endfor

stop

end
