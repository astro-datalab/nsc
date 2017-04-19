pro make_bok90prime_list,all

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/'

; Load all of the instcal exposures
if n_elements(all) eq 0 then begin
  all = mrdfits(dir+"bok90prime_archive_data.fits.gz",1)
  all.dtnsanam = strtrim(all.dtnsanam,2)
  all.dtacqnam = strtrim(all.dtacqnam,2)
  all.proctype = strtrim(all.proctype,2)
  all.prodtype2 = strtrim(all.prodtype2,2)
  all.date_obs = strtrim(all.date_obs,2)
  all.plver = strtrim(all.plver,2)
  ; Fix URI, ALREADY FIXED!!!
  ;uri = repstr(all.uri, 'irods:///noao-tuc-z1/', '/net/mss1/archive/')
  ;uri = strtrim(uri,2)
  ;all.uri = uri
  all.uri = strtrim(all.uri,2)
endif

; Get just the images
gdim = where(all.proctype eq 'InstCal' and all.prodtype2 eq 'image',ngim)
imstr = all[gdim]

; Get unique IDs
;  DTACQNAM        STRING    '/data_local/images/DTS/2013A-0609/DECam_00178042.fits.fz'
rawname = file_basename(imstr.dtacqnam)
uirname = uniq(rawname,sort(rawname))
urname = rawname[uirname]
nrname = n_elements(urname)

; Get latest version for each ID
; and use PLVER to get the most recent one
;  this also demotes blank PLVER entries
print,'Dealing with duplicates'
;  this could be faster with better index management
;alldbl = doubles(rawname,/all)
dbl = doubles(rawname,count=ndbl)
undefine,torem
for i=0,ndbl-1 do begin
  MATCH,rawname[alldbl],rawname[dbl[i]],ind1,ind2,/sort,count=nmatch
  dblind1 = alldbl[ind1]
  plver = imstr[dblind1].plver
  bestind = first_el(maxloc(plver))
  left = dblind1
  remove,bestind,left
  push,torem,left  
endfor
; removing duplicate entires
if n_elements(torem) gt 0 then remove,torem,imstr

; Make new structure with minimal columns
str = replicate({sp_id:'',sb_recno:0L,dtnsanam:'',dtacqnam:'',rawname:'',expnum:'',data_product_id:0L,uri:'',prop_id:'',ra:0.0d0,dec:0.0d0,$
                 exposure:0.0,release_date:'',date_obs:'',filter:'',mjd_obs:0.0d0,obstype:'',plver:'',proctype:'',prodtype:'',$
                 filename:'',pldname:'',fluxfile:'',maskfile:'',wtfile:''},nrname)
STRUCT_ASSIGN,imstr,str
str.rawname = file_basename(str.dtacqnam)
str.fluxfile = strtrim(str.uri,2)

; Get exposure number, d7429.0138.fits.fz
dum = strsplitter(str.rawname,'.',/extract)
expnum = strmid(reform(dum[0,*]),1)+reform(dum[1,*])  ; day number + expnum for that day
str.expnum = expnum

; Get mask and wtmap files
;-------------------------
print,'Getting mask and wtmap filenames'

; Deal with the c4d types first
gdnew = where(strmid(file_basename(str.fluxfile,'.fits.fz'),0,3) eq 'ksb',ngdnew)
base = file_basename(str[gdnew].fluxfile,'.fits.fz')
allbase = file_basename(all.uri,'.fits.fz')
; mask file
mbase = repstr(base,'oi_','od_')
MATCH,allbase,mbase,ind1,ind2,/sort
str[gdnew[ind2]].maskfile = all[ind1].uri
; wtmap file
wbase = repstr(base,'oi_','ow_')
MATCH,allbase,wbase,ind1,ind2,/sort
str[gdnew[ind2]].wtfile = all[ind1].uri

; Deal with the rest
;  use RAWNAME and PLVER and PRODTYPE
;bd = where(str.maskfile eq '' or str.wtfile eq '',nbd)
;allraw = file_basename(all.dtacqnam)
;; mask file
;strid = str[bd].rawname+'-'+str[bd].plver+'-InstCal-dqmask'
;allid = allraw+'-'+all.plver+'-'+all.proctype+'-'+all.prodtype2
;MATCH,allid,strid,ind1,ind2,/sort,count=nmatch
;str[bd[ind2]].maskfile = all[ind1].uri
;; wtmap file
;strid = str[bd].rawname+'-'+str[bd].plver+'-InstCal-wtmap'
;MATCH,allid,strid,ind1,ind2,/sort,count=nmatch
;str[bd[ind2]].wtfile = all[ind1].uri

; Only keeping ones with mask/weight files
gd3 = where(str.fluxfile ne '' and str.maskfile ne '' and str.wtfile ne '',ngd3,comp=bd3,ncomp=nbd3)
print,strtrim(nbd3,2),' exposures do NOT have flux/mask/wtmap files'
str = str[gd3]
print,strtrim(ngd3,2),' final exposures will the information we need'

;MWRFITS,str,dir+'bok90prime_instcal_list.fits',/create

stop

end
