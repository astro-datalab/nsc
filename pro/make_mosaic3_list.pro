pro make_mosaic3_list,all,file=file,version=version

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/'

if n_elements(version) eq 0 then version = 'v3' ;'v2'
if n_elements(file) eq 0 then file = dir+'instcal/'+version+'/lists/mosaic3_archive_info.fits.gz'

; Load all of the instcal exposures
if n_elements(all) eq 0 then begin
  all = mrdfits(file,1)
  all.dtnsanam = strtrim(all.dtnsanam,2)
  all.dtacqnam = strtrim(all.dtacqnam,2)
  all.proctype = strtrim(all.proctype,2)
  ;all.prodtype2 = strtrim(all.prodtype2,2)
  all.prodtype = strtrim(all.prodtype,2)
  all.date_obs = strtrim(all.date_obs,2)
  all.plver = strtrim(all.plver,2)
  ;; Fix URI, ALREADY FIXED!!!
  uri = repstr(all.uri, 'irods:///noao-tuc-z1/', '/net/mss1/archive/')
  uri = strtrim(uri,2)
  all.uri = uri
  all.uri = strtrim(all.uri,2)
endif

; Get just the images
;gdim = where(all.proctype eq 'InstCal' and all.prodtype2 eq 'image',ngim)
gdim = where(all.proctype eq 'InstCal' and all.prodtype eq 'image',ngim)
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
alldbl = doubles(rawname,/all)
dbl = doubles(rawname)
ndbl = n_elements(dbl)
print,strtrim(ndbl,2),' duplicates'
undefine,torem
for i=0,ndbl-1 do begin
  if i mod 1000 eq 0 then print,i
  MATCH,rawname[alldbl],rawname[dbl[i]],ind1,ind2,/sort,count=nmatch
  dblind1 = alldbl[ind1]
  plver = imstr[dblind1].plver
  bestind = first_el(maxloc(plver))
  left = dblind1
  remove,bestind,left
  push,torem,left  
endfor
; removing duplicate entires
remove,torem,imstr

; Make new structure with minimal columns
str = replicate({instrument:'',sp_id:'',sb_recno:0L,dtnsanam:'',dtacqnam:'',rawname:'',expnum:'',data_product_id:0L,uri:'',prop_id:'',ra:0.0d0,dec:0.0d0,$
                 exposure:0.0,release_date:'',date_obs:'',filter:'',mjd_obs:0.0d0,obstype:'',plver:'',proctype:'',prodtype:'',$
                 filename:'',pldname:'',fluxfile:'',maskfile:'',wtfile:''},nrname)
STRUCT_ASSIGN,imstr,str
str.instrument = 'k4m'  ; KPNO-4m+Mosaic3
str.rawname = file_basename(str.dtacqnam)
str.fluxfile = strtrim(str.uri,2)

; Get exposure number
dum = strsplitter(str.rawname,'.',/extract)
;expnum = strmid(file_basename(str.rawname,'.fits.fz'),5)
expnum = reform(dum[1,*])
; Problem exposures
bd = where(valid_num(expnum) eq 0 and strmid(str.rawname,0,4) eq 'mos3',nbd)  ; mos3107016.fits
if nbd gt 0 then expnum[bd]=strmid(reform(dum[0,bd]),4)
bd = where(valid_num(expnum) eq 0 and strmid(str.rawname,0,3) eq 'spw',nbd)  ; spw95892.fits
if nbd gt 0 then expnum[bd]=strmid(reform(dum[0,bd]),3)
bd = where(valid_num(expnum) eq 0 and strmid(str.rawname,0,2) eq 'ic',nbd)  ; ic101810.fits
if nbd gt 0 then expnum[bd]=strmid(reform(dum[0,bd]),2)
bd = where(valid_num(expnum) eq 0 and stregex(strmid(str.rawname,0,3),'n?_',/boolean) eq 1,nbd)  ; n1_102064.fits
if nbd gt 0 then expnum[bd]=strmid(reform(dum[0,bd]),3)
str.expnum = expnum
; there are 6 left with messed up names, get rid of them
bd = where(valid_num(str.expnum) eq 0,nbd)
if nbd gt 0 then begin
  print,'Getting rid of ',strtrim(nbd,2),' exposures with messed up raw names'
  remove,bd,str
endif

; Get mask and wtmap files
;-------------------------
print,'Getting mask and wtmap filenames'

; Deal with the c4d types first
gdnew = where(strmid(file_basename(str.fluxfile,'.fits.fz'),0,3) eq 'k4m',ngdnew)
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
bd = where(str.maskfile eq '' or str.wtfile eq '',nbd)
allraw = file_basename(all.dtacqnam)
; mask file
strid = str[bd].rawname+'-'+str[bd].plver+'-InstCal-dqmask'
;allid = allraw+'-'+all.plver+'-'+all.proctype+'-'+all.prodtype2
allid = allraw+'-'+all.plver+'-'+all.proctype+'-'+all.prodtype
MATCH,allid,strid,ind1,ind2,/sort,count=nmatch
str[bd[ind2]].maskfile = all[ind1].uri
; wtmap file
strid = str[bd].rawname+'-'+str[bd].plver+'-InstCal-wtmap'
MATCH,allid,strid,ind1,ind2,/sort,count=nmatch
str[bd[ind2]].wtfile = all[ind1].uri

; Only keeping ones with mask/weight files
gd3 = where(str.fluxfile ne '' and str.maskfile ne '' and str.wtfile ne '',ngd3,comp=bd3,ncomp=nbd3)
print,strtrim(nbd3,2),' exposures do NOT have flux/mask/wtmap files'
str = str[gd3]
print,strtrim(ngd3,2),' final exposures with the information we need'

; APPLY RELEASE DATE CUTS
release_date = strtrim(str.release_date,2)
release_year = long(strmid(release_date,0,4))
release_month = long(strmid(release_date,5,2))
release_day = long(strmid(release_date,8,2))
release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
;release_cutoff = [2017,4,24]  ; v1 - April 24, 2017
;release_cutoff = [2017,10,11]  ; v2 - Oct 11, 2017
release_cutoff = [2018,10,9]  ; v3 - Oct 9, 2018
release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
gdrelease = where(release_mjd le release_cutoff_mjd,ngdrelease,comp=bdrelease,ncomp=nbdrelease)
print,strtrim(ngdrelease,2),' exposures are PUBLIC'
str = str[gdrelease]  ; impose the public data cut

;MWRFITS,str,dir+'instcal/'+version+'/lists/mosaic3_instcal_list.fits',/create

stop

end
