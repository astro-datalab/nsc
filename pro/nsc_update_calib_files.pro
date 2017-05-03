pro nsc_update_calib_files,expdir
  
; Update the calibrated catalogs with sourceIDs, filter, MJD and depth (??)
; add depth to the metadata file

NSC_ROOTDIRS,dldir,mssdir,localdir

print,'Updating calibration catalog and metadata file for ',expdir

base = file_basename(expdir)
catfile = expdir+'/'+base+'_cat.fits'
metafile = expdir+'/'+base+'_meta.fits'
if file_test(catfile) eq 0 or file_test(metafile) eq 0 then begin
  print,'Catalog or metadata file NOT FOUND'
  return
endif

; What instrument is this?
instrument = 'c4d'  ; by default
if stregex(expdir,'/k4m/',/boolean) eq 1 then instrument='k4m'
if stregex(expdir,'/ksb/',/boolean) eq 1 then instrument='ksb'
print,'This is a '+instrument+' exposure'

; Load catalog
cat = MRDFITS(catfile,1,/silent,status=status)  
if status lt 0 then goto,BOMB
ncat = n_elements(cat)

if tag_exist(cat,'SOURCEID') eq 1 then begin
  print,'This file has already been updated'
  return
endif

; Load metadata
expstr = MRDFITS(metafile,1,/silent,status=status1)
if status1 lt 0 then goto,BOMB
nexpstr = n_elements(expstr)
chstr = MRDFITS(metafile,2,/silent,status=status2)
if status2 lt 0 then goto,BOMB
nchstr = n_elements(chstr)
chstr.filename = strtrim(chstr.filename,2)

; Load the logfile and get absolute flux filename
READLINE,expdir+'/'+base+'.log',loglines
ind = where(stregex(loglines,'Step #2: Copying InstCal images from mass store archive',/boolean) eq 1,nind)
line = loglines[ind[0]+1]
lo = strpos(line,'/archive')
; make sure the mss1 directory is correct for this server
fluxfile = mssdir+strtrim(strmid(line,lo+1),2)
head = headfits(fluxfile,exten=0)
expnum = sxpar(head,'expnum')
if instrument eq 'ksb' then begin  ; Bok doesn't have expnum
  ;DTACQNAM='/data1/batch/bok/20160102/d7390.0049.fits.fz'                                                                                                                          
  dtacqnam = sxpar(head,'DTACQNAM',count=ndtacqnam)
  if ndtacqnam eq 0 then begin
    print,'I cannot create an EXPNUM for this Bok exposure'
    return
  endif
  bokbase = file_basename(dtacqnam)
  dum = strsplit(bokbase,'.',/extract)
  boknight = strmid(dum[0],1)
  boknum = dum[1]
  expnum = boknight+boknum  ; concatenate the two numbers
endif

; Update Catalog schema
schema_cat = cat[0]
struct_assign,{dum:''},schema_cat
schema_cat = create_struct(schema_cat,'sourceid','','filter','','mjd',0.0d0)
newcat = replicate(schema_cat,ncat)
struct_assign,cat,newcat,/nozero
newcat.sourceid = instrument+'.'+strtrim(expnum,2)+'.'+strtrim(newcat.ccdnum,2)+'.'+strtrim(newcat.number,2)
newcat.filter = expstr.filter
newcat.mjd = expstr.mjd

; Update EXPSTR and CHSTR schema
newchstr = replicate({expdir:'',instrument:'',filename:'',ccdnum:0L,nsources:0L,cenra:999999.0d0,cendec:999999.0d0,$
                      gaianmatch:0L,rarms:999999.0,racoef:dblarr(4),decrms:999999.0,$
                      deccoef:dblarr(4),vra:dblarr(4),vdec:dblarr(4),zpterm:999999.0,$
                      zptermerr:999999.0,nrefmatch:0L,depth95:99.99,depth10sig:99.99},nchstr)
struct_assign,chstr,newchstr,/nozero
newchstr.expdir = expdir
newchstr.instrument = instrument

newexpstr = {file:'',instrument:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d0,dateobs:'',mjd:0.0d,filter:'',exptime:0.0,$
             airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,zpterm:999999.0,zptermerr:99999.0,$
             zptermsig:999999.0,zpspatialvar_rms:999999.0,zpspatialvar_range:999999.0,zpspatialvar_nccd:0,nrefmatch:0L,depth95:99.99,depth10sig:99.99}
struct_assign,expstr,newexpstr,/nozero
newexpstr.instrument = instrument

; Need good photometry
gdmag = where(cat.cmag lt 50,ngdmag)
if ngdmag eq 0 then goto,BOMB

; Measure the depth
; Get 95% percentile depth
cmag = cat[gdmag].cmag
si = sort(cmag)
cmag = cmag[si]
depth95 = cmag[round(0.95*ngdmag)-1]
print,'95% percentile depth = ',stringize(depth95,ndec=2),' mag'
; Get 10 sigma depth
;  S/N = 1.087/err
;  so S/N=5 is for err=1.087/5=0.2174
;  S/N=10 is for err=1.087/10=0.1087
depth10sig = 99.99
depind = where(cat.cmag lt 50 and cat.cmag gt depth95-3.0 and cat.cerr ge 0.0987 and cat.cerr le 0.1187,ndepind)
if ndepind lt 5 then depind = where(cat.cmag lt 50 and cat.cmag gt depth95-3.0 and cat.cerr ge 0.0787 and cat.cerr le 0.1387,ndepind)
if ndepind gt 5 then begin
  depth10sig = median([cat[depind].cmag])
endif else begin
  depind = where(cat.cmag lt 50,ndepind)
  if ndepind gt 0 then depth10sig=max([cat[depind].cmag])
endelse
print,'10sigma depth = ',stringize(depth10sig,ndec=2),' mag'

newchstr.depth95 = depth95
newchstr.depth10sig = depth10sig
newexpstr.depth95 = depth95
newexpstr.depth10sig = depth10sig

; Temporarily move the old metadata and catalogs files aside for safe
; keeping
file_move,catfile,catfile+'.bak',/allow,/over
file_move,metafile,metafile+'.bak',/allow,/over

; Write new versions of files
print,'Writing new catalog to ',catfile
MWRFITS,newcat,catfile,/create
print,'Writing new metadata to ',metafile
MWRFITS,newexpstr,metafile,/create
MWRFITS,newchstr,metafile,/silent

; Delete old files
file_delete,catfile+'.bak',/allow
file_delete,metafile+'.bak',/allow

BOMB:

end
