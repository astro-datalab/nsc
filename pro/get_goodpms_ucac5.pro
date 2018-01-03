pro get_goodpms_ucac5,ra,dec,radius,refname=refname,version=version,redo=redo

; Get the good PM values and crossmatch with reference catalog

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
cmbdir = dir+'combine/'
if n_elements(refname) eq 0 then refname='UCAC5'
radeg = 180.0d0 / !dpi

if n_elements(ra) eq 0 or n_elements(dec) eq 0 then begin
  print,'Syntax - get_goodpms_ucac5,ra,dec,radius,refname=refname,version=version,redo=redo'
  return
endif
if n_elements(radius) eq 0 then radius=1.0

if dec lt 0 then decsgn='-' else decsgn='+'
outfile = cmbdir+'ucac5pm/'+stringize(ra,ndec=5)+decsgn+stringize(abs(dec),ndec=5)+'-'+stringize(radius,ndec=2)+'.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,'Output file exists and /redo NOT set'
  return
endif

print,'Getting good proper motions for'
print,'RA  = ',stringize(ra,ndec=5)
print,'DEC = ',stringize(dec,ndec=5)
print,'RAD = ',stringize(radius,ndec=5)

; Get all the pixels within the radius
theta = (90-dec)/radeg
phi = ra/radeg
ANG2VEC,theta,phi,vec
QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive
undefine,obj
nobj = 0LL
for i=0,nlistpix-1 do begin
  ;print,strtrim(i+1,2),'/',strtrim(nlistpix,2),' ',strtrim(listpix[i],2)
  pixfile = cmbdir+strtrim(long(listpix[i])/1000,2)+'/'+strtrim(listpix[i],2)+'.fits.gz'
  if file_test(pixfile) eq 1 then begin
    head = headfits(pixfile,exten=2)
    nobj += sxpar(head,'NAXIS2')
  endif
endfor
cnt = 0LL
for i=0,nlistpix-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nlistpix,2),' ',strtrim(listpix[i],2)
  pixfile = cmbdir+strtrim(long(listpix[i])/1000,2)+'/'+strtrim(listpix[i],2)+'.fits.gz'
  if file_test(pixfile) eq 1 then begin
    obj1 = MRDFITS(pixfile,2,/silent)
    nobj1 = n_elements(obj1)
    if n_elements(obj) eq 0 then begin
      schema = obj1[0]
      struct_assign,{dum:''},schema
      obj = replicate(schema,nobj)
    endif
    obj[cnt:cnt+nobj1-1] = obj1
    cnt += nobj1
  endif
endfor
nobj = n_elements(obj)
if nobj eq 0 then begin
  print,'No objects found'
  return
endif
print,strtrim(nobj,2),' objects'

; Load the NSC data
;tempfile = MKTEMP('nscpm')
;cmd = "psql -h gp04.datalab.noao.edu -U datalab -d tapdb -w --pset footer -c 'SELECT * FROM nsc_dr1.object "+$
;      " WHERE q3c_radial_query(ra,dec,"+stringize(ra,ndec=4,/nocomma)+","+stringize(dec,ndec=4,/nocomma)+$
;      ","+stringize(radius,ndec=3)+")' > "+tempfile
;file_delete,tempfile,/allow
;spawn,cmd,out,outerr
;; Check for empty query
;READLINE,tempfile,tlines,nlineread=4
;if n_elements(tlines) lt 4 then begin
;  print,'No Results'
;  return
;endif                        
;;  Load ASCII file and create the FITS file
;obj = importascii(tempfile,/header,delim='|',skipline=2,/silent)
;file_delete,tempfile,/allow

; keep the good ones, only stars
;gdobj = where(obj.deltamjd ge 1000 and obj.ndet ge 5 and obj.gmag le 20 and obj.fwhm lt 1.5,ngdobj)
gdobj = where(finite(obj.pmra) eq 1 and finite(obj.pmdec) eq 1 and $
              (abs(obj.pmra)/obj.pmraerr ge 3 or obj.pmraerr lt 3) and $
              (abs(obj.pmdec)/obj.pmdecerr ge 3 or obj.pmdecerr lt 3) and $
              obj.deltamjd ge 100 and obj.ndet ge 3 and obj.gmag le 20 and obj.fwhm lt 1.5,ngdobj)
if ngdobj eq 0 then begin
  print,'No good ones to keep'
  return
endif
print,strtrim(ngdobj,2),' objects with good proper motions and stellar morphology kept'
obj = obj[gdobj]
nobj = n_elements(obj)

; Load the reference catalog
tempfile = MKTEMP('ucac5')
;cmd = "stilts tapquery tapurl='http://cds.vizier/tap' "
;cmd = "stilts tapquery tapurl='http://vizier.cfa.harvard.edu/tap' "
;cmd = "stilts tapquery tapurl='http://cds.vizier/tap' "
cmd = "stilts tapquery tapurl='http://tapvizier.u-strasbg.fr/TAPVizieR/tap' "
cmd = cmd+'adql="SELECT * FROM ucac5 WHERE '
;cmd = cmd+'1='+"CONTAINS(POINT('ICRS', RAJ2000, DEJ2000),CIRCLE('ICRS', "+strtrim(ra,2)+","+strtrim(dec,2)+","+strtrim(radius,2)+' ))" out='+tempfile+' ofmt=fits'
cmd = cmd+'1='+"CONTAINS(POINT('ICRS', ira, idc),CIRCLE('ICRS', "+strtrim(ra,2)+","+strtrim(dec,2)+","+strtrim(radius,2)+' ))" out='+tempfile+' ofmt=fits'
; Execute the command
print,'Getting UCAC5 data'
SPAWN,cmd,out,errout
; Check that the file is there
if file_test(tempfile) eq 0 then begin
  print,'Temporary UCAC5 file not found'
  return
endif
; Read in the UCAC5 data
ucac = MRDFITS(tempfile,1,/silent)
FILE_DELETE,tempfile,/allow

; Crossmatch
dcr = 2.0  ; 1.0
SRCMATCH,obj.ra,obj.dec,ucac.ira,ucac.idc,dcr,ind1,ind2,/sph,count=nmatch
if nmatch eq 0 then begin
  print,'No NSC and UCAC5 matches'
  return
endif
print,strtrim(nmatch,2),' matches'
obj2 = obj[ind1]
ucac2 = ucac[ind2]

; Add columns to obj table
objtags = tag_names(obj)
nobjtags = n_elements(objtags)
objtypes = lonarr(nobjtags)
for i=0,nobjtags-1 do objtypes[i]=size(obj[0].(i),/type)
newtags = [objtags,'ucac_dist','ucac_ra','ucac_dec','ucac_gra','ucac_gdec','ucac_pmra','ucac_pmdec','ucac_epmra','ucac_epmdec','ucac_gmag','ucac_umag','ucac_rmag','ucac_jmag','ucac_hmag','ucac_kmag']
newtypes = [objtypes,4,5,5,5,5,4,4,4,4,4,4,4,4,4,4]
schema = create_struct(newtags[0],fix('',type=newtypes[0]))
for i=1,n_elements(newtags)-1 do schema=create_struct(schema,newtags[i],fix('',type=newtypes[i]))
newobj = replicate(schema,nmatch)
STRUCT_ASSIGN,obj2,newobj,/nozero
newobj.ucac_ra = ucac2.ira
newobj.ucac_dec = ucac2.idc
newobj.ucac_gra = ucac2.rag
newobj.ucac_gdec = ucac2.dcg
newobj.ucac_pmra = ucac2.pmur
newobj.ucac_pmdec = ucac2.pmer
newobj.ucac_epmra = ucac2.pmud
newobj.ucac_epmdec = ucac2.pmed
newobj.ucac_gmag = ucac2.gmag
newobj.ucac_umag = ucac2.umag
newobj.ucac_jmag = ucac2.jmag
newobj.ucac_hmag = ucac2.hmag
newobj.ucac_kmag = ucac2.kmag
newobj.ucac_dist = sphdist(obj2.ra,obj2.dec,ucac2.ira,ucac2.idc,/deg)*3600

; Save the final file
outdir = file_dirname(outfile)
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
print,'Writing results to ',outfile
MWRFITS,newobj,outfile,/create

;stop

end
