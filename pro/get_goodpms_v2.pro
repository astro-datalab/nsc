pro get_goodpms_v2,pix,refname=refname,version=version,redo=redo

; Get the good PM values and crossmatch with reference catalog

if n_elements(pix) eq 0 then begin
  print,'Syntax - get_goodpms,pix,refname=refname,version=version,redo=redo'
  return
endif

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
cmbdir = dir+'combine/'
if n_elements(refname) eq 0 then refname='HSOY'

outfile = cmbdir+'hsoypm_v2/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,'Output file exists and /redo NOT set'
  return
endif

print,'Getting good proper motions for pixel = ',strtrim(pix,2)

pixfile = cmbdir+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits.gz'
meta = MRDFITS(pixfile,1,/silent)

; Check deltaMJD and Nexposures
;if range(meta.mjd) lt 1000 or n_elements(meta) lt 5 then begin
if range(meta.mjd) lt 100 or n_elements(meta) lt 3 then begin
  print,'deltaMJD or Nexposures not large enough'
  return
endif

; Load the object information
obj = MRDFITS(pixfile,2,/silent)
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

; Get the RA and DEC ranges
cendec = mean(minmax(obj.dec))
rdec = range(obj.dec)
cenra = mean(minmax(obj.ra))
rra = range(obj.ra)*cos(cendec/!radeg)
if range(obj.ra) gt 100 then begin
  ra = obj.ra
  bdra = where(ra gt 180,nbdra)
  if nbdra gt 0 then ra[bdra]-=360
  bdra2 = where(ra lt -180,nbdra2)
  if nbdra2 gt 0 then ra[bdra2]+=360
  cenra = mean(minmax(ra))
  if cenra lt 0 then cenra+=360
  rra = range(ra)
endif

; Load the reference catalog
rad = (rra > rdec)*0.7
tempfile = MKTEMP('hsoy')
cmd = "stilts tapquery tapurl='http://dc.g-vo.org/tap' "
cmd = cmd+'adql="SELECT TOP 30000000 * FROM hsoy.main WHERE '
cmd = cmd+'1='+"CONTAINS(POINT('ICRS', raj2000, dej2000),CIRCLE('ICRS', "+strtrim(cenra,2)+","+strtrim(cendec,2)+","+strtrim(rad,2)+' ))" out='+tempfile+' ofmt=fits'
; Execute the command
print,'Getting HSOY data'
SPAWN,cmd,out,errout
; Check that the file is there
if file_test(tempfile) eq 0 then begin
  print,'Temporary HSOY file not found'
  return
endif
; Read in the HSOY data
hsoy = MRDFITS(tempfile,1,/silent)
FILE_DELETE,tempfile,/allow

; Crossmatch
SRCMATCH,obj.ra,obj.dec,hsoy.raj2000,hsoy.dej2000,1.0,ind1,ind2,/sph,count=nmatch
if nmatch eq 0 then begin
  print,'No NSC and HSOY matches'
  return
endif
print,strtrim(nmatch,2),' matches'
obj2 = obj[ind1]
hsoy2 = hsoy[ind2]

; Add columns to obj table
objtags = tag_names(obj)
nobjtags = n_elements(objtags)
objtypes = lonarr(nobjtags)
for i=0,nobjtags-1 do objtypes[i]=size(obj[0].(i),/type)
newtags = [objtags,'hsoy_pmra','hsoy_pmdec','hsoy_epmra','hsoy_epmdec','gaiaid']
newtypes = [objtypes,4,4,4,4,14]
schema = create_struct(newtags[0],fix('',type=newtypes[0]))
for i=1,n_elements(newtags)-1 do schema=create_struct(schema,newtags[i],fix('',type=newtypes[i]))
newobj = replicate(schema,nmatch)
STRUCT_ASSIGN,obj2,newobj,/nozero
newobj.hsoy_pmra = hsoy2.pmra
newobj.hsoy_pmdec = hsoy2.pmde
newobj.hsoy_epmra = hsoy2.e_pmra
newobj.hsoy_epmdec = hsoy2.e_pmde
newobj.gaiaid = hsoy2.gaia_id

; Save the final file
outdir = file_dirname(outfile)
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
print,'Writing results to ',outfile
MWRFITS,newobj,outfile,/create

;stop

end
