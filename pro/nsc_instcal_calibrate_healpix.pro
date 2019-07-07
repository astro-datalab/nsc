pro nsc_instcal_calibrate_healpix,pix,nside=nside,version=version,redo=redo

; This program is a wrapper around NSC_INSTCAL_CALIBRATE
; for all exposures in the same region of the sky.
; It retrieves the reference catalogs ONCE which saves time.

if n_elements(nside) eq 0 then nside=64
radeg = 180.0d0 / !dpi

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v3'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(dir,/directory) eq 0 then file_mkdir,dir+'logs/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir

t00 = systime(1)

; Setting pool thread values
if n_elements(ncpu) eq 0 then ncpu=1
CPU, TPOOL_NTHREADS = ncpu

; Not enough inputs
if n_elements(pix) eq 0 then begin
  print,'Syntax - nsc_instcal_calibrate_healpix,pix,nside=nside,version=version,redo=redo'
  return
endif

; Load the list of exposures
listfile = dir+'/lists/nsc_calibrate_healpix_list.fits'
if file_test(listfile) eq 0 then begin
  print,listfile,' NOT FOUND'
  return
endif
list = MRDFITS(listfile,1,/silent)

; Get the exposures for this healpix
print,'Calibrating InstCal SExtractor catalogs for Healpix pixel = ',strtrim(pix,2)
MATCH,list.pix,pix,ind,ind2,/sort,count=nind
if nind eq 0 then begin
  print,'No exposures'
  return
endif
print,'NEXPOSURES = ',strtrim(nind,2)
list1 = list[ind]
list1.expdir = strtrim(list1.expdir,2)
list1.instrument = strtrim(list1.instrument,2)
list1.filter = strtrim(list1.filter,2)

; Central coordinates
PIX2ANG_RING,nside,pix,centh,cenphi
cenra = cenphi[0]*radeg
cendec = 90-centh[0]*radeg
print,'RA  = ',stringize(cenra,ndec=6)
print,'DEC = ',stringize(cendec,ndec=6)
glactc,cenra,cendec,2000.0,glon,glat,1,/deg
print,'L = ',stringize(glon,ndec=6)
print,'B = ',stringize(glat,ndec=6)

; List of instrument-filters
filters = strtrim(list1.instrument,2)+'-'+strtrim(strmid(list1.filter,0,2),2)
ui = uniq(filters,sort(filters))
filters = filters[ui]

; Get required radius
;  DECam      needs 1.1 deg
;  Mosaic3    needs 0.43 deg
;  Bok90prime needs 0.75 deg
dum = where(list1.instrument eq 'c4d',nc4d)
dum = where(list1.instrument eq 'ksb',nksb)
minradius = 0.43
if nksb gt 0 then minradius = minradius > 0.75
if nc4d gt 0 then minradius = minradius > 1.1
; Add extra area for the size of the healpix pixel
;   nside=128 is roughly 27' across
;   nside=64 is roughly 55' across
;   nside=32 is roughly 110' across
radius = minradius + 0.5

; Get all of the reference data that we need
print,''
ref = GETREFDATA(filters,cenra,cendec,radius)

; Loop over the exposures
for i=0,nind-1 do begin
  print,''
  print,'---- EXPOSURE ',strtrim(i+1,2),' OF ',strtrim(nind,2),' ----'
  print,''
  expdir = list1[i].expdir
  lo = strpos(expdir,'/dl1')
  expdir = dldir + strmid(expdir,lo+5)
  NSC_INSTCAL_CALIBRATE,expdir,ref,redo=redo
endfor

print,''
print,'Total time = ',stringize(systime(1)-t00,ndec=2),' sec'

;stop

end
