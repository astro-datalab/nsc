pro create_exposure_table

; Create the NSC exposure table

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
listdir = dir+'lists/'
cmbdir = dir +'combine/'

; exposures, EXPOSURE_mets.fits[1], nsc_instcal_calibrate.fits
;     -use the nsc_healpix_list.fits file to figure out which exposures passed the cut
;     -need to figure out how many sources were actually used, I think I need to actually load the object.fits.gz[1] meta-data info to get this

; Load the calibration summary file
sum = MRDFITS(listdir+'nsc_instcal_calibrate.fits',1)
sum.expdir = strtrim(sum.expdir,2)
sum.base = strtrim(sum.base,2)
sexpdir = sum.expdir
b = where(strmid(sexpdir,0,4) ne '/net',nb)
if nb gt 0 then sexpdir[b]='/net'+sexpdir[b]

; Load the final list of combine exposures
healstr = MRDFITS(listdir+'nsc_healpix_list.fits',1)
healstr.file = strtrim(healstr.file,2)
healstr.base = strtrim(healstr.base,2)
index = MRDFITS(listdir+'nsc_healpix_list.fits',2)
ui = uniq(healstr.file,sort(healstr.file))
uhstr = healstr[ui]
hexpdir = file_basename(uhstr.file)
b = where(strmid(hexpdir,0,4) ne '/net',nb)
if nb gt 0 then hexpdir[b]='/net'+hexpdir[b]

; Match HEALSTR exposures to SUM to get the combined exposures list
MATCH,sexpdir,hexpdir,ind1,ind2,/sort,count=nmatch
fsum = sum[ind1]  ; final list of exposures

; Final columns
tags = ['instrument', 'exposure','expnum','ra','dec','dateobs','mjd','filter','exptime','airmass','nsources','fwhm','nchips','rarms','rastderr',$
        'decrms','decstderr','ebv','ngaiamatch','zptype','zpterm','zptermerr','zptermsig','nrefmatch','depth95','depth10sig']
types = [7,7,3,5,5,7,5,7,4,4,3,4,2,4,4,4,4,4,3,1,4,4,4,3,4,4]
; ngoodgaiamatch, ngoodrefmatch
; zpspatialvar_XXX
schema = create_struct(tags[0],fix('',type=types[0]))
for i=1,n_elements(tags)-1 do schema=create_struct(schema,tags[i],fix('',type=types[i]))
nexp = n_elements(healstr)
expstr = REPLICATE(schema,nexp)
STRUCT_ASSIGN,fsum,expstr  ; copy over information
expstr.nsources = 0

; Loop over the exposures and get the final list of sources
for i=0,nexp-1 do begin
  stop  

endfor

; Cut out any exposures with NO sources
bdexp = where(expstr.nsources eq 0,nbdexp)
print,strtrim(nbdexp,2),' exposures with NO sources'
;REMOVE,bdexp,expstr


stop

end
