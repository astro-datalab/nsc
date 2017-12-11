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

; Load the final list of combine exposures
healstr = MRDFITS(listdir+'nsc_healpix_list.fits',1)
index = MRDFITS(listdir+'nsc_healpix_list.fits',2)


stop

end
