pro nsc_instcal_calibrate_coverage_main_combinesubs

;; Combine partial runs of nsc_instcal_calibrate_coverage_main.pro

version = 'v3'
subfiles = file_search('/dl1/users/dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_calibrate_coverage_sub*.fits',count=nsubfiles)
print,strtrim(nsubfiles,2),' sub files to combine'

nside = 512L
radeg = 180.0d0/!dpi

nmeas = 0L
for s=0,nsubfiles-1 do begin
  print,strtrim(s+1,2),' ',subfiles[s]
  covstr1 = mrdfits(subfiles[s],1)
  nmeas += total(covstr1.nmeas>0)
  if s eq 0 then begin
    covstr = covstr1
  endif else begin
    ;; combine
    covstr.nobj >= covstr1.nobj
    covstr.nmeas += covstr1.nmeas
    covstr.depth >= covstr1.depth
    covstr.nexp += covstr1.nexp
    covstr.unexp += covstr1.unexp
    covstr.udepth >= covstr.udepth
    covstr.gnexp += covstr1.gnexp
    covstr.gdepth >= covstr.gdepth
    covstr.rnexp += covstr1.rnexp
    covstr.rdepth >= covstr.rdepth
    covstr.inexp += covstr1.inexp
    covstr.idepth >= covstr.idepth
    covstr.znexp += covstr1.znexp
    covstr.zdepth >= covstr.zdepth
    covstr.ynexp += covstr1.ynexp
    covstr.ydepth >= covstr.ydepth
    covstr.vrnexp += covstr1.vrnexp
    covstr.vrdepth >= covstr.vrdepth

    ;stop
  endelse
endfor
print,nmeas

outfile = '/dl1/users/dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_calibrate_coverage.fits'
print,'Writing calib coverage information to ',outfile
MWRFITS,covstr,outfile,/create

stop

end
