pro calt3a

;; Calibrate exposure files for Healpix 148487 to t3a

meta = mrdfits('/dl1/users/dnidever/nsc/instcal/v2/combine/148/148487.fits.gz',1)
meta.base = strtrim(meta.base,2)
expstr = mrdfits('/dl1/users/dnidever/nsc/instcal/v2/lists/nscdr1_exposures.fits',1)
expstr.base = strtrim(expstr.base,2)
MATCH,expstr.base,meta.base,ind1,ind2,/sort
meta = meta[ind1]
expstr = expstr[ind2]
expstr.expdir = strtrim(expstr.expdir,2)

for i=0,n_elements(meta)-1 do begin

  expdir = repstr(expstr[i].expdir,'/v2/','/t3a/')
  nsc_instcal_calibrate_v3,expdir,/redo

  ;stop

endfor


stop

end
