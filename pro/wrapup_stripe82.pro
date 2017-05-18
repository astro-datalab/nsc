pro wrapup_stripe82

; Wrap-up the object and source catalogs
; for RA=220, DEC=0 region

str = mrdfits('/dl1/users/dnidever/nsc/instcal/combine/nsc_instcal_combine.fits',1)
g = where(str.ra ge 218 and str.ra lt 223 and str.dec gt -3 and str.dec lt 3 and str.nexposures gt 800,ng)
print,strtrim(ng,2),' healpix pixels'

undefine,allobj,allmeta,allsrc
for i=0,ng-1 do begin
  print,strtrim(i+1,2),' ',str[g[i]].pix
  objfile = '/dl1/users/dnidever/nsc/instcal/combine/'+strtrim(str[g[i]].pix,2)+'.fits.gz'
  obj = mrdfits(objfile,2)
  push,allobj,obj
  meta = mrdfits(objfile,1)
  push,allmeta,meta
  srcfile = '/dl1/users/dnidever/nsc/instcal/combine/source/'+strtrim(str[g[i]].pix,2)+'_source.fits'
  src = mrdfits(srcfile,1)
  push,allsrc,src
endfor

MWRFITS,allobj,'/dl1/users/dnidever/nsc/instcal/combine/stripe82_ra220dec0_object.fits',/create
MWRFITS,allmeta,'/dl1/users/dnidever/nsc/instcal/combine/stripe82_ra220dec0_expmeta.fits',/create
MWRFITS,allsrc,'/dl1/users/dnidever/nsc/instcal/combine/stripe82_ra220dec0_source.fits',/create

stop

end
