pro get_lmc_data,pix,redo=redo

;;  Make density map of LMC periphery with NSC data

dir = '/dl1/users/dnidever/nsc/instcal/v2/'
outdir = '/dl1/users/dnidever/nsc/lmc/'

if n_elements(pix) eq 0 then begin
  print,'Syntax - get_lmc_data,pix,redo=redo'
  return
endif

objfile = dir+'combine/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits.gz'

outfile = outdir+strtrim(pix,2)+'_lmc.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' already EXISTS and /redo NOT set'
  return
endif

if file_test(objfile) eq 0 then begin
  print,objfile,' NOT FOUND'
  return
endif
obj = mrdfits(objfile,2,/silent)
gg = where(obj.nphotg gt 0 and obj.nphotr gt 0 and obj.fwhm lt 1.5,ngg)
print,strtrim(pix,2),' ',strtrim(ngg,2)
;schema = {id:'',ra:0.0d0,dec:0.0d0,pa:0.0d0,rad:0.0d0,nphotg:0,gmag:99.99,gerr:9.99,$
schema = {id:'',ra:0.0d0,dec:0.0d0,mlon:0.0d0,mlat:0.0d0,nphotg:0,gmag:99.99,gerr:9.99,$
          nphotr:0,rmag:99.99,rerr:9.99,fwhm:99.99,class_star:99.9,ebv:99.99}
new = replicate(schema,ngg)
struct_assign,obj[gg],new
new.id = strtrim(new.id,2)
cel2lmc,new.ra,new.dec,pa,rad
glactc,new.ra,new.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
new.mlon = mlon
new.mlat = mlat
; save
mwrfits,new,outfile,/create

;stop

end
