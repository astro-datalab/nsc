pro nsc_combine_healpix,pix,obj,outfile=outfile,stp=stp

undefine,obj

; Combine multiple NSC healpix catalogs
npix = n_elements(pix)
if npix eq 0 then begin
  print,'Syntax - nsc_combine_healpix,pix,obj,outfile=outfile,stp=stp'
  return
endif

NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/decamcatalog/instcal/'

; Combine
print,'Combining ',strtrim(npix,2),' catalogs'
schema_obj = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,ndet:0L,$
              ndetu:0,nphotu:0,umag:0.0,uerr:0.0,uasemi:0.0,ubsemi:0.0,utheta:0.0,$
              ndetg:0,nphotg:0,gmag:0.0,gerr:0.0,gasemi:0.0,gbsemi:0.0,gtheta:0.0,$
              ndetr:0,nphotr:0,rmag:0.0,rerr:0.0,rasemi:0.0,rbsemi:0.0,rtheta:0.0,$
              ndeti:0,nphoti:0,imag:99.9,ierr:0.0,iasemi:0.0,ibsemi:0.0,itheta:0.0,$
              ndetz:0,nphotz:0,zmag:0.0,zerr:0.0,zasemi:0.0,zbsemi:0.0,ztheta:0.0,$
              ndety:0,nphoty:0,ymag:0.0,yerr:0.0,yasemi:0.0,ybsemi:0.0,ytheta:0.0,$
              ndetvr:0,nphotvr:0,vrmag:0.0,vrerr:0.0,vrasemi:0.0,vrbsemi:0.0,vrtheta:0.0,$
              x2:0.0,x2err:0.0,y2:0.0,y2err:0.0,xy:0.0,xyerr:0.0,cxx:0.0,cxxerr:0.0,$
              cxy:0.0,cxyerr:0.0,cyy:0.0,cyyerr:0.0,asemi:0.0,asemierr:0.0,bsemi:0.0,$
              bsemierr:0.0,theta:0.0,thetaerr:0.0,elongation:0.0,$
              ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
obj = replicate(schema_obj,1e7)
nobj = n_elements(obj)
cnt = 0LL
for i=0,npix-1 do begin
  print,strtrim(i+1,2),' ',pix[i]
  file = dir+'combine/'+strtrim(pix[i],2)+'.fits'
  if file_test(file) eq 1 then begin
    cat1 = MRDFITS(file,1,/silent)
    ncat1 = n_elements(cat1)
    print,'  ',strtrim(ncat1,2),' sources'
    cat1.id = strtrim(pix[i],2)+'.'+strtrim(cat1.id,2)

    ; Add new elements
    if cnt+ncat1 gt nobj then begin
      print,'Adding more elements'
      orig = obj
      obj = replicate(schema_obj,nobj+5e6)
      obj[0:nobj-1] = orig
      nobj = n_elements(obj)
      undefine,orig
    endif

    ; Stuff information in
    obj[cnt:cnt+ncat1-1] = cat1
    cnt += ncat1
    ;push,obj,cat1
  endif else begin
    print,file,' NOT FOUND'
  endelse
  ;stop
endfor
; Trim excess elements
obj = obj[0:cnt-1]
nobj = n_elements(obj)
print,strtrim(nobj,2),' objects'

; Write combined catalog
if n_elements(outfile) gt 0 then begin
  print,'Writing combined catalog to ',outfile
  MWRFITS,obj,outfile,/create

  ; Make figure
  setdisp
  !p.font = 0
  ; On sky distribution
  file = file_dirname(outfile)+'/'+file_basename(outfile,'.fits')+'_skycounts'
  ps_open,file,/color,thick=4,/encap
  hess,obj.ra,obj.dec,dum,dx=0.1,dy=0.1,/log,xtit='RA',ytit='DEC',tit='NOAO Source Catalog on sky distribution'
  ps_close
  ps2png,file+'.eps',/eps

  ; CMD
  file = file_dirname(outfile)+'/'+file_basename(outfile,'.fits')+'_gihess'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=7.5,ysize=9.5
  gd = where(obj.gmag lt 50 and obj.imag lt 50,ngd)
  hess,obj[gd].gmag-obj[gd].imag,obj[gd].gmag,dum,im,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[24,10],/log,/noplot,xarr=xarr,yarr=yarr
  displayc,im,xarr,yarr,xtit='g-i',ytit='i',tit='NOAO Source Catalog CMD',charsize=1.3,/yflip,/log
  ps_close
  ps2png,file+'.eps',/eps
endif


if keyword_set(stp) then stop

end
