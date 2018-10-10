pro make_astromrms_figure_multiexposure,src

  ;; Make figures of astrometric RMS (over multiple detection) of
  ;; stars vs. S/N for one Healpix

  ;; ran this on gp09
  dir = '/dl1/users/dnidever/nsc/instcal/v2/'
  
  if n_elements(src) eq 0 then begin
     pix = 98233  ;148487
     src = mrdfits(dir+'combine/source/'+strtrim(pix,2)+'_source.fits.gz',1)
     src.objectid = strtrim(src.objectid,2)
  endif

  index = create_index(src.objectid)
  n = n_elements(index.num)
  rms = dblarr(n)
  decsig = dblarr(n)
  rasig = dblarr(n)
  snr = fltarr(n)
  mag = fltarr(n)
  fwhm = fltarr(n)
  decerr = fltarr(n)

  for i=0,n-1 do begin
    ind = index.index[index.lo[i]:index.hi[i]]
    rms[i] = stddev(src[ind].dec)*3600
    decsig[i] = mad(src[ind].dec)*3600
    rasig[i] = mad(src[ind].ra)*3600.
    snr[i] = median([1.087/src[ind].cerr])
    mag[i] = median([src[ind].cmag])
    fwhm[i] = median([src[ind].fwhm_world*3600])
    ;decerr[i] = median([src[ind].decerr])
    decerr[i] = median([0.644*src[ind].fwhm_world*3600.*src[ind].cerr/1.087])
  endfor

  setdisp
  !p.font = 0
  
  ;; RMS vs. S/N
  file = 'decrms_snr_'+strtrim(pix,2)
  ps_open,file,/color,thick=4,/encap
  gd = where(index.num gt 100 and decsig gt 0.0 and fwhm lt 1.5,ngd)
  plot,snr[gd],decsig[gd]*1e3,ps=8,sym=0.5,xr=[4,400],yr=[0,100],xs=1,ys=1,/xlog,xtit='S/N',ytit='RMS in Declination (mas)',tit='RMS vs. S/N',charsize=1.4
  ;g = where(decsig lt 0.10 and decsig gt 0.0 and fwhm lt 1.1,ng)
  bindata,alog10(snr[gd]),decsig[gd]*1e3,xbin,ybin,binsize=0.1,/med,min=0.7,max=2.44
  oplot,10^xbin,ybin,co=250,thick=8
  al_legend,['HEALPix='+strtrim(pix,2),'N!dexposures!n > 100'],/top,/right,box=0,charsize=1.5
  ps_close
  ps2png,file+'.eps',/eps
  
  stop

  end
