pro make_nsc_photscatter_mag

  ;; Make histogram of the photometric scatter for one Healpix and filter
  ; I ran this on gp09
  
  pix = 181040
  objfile = '/dl1/users/dnidever/nsc/instcal/v2/combine/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits.gz'

  setdisp
  !p.font = 0

  meta = mrdfits(objfile,1,/silent)
  allobj = mrdfits(objfile,2,/silent)
  gdobj = where(allobj.fwhm lt 1.5 and allobj.class_star gt 0.5,ngdobj)
  allobj1 = allobj[gdobj]
  tags = tag_names(allobj)

  ;filters = ['u','g','r','i','z','y','vr']
  filters = 'g'
  nfilters = n_elements(filters)

  ; Filter loop
  for j=0,nfilters-1 do begin
    minmed = 99.99
    maglow = 99.99  
    magind = where(tags eq strupcase(filters[j])+'MAG',nmagind)
    scatind = where(tags eq strupcase(filters[j])+'RMS',nscatind)
    detind = where(tags eq 'NPHOT'+strupcase(filters[j]),ndetind)
    ;if str.nexp[j] lt 3 then goto,BOMB
    ; Need at least 3 detections to get a good RMS
    gdndet = where(allobj1.(detind) gt 2,ngdndet)
    if ngdndet lt 5 then goto,BOmB
    allobj2 = allobj1[gdndet]
    ;str.nexp[j] = median([allobj2.(detind)])
    mag = allobj2.(magind)
    scatter = allobj2.(scatind)
    bindata,mag,scatter,bin=0.2,xbin,med,min=15,max=20,/med,gdind=gdind
    bindata,mag,scatter,bin=0.2,xbin,num,min=15,max=20,/hist
    xbin = xbin[gdind]
    med = med[gdind]
    ; sometimes the values are still 99.99
    gdind2 = where(med lt 10,ngdind2)
    xbin = xbin[gdind2]
    med = med[gdind2]
    ; smooth
    smmed = gsmooth(med,2)
    ; bright star, low scatter region
    minmed0 = min(med)
    minind0 = first_el(minloc(med))
    gdlow = where(med lt 2*minmed0 and abs(xbin-xbin[minind0]) lt 2.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 2*minmed0 and abs(xbin-xbin[minind0]) lt 3.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 4*minmed0 and abs(xbin-xbin[minind0]) lt 4.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 5*minmed0 and abs(xbin-xbin[minind0]) lt 5.0,ngdlow)
    if ngdlow lt 2 then gdlow = where(med lt 6*minmed0 and abs(xbin-xbin[minind0]) lt 6.0,ngdlow)
    if ngdlow lt 2 then begin
      gdlow1 = where(abs(xbin-xbin[minind0]) lt 3.0,ngdlow1)
      si = sort(med[gdlow1])
      gdlow = gdlow1[si[0:2]]  ; take the lowest 3 points
    endif
    undefine,coef0,coef1,coef2
    coef0 = median([med[gdlow]])
    coef1 = robust_poly_fitq(xbin[gdlow],med[gdlow],1)
    if ngdlow gt 2 then begin
      coef2 = robust_poly_fitq(xbin[gdlow],med[gdlow],2)
      coef = coef2
    endif else coef=coef1
    minmed = min(med[gdlow])
    minind = gdlow[first_el(minloc(med[gdlow]))]
    maglow = xbin[minind]
    minfitmed = min(poly(xbin[gdlow],coef))
    ;str.minscatter[j] = minmed
    ;str.maglow[j] = maglow
    ;str.minfitscatter[j] = minfitmed
    print,' ',filters[j],'  ',minmed,' ',maglow

    file = 'nsc_photscatter_mag_'+strtrim(pix,2)+'_'+filters[j]
    ps_open,file,/color,thick=4,/encap

    ; select only stars with NGDET>2
    plot,allobj2.gmag,allobj2.grms,ps=8,sym=0.5,xr=[14,25],yr=[0,0.1],xs=1,ys=1,xtit='g',ytit='g RMS (mag)',$
         tit='RMS in g for HEALPix='+strtrim(pix,2),charsize=1.8
    bindata,allobj2.gmag,allobj2.grms,xbin,ybin,bin=0.5,/med,min=15,max=22.5
    oplot,xbin,ybin,co=250,thick=5
    minbin = min(ybin)
    oplot,[0,30],[0,0]+minbin,linestyle=2,co=80
    xyouts,24,0.007,stringize(minbin,ndec=4),align=0.5,charsize=1.6,co=80
    ;xyouts,24,0.007,stringize(minbin,ndec=4)+' mag',align=0.5,charsize=1.2,co=80
    
    ps_close
    ps2png,file+'.eps',/eps
    BOMB:
  endfor


  

  stop

  end
