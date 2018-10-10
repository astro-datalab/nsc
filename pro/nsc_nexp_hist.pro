pro nsc_nexp_hist,str

  ;; Histogram of number of exposure per object
  ;; I ran this on gp09

  ;; Use reverse indices
  dir = '/net/dl1/users/dnidever/nsc/instcal/v2/combine/coverage/'
  
  if n_elements(str) eq 0 then str = mrdfits(dir+'nsc_instcal_coverage.fits.gz',1)
  gd = where(str.nexp gt 0,ngd)
  nobj = str[gd].nobj
  nexp = str[gd].nexp
  lognexp = alog10(nexp)

  bin = 0.1
  pixhist = histogram(lognexp,bin=bin,min=0,max=3.3,locations=xhist,reverse_indices=rev,/l64)
  areahist = pixhist * 0.0002049  ; deg^2
  xhist += 0.5*bin
  nhist = n_elements(pixhist)
  
  ;; Loop over histogram elements and fill in the nobj values
  objhist = pixhist*0
  for i=0,nhist-1 do begin
    ;; Set all elements of A that are in the ith bin of H to 0.
    if rev[i] ne rev[i+1] then begin
       ind = rev[rev[i]:rev[i+1]-1]
       objhist[i] = total(nobj[ind])
    endif
  endfor

  ;; Extend to Nexposures=0
  xhist = [-2.0,xhist]
  areahist = [areahist[0],areahist]
  objhist = [objhist[0],objhist]
  
  setdisp
  !p.font = 0
  
  ;; Make the plot
  file = dir+'plots/nsc_nexp_hist'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=12,ysize=6.5
  ;; Area
  plot,[0],/nodata,/xlog,/ylog,xr=[1,2000],yr=[1,5000],xs=1,ys=1,charsize=1.3,$
       position=[0.08,0.10,0.49,0.98],xtit='Number of Exposures',ytit='Area (deg!u2!n)'
  oplot,10^xhist,areahist,ps=10
  xyouts,300,3000,'Area',charsize=1.7,align=0.5
  ;plot,xhist,areahist,ps=10,xtit='Number of Exposures',ytit='Area (deg!u2!n)',xr=[0,3.3],yr=[1e6,max(hist)*1.5],/ylog,xs=1,ys=1,charsize=1.3

  ;; Number of objects
  plot,[0],/nodata,/xlog,/ylog,xr=[1,2000],yr=[1e5,max(objhist)*1.5],xs=1,ys=1,charsize=1.3,$
       position=[0.58,0.10,0.99,0.98],/noerase,xtit='Number of Exposures',ytit='Number of Objects'
  oplot,10^xhist,objhist,ps=10
  xyouts,1500,max(objhist)*0.93,'Number of Objects',charsize=1.7,align=1
  ;plot,xhist,objhist,ps=10,xtit='Number of Exposures',ytit='Number of Objects',xr=[0,3.3],yr=[1e6,max(hist)*1.5],/ylog,xs=1,ys=1,charsize=1.3
  ; use logarithmic scale
  ps_close
  ps2png,file+'.eps',/eps

  ;; Cumulative versions

  ;; Make the plot
  file = dir+'plots/nsc_nexp_cumhist'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=12,ysize=6.5
  ;; Area
  cumareahist = reverse(total(reverse(areahist),/cum))
  plot,[0],/nodata,/xlog,/ylog,xr=[1,2000],yr=[1,40000],xs=1,ys=1,charsize=1.8,$
       position=[0.12,0.12,0.49,0.98],xtit='Number of Exposures',ytit='Cumulative Area (deg!u2!n)'
  oplot,10^xhist,cumareahist,ps=10
  xyouts,300,max(cumareahist)*0.70,'Area',charsize=1.8,align=0.5
  ;plot,xhist,areahist,ps=10,xtit='Number of Exposures',ytit='Area (deg!u2!n)',xr=[0,3.3],yr=[1e6,max(hist)*1.5],/ylog,xs=1,ys=1,charsize=1.3

  ;; Number of objects
  cumobjhist = reverse(total(reverse(objhist),/cum))
  plot,[0],/nodata,/xlog,/ylog,xr=[1,2000],yr=[1e5,max(cumobjhist)*1.3],xs=1,ys=1,charsize=1.8,$
       position=[0.62,0.12,0.99,0.98],/noerase,xtit='Number of Exposures',ytit='Cumulative Number of Objects'
  oplot,10^xhist,cumobjhist,ps=10
  xyouts,1500,max(cumobjhist)*0.70,'Number of Objects',charsize=1.8,align=1
  ;plot,xhist,objhist,ps=10,xtit='Number of Exposures',ytit='Number of Objects',xr=[0,3.3],yr=[1e6,max(hist)*1.5],/ylog,xs=1,ys=1,charsize=1.3
  ; use logarithmic scale
  ps_close
  ps2png,file+'.eps',/eps

  stop

  end
