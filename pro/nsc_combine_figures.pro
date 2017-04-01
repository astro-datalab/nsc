pro oplot_grid,nlon=nlon,nlat=nlat,thick=thick,linestyle=linestyle,color=color,aitoff=aitoff
if n_elements(nlon) eq 0 then nlon = 7 ;13  ; 30 deg
if n_elements(nlat) eq 0 then nlat = 7   ; 30 deg
dlon = 360/(nlon-1)
dlat = 180/(nlat-1)
if n_elements(thick) eq 0 then thick=2
if n_elements(linestyle) eq 0 then linestyle=1
if n_elements(color) eq 0 then if !d.name eq 'X' then color=255 else color=0

; Constant Longitude lines
for i=0,nlon-1 do begin
  if keyword_set(aitoff) then begin
    aitoff,fltarr(180)+i*dlon-180,findgen(180)-90,x,y
    oplot,x,y,thick=thick,linestyle=linestyle,color=color
  endif else begin
    oplot,fltarr(180)+i*dlon-180,findgen(180)-90,thick=thick,linestyle=linestyle,color=color
  endelse
endfor
; Constant Latitude lines
for i=0,nlat-1 do begin
  if keyword_set(aitoff) then begin
    aitoff,findgen(360)-180,fltarr(360)+i*dlat-90,x,y
    oplot,x,y,thick=thick,linestyle=linestyle,color=color
  endif else begin
    oplot,findgen(360)-180,fltarr(360)+i*dlat-90,thick=thick,linestyle=linestyle,color=color
  endelse
endfor

end

;---------

pro nsc_combine_figures

; Make summary figures for the NSC data

psdir = '/dl1/users/dnidever/nsc/instcal/combine/plots/'

sumstr = mrdfits('/dl1/users/dnidever/nsc/instcal/combine/nsc_instcal_combine.fits',1)

gd = where(sumstr.exists eq 1 and sumstr.success eq 1,ngd)
sumstr2 = sumstr[gd]
glactc,sumstr2.ra,sumstr2.dec,2000.0,glon,glat,1,/deg
aitoff,sumstr2.ra,sumstr2.dec,xequi,yequi
aitoff,glon,glat,xgal,ygal

setdisp
!p.font = 0

; Map of ALL objects

;file = psdir+'nsc_combine_equimap_nobjects'
;ps_open,file,/color,thick=4,/encap
;device,/inches,xsize=14.5,ysize=8.5
;plotc,xequi,yequi,sumstr2.nobjects,ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
;      xs=5,ys=5,tit='All Objects',/log,max=1e6,format='(F8.0)'
;xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
;xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
;oplot_grid,/aitoff
;ps_close
;ps2png,file+'.eps',/eps
;;spawn,['epstopdf',file+'.eps'],/noshell

;file = psdir+'nsc_combine_galmap_nobjects'
;ps_open,file,/color,thick=4,/encap
;device,/inches,xsize=14.5,ysize=8.5
;plotc,xgal,ygal,sumstr2.nobjects,ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
;      xs=5,ys=5,tit='All Objects',/log,max=1e6,format='(F8.0)'
;xyouts,0.52,0.03,'Galactic Longitude',align=0.5,/norm,color=0,charsize=1.4
;xyouts,0.05,0.45,'Galactic Latitude',align=0.5,orientation=90,/norm,color=0,charsize=1.4
;oplot_grid,/aitoff
;ps_close
;ps2png,file+'.eps',/eps
;;spawn,['epstopdf',file+'.eps'],/noshell

; Map of objects in ONE filter
filters = ['u','g','r','i','z','Y','VR']
nfilters = n_elements(filters)
for i=0,nfilters-1 do begin

  maxnobj = max(sumstr2.nobjfilt[i])
  maxnobj_pow = floor(alog10(maxnobj))
  maxnobj = round(maxnobj/10.^(maxnobj_pow-1)) * 10.^maxnobj_pow

;  ; All objects in this filter
;  gd = where(sumstr2.nobjfilt[i] gt 0,ngd)
;  file = psdir+'nsc_combine_equimap_'+filters[i]+'nobjects'
;  ps_open,file,/color,thick=4,/encap
;  device,/inches,xsize=14.5,ysize=8.5
;  plotc,xequi[gd],yequi[gd],sumstr2[gd].nobjfilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
;        xs=5,ys=5,tit=filters[i]+'-band Objects',/log,max=maxnobj,format='(F8.0)'
;  xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
;  xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
;  oplot_grid,/aitoff
;  ps_close
;  ps2png,file+'.eps',/eps
;  ;spawn,['epstopdf',file+'.eps'],/noshell

;  file = psdir+'nsc_combine_galmap_'+filters[i]+'nobjects'
;  ps_open,file,/color,thick=4,/encap
;  device,/inches,xsize=14.5,ysize=8.5
;  plotc,xgal[gd],ygal[gd],sumstr2[gd].nobjfilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
;        xs=5,ys=5,tit=filters[i]+'-band Objects',/log,max=maxnobj,format='(F8.0)'
;  xyouts,0.52,0.03,'Galactic Longitude',align=0.5,/norm,color=0,charsize=1.4
;  xyouts,0.05,0.45,'Galactic Latitude',align=0.5,orientation=90,/norm,color=0,charsize=1.4
;  oplot_grid,/aitoff
;  ps_close
;  ps2png,file+'.eps',/eps
;  ;spawn,['epstopdf',file+'.eps'],/noshell

  mindepth = 17.5
  maxdepth = 23.0

  ; Depth in this filter
  gd = where(sumstr2.depthfilt[i] lt 50,ngd)
  si = sort(sumstr2[gd].depthfilt[i])
  ;mindepth = round(10*sumstr2[gd[si[0.05*ngd]]].depthfilt[i])/10.
  ;maxdepth = round(10*sumstr2[gd[si[0.95*ngd]]].depthfilt[i])/10.
  file = psdir+'nsc_combine_equimap_'+filters[i]+'depth'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=14.5,ysize=8.5
  plotc,xequi[gd],yequi[gd],sumstr2[gd].depthfilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
        xs=5,ys=5,tit=filters[i]+'-band 10 sigma Depth',min=mindepth,max=maxdepth
  xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
  xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
  oplot_grid,/aitoff
  ps_close
  ps2png,file+'.eps',/eps
  ;spawn,['epstopdf',file+'.eps'],/noshell

  file = psdir+'nsc_combine_galmap_'+filters[i]+'depth'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=14.5,ysize=8.5
  plotc,xgal[gd],ygal[gd],sumstr2[gd].depthfilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
        xs=5,ys=5,tit=filters[i]+'-band 10 sigma Depth',min=mindepth,max=maxdepth
  xyouts,0.52,0.03,'Galactic Longitude',align=0.5,/norm,color=0,charsize=1.4
  xyouts,0.05,0.45,'Galactic Latitude',align=0.5,orientation=90,/norm,color=0,charsize=1.4
  oplot_grid,/aitoff
  ps_close
  ps2png,file+'.eps',/eps
  ;spawn,['epstopdf',file+'.eps'],/noshell

  ;; Total exptime in this filter
  ;file = psdir+'nsc_combine_equimap_'+filters[i]+'nobjects'
  ;ps_open,file,/color,thick=4,/encap
  ;device,/inches,xsize=14.5,ysize=8.5
  ;plotc,xequi,yequi,sumstr2.exptimefilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
  ;      xs=5,ys=5,tit=filters[i]+'-band Total Exposure Time',/log,max=1e6,format='(F8.0)'
  ;xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
  ;xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
  ;oplot_grid,/aitoff
  ;ps_close
  ;ps2png,file+'.eps',/eps
  ;spawn,['epstopdf',file+'.eps'],/noshell
  ;
  ;file = psdir+'nsc_combine_galmap_'+filters[i]+'nobjects'
  ;ps_open,file,/color,thick=4,/encap
  ;device,/inches,xsize=14.5,ysize=8.5
  ;plotc,xgal,ygal,sumstr2.nobjfilt[i],ps=8,sym=0.3,xr=[180,-180],yr=[-90,90],$
  ;      xs=5,ys=5,tit=filters[i]+'-band Total Exposure Time',/log,max=1e6,format='(F8.0)'
  ;xyouts,0.52,0.03,'Galactic Longitude',align=0.5,/norm,color=0,charsize=1.4
  ;xyouts,0.05,0.45,'Galactic Latitude',align=0.5,orientation=90,/norm,color=0,charsize=1.4
  ;oplot_grid,/aitoff
  ;ps_close
  ;ps2png,file+'.eps',/eps
  ;spawn,['epstopdf',file+'.eps'],/noshell

  ;stop

endfor

stop

; All galaxies plot
allgal = mrdfits('/dl1/users/dnidever/nsc/instcal/combine/gal/nsc_combine_gal_allcat.fits',1)
glactc,allgal.ra,allgal.dec,2000.0,glon,glat,1,/deg
aitoff,allgal.ra,allgal.dec,xequi,yequi
aitoff,glon,glat,xgal,ygal

file = psdir+'nsc_combine_equimap_galsmag'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=14.5,ysize=8.5
plotc,xequi,yequi,allgal.mag,ps=3,xr=[180,-180],yr=[-90,90],$
      xs=5,ys=5,tit='All Galaxies (color-coded by magnitude)'
xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
oplot_grid,/aitoff
ps_close
ps2png,file+'.eps',/eps
;spawn,['epstopdf',file+'.eps'],/noshell

file = psdir+'nsc_combine_galmap_galsmag'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=14.5,ysize=8.5
plotc,xgal,ygal,allgal.mag,ps=3,xr=[180,-180],yr=[-90,90],$
      xs=5,ys=5,tit='All Galaxies (color-coed by magnitude)'
xyouts,0.52,0.03,'Galactic Longitude',align=0.5,/norm,color=0,charsize=1.4
xyouts,0.05,0.45,'Galactic Latitude',align=0.5,orientation=90,/norm,color=0,charsize=1.4
oplot_grid,/aitoff
ps_close
ps2png,file+'.eps',/eps
;spawn,['epstopdf',file+'.eps'],/noshell

; Galaxies, FWHM

file = psdir+'nsc_combine_equimap_galsfwhm'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=14.5,ysize=8.5
plotc,xequi,yequi,allgal.fwhm,ps=3,xr=[180,-180],yr=[-90,90],$
      xs=5,ys=5,tit='All Galaxies (color-coded by FWHM)'
xyouts,0.52,0.03,'RA',align=0.5,/norm,color=0,charsize=1.4
xyouts,0.05,0.45,'DEC',align=0.5,orientation=90,/norm,color=0,charsize=1.4
oplot_grid,/aitoff
ps_close
ps2png,file+'.eps',/eps
;spawn,['epstopdf',file+'.eps'],/noshell


stop

end
