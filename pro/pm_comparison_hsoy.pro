pro pm_comparison_hsoy

; Compare the NSC proper motions to the HSOY ones

;pmdir = '/dl1/users/dnidever/nsc/instcal/v2/combine/hsoypm/'
pmdir = '/dl1/users/dnidever/nsc/instcal/v2/combine/hsoypm_v2/'
str = mrdfits(pmdir+'hsoypm.fits',1)
; the HSOY proper motion values are in deg/yr
str.hsoy_pmra *= 3600*1e3
str.hsoy_pmdec *= 3600*1e3
str.hsoy_epmra *= 3600*1e3
str.hsoy_epmdec *= 3600*1e3

gd = where(str.gmag lt 19,ngd)
str2 = str[gd]
glactc,str2.ra,str2.dec,2000.0,glon,glat,1,/deg
;gd2 = where(abs(glat) gt 20 and str2.pmraerr lt 10 and str2.pmdecerr lt 10 and str2.ndet gt 10,ngd2)
;gd2 = where(abs(str2.pmra)/str2.pmraerr gt 3 and abs(str2.pmdec)/str2.pmdecerr gt 3 and str2.ndet ge 10,ngd2)
;gd2 = where((abs(str2.pmra)/str2.pmraerr gt 3 or str2.pmraerr lt 2) and (abs(str2.pmdec)/str2.pmdecerr gt 3 or str2.pmdecerr lt 2) and str2.ndet ge 10,ngd2)
gd2 = where((abs(str2.pmra)/str2.pmraerr gt 3 or str2.pmraerr lt 3) and (abs(str2.pmdec)/str2.pmdecerr gt 3 or str2.pmdecerr lt 3) and str2.ndet ge 3,ngd2)
str3 = str2[gd2]
;gd = where((abs(str.pmra)/str.pmraerr gt 3 or str.pmraerr lt 3) and (abs(str.pmdec)/str.pmdecerr gt 3 or str.pmdecerr lt 3) and str.ndet ge 10,ngd)
;gd = where(abs(str.pmra)/str.pmraerr gt 3 and abs(str.pmdec)/str.pmdecerr gt 3 and str.ndet ge 10,ngd)
;str3 = str[gd]

; One-to-one
hess,str3.hsoy_pmra,str3.pmra,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log
oplot,[-1000,1000],[-1000,1000],linestyle=2

; Residuals
hess,str3.hsoy_pmra,str3.pmra-str3.hsoy_pmra,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log
oplot,[-1000,1000],[0,0],linestyle=2

; Hmm quite a few HSOY with PM>+/-90 but NO NSC values that high
; could this be from the 0.5" matching radius.
; The NSC PMRA/PMDEC values max out at ~150-200 mas/yr

setdisp
!p.font = 0
;file = '/dl1/users/dnidever/nsc/instcal/v2/combine/hsoypm/plots/pm_hsoy_comparison'
file = pmdir+'plots/pm_hsoy_comparison'
ps_open,file,/color,thick=4,/encap
;device,/inches,xsize=12,ysize=7.0
device,/inches,xsize=8.5,ysize=9.0
loadcol,3
black = fsc_color('black',0)
undefine,dum,im
hess,str3.hsoy_pmra,str3.pmra,dum,im,dx=2,dy=2,xr=[-150,150],yr=[-150,150],xarr=xarr,yarr=yarr,/noplot
displayc,im,xarr,yarr,/log,xtit='HSOY PMRA (mas/yr)',ytit='NSC PMRA (mas/yr)',$
     tit='NSC-HSOY Proper Motion Comparison',charsize=1.3,posim=[0.12,0.08,0.98,0.88],poscol=[0.12,0.95,0.98,0.97]

;     position=[0.0,0.0,0.5,1.0],charsize=1.3
setdisp
oplot,[-1000,1000],[-1000,1000],linestyle=2,co=250

; This looks nearly identical
;hess,str3.hsoy_pmdec,str3.pmdec,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log,xtit='HSOY PMDEC (mas/yr)',ytit='NSC PMDEC (mas/yr)',$
;     position=[0.5,0.0,1.0,1.0],charsize=1.3,/noerase
;oplot,[-1000,1000],[-1000,1000],linestyle=2

ps_close
ps2png,file+'.eps',/eps

stop

end
