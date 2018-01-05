pro pm_comparison_ucac5

; Compare the NSC proper motions to the HSOY ones

;pmdir = '/dl1/users/dnidever/nsc/instcal/v2/combine/hsoypm/'
pmdir = '/dl1/users/dnidever/nsc/instcal/v2/combine/ucac5pm/'
str = mrdfits(pmdir+'ucac5pm.fits',1)

; These two got messed up
;newobj.ucac_pmdec = ucac2.pmer
;newobj.ucac_epmra = ucac2.pmud
ucac_epmra = str.ucac_pmdec
ucac_pmdec = str.ucac_epmra
str.ucac_epmra = ucac_epmra
str.ucac_pmdec = ucac_pmdec

gd = where((abs(str.pmra)/str.pmraerr gt 3 or str.pmraerr lt 3) and (abs(str.pmdec)/str.pmdecerr gt 3 or str.pmdecerr lt 3) and str.ndet ge 3 and $
           (abs(str.ucac_pmra)/str.ucac_epmra gt 3 or str.ucac_epmra lt 3) and (abs(str.ucac_pmdec)/str.ucac_epmdec gt 3 or str.ucac_epmdec lt 3),ngd)            
str3 = str[gd]

; One-to-one
hess,str3.ucac_pmra,str3.pmra,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log
oplot,[-1000,1000],[-1000,1000],linestyle=2

; Residuals
hess,str3.ucac_pmra,str3.pmra-str3.ucac_pmra,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log
oplot,[-1000,1000],[0,0],linestyle=2

; Hmm quite a few UCAC with PM>+/-90 but NO NSC values that high
; could this be from the 0.5" matching radius.
; The NSC PMRA/PMDEC values max out at ~150-200 mas/yr

setdisp
!p.font = 0
file = pmdir+'plots/pm_ucac_comparison'
ps_open,file,/color,thick=4,/encap
;device,/inches,xsize=12,ysize=7.0
device,/inches,xsize=8.5,ysize=9.0
loadcol,3
black = fsc_color('black',0)
undefine,dum,im
hess,str3.ucac_pmra,str3.pmra,dum,im,dx=2,dy=2,xr=[-150,150],yr=[-150,150],xarr=xarr,yarr=yarr,/noplot
displayc,im,xarr,yarr,/log,xtit='UCAC5 PMRA (mas/yr)',ytit='NSC PMRA (mas/yr)',$
     tit='NSC-UCAC5 Proper Motion Comparison',charsize=1.3,posim=[0.12,0.08,0.98,0.88],poscol=[0.12,0.95,0.98,0.97]

;     position=[0.0,0.0,0.5,1.0],charsize=1.3
setdisp
oplot,[-1000,1000],[-1000,1000],linestyle=2,co=250

; This looks nearly identical
;hess,str3.ucac_pmdec,str3.pmdec,dx=2,dy=2,xr=[-150,150],yr=[-150,150],/log,xtit='UCAC PMDEC (mas/yr)',ytit='NSC PMDEC (mas/yr)',$
;     position=[0.5,0.0,1.0,1.0],charsize=1.3,/noerase
;oplot,[-1000,1000],[-1000,1000],linestyle=2

ps_close
ps2png,file+'.eps',/eps

stop

end
