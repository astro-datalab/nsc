function superfit,x,par,_extra=fa

; Superposition model of magnitudes/EBV and a constant term
; FA structure must have X1, X2, X3 .. vectors
; of the magnitudes and EBV.
tags = tag_names(fa)
npar = n_elements(par)
nvec = npar-1

model = x*0
for i=0,nvec-1 do begin
  ind = where(tags eq 'X'+strtrim(i+1,2),nind)
  model += par[i] * fa.(ind)
endfor
; Add the constant term
model += par[npar-1]

return,model
end

;--------

pro nsc_color_relations_smashuband_v3

; Make color-color relations plots for u-band using SMASH

;restore,'~/datalab/nsc/smash_matched_catalog.dat'
restore,'/dl1/users/dnidever/nsc/smash_matched_catalog_v3.dat'

gd1 = where(mobj.tmass_qflg eq 'AAA' and mobj.gaia_gmag lt 18 and mobj.tmass_jmag lt 15 and $
            mobj.uerr lt 0.1 and mobj.galex_nuverr lt 0.1,ngd1)
mobj1 = mobj[gd1]

;gd = where(mobj.tmass_qflg eq 'AAA' and mobj.gaia_gmag lt 15 and mobj.tmass_jmag lt 14 and $
;           mobj.uerr lt 0.1 and mobj.galex_nuverr lt 0.1 and mobj.ebv lt 0.05,ngd)
gd = where(mobj.tmass_qflg eq 'AAA' and mobj.gaia_gmag lt 15 and mobj.tmass_jmag lt 14 and $
           mobj.uerr lt 0.1 and mobj.galex_nuverr lt 0.1 and mobj.ebv lt 0.1,ngd)
mobj2 = mobj[gd]

setdisp
!p.font = 0
plotdir = '/dl1/users/dnidever/nsc/instcal/t3b/plots/'
if file_test(plotdir,/directory) eq 0 then file_mkdir,plotdir

;; G-J vs. EBV
file = plotdir+'nsc_color_relations_gj_ebv_v3'
ps_open,file,/color,thick=4,/encap
hess,mobj1.ebv,mobj1.gaia_gmag-mobj1.tmass_jmag,dx=0.02,dy=0.05,xr=[0.0,0.7],yr=[-1,3],xtit='EBV',ytit='GMAG-JMAG',/log
bindata,mobj1.ebv,mobj1.gaia_gmag-mobj1.tmass_jmag,xbin,ybin,binsize=0.05,/med
oplot,xbin,ybin,ps=-1,sym=2,co=255
;err = sqrt(mobj1.gaia_gerr^2 + mobj1.tmass_jerr^2) > 0.02
;gg = where(mobj1.ebv gt 0.15 and mobj1.ebv lt 1.0 and abs(mobj1.gaia_gmag-mobj1.tmass_jmag) lt 3,ngg)
;coef = dln_poly_fit(mobj1[gg].ebv,mobj1[gg].gaia_gmag-mobj2[gg].tmass_jmag,1,$
;                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'G-J vs. EBV'
;print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
coef = [0.98,1.12]
x = scale_vector(findgen(100),0,1)
oplot,x,poly(x,coef),co=250
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag',$
;           'Slope = '+stringize(coef[1],ndec=4)+' +/- '+stringize(coeferr[1],ndec=4)],textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; NUV-u vs. G-J
file = plotdir+'nsc_color_relations_nuvu_gj_v3'
ps_open,file,/color,thick=4,/encap
hess,mobj2.gaia_gmag-mobj2.tmass_jmag,mobj2.galex_nuv-mobj2.u,dx=0.02,dy=0.05,xr=[0.0,2.5],yr=[0.2,7],xtit='G-J',ytit='NUV-u',/log
err = sqrt(mobj2.galex_nuverr^2 + mobj2.uerr^2) > 0.02
gg = where(mobj2.gaia_gmag-mobj2.tmass_jmag ge 0.80 and mobj2.gaia_gmag-mobj2.tmass_jmag le 1.2,ngg)
coef = dln_poly_fit(mobj2[gg].gaia_gmag-mobj2[gg].tmass_jmag,mobj2[gg].galex_nuv-mobj2[gg].u,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'NUV-u vs. G-J'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.8,1.2)
oplot,x,poly(x,coef),co=0
al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag',$
           'Slope = '+stringize(coef[1],ndec=4)+' +/- '+stringize(coeferr[1],ndec=4)],textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Definitely sensitive to REDDENING


; <NUV+G>-u vs. G-J
file = plotdir+'nsc_color_relations_nuvgu_gj_v3'
ps_open,file,/color,thick=4,/encap
f = 0.34
hess,mobj2.gaia_gmag-mobj2.tmass_jmag,(f*mobj2.galex_nuv+(1-f)*mobj2.gaia_gmag)-mobj2.u,dx=0.02,dy=0.02,$
     xr=[0.0,2.5],yr=[-1,1],xtit='G-J',ytit='('+stringize(f,ndec=2)+'*NUV+'+stringize(1-f,ndec=2)+'*G)-u',/log
err = sqrt(mobj2.galex_nuverr^2 + mobj2.gaia_gerr^2 + mobj2.uerr^2) > 0.02
gg = where(mobj2.gaia_gmag-mobj2.tmass_jmag ge 0.70 and mobj2.gaia_gmag-mobj2.tmass_jmag le 1.1,ngg)
coef = dln_poly_fit(mobj2[gg].gaia_gmag-mobj2[gg].tmass_jmag,(f*mobj2[gg].galex_nuv+(1-f)*mobj2[gg].gaia_gmag)-mobj2[gg].u,0,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'('+stringize(f,ndec=2)+'*NUV+'+stringize(1-f,ndec=2)+'*G)-u vs. G-J'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.8,1.1)
oplot,x,poly(x,coef),co=0
al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag'],textcolor=[250],/top,/left,charsize=1.2
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag',$
;           'Slope = '+stringize(coef[1],ndec=4)+' +/- '+stringize(coeferr[1],ndec=4)],textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; (G-J)o = G-J-1.12*EBV
; (J-Ks)o = J-Ks-0.17*EBV
; this seems to work well by fitting line to main locus of points in
; G-J vs. EBV and J-Ks vs. EBV

; Solve linear equation, Ax=b
gg = where(mobj2.gaia_gmag-mobj2.tmass_jmag ge 0.70 and mobj2.gaia_gmag-mobj2.tmass_jmag le 1.1,ngg)
a = dblarr(2,ngg)
a[0,*] = mobj2[gg].galex_nuv
a[1,*] = mobj2[gg].gaia_gmag
b = mobj2[gg].u
SVDC, A, W, U, V
factor = SVSOL(U, W, V, B)
print,factor
;  0.30964776      0.70183440
; NOW, with extinction
gg1 = where(mobj1.gaia_gmag-mobj1.tmass_jmag ge 0.70 and mobj.gaia_gmag-mobj1.tmass_jmag le 1.1,ngg1)
a = dblarr(3,ngg1)
a[0,*] = mobj1[gg1].galex_nuv
a[1,*] = mobj1[gg1].gaia_gmag
a[2,*] = mobj1[gg1].ebv
b = mobj1[gg1].u
SVDC, A, W, U, V
factor = SVSOL(U, W, V, B)
print,factor
;    0.30436899      0.70834828      0.38482837
; very similar to what I got by trial-and-error
; should remove outliers/high-error

; <NUV+G>-u vs. G-J with E(B-V) correction
;  BUT WE NEED TO ADD THE EXTINCTION BACK IN FOR THE ZERO-POINT OR DO WE!!??
;  NO, NO, NO.  WE can reproduce the CALIBRATED u-band photometry
;  with this scheme that includes EBV.

file = plotdir+'nsc_color_relations_nuvgu_gj_ebv_v3'
ps_open,file,/color,thick=4,/encap
gd1 = where(mobj.tmass_qflg eq 'AAA' and mobj.gaia_gmag lt 15 and mobj.tmass_jmag lt 14 and $
            mobj.uerr lt 0.1 and mobj.galex_nuverr lt 0.1,ngd1)
mobj1 = mobj[gd]
f = 0.34
hess,mobj1.gaia_gmag-mobj1.tmass_jmag,(f*mobj1.galex_nuv+(1-f)*mobj1.gaia_gmag+0.4*mobj1.ebv)-mobj1.u,dx=0.02,dy=0.02,$
     xr=[0.0,2.5],yr=[-1,1],xtit='G-J',ytit='('+stringize(f,ndec=2)+'*NUV+'+stringize(1-f,ndec=2)+'*G)-u',/log
err = sqrt(mobj1.galex_nuverr^2 + mobj1.gaia_gerr^2 + mobj1.uerr^2) > 0.02
gg = where(mobj1.gaia_gmag-mobj1.tmass_jmag ge 0.70 and mobj1.gaia_gmag-mobj1.tmass_jmag le 1.1,ngg)
coef = dln_poly_fit(mobj1[gg].gaia_gmag-mobj1[gg].tmass_jmag,(f*mobj1[gg].galex_nuv+(1-f)*mobj1[gg].gaia_gmag+0.5*mobj1[gg].ebv)-mobj1[gg].u,0,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'('+stringize(f,ndec=2)+'*NUV+'+stringize(1-f,ndec=2)+'*G)-u vs. G-J'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.8,1.1)
oplot,x,poly(x,coef),co=0
al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag'],textcolor=[250],/top,/left,charsize=1.2
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag',$
;           'Slope = '+stringize(coef[1],ndec=4)+' +/- '+stringize(coeferr[1],ndec=4)],textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; NUV-u vs. J-Ks
file = plotdir+'nsc_color_relations_nuvu_jk_v3'
ps_open,file,/color,thick=4,/encap
hess,mobj2.tmass_jmag-mobj2.tmass_kmag,mobj2.galex_nuv-mobj2.u,dx=0.02,dy=0.05,xr=[-0.5,2.0],yr=[0.2,7],xtit='J-Ks',ytit='NUV-u',/log
err = sqrt(mobj2.galex_nuverr^2 + mobj2.uerr^2) > 0.02
gg = where(mobj2.tmass_jmag-mobj2.tmass_kmag ge 0.20 and mobj2.tmass_jmag-mobj2.tmass_kmag le 0.8,ngg)
coef = dln_poly_fit(mobj2[gg].tmass_jmag-mobj2[gg].tmass_kmag,mobj2[gg].galex_nuv-mobj2[gg].u,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'NUV-u vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.2,0.8)
oplot,x,poly(x,coef),co=0
al_legend,['Zero-point = '+stringize(coef[0],ndec=4)+' +/- '+stringize(coeferr[0],ndec=4)+' mag',$
           'Slope = '+stringize(coef[1],ndec=4)+' +/- '+stringize(coeferr[1],ndec=4)],textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; NOPE, too vertical

; More general case
;-------------------
gd1 = where(mobj.tmass_qflg eq 'AAA' and mobj.gaia_gmag lt 15 and mobj.tmass_jmag lt 14 and $
            mobj.uerr lt 0.1 and mobj.galex_nuverr lt 0.1,ngd1)
mobj1 = mobj[gd1]
; (G-J)o = G-J-1.12*EBV
gj0 = mobj1.gaia_gmag-mobj1.tmass_jmag-1.12*mobj1.ebv
gg = where(gj0 ge 0.7 and gj0 le 1.1 and mobj1.uerr gt 0,ngg)
;initpar = [0.0, 0.0, 0.0, 0.0, 0.0]
initpar = [0.2469, 0.7501, 0.5462, 0.35, 0.0052]
;initpar = [0.2469, 0.7501, 0.5462, 0.6809, 0.0052]   ; v2
;initpar = [0.23882, 0.754364, 0.407116, 0.489842, 0.126978]
;initpar = [0.0, 0.0, 0.0, 0.0]
;fa = {x1:mobj1[gg].galex_nuv,x2:gj0[gg],x3:mobj1[gg].ebv}
;fa = {x1:mobj1[gg].galex_nuv,x2:mobj1[gg].gaia_gmag,x3:mobj1[gg].tmass_kmag,x4:mobj1[gg].ebv}
;fa = {x1:mobj1[gg].galex_nuv,x2:mobj1[gg].gaia_gmag,x3:mobj1[gg].ebv}
fa = {x1:mobj1[gg].galex_nuv,x2:mobj1[gg].gaia_gmag,x3:gj0[gg],x4:mobj1[gg].ebv}
x = findgen(ngg)
y = mobj1[gg].u
err = mobj1[gg].uerr > 0.01
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
parinfo[3].fixed = 1  ; fix the extinction
;initpar[0] = 1.0
;parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
;initpar[2] = 0.0  ;0.04
;parinfo[2].fixed = 1
;initpar[3] = 0.223
;parinfo[3].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
print,'u-band:'
print,par
;     0.308744     0.695461     0.424305    0.0930483
;faall = {x1:mobj1.galex_nuv,x2:gj0,x3:mobj1.ebv}
;faall = {x1:mobj1.galex_nuv,x2:mobj1.tmass_jmag,x3:mobj1.tmass_kmag,x4:mobj1.ebv}
;faall = {x1:mobj1.galex_nuv,x2:mobj1.gaia_gmag,x3:mobj1.ebv}
faall = {x1:mobj1.galex_nuv,x2:mobj1.gaia_gmag,x3:gj0,x4:mobj1.ebv}
yfitall = superfit(mobj1.u*0,par,_extra=faall)

plotc,mobj1.ebv,yfitall-mobj1.u,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250


; Extinction plot
file = plotdir+'nsc_color_relations_super_u_ebv_gj_v3'
ps_open,file,/color,thick=4,/encap
plotc,mobj1[gg].ebv,yfitall[gg]-mobj1[gg].u,gj0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='u-band (color-coded by [G-J]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Scatter plot
file = plotdir+'nsc_color_relations_super_u_gj_scatter_v3'
ps_open,file,/color,thick=4,/encap
plotc,gj0,yfitall-mobj1.u,mobj1.ebv,ps=1,sym=0.5,xr=[-0.5,2.0],yr=[-0.5,0.5],$
      xtit='(G-J)!d0!n',ytit='Residuals',tit='u-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.7,1.1],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*NUV+'+stringize(par[1],ndec=3)+'*G+'+$
          stringize(par[2],ndec=3)+'*GJ0+'+stringize(par[3],ndec=3)+'*E(B-V)+'+$
          stringize(par[4],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Density plot
file = plotdir+'nsc_color_relations_super_u_gj_density_v3'
ps_open,file,/color,thick=4,/encap
hess,gj0,yfitall-mobj1.u,mobj1.ebv,dx=0.02,dy=0.02,xr=[-0.5,2.0],yr=[-0.5,0.5],$
      xtit='(G-J)!d0!n',ytit='Residuals',tit='u-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.7,1.1],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*NUV+'+stringize(par[1],ndec=3)+'*G+'+$
          stringize(par[2],ndec=3)+'*GJ0+'+stringize(par[3],ndec=3)+'*E(B-V)+'+$
          stringize(par[4],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-mobj1[gg].u)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

pdfcombine,plots+'.pdf',plotdir+'smashuband_combine.pdf',/clobber


stop

end
