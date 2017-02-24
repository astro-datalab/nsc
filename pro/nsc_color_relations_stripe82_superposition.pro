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

pro nsc_color_relations_stripe82_superposition

; Derive color-color relations for calibrating

;dir = '/datalab/users/dnidever/smash/cp/red/photred/catalogs/final/v4/'
;obj0 = mrdfits(dir+'Field100_combined_allobj.fits.gz',1)
;xmatch0 = mrdfits(dir+'Field100_combined_allobj_xmatch.fits.gz',1)
;obj0.id = strtrim(obj0.id,2)
;xmatch0.id = strtrim(xmatch0.id,2)
;match,obj0.id,xmatch0.id,ind1,ind2,/sort
;obj = obj0[ind1]
;xmatch = xmatch0[ind2]

str = mrdfits('~/datalab/nsc/Stripe82.fits.gz',1)
gd1 = where(str.ps1_gmag lt 20.0 and str.ps1_ng gt 5 and str.ps1_nr gt 5 and str.ps1_ni gt 5 and str.ps1_nz gt 5 and str.ps1_ny gt 5 and $
            str.ps1_gmag gt 12.0 and str.gaia_gmag ge 14 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 15.0 and $
            str.tmass_kmag lt 15.0,ngd1)
str1 = str[gd1]
gd2 = where(str.ps1_gmag lt 20.0 and str.ps1_ng gt 5 and str.ps1_nr gt 5 and str.ps1_gmag gt 12.0 and $
            str.gaia_gmag ge 14 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 15.0 and $
            str.tmass_kmag lt 15.0 and str.apass_gmag ge 14 and str.apass_gmag lt 18.0 and finite(str.apass_gmag) eq 1 and $
            str.apass_rmag lt 18.0 and finite(str.apass_rmag) eq 1,ngd2)
str2 = str[gd2]
; some surveys have problems on bright end, APASS?


; (G-J)o = G-J-1.12*EBV
; (J-Ks)o = J-Ks-0.17*EBV
; this seems to work well by fitting line to main locus of points in
; G-J vs. EBV and J-Ks vs. EBV


setdisp
!p.font = 0

; APASS_g-g vs. J-Ks
;-------------------
jk0 = str2.tmass_jmag-str2.tmass_kmag-0.17*str2.ebv
gg = where(jk0 ge 0.3 and jk0 le 0.7 and str2.ps1_gerr gt 0,ngg)
;gg = where(jk0 ge 0.3 and jk0 le 0.7 and str2.ps1_gerr gt 0 and str2.ebv gt 0.15,ngg)
initpar = [0.0, 0.0, 0.0, 0.0]
fa = {x1:str2[gg].apass_gmag,x2:jk0[gg],x3:str2[gg].ebv}
;fa = {x1:str2[gg].apass_gmag,x2:str2[gg].tmass_jmag,x3:str2[gg].tmass_kmag,x4:str2[gg].ebv}
x = findgen(ngg)
y = str2[gg].ps1_gmag
err = str2[gg].ps1_gerr
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
initpar[0] = 1.0
parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
;  this term depends a bit on what the color term is
initpar[2] = -0.05 ;-0.083 ; -0.074
parinfo[2].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
; Correct constant offset
medoff = median(yfit-y)
par[3] -= medoff
print,'g-band:'
print,par
;    1.00000    -0.143343   -0.0500000   -0.0138004
faall = {x1:str2.apass_gmag,x2:jk0,x3:str2.ebv}
;faall = {x1:str2.apass_gmag,x2:str2.tmass_jmag,x3:str2.tmass_kmag,x4:str2.ebv}
yfitall = superfit(str2.ps1_gmag*0,par,_extra=faall)

plotc,str2.ebv,yfitall-str2.ps1_gmag,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250

; Scatter plot
file = 'plots/nsc_color_relations_stripe82_super_g_jk_scatter'
ps_open,file,/color,thick=4,/encap
plotc,jk0,yfitall-str2.ps1_gmag,str2.ebv,ps=1,sym=0.5,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='g-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*APASS_g'+stringize(par[1],ndec=3)+'*JK0'+$
          stringize(par[2],ndec=3)+'*E(B-V)'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
;push,plots,file

; Density plot
file = 'plots/nsc_color_relations_stripe82_super_g_jk_density'
ps_open,file,/color,thick=4,/encap
hess,jk0,yfitall-str2.ps1_gmag,str2.ebv,dx=0.02,dy=0.02,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='g-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*APASS_g'+stringize(par[1],ndec=3)+'*JK0'+$
          stringize(par[2],ndec=3)+'*E(B-V)'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-str2[gg].ps1_gmag)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Extinction plot
file = 'plots/nsc_color_relations_stripe82_super_g_ebv_jk'
ps_open,file,/color,thick=4,/encap
plotc,str2[gg].ebv,yfitall[gg]-str2[gg].ps1_gmag,jk0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='g-band (color-coded by [J-Ks]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; APASS_r-r vs. J-Ks
;-------------------
jk0 = str2.tmass_jmag-str2.tmass_kmag-0.17*str2.ebv
gg = where(jk0 ge 0.3 and jk0 le 0.7 and str2.ps1_rerr gt 0,ngg)
initpar = [0.0, 0.0, 0.0, 0.0]
fa = {x1:str2[gg].apass_rmag,x2:jk0[gg],x3:str2[gg].ebv}
;fa = {x1:str2[gg].apass_rmag,x2:str2[gg].tmass_jmag,x3:str2[gg].tmass_kmag,x4:str2[gg].ebv}
x = findgen(ngg)
y = str2[gg].ps1_rmag
err = str2[gg].ps1_rerr
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
initpar[0] = 1.0
parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
initpar[2] = 0.0  ;0.04
parinfo[2].fixed = 1
;initpar[3] = 0.223
;parinfo[3].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
print,'r-band:'
print,par
;  1.00000   0.00740162      0.00000  0.000527566
faall = {x1:str2.apass_rmag,x2:jk0,x3:str2.ebv}
;faall = {x1:str2.apass_rmag,x2:str2.tmass_jmag,x3:str2.tmass_kmag,x4:str2.ebv}
yfitall = superfit(str2.ps1_rmag*0,par,_extra=faall)

plotc,str2.ebv,yfitall-str2.ps1_rmag,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250

; Scatter plot
file = 'plots/nsc_color_relations_stripe82_super_r_jk_scatter'
ps_open,file,/color,thick=4,/encap
plotc,jk0,yfitall-str2.ps1_rmag,str2.ebv,ps=1,sym=0.5,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='r-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*APASS_r+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
;push,plots,file

; Density plot
file = 'plots/nsc_color_relations_stripe82_super_r_jk_density'
ps_open,file,/color,thick=4,/encap
hess,jk0,yfitall-str2.ps1_rmag,str2.ebv,dx=0.02,dy=0.02,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='r-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*APASS_r+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-str2[gg].ps1_rmag)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Extinction plot
file = 'plots/nsc_color_relations_stripe82_super_r_ebv_jk'
ps_open,file,/color,thick=4,/encap
plotc,str2[gg].ebv,yfitall[gg]-str2[gg].ps1_rmag,jk0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='r-band (color-coded by [J-Ks]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; G-i vs. J-Ks
;--------------
jk0 = str1.tmass_jmag-str1.tmass_kmag-0.17*str1.ebv
gg = where(jk0 ge 0.3 and jk0 le 0.7 and str1.ps1_ierr gt 0,ngg)
initpar = [0.0, 0.0, 0.0, 0.0]
fa = {x1:str1[gg].gaia_gmag,x2:jk0[gg],x3:str1[gg].ebv}
;fa = {x1:str1[gg].gaia_gmag,x2:str1[gg].tmass_jmag,x3:str1[gg].tmass_kmag,x4:str1[gg].ebv}
;fa = {x1:str1[gg].gaia_gmag,x2:str1[gg].tmass_jmag,x3:str1[gg].ebv}
x = findgen(ngg)
y = str1[gg].ps1_imag
err = str1[gg].ps1_ierr > 0.01
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
initpar[0] = 1.0
parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
;initpar[2] = 0.0  ;0.04
;parinfo[2].fixed = 1
;initpar[3] = 0.223
;parinfo[3].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
print,'i-band:'
print,par
;  1.00000    -0.458742    -0.275646    0.0967721
faall = {x1:str1.gaia_gmag,x2:jk0,x3:str1.ebv}
;faall = {x1:str1.gaia_gmag,x2:str1.tmass_jmag,x3:str1.tmass_kmag,x4:str1.ebv}
;faall = {x1:str1.gaia_gmag,x2:str1.tmass_jmag,x3:str1.ebv}
yfitall = superfit(str1.ps1_imag*0,par,_extra=faall)

plotc,str1.ebv,yfitall-str1.ps1_imag,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250

; Scatter plot
file = 'plots/nsc_color_relations_stripe82_super_i_jk_scatter'
ps_open,file,/color,thick=4,/encap
plotc,jk0,yfitall-str1.ps1_imag,str1.ebv,ps=1,sym=0.5,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='i-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*G'+stringize(par[1],ndec=3)+'*JK0'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
;push,plots,file

; Density plot
file = 'plots/nsc_color_relations_stripe82_super_i_jk_density'
ps_open,file,/color,thick=4,/encap
hess,jk0,yfitall-str1.ps1_imag,str1.ebv,dx=0.02,dy=0.02,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='i-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*G'+stringize(par[1],ndec=3)+'*JK0'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-str1[gg].ps1_imag)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Extinction plot
file = 'plots/nsc_color_relations_stripe82_super_i_ebv_jk'
ps_open,file,/color,thick=4,/encap
plotc,str1[gg].ebv,yfitall[gg]-str1[gg].ps1_imag,jk0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='i-band (color-coded by [J-Ks]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; J-z vs. J-K
;--------------
gd1 = where(str.ps1_gmag lt 20.0 and str.ps1_ng gt 5 and str.ps1_nr gt 5 and str.ps1_ni gt 5 and str.ps1_nz gt 5 and str.ps1_ny gt 5 and $
            str.ps1_gmag gt 12.0 and str.gaia_gmag ge 14 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 14.0 and $
            str.tmass_kmag lt 15.0,ngd1)
str1 = str[gd1]
jk0 = str1.tmass_jmag-str1.tmass_kmag-0.17*str1.ebv
gg = where(jk0 ge 0.3 and jk0 le 0.7 and str1.ps1_zerr gt 0,ngg)
;gg = where(jk0 ge 0.3 and jk0 le 0.7 and str1.ps1_zerr gt 0 and str1.ebv gt 0.15,ngg)
initpar = [0.0, 0.0, 0.0, 0.0]
fa = {x1:str1[gg].tmass_jmag,x2:jk0[gg],x3:str1[gg].ebv}
x = findgen(ngg)
y = str1[gg].ps1_zmag
err = str1[gg].ps1_zerr > 0.01
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
initpar[0] = 1.0
parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
initpar[2] = 0.40  ;0.45 ;0.04
parinfo[2].fixed = 1
;initpar[3] = 0.223
;parinfo[3].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
print,'z-band:'
print,par
;     1.00000     0.765720     0.400000     0.605658
faall = {x1:str1.tmass_jmag,x2:jk0,x3:str1.ebv}
yfitall = superfit(str1.ps1_zmag*0,par,_extra=faall)

plotc,str1.ebv,yfitall-str1.ps1_zmag,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250

; Scatter plot
file = 'plots/nsc_color_relations_stripe82_super_z_jk_scatter'
ps_open,file,/color,thick=4,/encap
plotc,jk0,yfitall-str1.ps1_zmag,str1.ebv,ps=1,sym=0.5,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='z-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*J+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
;push,plots,file

; Density plot
file = 'plots/nsc_color_relations_stripe82_super_z_jk_density'
ps_open,file,/color,thick=4,/encap
hess,jk0,yfitall-str1.ps1_zmag,str1.ebv,dx=0.02,dy=0.02,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='z-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*J+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-str1[gg].ps1_zmag)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Extinction plot
file = 'plots/nsc_color_relations_stripe82_super_z_ebv_jk'
ps_open,file,/color,thick=4,/encap
plotc,str1[gg].ebv,yfitall[gg]-str1[gg].ps1_zmag,jk0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='z-band (color-coded by [J-Ks]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; J-Y vs. J-K
;--------------
gd1 = where(str.ps1_gmag lt 20.0 and str.ps1_ng gt 5 and str.ps1_nr gt 5 and str.ps1_ni gt 5 and str.ps1_nz gt 5 and str.ps1_ny gt 5 and $
            str.ps1_gmag gt 12.0 and str.gaia_gmag ge 14 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 14.0 and $
            str.tmass_kmag lt 15.0,ngd1)
str1 = str[gd1]
jk0 = str1.tmass_jmag-str1.tmass_kmag-0.17*str1.ebv
gg = where(jk0 ge 0.3 and jk0 le 0.7 and str1.ps1_yerr gt 0,ngg)
;gg = where(jk0 ge 0.3 and jk0 le 0.7 and str1.ps1_yerr gt 0 and str1.ebv gt 0.15,ngg)
initpar = [0.0, 0.0, 0.0, 0.0]
fa = {x1:str1[gg].tmass_jmag,x2:jk0[gg],x3:str1[gg].ebv}
x = findgen(ngg)
y = str1[gg].ps1_ymag
err = str1[gg].ps1_yerr > 0.01
parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
initpar[0] = 1.0
parinfo[0].fixed = 1
; FIXING THE EXTINCTION TERM
initpar[2] = 0.20  ;0.45 ;0.04
parinfo[2].fixed = 1
;initpar[3] = 0.223
;parinfo[3].fixed = 1
par = mpfitfun('superfit',x,y,err,initpar,functargs=fa,parinfo=parinfo,status=status,yfit=yfit,/quiet)
print,'Y-band:'
print,par
;   1.00000     0.544821     0.200000     0.663380
faall = {x1:str1.tmass_jmag,x2:jk0,x3:str1.ebv}
yfitall = superfit(str1.ps1_ymag*0,par,_extra=faall)

plotc,str1.ebv,yfitall-str1.ps1_ymag,ps=1,sym=0.5,xr=[0,1],yr=[-0.5,0.5] 
oplot,[-1,4],[0,0],linestyle=2,co=250

; Scatter plot
file = 'plots/nsc_color_relations_stripe82_super_y_jk_scatter'
ps_open,file,/color,thick=4,/encap
plotc,jk0,yfitall-str1.ps1_ymag,str1.ebv,ps=1,sym=0.5,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='Y-band (color-coded by E[B-V])'
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*J+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
;push,plots,file

; Density plot
file = 'plots/nsc_color_relations_stripe82_super_y_jk_density'
ps_open,file,/color,thick=4,/encap
hess,jk0,yfitall-str1.ps1_ymag,str1.ebv,dx=0.02,dy=0.02,xr=[-0.5,1.5],yr=[-0.5,0.5],$
      xtit='(J-Ks)!d0!n',ytit='Residuals',tit='Y-band fit',/log
oplot,[-1,3],[0,0],co=255
oplot,[0.3,0.7],[0,0],co=0
al_legend,stringize(par[0],ndec=3)+'*J+'+stringize(par[1],ndec=3)+'*JK0+'+$
          stringize(par[2],ndec=3)+'*E(B-V)+'+stringize(par[3],ndec=3),textcolor=250,/top,/left,charsize=1.2
rms = mad(yfitall[gg]-str1[gg].ps1_ymag)
al_legend,['RMS='+stringize(rms,ndec=3)+' mag'],textcolor=250,/bottom,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file

; Extinction plot
file = 'plots/nsc_color_relations_stripe82_super_y_ebv_jk'
ps_open,file,/color,thick=4,/encap
plotc,str1[gg].ebv,yfitall[gg]-str1[gg].ps1_ymag,jk0[gg],ps=1,sym=0.5,xr=[0,0.8],yr=[-0.5,0.5],xs=1,ys=1,$
      xtit='E(B-V)',ytit='Residuals',tit='Y-band (color-coded by [J-Ks]!d0!n)'
oplot,[-1,3],[0,0],co=250
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
push,plots,file


; Combine all of the figures
;plots = file_search('nsc_color_relations_stripe82_super_*_*.pdf')
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=plots/nsc_color_relations_stripe82_super.pdf '+strjoin(plots+'.pdf',' ')


stop

end
