pro nsc_color_relations_stripe82

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
gd1 = where(str.ps1_gmag lt 20.0 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 15.5,ngd1)
str1 = str[gd1]
gd2 = where(str.ps1_gmag lt 20.0 and str.tmass_match eq 1 and str.tmass_phqual eq 'AAA' and str.tmass_jmag lt 15.5 and $
            str.apass_gmag lt 19.0 and finite(str.apass_gmag) eq 1 and str.apass_rmag lt 19.0 and finite(str.apass_rmag) eq 1,ngd2)
str2 = str[gd2]

setdisp
!p.font = 0

; APASS_g-r vs. J-Ks
;--------------
file = 'nsc_color_relations_stripe82_agg_jk'
ps_open,file,/color,thick=4,/encap
hess,str2.tmass_jmag-str2.tmass_kmag,str2.apass_gmag-str2.ps1_gmag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-0.5,0.5],/log,xtit='J-Ks',ytit='APASS_g-g'
err = sqrt(str2.ps1_gerr^2 + str2.apass_gerr^2) > 0.01
gg = where(str2.tmass_jmag-str2.tmass_kmag ge 0.30 and str2.tmass_jmag-str2.tmass_kmag le 0.7,ngg)
coef = dln_poly_fit(str2[gg].tmass_jmag-str2[gg].tmass_kmag,str2[gg].apass_gmag-str2[gg].ps1_gmag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'APASS_g-g vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.30,0.70)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; APASS_r-r vs. J-Ks
;--------------
file = 'nsc_color_relations_stripe82_arr_jk'
ps_open,file,/color,thick=4,/encap
hess,str2.tmass_jmag-str2.tmass_kmag,str2.apass_rmag-str2.ps1_rmag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-0.6,0.5],/log,xtit='J-Ks',ytit='APASS_r-r'
err = sqrt(str2.ps1_rerr^2 + str2.apass_rerr^2) > 0.01
gg = where(str2.tmass_jmag-str2.tmass_kmag ge 0.30 and str2.tmass_jmag-str2.tmass_kmag le 0.77,ngg)
coef = dln_poly_fit(str2[gg].tmass_jmag-str2[gg].tmass_kmag,str2[gg].apass_rmag-str2[gg].ps1_rmag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'APASS_r-r vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.30,0.77)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-r vs. J-Ks
;--------------
file = 'nsc_color_relations_stripe82_gr_jk'
ps_open,file,/color,thick=4,/encap
hess,str1.tmass_jmag-str1.tmass_kmag,str1.gaia_gmag-str1.ps1_rmag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-1,0.5],/log,xtit='J-Ks',ytit='G-r'
err = sqrt(str1.ps1_rerr^2 + str1.gaia_gerr^2)
gg = where(str1.tmass_jmag-str1.tmass_kmag ge 0.25 and str1.tmass_jmag-str1.tmass_kmag le 0.68,ngg)
coef = dln_poly_fit(str1[gg].tmass_jmag-str1[gg].tmass_kmag,str1[gg].gaia_gmag-str1[gg].ps1_rmag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'G-r vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.25,0.68)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-r vs. G-J
;--------------
file = 'nsc_color_relations_stripe82_gr_gj'
ps_open,file,/color,thick=4,/encap
hess,str1.gaia_gmag-str1.tmass_jmag,str1.gaia_gmag-str1.ps1_rmag,dx=0.02,dy=0.02,xr=[-0.3,3.5],yr=[-1,0.5],/log,xtit='G-J',ytit='G-r'
err = sqrt(str1.ps1_rerr^2 + str1.gaia_gerr^2)
gg = where(str1.gaia_gmag-str1.tmass_jmag ge 0.84 and str1.gaia_gmag-str1.tmass_jmag le 1.5,ngg)
coef = dln_poly_fit(str1[gg].gaia_gmag-str1[gg].tmass_jmag,str1[gg].gaia_gmag-str1[gg].ps1_rmag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'G-r vs. G-J'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.92,1.6)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-i vs. J-Ks
;--------------
file = 'nsc_color_relations_stripe82_gi_jk'
ps_open,file,/color,thick=4,/encap
hess,str1.tmass_jmag-str1.tmass_kmag,str1.gaia_gmag-str1.ps1_imag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-0.6,1],/log,xtit='J-Ks',ytit='G-i'
err = sqrt(str1.ps1_ierr^2 + str1.gaia_gerr^2)
gg = where(str1.tmass_jmag-str1.tmass_kmag ge 0.3 and str1.tmass_jmag-str1.tmass_kmag le 0.72,ngg)
coef = dln_poly_fit(str1[gg].tmass_jmag-str1[gg].tmass_kmag,str1[gg].gaia_gmag-str1[gg].ps1_imag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'G-i vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.3,0.72)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-i vs. G-J
;--------------
file = 'nsc_color_relations_stripe82_gi_gj'
ps_open,file,/color,thick=4,/encap
hess,str1.gaia_gmag-str1.tmass_jmag,str1.gaia_gmag-str1.ps1_imag,dx=0.02,dy=0.02,xr=[-0.3,3.5],yr=[-0.6,1],/log,xtit='G-J',ytit='G-i'
err = sqrt(str1.ps1_ierr^2 + str1.gaia_gerr^2)
gg = where(str1.gaia_gmag-str1.tmass_jmag ge 0.92 and str1.gaia_gmag-str1.tmass_jmag le 1.6,ngg)
coef = dln_poly_fit(str1[gg].gaia_gmag-str1[gg].tmass_jmag,str1[gg].gaia_gmag-str1[gg].ps1_imag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'G-i vs. G-J'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.92,1.6)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; J-z vs. J-K
;--------------
file = 'nsc_color_relations_stripe82_jz_jk'
ps_open,file,/color,thick=4,/encap
hess,str1.tmass_jmag-str1.tmass_kmag,str1.tmass_jmag-str1.ps1_zmag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-2,-0.2],/log,xtit='J-Ks',ytit='J-z'
err = sqrt(str1.ps1_zerr^2 + str1.tmass_jerr^2)
gg = where(str1.tmass_jmag-str1.tmass_kmag ge 0.3 and str1.tmass_jmag-str1.tmass_kmag le 0.9,ngg)
coef = dln_poly_fit(str1[gg].tmass_jmag-str1[gg].tmass_kmag,str1[gg].tmass_jmag-str1[gg].ps1_zmag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'J-z vs. J-Ks'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.3,0.9)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; J-Y vs. J-K
;--------------
file = 'nsc_color_relations_stripe82_jy_jk'
ps_open,file,/color,thick=4,/encap
hess,str1.tmass_jmag-str1.tmass_kmag,str1.tmass_jmag-str1.ps1_ymag,dx=0.02,dy=0.02,xr=[-0.3,1.5],yr=[-2,-0.2],/log,xtit='J-Ks',ytit='J-Y'
err = sqrt(str1.ps1_yerr^2 + str1.tmass_jerr^2)
gg = where(str1.tmass_jmag-str1.tmass_kmag ge 0.3 and str1.tmass_jmag-str1.tmass_kmag le 0.9,ngg)
coef = dln_poly_fit(str1[gg].tmass_jmag-str1[gg].tmass_kmag,str1[gg].tmass_jmag-str1[gg].ps1_ymag,1,$
                    measure_errors=err[gg],sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
print,'J-Y vs. J-K'
print,coef
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=255
x = scale_vector(findgen(100),0.3,0.9)
oplot,x,poly(x,coef),co=0
;al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; Is there some dependence on extinction??

stop

end
