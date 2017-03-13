pro nsc_color_relations

; Derive color-color relations for calibrating

dir = '/datalab/users/dnidever/smash/cp/red/photred/catalogs/final/v4/'
obj0 = mrdfits(dir+'Field100_combined_allobj.fits.gz',1)
xmatch0 = mrdfits(dir+'Field100_combined_allobj_xmatch.fits.gz',1)
obj0.id = strtrim(obj0.id,2)
xmatch0.id = strtrim(xmatch0.id,2)
match,obj0.id,xmatch0.id,ind1,ind2,/sort
obj = obj0[ind1]
xmatch = xmatch0[ind2]
setdisp
!p.font = 0

; J-z vs. J-Ks
file = 'nsc_color_relations_jz_jk'
ps_open,file,/color,thick=4,/encap
gd = where(xmatch.tmass_match eq 1 and xmatch.tmass_jmag lt 14.5 and obj.z lt 50 and finite(obj.z) eq 1,ngd)
xmatch2 = xmatch[gd]
obj2 = obj[gd]
plot,xmatch2.tmass_jmag-xmatch2.tmass_kmag,xmatch2.tmass_jmag-obj2.z,ps=8,sym=0.5,xr=[-0.2,1.5],yr=[-2,-0.5],$
     xs=1,ys=1,xtit='J-Ks',ytit='J-z',tit='Color-color relations for J-z'
;coef = robust_poly_fit(xmatch2.tmass_jmag-xmatch2.tmass_kmag,xmatch2.tmass_jmag-obj2.z,1)
err = sqrt(xmatch2.tmass_jerr^2 + obj2.zerr^2)
coef = dln_poly_fit(xmatch2.tmass_jmag-xmatch2.tmass_kmag,xmatch2.tmass_jmag-obj2.z,1,$
                    measure_errors=err,sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
print,coef
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=250
al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/right,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-i vs. G-J
file = 'nsc_color_relations_gi_gj'
ps_open,file,/color,thick=4,/encap
gd = where(xmatch.tmass_match eq 1 and xmatch.tmass_jmag lt 15.0 and xmatch.gaia_match eq 1 and obj.i lt 50 and finite(obj.i) eq 1,ngd)
xmatch2 = xmatch[gd]
obj2 = obj[gd]
plot,xmatch2.gaia_gmag-xmatch2.tmass_jmag,xmatch2.gaia_gmag-obj2.i,ps=8,sym=0.5,xr=[0.4,3.2],yr=[-0.2,0.8],$
     xs=1,ys=1,xtit='G-J',ytit='G-i',tit='Color-color relations for G-i'
;coef = robust_poly_fit(xmatch2.tmass_jmag-xmatch2.tmass_kmag,xmatch2.tmass_jmag-obj2.z,1)
err = sqrt(xmatch2.gaia_gerr^2 + obj2.ierr^2)
coef = dln_poly_fit(xmatch2.gaia_gmag-xmatch2.tmass_jmag,xmatch2.gaia_gmag-obj2.i,1,$
                    measure_errors=err,sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
gd2 = where(xmatch2.gaia_gmag-xmatch2.tmass_jmag ge 1.0 and xmatch2.gaia_gmag-xmatch2.tmass_jmag lt 2.4,ngd2)
coef2 = dln_poly_fit(xmatch2[gd2].gaia_gmag-xmatch2[gd2].tmass_jmag,xmatch2[gd2].gaia_gmag-obj2[gd2].i,1,$
                    measure_errors=err[gd2],sigma=coeferr2,yerror=yerror2,status=status2,yfit=yfit2,/bootstrap)
print,coef2
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=250
al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps

; G-z vs. G-J
file = 'nsc_color_relations_gz_gj'
ps_open,file,/color,thick=4,/encap
gd = where(xmatch.tmass_match eq 1 and xmatch.tmass_jmag lt 15.0 and xmatch.gaia_match eq 1 and obj.z lt 50 and finite(obj.z) eq 1,ngd)
xmatch2 = xmatch[gd]
obj2 = obj[gd]
plot,xmatch2.gaia_gmag-xmatch2.tmass_jmag,xmatch2.gaia_gmag-obj2.z,ps=8,sym=0.5,xr=[0.4,3.2],yr=[-0.2,1.5],$
     xs=1,ys=1,xtit='G-J',ytit='G-z',tit='Color-color relations for G-i'
;coef = robust_poly_fit(xmatch2.tmass_jmag-xmatch2.tmass_kmag,xmatch2.tmass_jmag-obj2.z,1)
err = sqrt(xmatch2.gaia_gerr^2 + obj2.zerr^2)
coef = dln_poly_fit(xmatch2.gaia_gmag-xmatch2.tmass_jmag,xmatch2.gaia_gmag-obj2.z,1,$
                    measure_errors=err,sigma=coeferr,yerror=yerror,status=status,yfit=yfit,/bootstrap)
;coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
gd2 = where(xmatch2.gaia_gmag-xmatch2.tmass_jmag ge 1.0 and xmatch2.gaia_gmag-xmatch2.tmass_jmag lt 2.4,ngd2)
coef2 = dln_poly_fit(xmatch2[gd2].gaia_gmag-xmatch2[gd2].tmass_jmag,xmatch2[gd2].gaia_gmag-obj2[gd2].z,1,$
                    measure_errors=err[gd2],sigma=coeferr2,yerror=yerror2,status=status2,yfit=yfit2,/bootstrap)
print,coef2
x = scale_vector(findgen(100),-2,4)
oplot,x,poly(x,coef),co=250
al_legend,['Zero-point = '+stringize(coef[0],ndec=4),'Error = '+stringize(coeferr[0],ndec=4)]+' mag',textcolor=[250,250],/top,/left,charsize=1.2
ps_close
ps2png,file+'.eps',/eps



stop

end
