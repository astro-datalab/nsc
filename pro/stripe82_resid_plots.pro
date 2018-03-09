pro stripe82_resid_plots

; Plots of residuals for Stripe82 data
; same as in combine_stripe82_Xband_v2.pro but
; final version for paper

; Make the plot
!p.font = 0
setdisp

xr = [0.1,1.1]
yr = [-0.5,0.5]
charsize = 1.8 ;1.3
posim = [0.12, 0.11, 0.95, 0.86]
poscol = [0.12, 0.96, 0.95, 0.99]

dir = '/dl1/users/dnidever/nsc/instcal/combine/stripe82_v2/'

; ------- u-band ------

;restore,dir+'combine_stripe82_uband.dat'
;gj0 = allgaia.gmag - alltmass.jmag - 1.12*allcat.ebv
;model_mag = 0.2469*allgalex.nuv + 0.7501*allgaia.gmag + 0.5462*gj0 + 0.6809*allcat.ebv + 0.0052
;cmag = allcat.cmag
;gd = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0,ngd)
;save,model_mag,gj0,cmag,gd,file='stripe82_resid_plots_uband.dat'
restore,dir+'stripe82_resid_plots_uband.dat'

file = dir+'stripe82_uband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;gj0 = allgaia.gmag - alltmass.jmag - 1.12*allcat.ebv
; old version
;model_mag = 0.30874*allgalex.nuv + 0.6955*allgaia.gmag +
;0.424*allcat.ebv + 0.0930
; new version
;model_mag = 0.2469*allgalex.nuv + 0.7501*allgaia.gmag + 0.5462*gj0 + 0.6809*allcat.ebv + 0.0052
;gd = where(allcat.class_star ge 0.8 and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,gj0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=[0.6,1.4],yr=yr,xtit='(G-J)!d0!n',ytit='Model-Mag',tit='u-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,gj0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0.5,max=1.5
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.8 and xbin le 1.1,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;  0.16465859     -0.21008923
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.8,0.8],[-2,2],linestyle=2,co=250,thick=6
oplot,[1.1,1.1],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(G-J)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell


; ------- g-band ------

;restore,dir+'combine_stripe82_gband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = allapass.g_mag - 0.0421*jk0 - 0.05*allcat.ebv - 0.0620
;cmag = allcat.cmag
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_gband.dat'
restore,dir+'stripe82_resid_plots_gband.dat'

file = dir+'stripe82_gband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; old version
;model_mag = allapass.g_mag - 0.1433*jk0 - 0.05*allcat.ebv - 0.0138
; new version
;model_mag = allapass.g_mag - 0.0421*jk0 - 0.05*allcat.ebv - 0.0620
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='g-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.3 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;  0.0461191    -0.102567
; 0.0482447    -0.101153
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.3,0.3],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.7,0.7],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; -------- r-band ----------

;restore,dir+'combine_stripe82_rband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = allapass.r_mag - 0.0861884*jk0 + 0.0548607
;cmag = allcat.cmag
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_rband.dat'
restore,dir+'stripe82_resid_plots_rband.dat'

file = dir+'stripe82_rband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; old version
;model_mag = allapass.r_mag + 0.00740*jk0 + 0.000528
; new version
;model_mag = allapass.r_mag - 0.0861884*jk0 + 0.0548607
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='r-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.3 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;   -0.0543327    0.0935884
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.3,0.3],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.7,0.7],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; ----------- i-band -------------

;restore,dir+'combine_stripe82_iband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = allgaia.gmag - 0.4587*jk0 - 0.276*allcat.ebv + 0.0967721
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;cmag = allcat.cmag
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_iband.dat'
restore,dir+'stripe82_resid_plots_iband.dat'

file = dir+'stripe82_iband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; equation stayed the same
;model_mag = allgaia.gmag - 0.4587*jk0 - 0.276*allcat.ebv + 0.0967721
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='i-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.25 and xbin le 0.65,ngd)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;  -0.0143787    0.0222694
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.25,0.25],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.65,0.65],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; ----------- z-band ---------------

;restore,dir+'combine_stripe82_zband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = alltmass.jmag + 0.765720*jk0 + 0.40*allcat.ebv +  0.605658
;cmag = allcat.cmag
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_zband.dat'
restore,dir+'stripe82_resid_plots_zband.dat'

file = dir+'stripe82_zband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; equation stayed the same
;model_mag = alltmass.jmag + 0.765720*jk0 + 0.40*allcat.ebv +  0.605658
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='z-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.4 and xbin le 0.65,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
; -0.0371434    0.0832948
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.4,0.4],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.65,0.65],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

; --------- Y-band -------------

;restore,dir+'combine_stripe82_yband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = alltmass.jmag + 0.54482*jk0 + 0.20*allcat.ebv + 0.663380
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;cmag = allcat.cmag
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_yband.dat'
restore,dir+'stripe82_resid_plots_yband.dat'

file = dir+'stripe82_yband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; the equation stayed the same
;model_mag = alltmass.jmag + 0.54482*jk0 + 0.20*allcat.ebv + 0.663380
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='Y-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.4 and xbin le 0.7,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
;   0.00169113   0.00573150
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.4,0.4],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.7,0.7],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell


; ------------ VR-band -------------

;restore,dir+'combine_stripe82_vrband.dat'
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
;model_mag = allgaia.gmag
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
;cmag = allcat.cmag
;save,model_mag,jk0,cmag,gd,file='stripe82_resid_plots_vrband.dat'
restore,dir+'stripe82_resid_plots_vrband.dat'

file = dir+'stripe82_vrband_magdiff_color'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=10.5,ysize=9.5
;jk0 = alltmass.jmag-alltmass.kmag-0.17*allcat.ebv
; the equation stayed the same
;model_mag = allgaia.gmag
;gd = where(allcat.class_star gt 0.8 and alltmass.qflg eq 'AAA' and allcat.fwhm_world*3600 lt 2.0,ngd)
hess,jk0[gd],model_mag[gd]-cmag[gd],dx=0.02,dy=0.02,xr=xr,yr=yr,xtit='(J-Ks)!d0!n',ytit='Model-Mag',tit='VR-band',charsize=charsize,$
     posim=posim,poscol=poscol
bindata,jk0[gd],model_mag[gd]-cmag[gd],xbin,ybin,binsize=0.05,/med,min=0,max=1.2
oplot,xbin,ybin,ps=-1,co=255
gdbin = where(xbin ge 0.2 and xbin le 0.6,ngdbin)
coef = robust_poly_fitq(xbin[gdbin],ybin[gdbin],1)
xx = scale_vector(findgen(100),-1,3)
;oplot,xx,poly(xx,coef),co=250
oplot,[-1,3],[0,0],linestyle=2,co=255
oplot,[0.2,0.2],[-2,2],linestyle=2,co=250,thick=6
oplot,[0.6,0.6],[-2,2],linestyle=2,co=250,thick=6
;al_legend,[stringize(coef[1],ndec=3)+'*(J-Ks)!d0!n+'+stringize(coef[0],ndec=3)],textcolor=[250],/top,/left,charsize=1.4
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

stop

end
