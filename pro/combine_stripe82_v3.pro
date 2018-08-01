pro combine_stripe82_v3

;; Combine PS1, Gaia DR2, 2MASS, APASS, Galex and Skymapper catalogs for Stripe82
;; This is for NSC DR2 (v3)

dir = '/datalab/users/dnidever/nsc/'
dcr = 1.0

ps1 = mrdfits(dir+'Stripe82_PS1.fits.gz',1)
;gaia = mrdfits(dir+'Stripe82_Gaia.fits',1)
gaia = mrdfits(dir+'gaiadr2_stripe82.fits',1)
SRCMATCH,ps1.ramean,ps1.decmean,gaia.ra,gaia.dec,dcr,ind1,ind2,/sph,domains=10000L,count=nmatch
print,strtrim(nmatch,2),' matches between PS1 and Gaia'
; require at least PS1 and Gaia
schema = {ra:0.0d0,dec:0.0d0,ps1_objid:0L,ps1_ndetections:0,ps1_ng:0,ps1_nr:0,ps1_ni:0,ps1_nz:0,ps1_ny:0,$
          ps1_gmag:0.0,ps1_gerr:0.0,ps1_rmag:0.0,ps1_rerr:0.0,ps1_imag:0.0,ps1_ierr:0.0,ps1_zmag:0.0,ps1_zerr:0.0,$
          ps1_ymag:0.0,ps1_yerr:0.0,gaia_id:0LL,gaia_gmag:99.99,gaia_gerr:9.99,gaia_bp:99.99,gaia_bperr:9.99,gaia_rp:99.99,gaia_rperr:9.99,$
          tmass_match:0,tmass_designation:'',tmass_jmag:99.99,tmass_jerr:9.99,tmass_hmag:99.99,$
          tmass_herr:9.99,tmass_kmag:99.99,tmass_kerr:9.99,tmass_phqual:'',apass_match:0,apass_recno:0L,apass_nobs:0L,$
          apass_mobs:0L,apass_bv:99.99,apass_ebv:9.99,apass_vmag:99.99,apass_verr:9.99,apass_bmag:99.99,apass_berr:9.99,$
          apass_gmag:99.99,apass_gerr:9.99,apass_rmag:99.99,apass_rerr:9.99,apass_imag:99.99,apass_ierr:9.99,$
          galex_match:0,galex_fuv:99.99,galex_fuverr:9.99,galex_nuv:99.99,galex_nuverr:9.99,$
          sm_match:0,sm_id:0L,sm_umag:99.99,sm_uerr:9.99,sm_vmag:99.99,sm_verr:9.99,sm_gmag:99.99,sm_gerr:9.99,sm_rmag:99.99,sm_rerr:9.99,$
          sm_imag:99.99,sm_ierr:9.99,sm_zmag:99.99,sm_zerr:9.99,ebv:0.0}
all = replicate(schema,nmatch)
all.ps1_objid = ps1[ind1].objid
all.ps1_ndetections = ps1[ind1].ndetections
all.ps1_ng = ps1[ind1].ng
all.ps1_nr = ps1[ind1].nr
all.ps1_ni = ps1[ind1].ni
all.ps1_nz = ps1[ind1].nz
all.ps1_ny = ps1[ind1].ny
all.ps1_gmag = ps1[ind1].gmeanpsfmag
all.ps1_gerr = ps1[ind1].gmeanpsfmagerr
all.ps1_rmag = ps1[ind1].rmeanpsfmag
all.ps1_rerr = ps1[ind1].rmeanpsfmagerr
all.ps1_imag = ps1[ind1].imeanpsfmag
all.ps1_ierr = ps1[ind1].imeanpsfmagerr
all.ps1_zmag = ps1[ind1].zmeanpsfmag
all.ps1_zerr = ps1[ind1].zmeanpsfmagerr
all.ps1_ymag = ps1[ind1].ymeanpsfmag
all.ps1_yerr = ps1[ind1].ymeanpsfmagerr
; use Gaia coordinates
all.ra = gaia[ind2].ra
all.dec = gaia[ind2].dec
all.gaia_id  = gaia[ind2].source_id
all.gaia_gmag  = gaia[ind2].phot_g_mean_mag
all.gaia_gerr  = 2.5*alog10(1.0+gaia[ind2].phot_g_mean_flux_error/gaia[ind2].phot_g_mean_flux)
all.gaia_bp  = gaia[ind2].phot_bp_mean_mag
all.gaia_bperr = 2.5*alog10(1.0+gaia[ind2].phot_bp_mean_flux_error/gaia[ind2].phot_bp_mean_flux)
all.gaia_rp = gaia[ind2].phot_rp_mean_mag
all.gaia_rperr = 2.5*alog10(1.0+gaia[ind2].phot_rp_mean_flux_error/gaia[ind2].phot_rp_mean_flux)
;all.gaia_id  = gaia[ind2].id
;all.gaia_gmag  = gaia[ind2].gmag
;all.gaia_gerr  = gaia[ind2].gerr
;all.gaia_bp  = gaia[ind2].bp
;all.gaia_bperr = gaia[ind2].bperr
;all.gaia_rp = gaia[ind2].rp
;all.gaia_rperr = gaia[ind2].rperr

; Now match to 2MASS
tmass = mrdfits(dir+'Stripe82_2mass.fits',1)
SRCMATCH,all.ra,all.dec,tmass.ra,tmass.dec,dcr,tind1,tind2,/sph,domains=10000L,count=ntmatch
print,strtrim(ntmatch,2),' 2MASS matches'
all[tind1].tmass_match = 1
all[tind1].tmass_designation = tmass[tind2].designation
all[tind1].tmass_jmag = tmass[tind2].jmag
all[tind1].tmass_jerr = tmass[tind2].e_jmag
all[tind1].tmass_hmag = tmass[tind2].hmag
all[tind1].tmass_herr = tmass[tind2].e_hmag
all[tind1].tmass_kmag = tmass[tind2].kmag
all[tind1].tmass_kerr = tmass[tind2].e_kmag
all[tind1].tmass_phqual = tmass[tind2].ph_qual

; Match to APASS
apass = mrdfits(dir+'Stripe82_apass.fits.gz',1)
SRCMATCH,all.ra,all.dec,apass.raj2000,apass.dej2000,dcr,aind1,aind2,/sph,domains=10000L,count=namatch
print,strtrim(namatch,2),' APASS matches'
all[aind1].apass_match = 1
all[aind1].apass_recno = apass[aind2].recno
all[aind1].apass_nobs = apass[aind2].nobs
all[aind1].apass_mobs = apass[aind2].mobs
all[aind1].apass_bv = apass[aind2].b_v
all[aind1].apass_ebv = apass[aind2].e_b_v
all[aind1].apass_vmag = apass[aind2].vmag
all[aind1].apass_verr = apass[aind2].e_vmag
all[aind1].apass_bmag = apass[aind2].bmag
all[aind1].apass_berr = apass[aind2].e_bmag
all[aind1].apass_gmag = apass[aind2].g_mag
all[aind1].apass_gerr = apass[aind2].e_g_mag
all[aind1].apass_rmag = apass[aind2].r_mag
all[aind1].apass_rerr = apass[aind2].e_r_mag
all[aind1].apass_imag = apass[aind2].i_mag
all[aind1].apass_ierr = apass[aind2].e_i_mag

; Match to Galex
galex = mrdfits('/dl1/users/dnidever/nsc/instcal/galex_stripe82.fits',1)
SRCMATCH,all.ra,all.dec,galex.raj2000,galex.dej2000,dcr,gind1,gind2,/sph,domains=10000L,count=ngmatch
print,strtrim(ngmatch,2),' Galex matches'
all[gind1].galex_match = 1
all[gind1].galex_fuv = galex[gind2].fuv
all[gind1].galex_fuverr = galex[gind2].e_fuv
all[gind1].galex_nuv = galex[gind2].nuv
all[gind1].galex_nuverr = galex[gind2].e_nuv

; Match to Skymapper
sm = mrdfits('/dl1/users/dnidever/skymapper/skymapper_stripe82.fits.gz',1)
SRCMATCH,all.ra,all.dec,sm.raj2000,sm.dej2000,dcr,sind1,sind2,/sph,domains=10000L,count=nsmatch
print,strtrim(nsmatch,2),' Skymapper matches'
all[sind1].sm_match = 1
all[sind1].sm_id = sm[sind2].object_id
all[sind1].sm_umag = sm[sind2].u_psf
all[sind1].sm_uerr = sm[sind2].e_u_psf
all[sind1].sm_vmag = sm[sind2].v_psf
all[sind1].sm_verr = sm[sind2].e_v_psf
all[sind1].sm_gmag = sm[sind2].g_psf
all[sind1].sm_gerr = sm[sind2].e_g_psf
all[sind1].sm_rmag = sm[sind2].r_psf
all[sind1].sm_rerr = sm[sind2].e_r_psf
all[sind1].sm_imag = sm[sind2].i_psf
all[sind1].sm_ierr = sm[sind2].e_i_psf
all[sind1].sm_zmag = sm[sind2].z_psf
all[sind1].sm_zerr = sm[sind2].e_z_psf
;all[sind1].sm_id = sm[sind2].id 
;all[sind1].sm_umag = sm[sind2].umag
;all[sind1].sm_uerr = sm[sind2].uerr
;all[sind1].sm_vmag = sm[sind2].vmag
;all[sind1].sm_verr = sm[sind2].verr
;all[sind1].sm_gmag = sm[sind2].gmag
;all[sind1].sm_gerr = sm[sind2].gerr
;all[sind1].sm_rmag = sm[sind2].rmag
;all[sind1].sm_rerr = sm[sind2].rerr
;all[sind1].sm_imag = sm[sind2].imag
;all[sind1].sm_ierr = sm[sind2].ierr
;all[sind1].sm_zmag = sm[sind2].zmag
;all[sind1].sm_zerr = sm[sind2].zerr

; Add E(B-v)
glactc,all.ra,all.dec,2000.0,glon,glat,1,/deg
ebv = dust_getval(glon,glat,/noloop,/interp)
all.ebv = ebv

; Write final file
print,'Final file to ',dir+'Stripe82_v3.fits'
mwrfits,all,dir+'Stripe82_v3.fits',/create

stop

end
