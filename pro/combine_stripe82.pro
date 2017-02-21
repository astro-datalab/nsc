pro combine_stripe82

; Combine PS1, Gaia, 2MASS and APASS catalogs for Stripe82

dir = '/datalab/users/dnidever/nsc/'

ps1 = mrdfits(dir+'Stripe82_PS1.fits.gz',1)
gaia = mrdfits(dir+'Stripe82_Gaia.fits',1)
SRCMATCH,ps1.ramean,ps1.decmean,gaia.ra,gaia.dec,0.5,ind1,ind2,/sph,domains=10000L,count=nmatch
print,strtrim(nmatch,2),' matches between PS1 and Gaia'
; require at least PS1 and Gaia
schema = {ra:0.0d0,dec:0.0d0,ps1_objid:0L,ps1_ndetections:0,ps1_ng:0,ps1_nr:0,ps1_ni:0,ps1_nz:0,ps1_ny:0,$
          ps1_gmag:0.0,ps1_gerr:0.0,ps1_rmag:0.0,ps1_rerr:0.0,ps1_imag:0.0,ps1_ierr:0.0,ps1_zmag:0.0,ps1_zerr:0.0,$
          ps1_ymag:0.0,ps1_yerr:0.0,gaia_sourceid:0LL,gaia_gnobs:0L,gaia_gflux:0.0,gaia_gfluxerr:0.0,gaia_gmag:0.0,$
          gaia_gerr:0.0,tmass_match:0,tmass_designation:'',tmass_jmag:99.99,tmass_jerr:9.99,tmass_hmag:99.99,$
          tmass_herr:9.99,tmass_kmag:99.99,tmass_kerr:9.99,tmass_phqual:'',apass_match:0,apass_recno:0L,apass_nobs:0L,$
          apass_mobs:0L,apass_bv:99.99,apass_ebv:9.99,apass_vmag:99.99,apass_verr:9.99,apass_bmag:99.99,apass_berr:9.99,$
          apass_gmag:99.99,apass_gerr:9.99,apass_rmag:99.99,apass_rerr:9.99,apass_imag:99.99,apass_ierr:9.99,ebv:0.0}
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
all.gaia_sourceid  = gaia[ind2].source_id
all.gaia_gnobs  = gaia[ind2].phot_g_n_obs
all.gaia_gflux  = gaia[ind2].phot_g_mean_flux
all.gaia_gfluxerr  = gaia[ind2].phot_g_mean_flux_error
all.gaia_gmag = gaia[ind2].phot_g_mean_mag
all.gaia_gerr = 2.5*alog10(1.0+all.gaia_gfluxerr/all.gaia_gflux)
;all.gaia_gerr = 2.5*alog10(1.0+mgaia[gdphot].e__fg_/mgaia[gdphot]._fg_)

; Now match to 2MASS
tmass = mrdfits(dir+'Stripe82_2mass.fits',1)
SRCMATCH,all.ra,all.dec,tmass.ra,tmass.dec,0.5,tind1,tind2,/sph,domains=10000L,count=ntmatch
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
SRCMATCH,all.ra,all.dec,apass.raj2000,apass.dej2000,0.5,aind1,aind2,/sph,domains=10000L,count=namatch
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

; Add E(B-v)
glactc,all.ra,all.dec,2000.0,glon,glat,1,/deg
ebv = dust_getval(glon,glat,/noloop,/interp)
all.ebv = ebv

; Write final file
;mwrfits,all,dir+'Stripe82.fits',/create

stop

end
