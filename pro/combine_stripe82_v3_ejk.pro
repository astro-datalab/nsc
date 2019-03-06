pro combine_stripe82_v3_ejk

;; Add the new extinction information to Stripe82 file

; Stripe82_v3.fits has PS1, GaiaDR2, 2MASS, APASS, GALEX, SkymapperDR1, 2 million stars
;  need to add allwise |b|<16 or max(EBV)>0.2
; 55<ra<60.3 and |dec|<1.1
; 299.5<ra<311 and |dec|<1.1


ref = mrdfits('/dl1/users/dnidever/nsc/Stripe82_v3.fits',1)
old = ref
schema = ref[0]
struct_assign,{dum:''},schema
schema = create_struct(schema,'allwise_w1mag',99.99,'allwise_w1err',9.99,'allwise_w2mag',99.99,'allwise_w2err',9.99,$
                              'ejk',99.99,'e_ejk',99.99,'ext_type',0)
ref = replicate(schema,n_elements(old))
struct_assign,old,ref,/nozero
undefine,old

;; Get ALLWISE data
wise1=getrefcat(58.0,0.0,3.0,'ALLWISE')
wise2a=getrefcat(302,0.0,3.2,'ALLWISE')
wise2b=getrefcat(308.0,0.0,3.2,'ALLWISE')
srcmatch,wise2a.ra,wise2a.dec,wise2b.ra,wise2b.dec,1.0,ind1,ind2,/sph
left = wise2b
remove,ind2,left
wise2 = [wise2a,left]


srcmatch,ref.ra,ref.dec,wise1.ra,wise1.dec,1.0,ind1,ind2,/sph
ref[ind1].allwise_w1mag = wise1[ind2].w1mag
ref[ind1].allwise_w1err = wise1[ind2].e_w1mag
ref[ind1].allwise_w2mag = wise1[ind2].w2mag
ref[ind1].allwise_w2err = wise1[ind2].e_w2mag

srcmatch,ref.ra,ref.dec,wise2.ra,wise2.dec,1.0,ind1,ind2,/sph
ref[ind1].allwise_w1mag = wise2[ind2].w1mag
ref[ind1].allwise_w1err = wise2[ind2].e_w1mag
ref[ind1].allwise_w2mag = wise2[ind2].w2mag
ref[ind1].allwise_w2err = wise2[ind2].e_w2mag

;; Add new reddening information
ejk_sfd = 1.5*0.302*ref.ebv
ref.ejk = ejk_sfd
ref.e_ejk = 0.1
ref.ext_type = 1
gdwise = where(ref.tmass_jmag lt 50 and ref.tmass_hmag lt 50 and ref.tmass_kmag lt 50 and $
               ref.allwise_w2mag lt 50,ngdwise)
if ngdwise gt 0 then begin
  ejk = 1.5*0.918*(ref[gdwise].tmass_hmag-ref[gdwise].allwise_w2mag-0.05)
  e_ejk = 1.5*0.918*sqrt(ref[gdwise].tmass_herr^2+ref[gdwise].allwise_w2err^2)
  gdejk = where(ejk lt ejk_sfd[gdwise],ngdejk)
  if ngdejk gt 0 then begin
    ref[gdwise[gdejk]].ejk = ejk[gdejk] > 0.0
    ref[gdwise[gdejk]].e_ejk = e_ejk[gdejk]
    ref[gdwise[gdejk]].ext_type = 2
 endif
endif
;mwrfits,ref,'/dl1/users/dnidever/nsc/Stripe82_v3_ejk.fits',/create


stop

end
