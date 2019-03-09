pro combine_stripe82_v3_midplane

;; Get reference data for midplane fields

str = mrdfits('/dl1/users/dnidever/nsc/instcal/t3b/lists/nsc_instcal_calibrate.fits',1,/silent)
str.expdir = strtrim(str.expdir,2)
str.base = strtrim(str.base,2)

;;; I reran any field within 2 deg +/- of these coordinates
glactc,ra1,dec1,2000.0,5,1,2,/deg
glactc,ra2,dec2,2000.0,5,5,2,/deg
glactc,ra3,dec3,2000.0,5,10,2,/deg
glactc,ra4,dec4,2000.0,5,15,2,/deg
glactc,ra5,dec5,2000.0,10,20,2,/deg
glactc,ra6,dec6,2000.0,15,20,2,/deg
glactc,ra7,dec7,2000.0,30,25,2,/deg
rr = [ra1,ra2,ra3,ra4,ra5,ra6,ra7]
dd = [dec1,dec2,dec3,dec4,dec5,dec6,dec7]
radius = [1.0,1.0,1.0,1.0,2.0,2.0,2.0]
n = n_elements(rr)
;glon = [5.0,5.0,5.0,5.0,10.0,15.0,30.0]
;glat = [1.0,5.0,10.0,15.0,20.0,20.0,25.0]
;for i=0,n_elements(ra)-1 do begin
;  gd1 = where(str.ra ge ra[i]-2 and str.ra le ra[i]+2 and str.dec ge dec[i]-2 and str.dec le dec[i]+2,ngd1)
;  if ngd1 eq 0 then stop,'no exposures'
;  if ngd1 gt 0 then push,gd,gd1
;endfor
;ngd = n_elements(gd)
;
;;; Load the metadata to figure out the area covered
;for i=0,ngd-1 do begin
;  metafile = str[gd[i]].expdir+'/'+str[gd[i]].base+'_meta.fits'
;  meta = mrdfits(metafile,1,/silent)
;  chstr1 = mrdfits(metafile,2,/silent)
;  push,chstr,chstr1
;endfor
;mwrfits,chstr,'combine_stripe82_v3_midplane_chstr.fits',/create

;; Need 256<RA<

;rr = [258.43012,  257.87992, 254.67873, 255.85417, 258.43012, 260.63093,  262.80674, 264.93253, 267.05831, 269.484]
;dd = [8.872705,   -6.483075,  -10.434930,  -16.701443,  -18.169275,  -19.862927,  -21.725944,  -22.234040,  -23.927692, -24.943]
;radius = [3.0, 1.2, 3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
;n = n_elements(r)
;
;rr = [ 258.566626,  260.293983,  258.683339,  257.422835,  257.516206,  257.936374,  256.232359,  256.232359,  254.645057,  252.941042,  254.154861,$
;       254.364945,  256.115645,  256.115645,  257.866346,  257.843003,  259.430305,  258.986794,  261.040949,  260.434039,  262.838334,  264.472321,$
;       262.791649,  265.989595,  264.682405,  267.763638,  267.553554,  265.779511,  266.993330,  269.934506,  269.794450]
;dd = [  10.396394,    9.300361,    8.204327,    9.091593,    7.369254,   -6.461643,   -9.332207,  -11.419890,  -10.480433,  -11.524274,  -18.726779,$
;       -16.586904,  -15.699639,  -17.369785,  -16.430328,  -18.570203,  -17.630746,  -19.144316,  -18.622395,  -20.710078,  -21.492959,  -20.449117,$
;       -23.580641,  -21.284190,  -22.902144,  -22.484608,  -23.945986,  -23.737218,  -26.033669,  -23.424065,  -25.707364]
;n = n_elements(rr)
;radius = fltarr(n)+1.5
;radius = fltarr(n)+0.1

;plot,chstr.vra,chstr.vdec,ps=3
;pa=scale_vector(findgen(100),0,2*!dpi)
;for j=0,n-1 do oplot,rr[j]+radius[j]*sin(pa),dd[j]+radius[j]*cos(pa)   


;; Load the REF schema
schema = mrdfits('/dl1/users/dnidever/nsc/Stripe82_v3_ejk_schema.fits',1)
;; make all 99.99 by default    
for k=0,n_tags(schema)-1 do begin
  type = size(schema.(k),/type)
  if type eq 4 or type eq 5 then schema.(k)=99.99
endfor
ref = replicate(schema,1e7)
nref = n_elements(ref)
cnt = 0LL

;; Loop over the positions
For i=0,n-1 do begin
  cenra = rr[i]
  cendec = dd[i]
  rad = radius[i]
  print,''
  print,strtrim(i+1,2),' ',strtrim(cenra,2),' ',strtrim(cendec,2),' ',strtrim(rad,2)
  print,'-----------------------------'
  print,''

  filters='c4d-'+['u','g','r','i','z','Y','VR']
  ref0 = GETREFDATA(filters,cenra,cendec,rad)
  nref0 = n_elements(ref0)
  ;; gets ps, gaia, tmass, allwise, galex, ejk
  ;; put ref data in correct format
  ref1 = replicate(schema,nref0)
  ref1.ra = ref0.ra
  ref1.dec = ref0.dec
  ref1.ps1_gmag = ref0.ps_gmag
  ref1.ps1_rmag = ref0.ps_rmag
  ref1.ps1_imag = ref0.ps_imag
  ref1.ps1_zmag = ref0.ps_zmag
  ref1.ps1_ymag = ref0.ps_ymag
  ref1.gaia_gmag = ref0.gmag
  ref1.gaia_gerr = ref0.e_gmag
  ref1.gaia_bp = ref0.bp
  ref1.gaia_bperr = ref0.e_bp
  ref1.gaia_rp = ref0.rp
  ref1.gaia_rperr = ref0.e_rp
  ref1.tmass_jmag = ref0.jmag
  ref1.tmass_jerr = ref0.e_jmag
  ref1.tmass_hmag = ref0.hmag
  ref1.tmass_herr = ref0.e_hmag
  ref1.tmass_kmag = ref0.kmag
  ref1.tmass_kerr = ref0.e_kmag
  ref1.tmass_phqual = ref0.qflg
  ref1.galex_nuv = ref0.nuv
  ref1.galex_nuverr = ref0.e_nuv    
  ;; apass, sm, allwise
  if tag_exist(ref0,'w1mag') then begin  ;; if ALLWISE data was extracted
    ref1.allwise_w1mag = ref0.w1mag
    ref1.allwise_w1err = ref0.e_w1mag
    ref1.allwise_w2mag = ref0.w2mag
    ref1.allwise_w2err = ref0.e_w2mag
  endif
  ref1.ebv = ref0.ebv_sfd
  ref1.ejk = ref0.ejk
  ref1.e_ejk = ref0.e_ejk
  ref1.ext_type = ref0.ext_type

  ;; Add Skymapper data
  ;; need to get apass and skymapper
  ;apass0 = GETREFCAT(cenra,cendec,'APASS')
  sm = GETREFCAT(cenra,cendec,rad,'SKYMAPPER',count=nsm)
  if nsm gt 0 then begin
    SRCMATCH,ref1.ra,ref1.dec,sm.raj2000,sm.dej2000,1.0,ind1,ind2,/sph,coun=nmatch
    if nmatch gt 0 then begin
      ref1[ind1].sm_gmag = sm[ind2].sm_gmag
      ref1[ind1].sm_gerr = sm[ind2].e_sm_gmag
      ref1[ind1].sm_rmag = sm[ind2].sm_rmag
      ref1[ind1].sm_rerr = sm[ind2].e_sm_rmag
      ref1[ind1].sm_imag = sm[ind2].sm_imag
      ref1[ind1].sm_ierr = sm[ind2].e_sm_imag
      ref1[ind1].sm_zmag = sm[ind2].sm_zmag
      ref1[ind1].sm_zerr = sm[ind2].e_sm_zmag
    endif
  endif
  nref1 = n_elements(ref1)

  ;; Crossmatch
  if cnt gt 0 then begin
    SRCMATCH,ref[0:cnt-1].ra,ref[0:cnt-1].dec,ref1.ra,ref1.dec,1.0,ind1,ind2,/sph,count=nmatch
    ;; All matched, no new sources to add
    if nmatch eq nref1 then begin
      undefine,ref0,ref1
      print,'No new sources to add'
      goto,BOMB
    endif
    ;; Remove matches from ref1
    if nmatch gt 0 then REMOVE,ind2,ref1
    nref1 = n_elements(ref1)
  endif

  ;; Add new elements to REF
  if cnt+nref1 gt nref then begin
    print,'Adding new elements to REF'
    oldref = ref
    ref = replicate(schema,nref+1e7)
    ref[0:nref-1] = oldref
    undefine,oldred
    nref = n_elements(ref)
  endif

  ;; Stuff in new data
  ref[cnt:cnt+nref1-1] = ref1
  cnt += nref1

  ;stop
  BOMB:
Endfor
;; Trim the extra elements
ref = ref[0:cnt-1]
;; Save the final structure
;MWRFITS,ref,'/dl1/users/dnidever/nsc/Stripe82_v3_ejk_midplane.fits',/create


stop

end
