;+
;
; GETREDDENING
;
; This calculates E(J-Ks) reddening using reference catalog data for
; the NSC.
;
; INPUTS:
;  ref      The structure with the reference catalog data.
;  exttype  Extinction type:
;             1 - SFD
;             2 - RJCE ALLWISE
;             3 - RJCE GlIMPSE
;             4 - RJCE SAGE
;
; OUTPUTS:
;  The REF structure EJK column and EXT_TYPE columns are updated.
;
; USAGE:
;  IDL>getreddening,ref,ext_type
;
; By D. Nidever  Feb 2019
;-

pro getreddening,ref,ext_type
  
;; Not enough inputs
if n_elements(ref) eq 0 or ext_type eq 0 then begin
  print,'Syntax - getreddening,ref,ext_type'
  return
endif

;; Add SFD reddening
GLACTC,ref.ra,ref.dec,2000.0,glon,glat,1,/deg
ebv = dust_getval(glon,glat,/noloop,/interp)
ref.ebv_sfd = ebv

;; Start with SFD extinction for all 
ejk_sfd = 1.5*0.302*ref.ebv_sfd
ref.ejk = ejk_sfd
ref.e_ejk = 0.1   ; not sure what to put here
bd = where(ref.ebv_sfd gt 0.3,nbd)   ; E(B-V)>0.3 is roughly where things "break down"
if nbd gt 0 then ref[bd].e_ejk=1.0  
ref.ext_type = 1

;; RJCE extinction
if ext_type ge 2 then begin
   
  ;; RJCE GLIMPSE, type=3
  if ext_type eq 3 then begin
    gdglimpse = where(ref.jmag lt 50 and ref.hmag lt 50 and ref.kmag lt 50 and $
                      ref.gl_45mag lt 50,ngdglimpse)
    if ngdglimpse gt 0 then begin
      ejk = 1.5*0.918*(ref[gdglimpse].hmag-ref[gdglimpse].gl_45mag-0.08)
      e_ejk = 1.5*0.918*sqrt(ref[gdglimpse].e_hmag^2+ref[gdglimpse].e_gl_45mag^2)
      ;gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdglimpse] and e_ejk lt 0.2,ngdejk)
      gdejk = where(ejk lt ejk_sfd[gdglimpse],ngdejk)
      if ngdejk gt 0 then begin
        ref[gdglimpse[gdejk]].ejk = ejk[gdejk] > 0.0
        ref[gdglimpse[gdejk]].e_ejk = e_ejk[gdejk]
        ref[gdglimpse[gdejk]].ext_type = 3
      endif
    endif
  endif
  
  ;; RJCE SAGE, type=4
  if ext_type eq 4 then begin
    gdsage = where(ref.jmag lt 50 and ref.hmag lt 50 and ref.kmag lt 50 and $
                   ref.sage_45mag lt 50,ngdsage)
    if ngdsage gt 0 then begin
      ejk = 1.5*0.918*(ref[gdsage].hmag-ref[gdsage].sage_45mag-0.08)
      e_ejk = 1.5*0.918*sqrt(ref[gdsage].e_hmag^2+ref[gdsage].e_sage_45mag^2)
      ;gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdsage] and e_ejk lt 0.2,ngdejk)
      gdejk = where(ejk lt ejk_sfd[gdsage],ngdejk)
      if ngdejk gt 0 then begin
        ref[gdsage[gdejk]].ejk = ejk[gdejk] > 0.0
        ref[gdsage[gdejk]].e_ejk = e_ejk[gdejk]
        ref[gdsage[gdejk]].ext_type = 4
      endif
    endif
  endif
  
  ;; RJCE ALLWISE, type=2
  ;;   Use ALLWISE for those that don't have IRAC data
  gdwise = where(ref.jmag lt 50 and ref.hmag lt 50 and ref.kmag lt 50 and $
                 ref.w2mag lt 50 and ref.ext_type le 1,ngdwise)
  if ngdwise gt 0 then begin
    ejk = 1.5*0.918*(ref[gdwise].hmag-ref[gdwise].w2mag-0.05)
    e_ejk = 1.5*0.918*sqrt(ref[gdwise].e_hmag^2+ref[gdwise].e_w2mag^2)
    ;gdejk = where(ejk gt 0 and ejk lt ejk_sfd[gdwise] and e_ejk lt 0.2,ngdejk)
    gdejk = where(ejk lt ejk_sfd[gdwise],ngdejk)
    if ngdejk gt 0 then begin
      ref[gdwise[gdejk]].ejk = ejk[gdejk] > 0.0
      ref[gdwise[gdejk]].e_ejk = e_ejk[gdejk]
      ref[gdwise[gdejk]].ext_type = 2
    endif
  endif
endif

;; Fix NANs in E_EJK
bd = where(finite(ref.e_ejk) eq 0,nbd)
if nbd gt 0 then ref[bd].e_ejk = 9.99

end

