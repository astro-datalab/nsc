function getmodelmag,cat,filter

; This calculates the model magnitude for stars given the
; the magnitudes in reference catalogs
; NUV - Galex NUV magnitude
; GMAG - Gaia G magnitude
; JMAG - 2MASS J magnitude
; KMAG - 2MASS Ks magnitude
; APASS_GMAG - APASS g magnitue
; APASS_RMAG - APASS r magnitude
; EBV  - E(B-V) reddening

CASE instfilt of
; ---- DECam u-band ----
'u': begin
  ; (G-J)o = G-J-1.12*EBV
  col = cat.gmag - cat.jmag - 1.12*cat.ebv
  ; u = 0.2469*NUV + 0.7501*G + 0.5462*GJ0 + 0.6809*EBV + 0.0052  ; for 0.8<GJ0<1.1
  model_mag = 0.2469*cat.nuv + 0.7501*cat.gmag + 0.5462*col + 0.6809*cat.ebv + 0.0052
end

;---- DECam g-band ----
'g': begin
  ; (J-Ks)o = J-Ks-0.17*EBV
  col = cat.jmag-cat.kmag-0.17*cat.ebv
  ; g = APASS_G - 0.0421*JK0 - 0.05*EBV - 0.0620
  model_mag = cat.apass_gmag - 0.0421*col - 0.05*cat.ebv - 0.0620
end

; ---- DECam r-band ----
'r': begin
  ; (J-Ks)o = J-Ks-0.17*EBV
  col = cat.jmag-cat.kmag-0.17*cat.ebv
  ; r = APASS_r - 0.0861884*JK0 + 0.0*EBV + 0.0548607
  model_mag = cat.apass_rmag - 0.0861884*col + 0.0548607
end

; ---- DECam i-band ----
'i': begin
  ; (J-Ks)o = J-Ks-0.17*EBV
  col = cat.jmag-cat.kmag-0.17*cat.ebv
  ; i = G - 0.4587*JK0 - 0.276*EBV + 0.0967721
  model_mag = cat.gmag - 0.4587*col - 0.276*cat.ebv + 0.0967721
end

; ---- DECam z-band ----
'z': begin
  ; (J-Ks)o = J-Ks-0.17*EBV
  col = cat.jmag-cat.kmag-0.17*cat.ebv
  ; z = J + 0.765720*JK0 + 0.40*EBV +  0.605658
  model_mag = cat.jmag + 0.765720*col + 0.40*cat.ebv +  0.605658
end

; ---- DECam Y-band ----
'Y': begin
  ; (J-Ks)o = J-Ks-0.17*EBV
  col = cat.jmag-cat.kmag-0.17*cat.ebv
  ; Y = J + 0.54482*JK0 + 0.20*EBV + 0.663380
  model_mag = cat.jmag + 0.54482*col + 0.20*cat.ebv + 0.663380
end

; ---- DECam VR-band ----
'VR': begin
  model_mag = cat.gmag
end

return,model_mag

stop

end
