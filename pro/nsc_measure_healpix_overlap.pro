function nsc_measure_healpix_overlap,nside,pix,vra,vdec,dx=dx

; Measure the fractional overlap of a chip and a healpix pixel

;Calculating area of overlap of two colygons
;-probabaly fastest to create a grid of cartesian points and
; then use roi_cut on both polygons to get map of "inside" pixels
; and use union of both to get overlap
;-could do similar with high-res healpix and query_polygon, not sure
;  which one is faster.
;Figuring out which healpix are fully INSIDE the chip polygon
;-use ispointinpolygon.pro (part of dopolygonsoverlap.pro) on all the
;  healpix vertices   

radeg = 180.0d0 / !dpi

PIX2ANG_RING,nside,pix,htheta,hphi
hcenra = hphi*radeg
hcendec = 90-htheta*radeg
PIX2VEC_RING,nside,pix,hvec,hvertex
hvertex = transpose(reform(hvertex))   ; [1,3,4] -> [4,3]
VEC2ANG,hvertex,hdec,hra,/astro
    
; Transform to tangent plane system
ROTSPHCEN,hra,hdec,hcenra,hcendec,hlon,hlat,/gnomic
ROTSPHCEN,vra,vdec,hcenra,hcendec,vlon,vlat,/gnomic

; Transfor to pixel-based tangent plane system
rlon = minmax([vlon,hlon])
rlat = minmax([vlat,hlat])
lon0 = rlon[0]
lat0 = rlat[0]
if nside eq 4096 then begin
  dx = 1d-3  ; gives ~20x20 pixel                                                                                                                                                  
endif else begin
  dx = (range(hlon) > range(hlat))/20
  ; round 2nd decimal place
  pow = floor(alog10(abs(dx)))
  dx = round(dx*10.^(-pow+1))*10.^(pow-1)
endelse
vx = (vlon-lon0)/dx
vy = (vlat-lat0)/dx
hx = (hlon-lon0)/dx
hy = (hlat-lat0)/dx
nx = ceil((rlon[1]-lon0)/dx)+1
ny = ceil((rlat[1]-lat0)/dx)+1
; should give ~300x300 region total                                                                                                                                              

; Number of fine pixels inside the healpix
inim = bytarr(nx,ny)
in = polyfillv(hx,hy,nx,ny)
inim[in] = 1
ninside = n_elements(in)
; Pixels in the chip
in1 = polyfillv(vx,vy,nx,ny)
inim1 = bytarr(nx,ny)
inim1[in1] = 1
overlap = float(total(inim AND inim1))/ninside  ; fractional coverage

return,overlap

end
