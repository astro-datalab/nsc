pro nsc_instcal_calibrate_coverage,expdirs,redo=redo

;; Get coverage information from one calibrated exposure

nside = 512L  ;4096L
radeg = 180.0d0 / !dpi

;; Exposure loop
for e=0L,n_elements(expdirs)-1 do begin
  expdir = expdirs[e]
  base = file_basename(expdir)

  ;; Check if output file already exists
  outfile = expdir+'/'+base+'_hlpmeta.fits'
  if file_test(outfile) eq 1 and not keyword_set(redo) then begin
    print,outfile,' EXISTS and /redo NOT set'
    goto,expbomb
  endif

  ;; Load metadata information for this exposures
  metafile = expdir+'/'+base+'_meta.fits'
  expstr = mrdfits(metafile,1,/silent)
  chstr = mrdfits(metafile,2,/silent)
  nchstr = n_elements(chstr)

  print,'Creating coverage structure for exposure ',expdir
  schema = {pix:0L,nmeas:0L,depth95:0.0}
  undefine,all
  ;; Loop over chips
  for i=0,nchstr-1 do begin
    chstr1 = chstr[i]

    cenra = chstr1.cenra
    cendec = chstr1.cendec
    ROTSPHCEN,chstr1.vra,chstr1.vdec,cenra,cendec,vlon,vlat,/gnomic
    radius = max(sqrt(vlon^2+vlat^2))*1.5
    area = polygonarea(vlon,vlat)

    ;; Get the pixel numbers healpix that overlap this chip
    ANG2VEC,cendec,cenra,vector,/astro
    QUERY_DISC,nside,vector,radius,listpix,nlist,/deg,/inclusive
    PIX2ANG_RING,nside,listpix,theta,phi
    listra = phi*radeg 
    listdec = 90-theta*radeg

    ;; Loop over potentially overlapping healpix
    ;plot,[vlon,vlon[0]],[vlat,vlat[0]],ps=-1,xr=[-0.4,0.4],yr=[-0.4,0.4],xs=1,ys=1
    for j=0,nlist-1 do begin
      ;; Check that they actually overlap
      PIX2VEC_RING,nside,listpix[j],vec_out,vertex
      vertex = transpose(reform(vertex))  ; [1,3,4] -> [4,3] 
      VEC2ANG,vertex,bnddec,bndra,/astro
      ROTSPHCEN,bndra,bnddec,cenra,cendec,bndlon,bndlat,/gnomic
      olap = dopolygonsoverlap(vlon,vlat,bndlon,bndlat)
      co = 250
      if olap eq 1 then co=150
      ;oplot,[bndlon,bndlon[0]],[bndlat,bndlat[0]],ps=-4,co=co
      ;print,j,' ',olap
      ;; get polygon intersection
      if olap eq 1 then begin
        polygonintersection,vlon,vlat,bndlon,bndlat,intlon,intlat
        intarea = polygonarea(intlon,intlat)
        ;print,intarea
        ;oplot,[intlon,intlon[0]],[intlat,intlat[0]],co=80
        new = schema
        new.pix = listpix[j]
        new.nmeas = round(chstr1.nmeas*(intarea/area))
        new.depth95 = chstr1.depth95   ;; depth should be identical for all chips
        push,all,new
      endif
    endfor  ; healpix list loop
  endfor  ; chip loop

  ;; Deal with duplicate
  ui = uniq(all.pix,sort(all.pix))
  if n_elements(ui) lt n_elements(all) then begin
    orig = all
    index = create_index(orig.pix)
    all = orig[index.index[index.lo]]
    for i=0,n_elements(all)-1 do begin
      ind = index.index[index.lo[i]:index.hi[i]]
      all[i].nmeas = total(orig[ind].nmeas)
    endfor
    undefine,orig
  endif

  ;; save in exposure directory as hlpmeta.fits
  print,'Writing output to ',outfile
  MWRFITS,all,outfile,/create

  EXPBOMB:
Endfor ;; exposure loop

;stop

end
