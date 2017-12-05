pro create_combine_groups,list

; Divide the GALACTIC sky up into 5 equal GLON region and 2 pole regions,
; gp09   0,  -36<GLON<+36,  -55<GLAT<+40
; hulk  +1,  +36<GLON<+108, -55<GLAT<+40
; thing -1  -108<GLON<-36, -55<GLAT<+40
; gp05  -2  -180<GLO<-108, -55<GLAT<+40
; gp06  +2  +108<GLON<+180, -55<GLAT<+40
; gp07  all GLON, GLAT>+40
; gp08  all GLON, GLAT<-55
; randomize within each group

str = REPLICATE({glon:[0L,0L],glat:[0L,0L],server:''},7)
str[0].server = 'gp09'
;str[0].glon = [-38,38]    ; -36<GLON<+36
str[0].glon = [0,24]     ;  0<GLON<+24
str[0].glat = [-55,40]    ; -55<GLAT<+40
str[1].server = 'thing'
;str[1].glon = [34,110]    ; +36<GLON<+108
str[1].glon = [-48,0]     ; -48<GLON<0
str[1].glat = [-55,40]    ; -55<GLAT<+40 
str[2].server = 'hulk'
;str[2].glon = [250,326]  ; -108<GLON<-36
str[2].glon = [24,222]    ; 24<GLON<222
str[2].glat = [-55,40]    ; -55<GLAT<+40 
str[3].server = 'gp05'
;str[3].glon = [178,254] ;  -180<GLON<-108
str[3].glon = [222,278]   ; 222<GLON<278
str[3].glat = [-55,40]    ; -55<GLAT<+40 
str[4].server = 'gp06'
;str[4].glon = [106,182]   ; +108<GLON<+180
str[4].glon = [278,312]   ; 278<GLON<312
str[4].glat = [-55,40]    ; -55<GLAT<+40 
str[5].server = 'gp07'
str[5].glon = [0,360]  ; all GLON
str[5].glat = [40,100]    ; GLAT>+40
str[6].server = 'gp08'
str[6].glon = [0,360]  ; all GLON
str[6].glat = [-100,-55]  ; GLAT<-55
nstr = n_elements(str)

; Load the list
NSC_ROOTDIRS,dldir,mssdir,localdir
version='v2'
nside = 128
radeg = 180.0d0/!dpi
listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
if n_elements(list) eq 0 then begin
  healstr = MRDFITS(listfile,1,/silent)
  healstr.file = strtrim(healstr.file,2)
  healstr.base = strtrim(healstr.base,2)
  ; Get unique ones
  ui = uniq(healstr.file,sort(healstr.file))
  list = healstr[ui]
  add_tag,list,'ra',0.0d0,list
  add_tag,list,'dec',0.0d0,list
  add_tag,list,'glon',0.0d0,list
  add_tag,list,'glat',0.0d0,list
  PIX2ANG_RING,nside,list.pix,theta,phi
  ra = phi*radeg
  dec = 90-theta*radeg
  list.ra = ra
  list.dec = dec
  glactc,ra,dec,2000.0,glon,glat,1,/deg
  list.glon = glon
  list.glat = glat
endif

index = MRDFITS(listfile,2,/silent)
add_tag,index,'ra',0.0d0,index
add_tag,index,'dec',0.0d0,index
add_tag,index,'glon',0.0d0,index
add_tag,index,'glat',0.0d0,index
PIX2ANG_RING,nside,index.pix,ptheta,pphi
pra = pphi*radeg
pdec = 90-ptheta*radeg
glactc,pra,pdec,2000.0,pglon,pglat,1,/deg
index.ra = pra
index.dec = pdec
index.glon = pglon
index.glat = pglat
dt = predictcombtime(pglon,pglat,index.nexp)

; Loop over the servers
undefine,allpix
for i=0,nstr-1 do begin
  off = 2  ; 2 deg overlap/buffer
  glr2 = [glr[0]-off, glr[1]+off]
  gbr2 = [gbr[0]-off, gbr[1]+off]

  ; Select files on GLON/GLAT with overlap
  if glr2[0] lt 0 then begin
    ind = where((list.glon le glr2[1] or list.glon ge (glr2[0]+360)) and $
                 list.glat ge gbr2[0] and list.glat le gbr2[1],nind)
  endif else begin
    ind = where(list.glon ge glr2[0] and list.glon le glr2[1] and $
                list.glat ge gbr2[0] and list.glat le gbr2[1],nind)
  endelse
  ; Select healpix pixels using GLON/GLAT
  if glr[0] lt 0 then begin
    pind = where((index.glon lt glr[1] or index.glon ge (glr[0]+360)) and $
                  index.glat ge gbr[0] and index.glat lt gbr[1],npind)
  endif else begin
    pind = where(index.glon ge glr[0] and index.glon lt glr[1] and $
                 index.glat ge gbr[0] and index.glat lt gbr[1],npind)
  endelse
  tfrac = total(dt[pind])/total(dt)
  print,server,' ',strtrim(nind,2),' exposures selected.  Time fraction = ',stringize(tfrac,ndec=3)

  ; Catalog and metafile list
  catfile = list[ind].file
  metafile = repstr(catfile,'_cat','_meta')
  allfiles = [catfile,metafile]
  if server eq 'gp05' or server eq 'gp06' or server eq 'gp07' or server eq 'gp08' then allfiles='/net'+allfiles
  ; Pixel list
  allpix = index[pind].pix
  push,everypix,allpix

  ; Write out
  ;outfile = '/dl1/users/dnidever/nsc/instcal/v2/lists/combine_list_'+server+'.txt'
  ;WRITELINE,outfile,allfiles
  outpixfile = '/dl1/users/dnidever/nsc/instcal/v2/lists/combine_pix_'+server+'.txt'
  WRITELINE,outpixfile,allpix

  ;stop

endfor


stop

end
