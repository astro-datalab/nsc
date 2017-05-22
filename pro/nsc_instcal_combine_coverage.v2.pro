pro nsc_instcal_combine_coverage,redo=redo

; Make the NSC coverage/depth map

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/nsc/instcal/'
nside = 128
nmulti = 30
radeg = 180.0d0 / !dpi

; Restore the calibration summary file
temp = MRDFITS(dir+'nsc_instcal_calibrate.fits',1,/silent)
schema = temp[0]
struct_assign,{dum:''},schema
schema = create_struct(schema,'chipindx',-1LL)
str = replicate(schema,n_elements(temp))
struct_assign,temp,str,/nozero
str.expdir = strtrim(str.expdir,2)
str.instrument = strtrim(str.instrument,2)
str.metafile = strtrim(str.metafile,2)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)
str.filter = strtrim(str.filter,2)
gd = where(str.success eq 1,nstr)
str = str[gd]
si = sort(str.expdir)
str = str[si]
chstr = mrdfits(dir+'nsc_instcal_calibrate.fits',2,/silent)
chstr.expdir = strtrim(chstr.expdir,2)
chstr.instrument = strtrim(chstr.instrument,2)
nchstr = n_elements(chstr)
; Get indices for CHSTR
siexp = sort(chstr.expdir)
chstr = chstr[siexp]
expdir = chstr[siexp].expdir
brklo = where(expdir ne shift(expdir,1),nbrk)
brkhi = [brklo[1:nbrk-1]-1,n_elements(expdir)-1]
nchexp = brkhi-brklo+1
if nstr ne n_elements(brklo) then stop,'number of exposures in STR and CHSTR do not match'
str.chipindx = brklo
str.nchips = nchexp
; Fixing absolute paths of flux filename
file = str.file
g1 = where(stregex(file,'/net/mss1/',/boolean) eq 1,ng1)
file[g1] = strmid(file[g1],10)
g2 = where(stregex(file,'/mss1/',/boolean) eq 1,ng2)
file[g2] = strmid(file[g2],6)

; Fix instrument in STR and CHSTR
print,'FIXING INSTRUMENT IN STR AND CHSTR'
type = ['c4d','k4m','ksb']
for i=0,n_elements(type)-1 do begin
  gd = where(stregex(str.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
  if ngd gt 0 then str[gd].instrument=type[i]
  gd = where(stregex(chstr.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
  if ngd gt 0 then chstr[gd].instrument=type[i]
endfor

; APPLY RELEASE-DATE CUTS
list1 = MRDFITS(dir+'decam_instcal_list.fits',1)
list2 = MRDFITS(dir+'mosaic3_instcal_list.fits',1)
list3 = MRDFITS(dir+'bok90prime_instcal_list.fits',1)
list = [list1,list2,list3]
list.fluxfile = strtrim(list.fluxfile,2)
fluxfile = strmid(list.fluxfile,10)
MATCH,fluxfile,file,ind1,ind2,/sort,count=nmatch
; some don't match because they were from a previous version
;  of the input list
release_date = strarr(n_elements(str))+'2020-01-01 00:00:00'
release_date[ind2] = list[ind1].release_date
release_year = long(strmid(release_date,0,4))
release_month = long(strmid(release_date,5,2))
release_day = long(strmid(release_date,8,2))
release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
release_cutoff = [2017,4,24]  ; April 24, 2017 for now
release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
gdrelease = where(release_mjd le release_cutoff_mjd,ngdrelease,comp=bdrelease,ncomp=nbdrelease)
print,strtrim(ngdrelease,2),' exposures are PUBLIC'
str = str[gdrelease]  ; impose the public data cut

; Zero-point structure
zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},10)
zpstr[0:6].instrument = 'c4d'
zpstr[0:6].filter = ['u','g','r','i','z','Y','VR']
;ZP(u) = -0.38340071*X -1.6202186
zpstr[0].amcoef = [-1.6202186,-0.38340071]
;ZP(g) = -0.19290602*X + 0.25785509
zpstr[1].amcoef = [0.25785509,-0.19290602]
;ZP(r) = -0.17150109*X + 0.58451093
zpstr[2].amcoef = [0.58451093,-0.17150109]
;ZP(i) = -0.10241476*X + 0.39506315
zpstr[3].amcoef = [0.39506315,-0.10241476]
;ZP(z) = -0.07762681*X + 0.097165843
zpstr[4].amcoef = [0.097165843,-0.07762681]
;ZP(Y) = -0.03525268*X -1.0885863
zpstr[5].amcoef = [-1.0885863,-0.03525268]
;ZP(VR) = -0.091242237*X + 1.0202652
zpstr[6].amcoef = [1.0202652,-0.091242237]
; Mosiac3 z-band
zpstr[7].instrument = 'k4m'
zpstr[7].filter = 'z'
zpstr[7].amcoef = [-2.4692841,  -1.0928824]
; Bok 90Prime, g and r
zpstr[8].instrument = 'ksb'
zpstr[8].filter = 'g'
zpstr[8].amcoef = [-2.80864,  -1.46826]
zpstr[9].instrument = 'ksb'
zpstr[9].filter = 'r'
zpstr[9].amcoef = [-4.11, 0.0]
nzpstr = n_elements(zpstr)

; APPLY QA CUTS IN ZEROPOINT AND SEEING
print,'APPLY QA CUTS IN ZEROPOINT AND SEEING'
fwhmthresh = 3.0  ; arcsec
;filters = ['u','g','r','i','z','Y','VR']
;nfilters = n_elements(filters)
;zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
;zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
badmask = bytarr(n_elements(str)) + 1

for i=0,nzpstr-1 do begin
  ind = where(str.instrument eq zpstr[i].instrument and str.filter eq zpstr[i].filter and str.success eq 1,nind)
  print,zpstr[i].instrument,'-',zpstr[i].filter,' ',strtrim(nind,2),' exposures'
  if nind gt 0 then begin
    zpterm = str[ind].zpterm
    am = str[ind].airmass
    mjd = str[ind].mjd
    ; Correct for airmass extinction effect
    zpterm -= poly(am,zpstr[i].amcoef)
    ; Remove temporal variations
    ;gd1 = where(abs(zpterm) lt 2,ngd1)
    ;si = sort(mjd[gd1])
    ;BINDATA,mjd[gd1[si]],zpterm[gd1[si]],xbin,ybin,binsize=60,/med,gdind=gdind
    ;xbin = xbin[gdind] & ybin=ybin[gdind]
    ; this doesn't work so well

    ;sigzp = mad(str[ind].zpterm)
    ; We are using ADDITIVE zpterm 
    ;  calmag = instmag + zpterm
    ; if there are clouds then instmag is larger/fainter
    ;  and zpterm is smaller (more negative)
    ;bdind = where(str[ind].zpterm-medzp lt -zpthresh[i],nbdind)
    gdind = where(zpterm ge -zpstr[i].thresh and zpterm le zpstr[i].thresh,ngdind,comp=bdind,ncomp=nbdind)
    print,'  ',strtrim(nbdind,2),' exposures with ZPTERM below the threshold'    
    if nbdind gt 0 then badmask[ind[gdind]] = 0
  endif
endfor
; Get bad DECaLS and SMASH exposures
badexp = bytarr(n_elements(str))
READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/smash_badexposures.txt',smashexpnum,format='A',comment='#',/silent
MATCH,long(str.expnum),long(smashexpnum),ind1,ind2,/sort,count=nmatch
badexp[ind1] = 1
badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'c4d')   ; make sure they are DECam exposures
READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/decals_bad_expid.txt',decalsexpnum,format='A',comment='#',/silent
MATCH,long(str.expnum),long(decalsexpnum),ind1,ind2,/sort,count=nmatch
badexp[ind1] = 1
badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'c4d')   ; make sure they are DECam exposures
READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/mzls_bad_expid.txt',mzlsexpnum,format='A',comment='#',/silent
MATCH,long(str.expnum),long(mzlsexpnum),ind1,ind2,/sort,count=nmatch
badexp[ind1] = 1
badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'k4m')   ; make sure they are Mosaic3 exposures

; Final QA cuts
;  Many of the short u-band exposures have weird ZPTERMs, not sure why
;  There are a few exposures with BAD WCS, RA>360!
bdexp = where(str.fwhm gt fwhmthresh or str.ra gt 360 or badmask eq 1 or badexp eq 1 or $
              str.rarms gt 0.2 or str.decrms gt 0.2 or $
              (str.instrument eq 'c4d' and str.zpspatialvar_nccd gt 5 and str.zpspatialvar_rms gt 0.1),nbdexp)
; rarms/decrms, nrefmatch
print,'QA cuts remove ',strtrim(nbdexp,2),' exposures'
; Remove
torem = bytarr(nchstr)
for i=0,nbdexp-1 do torem[str[bdexp[i]].chipindx:str[bdexp[i]].chipindx+str[bdexp[i]].nchips-1]=1
bdchstr = where(torem eq 1,nbdchstr)
REMOVE,bdchstr,chstr
REMOVE,bdexp,str
; Get new CHIPINDEX values
;   make two arrays of old and new indices to transfer 
;   the new index values into an array with the size of
;   the old CHSTR
trimoldindex = lindgen(nchstr)                    ; index into original array, but "bad" ones removed/trimed
remove,bdchstr,trimoldindex
trimnewindex = lindgen(n_elements(trimoldindex))  ; new index of trimmed array
newindex = lonarr(nchstr)-1
newindex[trimoldindex] = trimnewindex             ; new index in original array
newchipindex = newindex[str.chipindx]
str.chipindx = newchipindex
nstr = n_elements(str)
nchstr = n_elements(chstr)

; KLUDGE!!!!!!
print,'KLUDGE!!!! CREATING FAKE DEPTH INFORMATION!!'
add_tag,str,'depth95',99.99,str
add_tag,str,'depth10sig',99.99,str
str.depth95 = 20.+2.5*alog(str.exptime)+randomn(seed,nstr)*0.2
str.depth10sig = str.depth95

; CREATE LIST OF HEALPIX AND OVERLAPPING CHIPS
nside2 = 4096
nhealstr = 100000000LL
healstr = replicate({pix:0L,chipindx:-1LL},nhealstr)
cnt = 0LL
for i=0,nchstr-1 do begin
  if i mod 5e3 eq 0 then print,i

  theta = (90-chstr[i].vdec)/radeg
  phi = chstr[i].vra/radeg
  ANG2VEC,theta,phi,vec
  QUERY_POLYGON,nside2,vec,listpix,nlistpix,/inclusive

  ; Add new elements to array
  if cnt+nlistpix gt nhealstr then begin
    print,'Adding new elements'
    old = healstr
    nnew = 10000000LL ; 1e7
    healstr = replicate({pix:0L,chipindx:-1LL},nhealstr+nnew)
    healstr[0:nhealstr-1] = old
    nhealstr += nnew
    undefine,old
  endif

  ; Add to the structure
  healstr[cnt:cnt+nlistpix-1].pix = listpix
  healstr[cnt:cnt+nlistpix-1].chipindx = i
  cnt += nlistpix

endfor
; Trim extra elements
healstr = healstr[0:cnt-1]
nhealstr = n_elements(healstr)

; Get uniq pixels
ui = uniq(healstr.pix,sort(healstr.pix))
upix = healstr[ui].pix
nupix = n_elements(upix)
print,strtrim(nupix,2),' Healpix pixels have overlapping chip data'

; Get start/stop indices for each pixel
idx = sort(healstr.pix)
healstr = healstr[idx]
q = healstr.pix
lo = where(q ne shift(q,1),nlo)
;hi = where(q ne shift(q,-1))
hi = [lo[1:nlo-1]-1,nhealstr-1]
nchips = hi-lo+1
index = replicate({pix:0L,lo:0L,hi:0L,nchips:0L},nupix)
index.pix = upix
index.lo = lo
index.hi = hi
index.nchips = nchips
nindex = n_elements(index)

stop

; Initialize the coverage structure
nallpix = nside2npix(nside2)
covstr = replicate({pix:0L,coverage:0.0,ucoverage:0.0,udepth:99.99,gcoverage:0.0,gdepth:99.99,$
                    rcoverage:0.0,rdepth:99.99,icoverage:0.0,idepth:99.99,zcoverage:0.0,$
                    zdepth:99.99,ycoverage:0.0,ydepth:99.99,vrcoverage:9.9,vrdepth:99.99},nallpix)
covtags = tag_names(covstr)

; Now loop over each healpix
for i=0,nindex-1 do begin
  pix1 = index[i].pix

  ; Get pixel coordinates and boundary
  PIX2ANG_RING,nside2,pix1,htheta,hphi
  hcenra = hphi*radeg
  hcendec = 90-htheta*radeg
  PIX2VEC_RING,nside2,pix1,hvec,hvertex
  hvertex = transpose(reform(hvertex))   ; [1,3,4] -> [4,3]
  VEC2ANG,hvertex,hdec,hra,/astro

  ; Get all of the chips that overlap this healpix
  chstr1 = chstr[healstr[index[i].lo:index[i].hi].chipindx]
  nchstr1 = n_elements(chstr1)

  ; Transform to tangent plane system
  ROTSPHCEN,hra,hdec,hcenra,hcendec,hlon,hlat,/gnomic
  ROTSPHCEN,chstr1.vra,chstr1.vdec,hcenra,hcendec,vlon,vlat,/gnomic

  ; Transfor to pixel-based tangent plane system
  rlon = minmax([ (vlon)(*), hlon])
  rlat = minmax([ (vlat)(*), hlat])
  lon0 = rlon[0]
  lat0 = rlat[0]
  dx = 1d-3  ; gives ~20x20 pixel
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
  
  ; Get unique filters
  uifilt = uniq(chstr1.filter,sort(chstr1.filter))
  ufilt = chstr1[uifilt].filter
  nfilt = n_elements(ufilt)

  ; Loop over filters
  for j=0,nfilters-1 do begin
    filtind = where(chstr1.filter eq ufilt[i],nfiltind)
    ; Make the filter coverage map
    filtim = bytarr(nx,ny,nfiltind)

    ; Loop over the chips
    overlap = fltarr(nfiltind)
    for k=0,nfiltind-1 do begin
      in1 = polyfillv(vx[filtind[k],*],vy[filtind[k],*],nx,ny)
      inim1 = bytarr(nx,ny)
      inim1[in1] = 1
      filtim[*,*,k] = (inim AND inim1)
      overlap[k] = float(total(filtim[*,*,k]))/ninside  ; fractional coverage
    endfor  ; chip loop


    stop

    ; Get columns for this filter
    fcovind = where(covtags eq strupcase(chstr1[j].filter)+'COVERAGE',nfcovind)
    fdepind = where(covtags eq strupcase(chstr1[j].filter)+'DEPTH',nfdepind)
    covstr[pix1].coverage >= overlap
    if 
    covstr[pix1].(fcovind) = 
    covstr[pix1].(fdepind) = 
    stop


  endfor  ; filter loop


  stop

endfor

;Calculating area of overlap of two colygons
;-probabaly fastest to create a grid of cartesian points and
; then use roi_cut on both polygons to get map of "inside" pixels
; and use union of both to get overlap
;-could do similar with high-res healpix and query_polygon, not sure
;  which one is faster.
;Figuring out which healpix are fully INSIDE the chip polygon
;-use ispointinpolygon.pro (part of dopolygonsoverlap.pro) on all the healpix
;   vertices



; TWO WAYS TO DO THIS
; 1) Create the chip/healpix index
;      -loop over each healpix, get overlapping chips (only need pix
;         and chipindex)
;      -do finer scale healpix and get depth for each
;      -use that information to calculate final depth and covering fraction
; 2) Loop through all chips
;    -find the healpix it overlaps
;    -calculate covering fraction
;    -Use the depth of this chip if it is the deepest
;      exposure for any of the pixels



stop

end
