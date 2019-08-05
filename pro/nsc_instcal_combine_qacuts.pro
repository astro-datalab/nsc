pro nsc_instcal_combine_qacuts,version,redo=redo

if n_elements(version) eq 0 then begin
  print,'Syntax - nsc_instcal_combine_qacuts,version,redo=redo'
  return
endif

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
if file_test(localdir+'dnidever/nsc/instcal/'+version+'/') eq 0 then file_mkdir,localdir+'dnidever/nsc/instcal/'+version+'/'
plotsdir = dir+'plots/'
if file_test(plotsdir,/directory) eq 0 then file_mkdir,plotsdir
radeg = 180.0d0 / !dpi


; Restore the calibration summary file
temp = MRDFITS(dir+'lists/nsc_calibrate_summary.fits.gz',1,/silent)
schema = temp[0]
struct_assign,{dum:''},schema
schema = create_struct(schema,'chipindx',-1LL,'NGOODCHIPWCS',0)
str = replicate(schema,n_elements(temp))
struct_assign,temp,str,/nozero
str.expdir = strtrim(str.expdir,2)
str.instrument = strtrim(str.instrument,2)
str.metafile = strtrim(str.metafile,2)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)
str.filter = strtrim(str.filter,2)
; Add WCSCAL and TELSTAT information
add_tag,str,'wcscal','',str
add_tag,str,'telstat','',str
coords = MRDFITS(dir+'lists/allcoords.fits',1)
coords.file = strtrim(coords.file,2)
coords.wcscal = strtrim(coords.wcscal,2)
coords.telstat = strtrim(coords.telstat,2)
fluxfile = str.file
g = where(strmid(fluxfile,0,4) eq '/net',ng)
if ng gt 0 then fluxfile[g]=strmid(fluxfile[g],4)
MATCH,fluxfile,coords.file,ind1,ind2,/sort
str[ind1].wcscal = coords[ind2].wcscal    ; Failed (3153), Poor (14), Successful (308190)
str[ind1].telstat = coords[ind2].telstat  ; NAN (68188), Not (1222), Track (241826), UNKNOWN (116), Unknown (5)
; the 2054 failed exposures did not match b/c no fluxfile info
; Only want exposures with successful SE processing
gd = where(str.success eq 1,nstr)
str = str[gd]
si = sort(str.expdir)
str = str[si]
chstr = mrdfits(dir+'lists/nsc_calibrate_summary.fits.gz',2,/silent)
chstr.expdir = strtrim(chstr.expdir,2)
chstr.instrument = strtrim(chstr.instrument,2)
nchstr = n_elements(chstr)
; Get indices for CHSTR
siexp = sort(chstr.expdir)
chstr = chstr[siexp]
expdir = chstr.expdir
brklo = where(expdir ne shift(expdir,1),nbrk)
brkhi = [brklo[1:nbrk-1]-1,n_elements(expdir)-1]
nchexp = brkhi-brklo+1
if nstr ne n_elements(brklo) then stop,'number of exposures in STR and CHSTR do not match'
str.chipindx = brklo
str.nchips = nchexp
; Getting number of good chip WCS for each exposures
for i=0,n_elements(str)-1 do str[i].ngoodchipwcs = total(chstr[brklo[i]:brkhi[i]].ngaiamatch gt 0)
; Fixing absolute paths of flux filename
file = str.file
g1 = where(stregex(file,'/net/mss1/',/boolean) eq 1,ng1)
if ng1 gt 0 then file[g1] = strmid(file[g1],10)
g2 = where(stregex(file,'/mss1/',/boolean) eq 1,ng2)
if ng2 gt 0 then file[g2] = strmid(file[g2],6)
; Fixing very negative RAs
print,'FIXING NEGATIVE RAs in STR and CHSTR'
;bdra = where(chstr.cenra lt -180,nbdra)
bdra = where(chstr.cenra lt -0,nbdra)
uibd = uniq(chstr[bdra].expdir,sort(chstr[bdra].expdir))
MATCH,str.expdir,chstr[bdra[uibd]].expdir,ind1,ind2,/sort,count=nmatch
for i=0,nmatch-1 do begin
  MATCH,chstr[bdra].expdir,str[ind1[i]].expdir,ind3,ind4,/sort
  ; Fix STR RA
  chra = chstr[bdra[ind3]].cenra
  bd1 = where(chra lt -180,nbd1)
  if nbd1 gt 0 then chra[bd1]+=360
  cenra = mean(minmax(chra))
  if cenra lt 0 then cenra+=360
  str[ind1[i]].ra = cenra
  ; Fix CHSTR CENRA
  bd2 = where(chra lt 0,nbd2)
  if nbd2 gt 0 then chra[bd2]+=360
  chstr[bdra[ind3]].cenra = chra
  ; Fix CHSTR VRA
  vra = chstr[bdra[ind3]].vra
  bd3 = where(vra lt 0,nbd3)
  if nbd3 gt 0 then vra[bd3]+=360
  chstr[bdra[ind3]].vra = vra
endfor

; Fix instrument in STR and CHSTR
print,'FIXING INSTRUMENT IN STR AND CHSTR'
type = ['c4d','k4m','ksb']
for i=0,n_elements(type)-1 do begin
  gd = where(stregex(str.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
  if ngd gt 0 then str[gd].instrument=type[i]
  gd = where(stregex(chstr.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
  if ngd gt 0 then chstr[gd].instrument=type[i]
endfor

;; Fix missing AIRMASS                                                                                                                                                           
;bdam = where(str.airmass lt 0.9,nbdam)
;for i=0,nbdam-1 do begin
;  type = ['c4d','k4m','ksb']
;  obs = ['ctio','kpno','kpno']
;  MATCH,str[bdam[i]].instrument,type,ind1,ind2,/sort
;  obsname = obs[ind2]
;  OBSERVATORY,obsname,obstr
;  lat = obstr.latitude
;  lon = obstr.longitude
;  jd = date2jd(str[bdam[i]].dateobs)
;  ra = str[bdam[i]].ra
;  dec = str[bdam[i]].dec
;  str[bdam[i]].airmass = AIRMASS(jd,ra,dec,lat,lon)
;endfor
; THIS IS STILL RETURNING -1, IS ONE OF THE VALUES WRONG??

; This is now done at the beginning when the lists are created
;; APPLY RELEASE-DATE CUTS
;list1 = MRDFITS(dir+'lists/decam_instcal_list.fits',1)
;list2 = MRDFITS(dir+'lists/mosaic3_instcal_list.fits',1)
;list3 = MRDFITS(dir+'lists/bok90prime_instcal_list.fits',1)
;list = [list1,list2,list3]
;list.fluxfile = strtrim(list.fluxfile,2)
;fluxfile = strmid(list.fluxfile,10)
;MATCH,fluxfile,file,ind1,ind2,/sort,count=nmatch
;; some don't match because they were from a previous version
;;  of the input list
;release_date = strarr(n_elements(str))+'2020-01-01 00:00:00'
;release_date[ind2] = list[ind1].release_date
;release_year = long(strmid(release_date,0,4))
;release_month = long(strmid(release_date,5,2))
;release_day = long(strmid(release_date,8,2))
;release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
;;release_cutoff = [2017,4,24]  ; v1 - April 24, 2017
;release_cutoff = [2017,10,11]  ; v2 - Oct 11, 2017
;release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
;gdrelease = where(release_mjd le release_cutoff_mjd,ngdrelease,comp=bdrelease,ncomp=nbdrelease)
;print,strtrim(ngdrelease,2),' exposures are PUBLIC'
;str = str[gdrelease]  ; impose the public data cut

; Zero-point structure
zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},10)
zpstr[0:6].instrument = 'c4d'
zpstr[0:6].filter = ['u','g','r','i','z','Y','VR']
zpstr[0].amcoef = [-1.60273, -0.375253]   ; c4d-u
zpstr[1].amcoef = [0.277124, -0.198037]   ; c4d-g
zpstr[2].amcoef = [0.516382, -0.115443]   ; c4d-r
zpstr[3].amcoef = [0.380338, -0.067439]   ; c4d-i
zpstr[4].amcoef = [0.123924, -0.096877]   ; c4d-z
zpstr[5].amcoef = [-1.06529, -0.051967]   ; c4d-Y
zpstr[6].amcoef = [1.004357, -0.081105]   ; c4d-VR
; Mosiac3 z-band
zpstr[7].instrument = 'k4m'
zpstr[7].filter = 'z'
zpstr[7].amcoef = [-2.687201, -0.73573]   ; k4m-z
; Bok 90Prime, g and r
zpstr[8].instrument = 'ksb'
zpstr[8].filter = 'g'
zpstr[8].amcoef = [-2.859646, -1.40837]   ; ksb-g
zpstr[9].instrument = 'ksb'
zpstr[9].filter = 'r'
zpstr[9].amcoef = [-4.008771, -0.25718]   ; ksb-r
nzpstr = n_elements(zpstr)

STOP,'DOUBLE-CHECK THESE ZERO-POINTS!!!'

; APPLY QA CUTS IN ZEROPOINT AND SEEING
If not keyword_set(nocuts) then begin
  print,'APPLYING QA CUTS'
  ;fwhmthresh = 3.0  ; arcsec, v1
  fwhmthresh = 2.0  ; arcsec, v2
  ;filters = ['u','g','r','i','z','Y','VR']
  ;nfilters = n_elements(filters)
  ;zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
  ;zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
  badzpmask = bytarr(n_elements(str)) + 1

  for i=0,nzpstr-1 do begin
    ind = where(str.instrument eq zpstr[i].instrument and str.filter eq zpstr[i].filter and str.success eq 1,nind)
    print,zpstr[i].instrument,'-',zpstr[i].filter,' ',strtrim(nind,2),' exposures'
    if nind gt 0 then begin
      str1 = str[ind]
      zpterm = str1.zpterm
      bdzp = where(finite(zpterm) eq 0,nbdzp)  ; fix Infinity/NAN
      if nbdzp gt 0 then zpterm[bdzp] = 999999.9
      am = str1.airmass
      mjd = str1.mjd
      bdam = where(am lt 0.9,nbdam)
      if nbdam gt 0 then am[bdam] = median(am)
      glactc,str1.ra,str1.dec,2000.0,glon,glat,1,/deg

      ; Measure airmass dependence
      gg0 = where(abs(zpterm) lt 50 and am lt 2.0,ngg0)
      coef0 = robust_poly_fitq(am[gg0],zpterm[gg0],1)
      zpf = poly(am,coef0)
      sig0 = mad(zpterm[gg0]-zpf[gg0])
      gg = where(abs(zpterm-zpf) lt (3.5*sig0 > 0.2),ngg)
      coef = robust_poly_fitq(am[gg],zpterm[gg],1)
      print,zpstr[i].instrument+'-'+zpstr[i].filter,' ',coef
      ; Trim out bad exposures to determine the correlations and make figures
      gg = where(abs(zpterm-zpf) lt (3.5*sig0 > 0.2) and str1.airmass lt 2.0 and str1.fwhm lt 2.0 and str1.rarms lt 0.15 and $
                 str1.decrms lt 0.15 and str1.success eq 1 and str1.wcscal eq 'Successful' and str1.zptermerr lt 0.05 and $
                 str1.zptermsig lt 0.08 and (str1.ngoodchipwcs eq str1.nchips) and $
                 (str1.instrument ne 'c4d' or str1.zpspatialvar_nccd le 5 or (str1.instrument eq 'c4d' and str1.zpspatialvar_nccd gt 5 and str1.zpspatialvar_rms lt 0.1)) and $
                 abs(glat) gt 10 and str1.nrefmatch gt 100 and str1.exptime ge 30,ngg)

      ; Zpterm with airmass dependence removed
      relzpterm = zpterm + 25   ; 25 to get "absolute" zpterm
      relzpterm -= zpstr[i].amcoef[1]*(am-1)

      ; CURRENTLY K4M/KSB HAVE EXPTIME-DEPENDENCE IN THE ZEROPOINTS!!
      if zpstr[i].instrument eq 'k4m' or zpstr[i].instrument eq 'ksb' then begin
        print,'REMOVING EXPTIME-DEPENDENCE IN K4M/KSB ZEROPOINTS!!!'
        relzpterm += 2.5*alog10(str1.exptime)
      endif

      ; Fit temporal variation in zpterm
      mjd0 = 56200L
      xx = str1[gg].mjd-mjd0
      yy = relzpterm[gg]
      invvar = 1.0/str1[gg].zptermerr^2
      nord = 3
      bkspace = 200 ;20
      sset1 = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=yfit1)
      sig1 = mad(yy-yfit1)
      gd = where(yy-yfit1 gt -3*sig1,ngd)      
      ; refit
      sset = bspline_iterfit(xx[gd],yy[gd],invvar=invvar[gd],nord=nord,bkspace=bkspace)
      yfit = bspline_valu(xx,sset)
      allzpfit = bspline_valu(str1.mjd-mjd0,sset)

      ; Make some figures
      setdisp,/silent
      !p.font = 0
      ; ZPterm vs. airmass
      file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_airmass'
      ps_open,file,/color,thick=4,/encap
      hess,am[gg],relzpterm[gg],dx=0.01,dy=0.02,xr=[0.9,2.5],yr=[-0.5,0.5]+median(relzpterm[gg]),xtit='Airmass',ytit='Zero-point',$
           tit=zpstr[i].instrument+'-'+zpstr[i].filter
      x = scale_vector(findgen(100),0.5,2.0)
      oplot,x,poly(x,coef),co=250
      ps_close
      ps2png,file+'.eps',/eps
      ; ZPterm vs. time (density)
      file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_time_density'
      ps_open,file,/color,thick=4,/encap
      hess,str1[gg].mjd-mjd0,relzpterm[gg],dx=2,dy=0.02,yr=[-0.5,0.5]+median(relzpterm[gg]),xtit='Time (days)',ytit='Zero-point',$
           tit=zpstr[i].instrument+'-'+zpstr[i].filter
      oplot,str1[gg].mjd-mjd0,allzpfit[gg],ps=1,sym=0.3,co=250
      xyouts,50,-0.45+median(relzpterm[gg]),'MJD!d0!n = '+strtrim(mjd0,2),align=0,charsize=1.2
      ps_close
      ps2png,file+'.eps',/eps
      ; ZPterm vs. time (points)
      file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_time'
      ps_open,file,/color,thick=4,/encap
      plot,str1[gg].mjd-mjd0,relzpterm[gg],ps=1,sym=0.5,yr=[-0.5,0.5]+median(relzpterm[gg]),xs=1,ys=1,xtit='Time (days)',ytit='Zero-point',$
           tit=zpstr[i].instrument+'-'+zpstr[i].filter,thick=1
      oplot,str1[gg].mjd-mjd0,allzpfit[gg],ps=1,sym=0.3,co=250
      xyouts,50,-0.45+median(relzpterm[gg]),'MJD!d0!n = '+strtrim(mjd0,2),align=0,charsize=1.2
      ps_close
      ps2png,file+'.eps',/eps

      ; Remove temporal variations to get residual values
      relzpterm -= allzpfit


      ; Find the GOOD exposures
      ;------------------------
      ; We are using ADDITIVE zpterm 
      ;  calmag = instmag + zpterm
      ; if there are clouds then instmag is larger/fainter
      ;  and zpterm is smaller (more negative)
      ;bdind = where(str[ind].zpterm-medzp lt -zpthresh[i],nbdind)
      gdind = where(relzpterm ge -zpstr[i].thresh and relzpterm le zpstr[i].thresh,ngdind,comp=bdind,ncomp=nbdind)
      print,'  ',strtrim(nbdind,2),' exposures with ZPTERM below the threshold'    
      if ngdind gt 0 then badzpmask[ind[gdind]] = 0

    endif
  endfor
  ; Get bad DECaLS and SMASH exposures
  badexp = bytarr(n_elements(str))
  READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/smash_badexposures.txt',smashexpnum,format='A',comment='#',/silent
  MATCH,long(str.expnum),long(smashexpnum),ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    badexp[ind1] = 1
    badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'c4d')   ; make sure they are DECam exposures
  endif
  READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/decals_bad_expid.txt',decalsexpnum,format='A',comment='#',/silent
  MATCH,long(str.expnum),long(decalsexpnum),ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    badexp[ind1] = 1
    badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'c4d')   ; make sure they are DECam exposures
  endif
  READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/mzls_bad_expid.txt',mzlsexpnum,format='A',comment='#',/silent
  MATCH,long(str.expnum),long(mzlsexpnum),ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    badexp[ind1] = 1
    badexp[ind1] = badexp[ind1] AND (str[ind1].instrument eq 'k4m')   ; make sure they are Mosaic3 exposures
  endif

  ; Final QA cuts
  ;  Many of the short u-band exposures have weird ZPTERMs, not sure why
  ;  There are a few exposures with BAD WCS, RA>360!
  bdexp = where(str.success eq 0 or $                              ; SE failure
                str.wcscal ne 'Successful' or $                    ; CP WCS failure
                str.fwhm gt fwhmthresh or $                        ; bad seeing
                str.ra gt 360 or $                                 ; bad WCS/coords
                str.rarms gt 0.15 or str.decrms gt 0.15 or $       ; bad WCS
                badzpmask eq 1 or $                                ; bad ZPTERM
                str.zptermerr gt 0.05 or $                         ; bad ZPTERMERR
                str.nrefmatch lt 5 or $                            ; few phot ref match
                badexp eq 1 or $                                   ; bad SMASH/LS exposure
                ;str.ngoodchipwcs lt str.nchips or $                ; not all chips astrom calibrated
                (str.instrument eq 'c4d' and str.zpspatialvar_nccd gt 5 and str.zpspatialvar_rms gt 0.1),nbdexp)  ; bad spatial zpterm
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

; SHOULD INCLUDE CUTS ON ZTERMERR OR NPHOTMATCH
STOP,'SHOULD INCLUDE CUTS ON ZTERMERR OR NPHOTMATCH'

Endif else print,'SKIPPING QA CUTS'

; CREATE LIST OF HEALPIX AND OVERLAPPING EXPOSURES
; Which healpix pixels have data
listfile = dir+'lists/nsc_instcal_combine_healpix_list.fits'
if file_test(listfile) eq 0 or keyword_set(redo) then begin
  print,'Finding the Healpix pixels with data'
  radius = 1.1
  healstr = replicate({file:'',base:'',pix:0L},1e5)
  nhealstr = n_elements(healstr)
  cnt = 0LL
  for i=0,nstr-1 do begin
    if i mod 1e3 eq 0 then print,i
    theta = (90-str[i].dec)/radeg
    phi = str[i].ra/radeg
    ANG2VEC,theta,phi,vec
    QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

    ; Use the chip corners to figure out which ones actually overlap
    chstr1 = chstr[str[i].chipindx:str[i].chipindx+str[i].nchips-1]
    ;  rotate to tangent plane so it can handle RA=0/360 and poles properly
    ROTSPHCEN,chstr1.vra,chstr1.vdec,str[i].ra,str[i].dec,vlon,vlat,/gnomic
    ;  loop over healpix
    overlap = bytarr(nlistpix)
    for j=0,nlistpix-1 do begin
      PIX2VEC_RING,nside,listpix[j],vec,vertex
      vertex = transpose(reform(vertex))  ; [1,3,4] -> [4,3]
      VEC2ANG,vertex,hdec,hra,/astro
      ROTSPHCEN,hra,hdec,str[i].ra,str[i].dec,hlon,hlat,/gnomic
      ;  loop over chips
      for k=0,str[i].nchips-1 do overlap[j] >= DOPOLYGONSOVERLAP(hlon,hlat,vlon[*,k],vlat[*,k])
    endfor
    ; Only keep the healpix with real overlaps
    gdlistpix = where(overlap eq 1,ngdlistpix)
    if ngdlistpix gt 0 then begin
      listpix = listpix[gdlistpix]
      nlistpix = ngdlistpix
    endif else begin
      undefine,listpix
      nlistpix = 0
    endelse

if nlistpix eq 0 then stop,'No healpix for this exposure.  Something is wrong!'

    ; Add new elements to array
    if cnt+nlistpix gt nhealstr then begin
      old = healstr
      healstr = replicate({file:'',base:'',pix:0L},nhealstr+1e4)
      healstr[0:nhealstr-1] = old
      nhealstr += 1e4
      undefine,old
    endif

    ; Add to the structure
    healstr[cnt:cnt+nlistpix-1].file = str[i].expdir+'/'+str[i].base+'_cat.fits'
    healstr[cnt:cnt+nlistpix-1].base = str[i].base
    healstr[cnt:cnt+nlistpix-1].pix = listpix
    cnt += nlistpix

  endfor
  ; Trim extra elements
  healstr = healstr[0:cnt-1]
  nhealstr = n_elements(healstr)

  ; Get uniq pixels
  ui = uniq(healstr.pix,sort(healstr.pix))
  upix = healstr[ui].pix
  nupix = n_elements(upix)
  print,strtrim(nupix,2),' Healpix pixels have overlapping data'

  ; Get start/stop indices for each pixel
  idx = sort(healstr.pix)
  healstr = healstr[idx]
  q = healstr.pix
  lo = where(q ne shift(q,1),nlo)
  ;hi = where(q ne shift(q,-1))
  hi = [lo[1:nlo-1]-1,nhealstr-1]
  nexp = hi-lo+1
  index = replicate({pix:0L,lo:0L,hi:0L,nexp:0L},nupix)
  index.pix = upix
  index.lo = lo
  index.hi = hi
  index.nexp = nexp
  npix = n_elements(index)

  ; Replace /net/dl1/ with /dl1/ so it will work on all machines
  healstr.file = repstr(healstr.file,'/net/dl1/','/dl1/')

  ; Write the full list plus an index
  print,'Writing list to ',listfile
  MWRFITS,healstr,listfile,/create
  MWRFITS,index,listfile,/silent
  if file_test(listfile+'.gz') eq 1 then file_delete,listfile+'.gz',/allow
  spawn,['gzip',listfile],/noshell
  ; PUT NSIDE IN HEADER!!

endif

stop

end
