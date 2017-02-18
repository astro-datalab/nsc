pro nsc_instcal_combine_main

; Combine all of the data

dir = '/datalab/users/dnidever/decamcatalog/instcal/'
;nside = 256
nside = 128
nmulti = 10
radeg = 180.0d0 / !dpi

; Restore the calibration summary file
str0 = mrdfits(dir+'nsc_instcal_calibrate.fits',1,/silent)
gd = where(str0.success eq 1,ngd)
str = str0[gd]
nstr = n_elements(str)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)

; APPLY QA CUTS IN ZEROPOINT AND SEEING
print,'APPLY QA CUTS IN ZEROPOINT AND SEEING'
fwhmthresh = 3.0  ; arcsec
filters = ['u','g','r','i','z','Y','VR']
nfilters = n_elements(filters)
zpthresh = [1.0,1.0,1.0,1.0,1.0,1.0,1.0]
badmask = lonarr(n_elements(str))
for i=0,nfilters-1 do begin
  ind = where(str.filter eq filters[i],nind)
  print,filters[i],' ',strtrim(nind,2),' exposures'
  if nind gt 0 then begin
    medzp = median(str[ind].zpterm)
    sigzp = mad(str[ind].zpterm)
    ; We are using ADDITIVE zpterm 
    ;  calmag = instmag + zpterm
    ; if there are clouds then instmag is larger/fainter
    ;  and zpterm is smaller (more negative)
    bdind = where(str[ind].zpterm-medzp lt -zpthresh[i],nbdind)
    print,'  ',strtrim(nbdind,2),' exposures with ZPTERM below the threshold'    
    if nbdind gt 0 then badmask[ind[bdind]] = 1
  endif
endfor
bdexp = where(str.fwhm gt fwhmthresh or badmask eq 1,nbdexp)
print,'QA cuts remove ',strtrim(nbdexp,2),' exposures'
REMOVE,bdexp,str]

; Which healpix pixels have data
print,'Finding the Healpix pixels with data'
radius = 1.1
healstr = replicate({file:'',base:'',pix:0L},1e5)
nhealstr = n_elements(healstr)
cnt = 0LL
for i=0,nstr-1 do begin
  if i mod 1e3 eq 0 then print,i
  ;head = headfits(str[i].file,exten=0)
  ;sra = sxpar(head,'ra')
  ;sdec = sxpar(head,'dec')
  ;ra = sexig2ten(sra)*15.0d0
  ;dec = sexig2ten(sdec)
  theta = (90-str[i].dec)/radeg
  phi = str[i].ra/radeg
  ANG2VEC,theta,phi,vec
  QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

  ; Add new elements to array
  if cnt+nlistpix gt nhealstr then begin
    old = healstr
    healstr = replicate({file:'',base:'',pix:0L},nhealstr+1e4)
    healstr[0:nhealstr-1] = old
    nhealstr += 1e4
    undefine,old
  endif

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

; Write the full list plus an index
print,'Writing list to ',dir+'combine/healpix_list.fits'
MWRFITS,healstr,dir+'combine/healpix_list.fits',/create
MWRFITS,index,dir+'combine/healpix_list.fits',/silent
; PUT NSIDE IN HEADER!!

;; Loop over each pixel and get list of overlapping exposures
;for i=0,nupix-1 do begin
;  healstr1 = healstr[idx[lo[i]:hi[i]]]
;  ; Create a file with this list
;  ;lines = healstr1.file+'  '+healstr1.base+'  '+healstr1.pix
;  WRITELINE,dir+'combine/lists/'+strtrim(upix[i],2)+'.lst',healstr1.file
;  ;MWRFITS,healstr1,dir+'combine/lists/'+strtrim(upix[i],2)+'.fits',/create
;endfor

; Now run the combination program on each healpix pixel
cmd = "nsc_instcal_combine,'"+strtrim(upix,2)+"',nside="+strtrim(nside,2)
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(nupix)+'/data0/dnidever/decamcatalog/tmp/'
;stop
;PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='nsccmb',jobs=jobs,nmulti=nmulti,wait=10

; Combine everything
schema_obj = {id:'',pix:0L,ra:0.0d0,dec:0.0d0,ndet:0L,umag:0.0,uerr:0.0,ndetu:0,gmag:0.0,$
              gerr:0.0,ndetg:0,rmag:0.0,rerr:0.0,ndetr:0,imag:99.9,ierr:0.0,ndeti:0,$
              zmag:0.0,zerr:0.0,ndetz:0,ymag:0.0,yerr:0.0,ndety:0,vrmag:0.0,vrerr:0.0,$
              ndetvr:0,x2:0.0,y2:0.0,xy:0.0,cxx:0.0,cxy:0.0,cyy:0.0,asemi:0.0,bsemi:0.0,theta:0.0,$
              elongation:0.0,ellipticity:0.0,fwhm:0.0,flags:0,class_star:0.0,ebv:0.0}
obj = replicate(schema_obj,1e7)
nobj = n_elements(obj)
cnt = 0LL
for i=0,nupix-1 do begin
  print,strtrim(i+1,2),' ',index[i].pix
  file = dir+'combine/'+strtrim(index[i].pix,2)+'.fits'
  if file_test(file) eq 1 then begin
    cat1 = MRDFITS(file,1,/silent)
    ncat1 = n_elements(cat1)
    print,'  ',strtrim(ncat1,2),' sources'
    cat1.id = strtrim(index[i].pix,2)+'.'+strtrim(cat1.id,2)

    ; Add new elements
    if cnt+ncat1 gt nobj then begin
      print,'Adding more elements'
      orig = obj
      obj = replicate(schema_obj,nobj+5e6)
      obj[0:nobj-1] = orig
      nobj = n_elements(obj)
      undefine,orig
    endif

    ; Stuff information in
    obj[cnt:cnt+ncat1-1] = cat1
    cnt += ncat1
    ;push,obj,cat1
  endif else begin
    print,file,' NOT FOUND'
  endelse
  ;stop
endfor
; Trim excess elements
obj = obj[0:cnt-1]
nobj = n_elements(obj)
print,strtrim(nobj,2),' objects'

; Write combined catalog
print,'Writing combined catalog to ',dir+'combine/nsc_combine.fits'
MWRFITS,obj,dir+'combine/nsc_combine.fits',/create

; Make figure
setdisp
!p.font = 0
; On sky distribution
file = 'nsc_skycounts'
ps_open,file,/color,thick=4,/encap
hess,obj.ra,obj.dec,dum,dx=0.1,dy=0.1,/log,xtit='RA',ytit='DEC',tit='NOAO Source Catalog on sky distribution'
ps_close
ps2png,file+'.eps',/eps

; CMD
file = 'nsc_gihess'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=7.5,ysize=9.5
gd = where(obj.gmag lt 50 and obj.imag lt 50,ngd)
hess,obj[gd].gmag-obj[gd].imag,obj[gd].gmag,dum,im,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[24,10],/log,/noplot,xarr=xarr,yarr=yarr
displayc,im,xarr,yarr,xtit='g-i',ytit='i',tit='NOAO Source Catalog CMD',charsize=1.3,/yflip,/log
ps_close
ps2png,file+'.eps',/eps

stop

end
