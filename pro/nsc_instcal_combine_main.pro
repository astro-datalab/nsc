pro nsc_instcal_combine_main,redo=redo,nmulti=nmulti

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir
dir = dldir+'users/dnidever/decamcatalog/instcal/'
;dir = '/datalab/users/dnidever/decamcatalog/instcal/'
;nside = 256
nside = 128
nmulti = 30
radeg = 180.0d0 / !dpi

; Restore the calibration summary file
str0 = mrdfits(dir+'nsc_instcal_calibrate.fits',1,/silent)
gd = where(str0.success eq 1,ngd)
str = str0[gd]
nstr = n_elements(str)
str.expdir = strtrim(str.expdir,2)
str.metafile = strtrim(str.metafile,2)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)
str.filter = strtrim(str.filter,2)

; APPLY RELEASE-DATE CUTS

; APPLY QA CUTS IN ZEROPOINT AND SEEING
print,'APPLY QA CUTS IN ZEROPOINT AND SEEING'
fwhmthresh = 3.0  ; arcsec
filters = ['u','g','r','i','z','Y','VR']
nfilters = n_elements(filters)
zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
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
REMOVE,bdexp,str
nstr = n_elements(str)

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

; Make the commands
cmd = "nsc_instcal_combine,"+strtrim(upix,2)+",nside="+strtrim(nside,2)
if keyword_set(redo) then cmd+=',/redo'
cmddir = strarr(nupix)+localdir+'dnidever/nsc/instcal/tmp/'
;dirs = strarr(nupix)+'/data0/dnidever/decamcatalog/tmp/'

; Check if the output file exists
if not keyword_set(redo) then begin
  outfiles = dir+'combine/'+strtrim(upix,2)+'.fits'
  test = file_test(outfiles)
  gd = where(test eq 0,ngd,comp=bd,ncomp=nbd)
  if nbd gt 0 then begin
    print,strtrim(nbd,2),' files already exist and /redo not set.'
  endif 
  if ngd eq 0 then begin
    print,'No files to process'
    return
  endif
  print,strtrim(ngd,2),' files left to process'
  cmd = cmd[gd]
  cmddir = cmddir[gd]
endif

stop

; Now run the combination program on each healpix pixel
PBS_DAEMON,cmd,cmddir,/hyperthread,/idle,prefix='nsccmb',jobs=jobs,nmulti=nmulti,wait=1

stop

end
