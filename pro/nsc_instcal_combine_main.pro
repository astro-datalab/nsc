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
;stop

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
stop
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='nsccmb',jobs=jobs,nmulti=nmulti,wait=10

stop

end
