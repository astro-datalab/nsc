pro nsc_instcal_combine_main

; Combine all of the data

dir = '/datalab/users/dnidever/decamcatalog/instcal/'
nside = 256
nmulti = 10
radeg = 180.0d0 / !dpi

; Restore the calibration summary file
str0 = mrdfits(dir+'nsc_instcal_calibrate.fits',1,/silent)
gd = where(str0.success eq 1,ngd)
str = str0[gd]
nstr = n_elements(str)
str.file = strtrim(str.file,2)
str.base = strtrim(str.base,2)

; Which healpix pixels have data
print,'Finding the Healpix pixels with data'
radius = 1.1
healstr = replicate({base:'',pix:0L},1e5)
nhealstr = n_elements(healstr)
cnt = 0LL
for i=0,nstr-1 do begin
  if i mod 1e3 eq 0 then print,i
  head = headfits(str[i].file,exten=0)
  sra = sxpar(head,'ra')
  sdec = sxpar(head,'dec')
  ra = sexig2ten(sra)*15.0d0
  dec = sexig2ten(sdec)
  theta = (90-dec)/radeg
  phi = ra/radeg
  ang2vec,theta,phi,vec
  query_disc,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

  ; Add new elements to array
  if cnt+nlistpix gt nhealstr then begin
    old = healstr
    healstr = replicate({base:'',pix:0L},nhealstr+1e4)
    healstr[0:nhealstr-1] = old
    nhealstr += 1e4
    undefine,old
  endif

  healstr[cnt:cnt+nlistpix-1].base = str[i].base
  healstr[cnt:cnt+nlistpix-1].pix = listpix
  cnt += nlistpix

  ;stop

endfor
; Trim extra elements
healstr = healstr[0:cnt-1]
nhealstr = n_elements(healstr)

; Get uniq pixels
ui = uniq(healstr.pix,sort(healstr.pix))
upix = healstr[ui].pix
nupix = n_elements(upix)

; Get start/stop indices for each pixel
idx = sort(healstr.pix)
q = healstr[idx].pix
lo = where(q ne shift(q,1),nlo)
;hi = where(q ne shift(q,-1))
hi = [lo[1:nlo-1]-1,nhealstr-1]
nexp = hi-lo+1

; Loop over each pixels and get list of overlapping exposures
for i=0,nupix-1 do begin
  healstr1 = healstr[idx[lo[i]:hi[i]]]
  ; Create a file with this list
  MWRFITS,healstr1,dir+'combine/lists/'+strtrim(upix[i],2)+'.fits',/create
endfor

; Now run the combination program on each healpix pixel
cmd = 'nsc_instcal_combine,"'+strtrim(upix,2)+'"'
dirs = strarr(nupix)+'/data0/dnidever/decamcatalog/tmp/'
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='nsccmb',jobs=jobs,nmulti=nmulti,wait=10

stop

end
