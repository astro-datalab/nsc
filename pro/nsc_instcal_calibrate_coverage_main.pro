pro nsc_instcal_calibrate_coverage_main,lohi,outfile

;; combine all of the individual exposure coverage files

version = 'v3'

;; list of "good" exposures
liststr = mrdfits('/dl1/users/dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_combine_healpix_list.fits.gz',1)
liststr.file = strtrim(liststr.file,1)
expdir = file_dirname(liststr.file)
; get unique ones
ui = uniq(expdir,sort(expdir))
expdir = expdir[ui]
; only do a fraction of them
if n_elements(lohi) gt 0 then expdir=expdir[lohi[0]:lohi[1]]
nexpdir = n_elements(expdir)
print,strtrim(nexpdir,2),' good exposures'

nside = 512L
radeg = 180.0d0/!dpi
covstr = replicate({pix:0L,ra:0.0d0,dec:0.0d0,nobj:0L,nmeas:0LL,depth:0.0,nexp:0L,unexp:0,udepth:-9999.0,gnexp:0,gdepth:-9999.0,$
                    rnexp:0,rdepth:0.0,inexp:0,idepth:-9999.0,znexp:0,zdepth:-9999.0,ynexp:0,ydepth:-9999.0,vrnexp:0,vrdepth:-9999.0},3145728L)
covstr.pix = lindgen(3145728L)
pix2ang_ring,nside,covstr.pix,theta,phi
covstr.dec = 90-theta*radeg
covstr.ra = phi*radeg

nmeas = 0L
for i=0,nexpdir-1 do begin
  if i mod 1000 eq 0 then print,i
  expdir1 = expdir[i]
  base = file_basename(expdir1)
  ;; Load exposure meta data information
  metafile = expdir1+'/'+base+'_meta.fits'
  if file_test(metafile) eq 0 then goto,BOMB
  metastr = mrdfits(metafile,1,/silent)
  ;; Load exposure coverage information
  cfile = expdir1+'/'+base+'_hlpmeta.fits'
  if file_test(cfile) eq 0 then goto,BOMB
  cstr = mrdfits(cfile,1,/silent)

  nmeas += total(cstr.nmeas>0)

  ;; nmeas
  covstr[cstr.pix].nmeas += cstr.nmeas
  ;; depth
  covstr[cstr.pix].depth >= cstr.depth95
  ;; nexp
  covstr[cstr.pix].nexp++
  ;; nobj
  covstr[cstr.pix].nobj >= cstr.nmeas
  ;; Add filter depth and nexp
  case metastr.filter of
  'u': begin
         covstr[cstr.pix].unexp++
         covstr[cstr.pix].udepth >= cstr.depth95
      end
  'g': begin
         covstr[cstr.pix].gnexp++
         covstr[cstr.pix].gdepth >= cstr.depth95
      end
  'r': begin
         covstr[cstr.pix].rnexp++
         covstr[cstr.pix].rdepth >= cstr.depth95
      end
  'i': begin
         covstr[cstr.pix].inexp++
         covstr[cstr.pix].idepth >= cstr.depth95
      end
  'z': begin
         covstr[cstr.pix].znexp++
         covstr[cstr.pix].zdepth >= cstr.depth95
      end
  'Y': begin
         covstr[cstr.pix].ynexp++
         covstr[cstr.pix].ydepth >= cstr.depth95
      end
  'VR': begin
         covstr[cstr.pix].vrnexp++
         covstr[cstr.pix].vrdepth >= cstr.depth95
       end
  else:
  endcase

  BOMB:
endfor

print,strtrim(nmeas,2),' total measurements'

;; Save the coverage file
if n_elements(outfile) eq 0 then $
  outfile = '/dl1/users/dnidever/nsc/instcal/'+version+'/lists/nsc_instcal_calibrate_coverage.fits'
print,'Writing calib coverage information to ',outfile
MWRFITS,covstr,outfile,/create

;stop

end
