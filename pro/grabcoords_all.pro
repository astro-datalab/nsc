pro grabcoords_all,version
  
; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v3'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(dir,/directory) eq 0 then file_mkdir,dir+'logs/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir
  
; Load the list of exposures
list1 = MRDFITS(dir+'/lists/decam_instcal_list.fits.gz',1)
list2 = MRDFITS(dir+'/lists/mosaic3_instcal_list.fits.gz',1)
list3 = MRDFITS(dir+'/lists/bok90prime_instcal_list.fits.gz',1)
str = [list1,list2,list3]
undefine,list1,list2,list3
nstr = n_elements(str)
str.fluxfile = strtrim(str.fluxfile,2)
str.maskfile = strtrim(str.maskfile,2)
str.wtfile = strtrim(str.wtfile,2)
print,strtrim(nstr,2),' InstCal images'

; Submit them 100 at a time
fluxfiles = strmid(str.fluxfile,4)
nbin = 300
nbatch = ceil(nstr/float(nbin))
cmd = strarr(nbatch)
dir = strarr(nbatch)+localdir+'dnidever/nsc/instcal/'+version+'/tmp/coords/'
for i=0L,nbatch-1 do begin
  lo = i*nbin
  hi = (lo+nbin-1) < (nstr-1)
  fils = fluxfiles[lo:hi]
  cmd[i] = 'grabcoords,["'+strjoin(fils,'","')+'"]'
endfor

stop

pbs_daemon,cmd,dir,jobs=jobs,/hyperthread,/idle,prefix='coords',wait=0.1,nmulti=20,verbose=0

stop

end
