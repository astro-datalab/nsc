pro update_metafile_nmeas_all

;; Add NMEAS to the exposure structure in the metadatafile

version = 'v3'
NSC_ROOTDIRS,dldir,mssdir,localdir,host
hostname = first_el(strsplit(host,'.',/extract))
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'

; Find all of the directories
print,'Getting the exposure directories'
c4d_expdirs = file_search(dir+'c4d/20??????/*',/test_directory,count=nc4d_expdirs)
if nc4d_expdirs gt 0 then push,expdirs,c4d_expdirs
k4m_expdirs = file_search(dir+'k4m/20??????/*',/test_directory,count=nk4m_expdirs)
if nk4m_expdirs gt 0 then push,expdirs,k4m_expdirs
ksb_expdirs = file_search(dir+'ksb/20??????/*',/test_directory,count=nksb_expdirs)
if nksb_expdirs gt 0 then push,expdirs,ksb_expdirs
expdirs = trailingslash(repstr(expdirs,'/net/dl1/','/dl1/'))
nexpdirs = n_elements(expdirs)
print,strtrim(nexpdirs,2),' exposures directories'

;; Do groups of 100
ngroups = ceil(nexpdirs/100)
cmd = strarr(ngroups)
dirs = tmpdir+strarr(ngroups)
for i=0,ngroups-1 do begin
  lo = i*100
  hi = (lo+99) < (nexpdirs-1)
  cmd[i] = 'update_metafile_nmeas,["'+strjoin(expdirs[lo:hi],'","')+'"]'
endfor

stop

pbs_daemon,cmd,dirs,/idle,/hyperthread,prefix='metaupdate',nmulti=20,wait=1


stop

end
