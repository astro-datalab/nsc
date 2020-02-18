pro remove_old_exposure_directories.pro

;; Remove old reduction exposure directories

;; DELETING FILES FOR OLD REDUCTION VERSIONS that I'm not using

print,'Getting the exposure directories'
dir = '/net/dl1/users/dnidever/nsc/instcal/v3/'
c4d_expdirs = file_search(dir+'c4d/20??????/*',/test_directory,count=nc4d_expdirs)
if nc4d_expdirs gt 0 then push,expdirs,c4d_expdirs
k4m_expdirs = file_search(dir+'k4m/20??????/*',/test_directory,count=nk4m_expdirs)
if nk4m_expdirs gt 0 then push,expdirs,k4m_expdirs
ksb_expdirs = file_search(dir+'ksb/20??????/*',/test_directory,count=nksb_expdirs)
if nksb_expdirs gt 0 then push,expdirs,ksb_expdirs
; 669143 exposure directories

sum = mrdfits(dir+'lists/nsc_calibrate_summary.fits.gz',1)
sum.expdir = strtrim(sum.expdir,2)
; 490623 exposure directories

MATCH,expdirs+'/','/net'+sum.expdir,ind1,ind2,/sort
; 490611 matches, 12 don't match don't have directories

torem = expdirs
remove,ind1,torem
; 178532 to remove, 131414 have v1 versions

;; group by 100
ntorem = n_elements(torem)
bins = 100L
ngroups = ceil(ntorem/float(bins))
cmd = strarr(ngroups)
for i=0,ngroups-1 do cmd[i]="remove_exposure_directories,['"+strjoin(torem[i*bins:((i+1)*bins-1)<(ntorem-1)],"','")+"']"
cmddir = '/data0/dnidever/nsc/instcal/v3/tmp/'+strarr(n_elements(cmd))
pbs_daemon,cmd,cmddir,jobs=jobs,/idle,/hyperthread,prefix='remdir',wait=1,nmulti=10
;; 1786 jobs, on hulk, started at 4:15pm

stop

end
