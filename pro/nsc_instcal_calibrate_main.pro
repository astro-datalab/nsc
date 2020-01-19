;+
;
; NSC_INSTCAL_CALIBRATE_MAIN
;
; This calibrates the NSC exposures.
; This is a wrapper around nsc_instcal_calibrate.pro which runs
; on one individual exposure.
;
; INPUTS:
;  version   The version name, e.g. 'v3'. 
;  =hosts    Array of hosts to run.  The default is gp09,hulk and thing.
;  =redo     Rerun on exposures that were previously processed.
;  =nmulti   The number of simultaneously jobs to run. Default is 15.
;
; OUTPUTS:
;  A log "journal" file is put in ROOTDIR+users/dnidever/nsc/instcal/v#/logs/
;  as well as a structure with information on the jobs that were run.
;  The individual catalogs are put in ROOTDIR+users/dnidever/nsc/instcal/v#/NIGHT/EXPOSURENAME/.
;
; USAGE:
;  IDL>nsc_instcal_calibrate_main,'v3'
;
; By D.Nidever  Feb 2017
;-

pro nsc_instcal_calibrate_main,version,hosts=hosts,nmulti=nmulti,redo=redo

if n_elements(version) eq 0 then begin
  print,'Syntax - nsc_instcal_calibrate_main,version,hosts=hosts,nmulti=nmulti,redo=redo'
  return
endif

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir,host
hostname = first_el(strsplit(host,'.',/extract))
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(dir,/directory) eq 0 then file_mkdir,dir+'logs/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir
if file_test(dir+'logs/',/directory) eq 0 then file_mkdir,dir+'logs/'
;; Hosts
if n_elements(hosts) eq 0 then hosts = ['gp09','hulk','thing']
if total(hosts eq hostname) eq 0 then begin
  print,'Current HOST='+host+' not in list of HOSTS = [ '+strjoin(hosts,', ')+' ] '
  return
endif
if n_elements(nmulti) eq 0 then nmulti = 15
wait = 1

; Log file
;------------------
; format is nds_main_laf.DATETIME.log
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logtime = smonth+sday+syear+shour+sminute+ssecond
logfile = dir+'logs/nsc_instcal_calibrate_main.'+logtime+'.log'
JOURNAL,logfile


print, "Calibrating DECam/Mosiac3/Bok InstCal SExtractor catalogs"

; Find all of the directories
print,'Getting the exposure directories'
c4d_expdirs = file_search(dir+'c4d/20??????/*',/test_directory,count=nc4d_expdirs)
if nc4d_expdirs gt 0 then push,expdirs,c4d_expdirs
k4m_expdirs = file_search(dir+'k4m/20??????/*',/test_directory,count=nk4m_expdirs)
if nk4m_expdirs gt 0 then push,expdirs,k4m_expdirs
ksb_expdirs = file_search(dir+'ksb/20??????/*',/test_directory,count=nksb_expdirs)
if nksb_expdirs gt 0 then push,expdirs,ksb_expdirs
expdirs = trailingslash(repstr(expdirs,'/net/dl1/','/dl1/'))
; Match these to our lists
list1 = MRDFITS(dir+'lists/decam_instcal_list.fits.gz',1)
list2 = MRDFITS(dir+'lists/mosaic3_instcal_list.fits.gz',1)
list3 = MRDFITS(dir+'lists/bok90prime_instcal_list.fits.gz',1)
list = [list1,list2,list3]
lbase = file_basename(strtrim(list.fluxfile,2),'.fits.fz')
MATCH,file_basename(expdirs),lbase,ind1,ind2,/sort,count=nmatch
expdirs = expdirs[ind1]
; GDL file_search doesn't have /test_directory
;expdirs = file_search(dir+'instcal/20??????/*',count=nexpdirs)
;gddirs = where(file_test(expdirs,/directory) eq 1,ngddirs)
;expdirs = expdirs[gddirs]
nexpdirs = n_elements(expdirs)
print,strtrim(nexpdirs,2),' exposure directories'

allcmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
if keyword_set(redo) then allcmd+=',/redo'
alldirs = strarr(nexpdirs)+tmpdir
;dirs = strarr(nexpdirs)+'/data0/dnidever/decamcatalog/instcal/tmp/'
nallcmd = n_elements(allcmd)

; RANDOMIZE
rnd = sort(randomu(0,nexpdirs))
allcmd = allcmd[rnd]
alldirs = alldirs[rnd]
expdirs = expdirs[rnd]

;; Only run failed exposures
;sumstr = mrdfits(dir+'lists/nsc_calibrate_summary.fits.gz',1)
;sumstr.expdir = strtrim(sumstr.expdir,2)
;sumstr.base = strtrim(sumstr.base,2)
;sumstr.filter = strtrim(sumstr.filter,2)
;sumstr.instrument = strtrim(sumstr.instrument,2)
;bd = where(sumstr.success eq 0 and sumstr.fwhm le 2 and sumstr.exptime ge 30 and sumstr.nsources ge 1000 and $
;           sumstr.filter ne 'u' and sumstr.filter ne 'N9' and sumstr.filter ne 'N6',nbd)
;failed_expdirs = sumstr[bd].expdir
;;sumstr = mrdfits(dir+'lists/nsc_instcal_calibrate_failures.fits',1)
;;bd = where(sumstr.nsources gt 100 and sumstr.fwhm le 2 and sumstr.exptime ge 30 and sumstr.meta_exists eq 0,nbd)
;;failed_expdirs = strtrim(sumstr[bd].expdir,2)
;;failed_expdirs = trailingslash(repstr(failed_expdirs,'/net/dl1/','/dl1/'))
;MATCH,expdirs,failed_expdirs,ind1,ind2,/sort,count=nmatch_failed
;print,strtrim(nmatch_failed,2),' failed exposures to run'
;allcmd = allcmd[ind1]
;alldirs = alldirs[ind1]
;expdirs = expdirs[ind1]
;nallcmd = n_elements(allcmd)

;; Only run on the ~180,000 "new" exposures and the ~15,000 exposures
;; with extension name problems
sub1 = mrdfits(dir+'lists/nsc_instcal_measure_main.gp06.datalab.noao.edu.122719165718_run.fits',1)
sub2 = mrdfits(dir+'lists/nsc_instcal_measure_main.gp07.datalab.noao.edu.122719165757_run.fits',1)
sub3 = mrdfits(dir+'lists/nsc_instcal_measure_main.gp08.datalab.noao.edu.122719165811_run.fits',1)
g1 = where(sub1.submitted eq 'T',ng1)
g2 = where(sub2.submitted eq 'T',ng2)
g3 = where(sub3.submitted eq 'T',ng3)
sub = [sub1[g1],sub2[g2],sub3[g3]]
extstr = mrdfits(dir+'lists/badexposures_extnameproblem2.fits',1)
fluxfiletorun = [sub.fluxfile,extstr.fluxfile]
fluxfiletorun = strtrim(fluxfiletorun,2)
MATCH,file_basename(expdirs),file_basename(fluxfiletorun,'.fits.fz'),ind1,ind2,/sort,count=nmatch
print,strtrim(nmatch,2),' exposures to rerun'
allcmd = allcmd[ind1]
alldirs = alldirs[ind1]
expdirs = expdirs[ind1]
nallcmd = n_elements(allcmd)


;; Parcel out the jobs
nhosts = n_elements(hosts)
torun = lindgen(nallcmd)
nperhost = nallcmd/nhosts
for i=0,nhosts-1 do $
  if stregex(host,hosts[i],/boolean) eq 1 then torun=torun[i*nperhost:(i+1)*nperhost-1]
ntorun = n_elements(torun)
cmd = allcmd[torun]
dirs = alldirs[torun]
expdirs = expdirs[torun]
print,'Running ',strtrim(n_elements(torun),2),' on ',hostname

; Saving the structure of jobs to run
runfile = dir+'lists/nsc_instcal_calibrate_main.'+hostname+'.'+logtime+'_run.fits'
print,'Writing running information to ',runfile
runstr = replicate({expdir:'',cmd:'',dir:''},n_elements(cmd))
runstr.expdir = expdirs
runstr.cmd = cmd
runstr.dir = dirs
MWRFITS,runstr,runfile,/create

stop

a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccalib',wait=wait,nmulti=nmulti


; End logfile
;------------
JOURNAL

stop

end
