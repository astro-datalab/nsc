;+
;
; NSC_INSTCAL_MEASURE_UPDATE_MAIN
;
; Driver for nsc_instcal_measure_update.pro to update the objectIDs
; in the measurement catalogs.
;
; INPUTS:
;  version   The version name, e.g. 'v3'.
;  =hosts    Array of hosts to run.  The default is gp06-gp09,hulk and thing.
;  =redo     Rerun on exposures that were previously processed.
;  =nmulti   The number of simultaneously jobs to run. Default is 30.
;
; OUTPUTS:
;  A log "journal" file is put in ROOTDIR+users/dnidever/nsc/instcal/v#/logs/
;  as well as a structure with information on the jobs that were run.
;  The individual catalogs are put in ROOTDIR+users/dnidever/nsc/instcal/v#/NIGHT/EXPOSURENAME/.
;
; USAGE:
;  IDL>nsc_instcal_measure_update_main,'v3'
;
; By D.Nidever  Aug 2019
;-

pro nsc_instcal_measure_update_main,version,hosts=hosts,redo=redo,nmulti=nmulti

if n_elements(version) eq 0 then begin
  print,'Syntax - nsc_instcal_measure_update_main,version,hosts=hosts,redo=redo,nmulti=nmulti'
  return
endif

; Main NOAO DECam source catalog
NSC_ROOTDIRS,dldir,mssdir,localdir,host
hostname = first_el(strsplit(host,'.',/extract))
if n_elements(maxjobs) eq 0 then maxjobs=48300L
if n_elements(nmulti) eq 0 then nmulti=30
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
tmpdir = localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then file_mkdir,tmpdir
subdirs = ['logs','c4d','k4m','ksb']
for i=0,n_elements(subdirs)-1 do if file_test(dir+subdirs[i],/directory) eq 0 then file_mkdir,dir+subdirs[i]
;; Hosts
if n_elements(hosts) eq 0 then hosts = ['gp06','gp07','gp08','gp09','hulk','thing']
if total(hosts eq hostname) eq 0 then begin
  print,'Current HOST='+hostname+' not in list of HOSTS = [ '+strjoin(hosts,', ')+' ] '
  return
endif

t0 = systime(1)

; Log file
;------------------
; format is nsc_main_laf.DATETIME.log
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
logfile = dir+'logs/nsc_instcal_measure_update_main.'+logtime+'.log'
JOURNAL,logfile

print, "Updating Object IDs in measurement catalogs"

; Loading the exposure list that combine was run ont
list = MRDFITS(dir+'/lists/nsc_instcal_combine_healpix_list.fits.gz',1)
list.file = strtrim(list.file,2)
list.base = strtrim(list.base,2)
expdir = file_dirname(list.file)
ui = uniq(expdir,sort(expdir))  ;; get unique exposures
expdir = expdir[ui]
nexpdir = n_elements(expdir)

;; Putting them in RANDOM but REPEATABLE order
seed = 1
print,'RANDOMIZING WITH SEED=1'
si = sort(randomu(seed,nexpdir))
expdir = expdir[si]

;; Make the jobs
allcmd = 'nsc_instcal_measure_update,"'+expdir+'"'
alldir = strarr(nexpdir)+tmpdir

;; Parcel out the jobs
nhosts = n_elements(hosts)
nperhost = nexpdir/nhosts
torun = lindgen(nexpdir)
for i=0,nhosts-1 do $
  if stregex(host,hosts[i],/boolean) eq 1 then torun=torun[i*nperhost:(i+1)*nperhost-1]
ntorun = n_elements(torun)

if ntorun eq 0 then begin
  print,'No exposures to process.'
  return
endif

cmd = allcmd[torun]
cmddir = alldir[torun]
expdir = expdir[torun]

; Saving the structure of jobs to run
runfile = dir+'lists/nsc_instcal_measure_update_main.'+hostname+'.'+logtime+'_run.fits'
print,'Writing running information to ',runfile
expstr = replicate({expdir:'',cmd:''},ntorun)
expstr.expdir = expdir
expstr.cmd = cmd
MWRFITS,expstr,runfile,/create

; Run PBS_DAEMON
stop
a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,cmddir,jobs=jobs,/idle,/hyperthread,prefix='nscmeasup',wait=5,nmulti=nmulti


print,'dt=',stringize(systime(1)-t0,ndec=2),' sec'

; End logfile
;------------
JOURNAL

stop

end
