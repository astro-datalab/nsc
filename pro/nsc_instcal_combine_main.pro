;+
;
; NSC_INSTCAL_COMBINE_MAIN
;
; This combines the NSC exposures in one region of the sky.
;
; INPUTS:
;  version   The version name, e.g. 'v3'.
;  =hosts    Array of hosts to run.  The default is gp09, hulk and thing.
;  =redo     Rerun on exposures that were previously processed.
;  =nmulti   The number of simultaneously jobs to run. Default is 30.
;
; OUTPUTS:
;  A log "journal" file is put in ROOTDIR+users/dnidever/nsc/instcal/v#/logs/ 
;  as well as a structure with information on the jobs that were run.
;  The individual catalogs are put in ROOTDIR+users/dnidever/nsc/instcal/v#/combine/HEALPIX1000/.
;
; USAGE:
;  IDL>nsc_instcal_combine_main,'v3'
;
; By D. Nidever 2017
;-

pro nsc_instcal_combine_main,version,hosts=hosts,redo=redo,nmulti=nmulti

if n_elements(version) eq 0 then begin
  print,'Syntax - nsc_instcal_combine_main,version,hosts=hosts,redo=redo,nmulti=nmulti'
  return
endif

; Combine all of the data
NSC_ROOTDIRS,dldir,mssdir,localdir,longhost
host = first_el(strsplit(longhost,'.',/extract))
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
if file_test(dir+'combine/',/directory) eq 0 then file_mkdir,dir+'combine/'
if file_test(dir+'combine/logs/',/directory) eq 0 then file_mkdir,dir+'combine/logs/'
if file_test(localdir+'dnidever/nsc/instcal/'+version+'/') eq 0 then file_mkdir,localdir+'dnidever/nsc/instcal/'+version+'/'
plotsdir = dir+'plots/'
if file_test(plotsdir,/directory) eq 0 then file_mkdir,plotsdir
;; Hosts
if n_elements(hosts) eq 0 then hosts = ['gp09','hulk','thing']
if total(hosts eq host) eq 0 then begin
  print,'Current HOST='+host+' not in list of HOSTS = [ '+strjoin(hosts,', ')+' ] '
  return
endif
nside = 128
nmulti = 30
radeg = 180.0d0 / !dpi

; Log file
;------------------
; format is nsc_combine_main.DATETIME.log
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
logfile = dir+'combine/logs/nsc_instcal_combine_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

print, "Combining NOAO InstCal catalogs"

; CREATE LIST OF HEALPIX AND OVERLAPPING EXPOSURES
; Which healpix pixels have data
listfile = dir+'lists/nsc_instcal_combine_healpix_list.fits.gz'
if file_test(listfile) eq 0 then begin
  print,listfile,' NOT FOUND.  Run nsc_instcal_combine_qacuts.pro'
  return
endif

print,'Reading list from ',listfile
healstr = MRDFITS(listfile,1)
index = MRDFITS(listfile,2)
upix = index.pix
npix = n_elements(index)
; Copy to local directory for faster reading speed
file_copy,listfile,localdir+'dnidever/nsc/instcal/'+version+'/',/over


;; Create the commands
allcmd = "nsc_instcal_combine,"+strtrim(upix,2)+",nside="+strtrim(nside,2)+",version='"+version+"'"
if keyword_set(redo) then allcmd+=',/redo'
alldirs = strarr(npix)+localdir+'dnidever/nsc/instcal/'+version+'/tmp/'
nallcmd = n_elements(allcmd)

; RANDOMIZE
rnd = sort(randomu(0,npix))
allcmd = allcmd[rnd]
alldirs = alldirs[rnd]

;; Parcel out the jobs
nhosts = n_elements(hosts)
torun = lindgen(nallcmd)
nperhost = nallcmd/nhosts
for i=0,nhosts-1 do $
  if stregex(host,hosts[i],/boolean) eq 1 then torun=torun[i*nperhost:(i+1)*nperhost-1]
ntorun = n_elements(torun)
cmd = allcmd[torun]
dirs = alldirs[torun]
print,'Running ',strtrim(n_elements(torun),2),' on ',hostname

; Saving the structure of jobs to run
runfile = dir+'lists/nsc_instcal_combine_main.'+hostname+'.'+logtime+'_run.fits'
print,'Writing running information to ',runfile
runstr = replicate({pix:0L},n_elements(cmd))
runstr.pix = pix
MWRFITS,runstr,runfile,/create

stop

a = '' & read,a,prompt='Press RETURN to start'
PBS_DAEMON,cmd,dirs,jobs=jobs,/hyperthread,/idle,prefix='nsccmb',wait=wait,nmulti=nmulti


stop

; Run nsc_combine_summary.pro when it's done


; End logfile
;------------
JOURNAL


stop

end
