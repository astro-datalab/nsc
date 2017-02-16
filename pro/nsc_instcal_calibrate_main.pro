pro nsc_instcal_calibrate_main,redo=redo

; Main NOAO DECam source catalog
dir = "/datalab/users/dnidever/decamcatalog/"

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
logfile = dir+'nsc_instcal_calibrate_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

print, "Calibrating DECam InstCal SExtractor catalogs"

; Find all of the directories
expdirs = file_search(dir+'instcal/20??????/*',/test_directory,count=nexpdirs)

cmd = 'nsc_instcal_calibrate,"'+expdirs+'"'
dir = strarr(nexpdirs)+'/data0/dnidever/decamcatalog/instcal/tmp/'
stop

; Run PBS_DAEMON
; this works
PBS_DAEMON,cmd,dir,/hyperthread,/idle,prefix='nsccalib',wait=10,nmulti=10
;PBS_DAEMON,cmd,dir,/hyperthread,prefix='nds',wait=10,nmulti=20

; End logfile
;------------
JOURNAL

;stop

end
