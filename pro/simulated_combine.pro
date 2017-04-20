pro simulated_combine

; Simulate how long it will take to run the combine stage
; using various ordering strategies

dir = '/dl1/users/dnidever/nsc/instcal/'

; Load the table healpix
str = mrdfits(dir+'combine/nsc_instcal_combine.fits',1)
; only keep successful ones
gd = where(str.success eq 1,ngd)
str = str[gd]

; Load the times
sum1 = mrdfits(dir+'nsccmb_summary_hulk.fits',1)
sum2 = mrdfits(dir+'nsccmb_summary_thing.fits',1)
sum3 = mrdfits(dir+'nsccmb_summary_gp09.fits',1)
sum = [sum1,sum2,sum3]
si = sort(sum.mtime)
sum = sum[si]
; only keep fairly recent ones
gd = where(sum.mtime gt 1.4897704e+09,ngd)
sum = sum[gd]

; Deal with duplicates
dbl = doubles(sum.pix,count=ndbl)
alldbl = doubles(sum.pix,/all,count=nalldbl)
torem = bytarr(nalldbl)
for i=0,ndbl-1 do begin
  MATCH,sum[alldbl].pix,sum[dbl[i]].pix,ind1,ind2,/sort,count=nmatch
  torem[ind1[0:nmatch-2]] = 1
endfor
bd=where(torem eq 1,nbd)
remove,alldbl[bd],sum

; Match them up and get DT
MATCH,str.pix,sum.pix,ind1,ind2,/sort,count=nmatch
add_tag,str,'dt',-1.0,str
str[ind1].dt = sum[ind2].dt
; only keep ones with DT values
gd = where(str.dt gt 0,ngd)
str = str[gd]
nstr = n_elements(str)

; Do the simulation
;-------------------

; 1) Current order for now
; Total time = 8063.44 hours
; Total time per cpu = 89.5937 hours
; Run statistics:
; Total time = 8063.44 hours
; Average time per cpu = 89.5938 hours
; Maximum time per cpu = 173.295 hours
; Some jobs take ~100 hours to run.

; 2) Inverse sort by DT, biggest jobs FIRST
;si = reverse(sort(str.dt))
;str = str[si]
;Total time = 8063.43 hours
;Total time per cpu = 89.5937 hours
;Run statistics:
;Total time = 8063.44 hours
;Average time per cpu = 89.5938 hours
;Maximum time per cpu = 106.881 hours

; 3) Random
;si = sort(randomu(seed,nstr))
;str = str[si]
;Total time = 8063.44 hours
;Total time per cpu = 89.5938 hours
;Run statistics:
;Total time = 8063.44 hours
;Average time per cpu = 89.5938 hours
;Maximum time per cpu = 153.046 hours

tottime = total(str.dt)/3600.
print,'Total time = ',strtrim(tottime,2),' hours'

; Initialize the time array
ncpu = 90
print,'Total time per cpu = ',strtrim(tottime/ncpu,2),' hours'
timearr = fltarr(ncpu)

; Loop through all healpix
runstr = replicate({pix:0L,dt:0.0,time:0.0,cpu:0},nstr)
runstr.pix = str.pix
runstr.dt = str.dt
For i=0,nstr-1 do begin
  ;if (i+1) mod 10000 eq 0 then print,i
  ; At each step find out which cpu has the lowest overall
  ; cumulative time and add the next job/pix to it
  ; If multiple have the same, take the first one
  ind = first_el(minloc(timearr))
  runstr[i].time = timearr[ind]
  runstr[i].cpu = ind
  timearr[ind] += runstr[i].dt
Endfor
print,'Run statistics:'
print,'Total time = ',strtrim(total(timearr)/3600,2),' hours'
print,'Average time per cpu = ',strtrim(mean(timearr)/3600,2),' hours'
print,'Maximum time per cpu = ',strtrim(max(timearr)/3600,2),' hours'


stop

end
