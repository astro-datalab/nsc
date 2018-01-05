pro get_goodpms_all_ucac5

; Roll-up the NSC-UCAC5 proper motion files

if n_elements(nside) eq 0 then nside = 128
NSC_ROOTDIRS,dldir,mssdir,localdir
if n_elements(version) eq 0 then version='v2'
dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
pmdir = dir+'combine/ucac5pm/'

files = file_search(pmdir+'*.fits',count=nfiles)
print,strtrim(nfiles,2),' PM files found'

; Loop over the fits files
for i=0,nfiles-1 do begin
  if i mod 100 eq 0 then print,i,' ',files[i]
  str1 = mrdfits(files[i],1,/silent,status=status)
  if status ne 0 then goto,BOMB
  nstr1 = n_elements(str1)

  ; Initialize the structure
  if i eq 0 then begin
    schema = str1[0]
    struct_assign,{dum:''},schema
    str = replicate(schema,1e7)
    nstr = n_elements(str)
    cnt = 0LL
  endif

  ; Add more elements
  if cnt+nstr1 gt nstr then begin
    old = str
    nnew = 1e6 > nstr1
    str = replicate(schema,nstr+nnew)
    str[0:nstr-1] = old
    nstr = n_elements(str)
    undefine,old
  endif

  ; Load in this pixel
  str[cnt:cnt+nstr1-1] = str1
  cnt += nstr1

  BOMB:

 ; stop
endfor
; Trim off extra elements
str = str[0:cnt-1]

; Make sure they are unique
str.id = strtrim(str.id,2)
ui = uniq(str.id,sort(str.id))
str = str[ui]
print,strtrim(n_elements(str),2),' unique objects'

; Save the file
MWRFITS,str,pmdir+'ucac5pm.fits',/create

stop

end
