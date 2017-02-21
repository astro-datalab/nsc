pro nsc_rootdirs,dldir,mssdir,localdir

; Return the root directories for various machines

undefine,dldir,mssdir,localdir

spawn,'hostname',out,errout,/noshell
host = (strsplit(out[0],' ',/extract))[0]

if stregex(host,'thing',/boolean) eq 1 or stregex(host,'hulk',/boolean) eq 1 then begin
  dldir = '/d0/'
  mssdir = '/mss1/'
  localdir = '/dl1/'
  return
endif
if stregex(host,'gp09',/boolean) eq 1 then begin
  dldir = '/net/dl1/'
  mssdir = '/net/mss1/'
  localdir = '/data0/'
  return
endif

; Not sure about this host
print,host,' UNKNOWN'

end
