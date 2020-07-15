pro nsc_rootdirs,dldir,mssdir,localdir,host

; Return the root directories for various machines

undefine,dldir,mssdir,localdir,host

spawn,'hostname',out,errout,/noshell
host = (strsplit(out[0],' ',/extract))[0]

if stregex(host,'thing',/boolean) eq 1 or stregex(host,'hulk',/boolean) eq 1 then begin
  ;dldir = '/dl1/users/'
  dldir = '/net/dl2/'
  mssdir = '/mss1/'
  localdir = '/d0/'
  return
endif
if stregex(host,'gp09',/boolean) eq 1 or stregex(host,'gp08',/boolean) eq 1 or stregex(host,'gp07',/boolean) eq 1 or $
   stregex(host,'gp06',/boolean) eq 1 or stregex(host,'gp05',/boolean) eq 1 then begin
  ;dldir = '/net/dl1/users/'
  dldir = '/net/dl2/'
  mssdir = '/net/mss1/'
  localdir = '/data0/'
  return
endif

; Not sure about this host
print,host,' UNKNOWN'

end
