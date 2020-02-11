pro remove_exposure_directories,expdir

n = n_elements(expdir)
print,strtrim(n,2),' exposure directories to delete'
for i=0,n-1 do begin
  files = file_search(expdir[i]+'/*.fits*',count=nfiles)
  print,strtrim(i+1,2),' ',expdir[i],' ',strtrim(nfiles,2),' fits files'
  if nfiles gt 0 then FILE_DELETE,files,/allow
  ;; leave the log and config files
endfor

;stop

end
