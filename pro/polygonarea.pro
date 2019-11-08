function polygonarea,x,y

; (X[i], Y[i]) are coordinates of i'th point. 
;double polygonArea(double X[], double Y[], int n) 
;
;https://www.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/

; Initialze area 
area = 0.0
  
; Calculate value of shoelace formula 
n = n_elements(x)
j = n - 1; 
for i=0,n-1 do begin
  area += (x[j] + x[i]) * (y[j] - y[i])
  j = i                                 ; j is previous vertex to i 
endfor
area = abs(area/2.0)  

; Return absolute value 
return,area

end
