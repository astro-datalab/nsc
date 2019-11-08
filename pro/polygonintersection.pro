function computeintersection,p1,p2,edge

;m = (y2-y1)/(x2-x1)
;b = y1 - m*x1

;; p1 and p2 define the first line
m1 = (p1[1]-p2[1])/(p1[0]-p2[0])
b1 = p1[1]-m1*p1[0]

;; edge [x1,x2, y1,y2], defines the second line
m2 = (edge[3]-edge[2])/(edge[1]-edge[0])
b2 = edge[2]-m2*edge[0]

;; Set them equal to each other
;; m1*x+b1 = m2*x+b2
;; x*(m1-m2)=(b2-b1)
;; x=(b2-b1)/(m1-m2)
xout = (b2-b1)/(m1-m2)
yout = m1*xout+b1

return,[xout,yout]

end

;;


pro polygonintersection,cx,cy,sx,sy,xout,yout

;; cx/cy is the clip polygon, while sx/sy is the subject polygon

;; https://en.wikipedia.org/wiki/Sutherland–Hodgman_algorithm
;; https://www.geeksforgeeks.org/polygon-clipping-sutherland-hodgman-algorithm-please-change-bmp-images-jpeg-png/

nc = n_elements(cx)
ns = n_elements(sx)

;; Are the clipper polygon vertices ordered clockwise or
;; counter-clockwise?  use cross-product of the first two segments
;; make the second point the origin
v1 = [cx[0]-cx[1],cy[0]-cy[1],0]
v2 = [cx[2]-cx[1],cy[2]-cy[1],0]
cp = crossp(v1,v2)
;; negative is counterclockwise
if cp[2] le 0 then begin
  ;; reverse the order of the clipping polgon vertices
  ;;   will revert at the end
  cx = reverse(cx)
  cy = reverse(cy)
  counterclockwise=1
endif else counterclockwise=0

;; Initialize xout 
xout = sx
yout = sy
  
;; Loop over edges in clip polygon
;  for (Edge clipEdge in clipPolygon) do
for i=0L,nc-1 do begin
  ;; Clip edge
  if i lt (nc-1) then begin
    clipx = cx[i:i+1]
    clipy = cy[i:i+1]
  endif else begin
    clipx = [cx[i],cx[0]]
    clipy = [cy[i],cy[0]]
  endelse

  ;; The current, temporary vertices
  x = xout
  y = yout
  ;; Initialize xout/yout
  xout = dblarr(n_elements(x)*2)
  yout = dblarr(n_elements(x)*2)
  nout = 0L
   
  ;; Loop over vertices in subject polygon
  for j=0L,n_elements(x)-1 do begin
    x1 = x[j]
    y1 = y[j]

    if j eq 0 then pind=n_elements(x)-1 else pind=j-1
    ;Point prev_point = inputList[(i + inputList.count - 1) % inputList.count];
    px = x[pind]
    py = y[pind]   

    ;Point Intersecting_point = ComputeIntersection(prev_point, current_point, clipEdge)
    intersect = computeintersection([px,py],[x1,y1],[clipx,clipy])
   
    ;; Calculate if the point is inside or outside the clip edge
    ;; requires the vertices of the clipper polygon to be in clockwise order
    ;p = (x2-x1)*(py-y1)-(y2-y1)*(px-x1)
    pos = (clipx[1]-clipx[0])*(y1-clipy[0])-(clipy[1]-clipy[0])*(x1-clipx[0])   ;; current point
    ppos = (clipx[1]-clipx[0])*(py-clipy[0])-(clipy[1]-clipy[0])*(px-clipx[0])  ;; previous point
    ;; p<0, the point is on the right side of the line, INSIDE
    ;; p=0, the point is on the edge
    ;; p>0, the point is on the ledge side of the line, OUTSIDE

    ;; If point inside the clip edge
    if pos le 0 then begin
      ;; If previous point not inside the clip edge
      if ppos gt 0 then begin
         ;outputList.add(Intersecting_point)
         xout[nout] = intersect[0]
         yout[nout] = intersect[1]
         nout++
      endif
      ;outputList.add(current_point)
      xout[nout] = x1
      yout[nout] = y1
      nout++

    ;; point outside the clip edge
    endif else begin
      ;; if previous point inside clip edge
      if ppos le 0 then begin
        ;outputList.add(Intersecting_point)
        xout[nout] = intersect[0]
        yout[nout] = intersect[1]
        nout++
      endif
    endelse
  endfor
  ;; Trim excess elements
  if nout lt n_elements(xout) then begin
    xout = xout[0:nout-1]
    yout = yout[0:nout-1]
  endif

endfor

;; Reverse clipping polygon vertices if counterclockwise
if counterclockwise eq 1 then begin
  cx = reverse(cx)
  cy = reverse(cy)
endif

end
