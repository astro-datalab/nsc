function ellipsecoords, rmax, rmin, xc, yc, pos_ang, $
                        pix2world=pix2world,world2pix=world2pix, $
                        head=head, NPOINTS = npoints
        
;+
; NAME:
;      ELLIPSECOORDS (MODIFIED TVELLIPSE)
;
; PURPOSE:
;      Draw an ellipse on the current graphics device and get the coordinates.
;
; CALLING SEQUENCE:
;      coords = ELLIPSECOORDS(rmax, rmin, xc, yc, pos_ang, NPOINTS=)
;
; INPUTS:
;       RMAX,RMIN - Scalars giving the major and minor axis of the ellipse
; OPTIONAL INPUTS:
;       XC,YC - Scalars giving the position on the TV of the ellipse center
;               If not supplied (or if XC, YC are negative and /DATA is not set), 
;               and an interactive graphics device (e.g. not postscript) is set,
;               then the user will be prompted for X,Y
;       POS_ANG - Position angle of the major axis, measured counter-clockwise
;                 from the X axis.  In degrees. Default is 0.  If
;                 /world2pix is set, then the position angle is
;                 assumed to be CW from North.
;       COLOR - Scalar  giving intensity level to draw ellipse.   The color
;               can be specified either with either this parameter or with the 
;               COLOR keyword.   Default is !P.COLOR
;
; OPTIONAL KEYWORD INPUT:
;        COLOR - Intensity value used to draw the circle, overrides parameter
;               value.  Default = !P.COLOR
;        /DATA - if this keyword is set and non-zero, then the ellipse radii and
;               X,Y position center are interpreted as being in DATA 
;               coordinates.   Note that the data coordinates must have been 
;               previously defined (with a PLOT or CONTOUR call).
;        THICK - Thickness of the drawn ellipse, default = !P.THICK
;        LINESTLYLE - Linestyle used to draw ellipse, default = !P.LINESTYLE
;        NPOINTS - Number of points to connect to draw ellipse, default = 120
;                  Increase this value to improve smoothness
; RESTRICTIONS:
;        TVELLIPSE does not check whether the ellipse is within the boundaries
;        of the window.
;
;        The ellipse is evaluated at NPOINTS (default = 120) points and 
;        connected by straight lines, rather than using the more sophisticated 
;        algorithm used by TVCIRCLE
;
;        TVELLIPSE does not accept normalized coordinates.
;
;        TVELLIPSE is not vectorized; it only draws one ellipse at a time
; EXAMPLE:
;        Draw an ellipse of major axis 50 pixels, minor axis 30 pixels, centered
;        on (250,100), with the major axis inclined 25 degrees counter-clockwise
;        from the X axis.   Use a double thickness line and device coordinates 
;        (default)
;
;	IDL> tvellipse,50,30,250,100,25,thick=2
; NOTES:
;        Note that the position angle for TVELLIPSE (counter-clockwise from the
;        X axis) differs from the astronomical position angle (counter-clockwise
;        from the Y axis). 
;
; REVISION HISTORY:
;        Written  W. Landsman STX          July, 1989            
;        Converted to use with a workstation.  M. Greason, STX, June 1990
;        LINESTYLE keyword, evaluate at 120 points,  W. Landsman HSTX Nov 1995
;        Added NPOINTS keyword, fixed /DATA keyword W. Landsman HSTX Jan 1996
;        Check for reversed /DATA coordinates  P. Mangiafico, W.Landsman May 1996
;        Converted to IDL V5.0   W. Landsman   September 1997
;        Work correctly when X & Y data scales are unequal  December 1998
;        Removed cursor input when -ve coords are entered with /data 
;        keyword set  P. Maxted, Keele, 2002
;-
 On_error,2                              ;Return to caller

 if N_params() lt 4 then begin
   print,'Syntax - coords= ELLIPSECOORDS(rmax, rmin, xc, yc pos_ang,'
   print,'                               NPOINTS =)'
   return,-1
 endif

 if N_params() LT 5 then pos_ang = 0.    ;Default position angle

 if not keyword_set(NPOINTS) then npoints = 120   ;Number of points to connect
 phi = 2*!pi*(findgen(npoints)/(npoints-1))       ;Divide circle into Npoints
 ang = pos_ang/!radeg                             ;Position angle in radians
 ; If world coordinates, assume E of N (CW)
 ;  convert to CCW from X/RA axis
 if keyword_set(world2pix) then ang=(90-pos_ang)/!radeg
 cosang = cos(ang)
 sinang = sin(ang)

 x =  rmax*cos(phi)              ;Parameterized equation of ellipse
 y =  rmin*sin(phi)

 xprime = xc + x*cosang - y*sinang   	;Rotate to desired position angle
 if keyword_set(world2pix) then xprime = xc + (x*cosang - y*sinang)/cos(yc/!radeg)  ; correct for cos(dec)
 yprime = yc + x*sinang + y*cosang

 ; Convert from RA/DEC to X/Y or vice versa
 if (keyword_set(pix2world) or keyword_set(world2pix)) and n_elements(head) gt 0 then begin

   ; Pixel to world coordinates (X/Y -> RA/DEC)
   if keyword_set(pix2world) then begin
     x = xprime
     y = yprime
     undefine,xprime,yprime
     head_xyad,head,x,y,xprime,yprime,/deg
      
   ; World to pixel coordinates (RA/DEC -> X/Y)
   endif else begin
     ra = xprime
     dec = yprime
     undefine,xprime,yprime
     head_adxy,head,ra,dec,xprime,yprime,/deg
   endelse
    
 endif
    
 coords = fltarr(2,npoints)
 coords(0,*) = xprime
 coords(1,*) = yprime

 return,coords
 
 end
