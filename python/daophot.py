#!/usr/bin/env python

import os
import sys
import numpy as np
import scipy
from astropy.io import fits
import photutils
from skimage import measure, morphology
from scipy.cluster import vq
#import gaps
import matplotlib.pyplot as plt
import pylab
from scipy.signal import argrelmin
import scipy.ndimage.filters as filters
import time
from dlnpyutils import utils as dln

def profile(ipstyp,dx,dy,par,dhdxc,dhdyc,term,ideriv):
    """ PSF analytical profile."""

    #C#######################################################################
    #C
    #      REAL FUNCTION  PROFIL  (IPSTYP, DX, DY, PAR, DHDXC, DHDYC, 
    #     .     TERM, IDERIV)
    #C
    #C Compute the value of an ANALYTIC prfile for a point DX,DY distant
    #C from the centroid.  Return both the computed value and its
    #C first derivatives with respect to x and y.  If IDERIV .NE. 0,
    #C return also the first derivatives with respect to all the parameters
    #C defining the profile.
    #C
    #      IMPLICIT NONE
    #      INTEGER MAXPAR, MAXPT
    #      PARAMETER (MAXPAR=6, MAXPT=4)
    #C
    #      REAL PAR(MAXPAR), TERM(MAXPAR)
    #      REAL D(MAXPT,MAXPT), W(MAXPT,MAXPT)
    #      REAL X(MAXPT), XSQ(MAXPT), P1XSQ(MAXPT)
    #C
    #      REAL EXP, DAOERF
    #C
    #      REAL DX, DY, DHDXC, DHDYC, WFSQ, Y, WT, WF, ONEMP3
    #      REAL RSQ, E, TALPHA, P1SQ, P2SQ, XY, DENOM
    #      REAL FUNC, YSQ, WP4FOD, P4FOD, F, P1P2, ERFX, DHDSX, ERFY
    #      REAL DEBY, DFBY, DBYX0, DBYY0
    #      REAL DHDSY, ALPHA, P2YSQ
    #      INTEGER I, IPSTYP, IDERIV, IX, IY, NPT
    #C
    #      DATA D / 0.00000000,  0.0,        0.0       , 0.0       ,
    #     .        -0.28867513,  0.28867513, 0.0       , 0.0       ,
    #     .        -0.38729833,  0.00000000, 0.38729833, 0.0       ,
    #     .        -0.43056816, -0.16999052, 0.16999052, 0.43056816/
    #      DATA W / 1.00000000,  0.0       , 0.0       , 0.0       ,
    #     .         0.50000000,  0.50000000, 0.0       , 0.0       ,
    #     .         0.27777778,  0.44444444, 0.27777778, 0.0       ,
    #     .         0.17392742,  0.32607258, 0.32607258, 0.17392742/
    #C
    #      PROFIL = 0.
    #      DHDXC = 0.
    #      DHDYC = 0.
    #C
    #      IF (IDERIV .GT. 0) THEN
    #         DO I=1,MAXPAR
    #            TERM(I) = 0.
    #         END DO
    #      END IF
    #C
    #      IF (IPSTYP .EQ. 1) THEN
    #C
    #C GAUSSIAN
    #C
    #C     F = ERFX * ERFY / (PAR(1) * PAR(2))
    #C
    #C PAR(1) is the HWHM in X; sigma(x) = 0.8493218 * HWHM
    #C PAR(2) is the HWHM in Y; ditto
    #C
    #         P1P2 = PAR(1)*PAR(2)
    #         ERFX = DAOERF(DX, 0., PAR(1), DHDXC, DHDSX)
    #         ERFY = DAOERF(DY, 0., PAR(2), DHDYC, DHDSY)
    #         PROFIL = ERFX*ERFY/P1P2
    #         DHDXC = DHDXC*ERFY/P1P2
    #         DHDYC = DHDYC*ERFX/P1P2
    #         IF (IDERIV .GT. 0) THEN
    #            TERM(1) = (DHDSX-ERFX/PAR(1))*ERFY/P1P2
    #            TERM(2) = (DHDSY-ERFY/PAR(2))*ERFX/P1P2
    #         END IF
    #C
    #      ELSE IF (IPSTYP .EQ. 2) THEN
    #C
    #C MOFFAT function   BETA = 1.5
    #C
    #C                            BETA-1
    #C F = --------------------------------------------------------
    #C      Ax * Ay * [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]**BETA
    #C
    #C PAR(1) is the HWHM in x at y = 0: 
    #C
    #C             1/2 = 1/[1 + (PAR(1)/Ax)**2]**BETA
    #C so
    #C             2**(1/BETA) - 1 = (PAR(1)/Ax)**2
    #C
    #C             Ax**2 = PAR(1)**2/[2**(1/BETA) - 1]
    #C
    #C When BETA = 1.5, Ax**2 = 1.7024144 * PAR(1)**2
    #C
    #C Hence, let us use
    #C
    #C                                  1
    #C F = ---------------------------------------------------------------
    #C     P(1)*P(2)*{1+0.5874011*[(X/P(1))**2+(Y/P(2))**2+(XY*P(3))]**1.5
    #C 
    #C neglecting a constant of proportionality.
    #C
    #         ALPHA = 0.5874011
    #         TALPHA = 1.1748021
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         P1P2 = PAR(1)*PAR(2)
    #         XY = DX*DY
    #C
    #         DENOM = 1. + ALPHA*(DX**2/P1SQ + DY**2/P2SQ + XY*PAR(3))
    #         IF (DENOM .GT. 5.E6) RETURN
    #         FUNC = 1. / (P1P2 * DENOM**PAR(4))
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = (PAR(4) - 1.) * FUNC
    #            P4FOD = PAR(4)*ALPHA*PROFIL/DENOM
    #            DHDXC = P4FOD*(2.*DX/P1SQ + DY*PAR(3))
    #            DHDYC = P4FOD*(2.*DY/P2SQ + DX*PAR(3))
    #            IF (IDERIV .GT. 0) THEN
    #               TERM(1) = (2.*P4FOD*DX**2/P1SQ-PROFIL)/PAR(1)
    #               TERM(2) = (2.*P4FOD*DY**2/P2SQ-PROFIL)/PAR(2)
    #               TERM(3) = - P4FOD*XY
    #C              TERM(4) = PROFIL*(1./(PAR(4)-1.)-ALOG(DENOM))
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            XSQ(IX) = X(IX)**2
    #            P1XSQ(IX) = XSQ(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            YSQ = Y**2
    #            P2YSQ = YSQ/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY=X(IX)*Y
    #               DENOM = 1. + ALPHA*(P1XSQ(IX) + P2YSQ + XY*PAR(3))
    #               FUNC = (PAR(4) - 1.) / (P1P2 * DENOM**PAR(4))
    #               P4FOD = PAR(4)*ALPHA*FUNC/DENOM
    #               WP4FOD = WT*P4FOD
    #               WF = WT*FUNC
    #               PROFIL = PROFIL + WF
    #               DHDXC = DHDXC + WP4FOD*(2.*X(IX)/P1SQ + Y*PAR(3))
    #               DHDYC = DHDYC + WP4FOD*(2.*Y/P2SQ + X(IX)*PAR(3))
    #               IF (IDERIV .GT. 0) THEN
    #                  TERM(1) = TERM(1) + 
    #     .                     (2.*WP4FOD*P1XSQ(IX)-WF)/PAR(1)
    #                  TERM(2) = TERM(2) + 
    #     .                     (2.*WP4FOD*P2YSQ-WF)/PAR(2)
    #                  TERM(3) = TERM(3) - WP4FOD*XY
    #C                 TERM(4) = TERM(4) + WF*(1./(PAR(4)-1.)-ALOG(DENOM))
    #               END IF
    #            END DO
    #         END DO
    #C
    #      ELSE IF (IPSTYP .EQ. 3) THEN
    #C
    #C MOFFAT function  BETA = 2.5
    #C
    #C                            BETA-1
    #C F = --------------------------------------------------------
    #C      Ax * Ay * [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]**BETA
    #C
    #C PAR(1) is the HWHM in x at y = 0: 
    #C
    #C             1/2 = 1/[1 + (PAR(1)/Ax)**2]**BETA
    #C so
    #C             2**(1/BETA) - 1 = (PAR(1)/Ax)**2
    #C
    #C             Ax**2 = PAR(1)**2/[2**(1/BETA) - 1]
    #C
    #C When BETA = 2.5, Ax**2 = 3.129813 * PAR(1)**2
    #C
    #C Hence, let us use
    #C
    #C                                  1
    #C F = ---------------------------------------------------------------
    #C     P(1)*P(2)*{1+0.3195079*[(X/P(1))**2+(Y/P(2))**2+(XY*P(3))]**2.5
    #C 
    #C neglecting a constant of proportionality.
    #C
    #         ALPHA = 0.3195079
    #         TALPHA = 0.6390158
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         P1P2 = PAR(1)*PAR(2)
    #         XY = DX*DY
    #C
    #         DENOM = 1. + ALPHA*(DX**2/P1SQ + DY**2/P2SQ + XY*PAR(3))
    #         IF (DENOM .GT. 1.E4) RETURN
    #         FUNC = 1. / (P1P2 * DENOM**PAR(4))
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = (PAR(4) - 1.) * FUNC
    #            P4FOD = PAR(4)*ALPHA*PROFIL/DENOM
    #            DHDXC = P4FOD*(2.*DX/P1SQ + DY*PAR(3))
    #            DHDYC = P4FOD*(2.*DY/P2SQ + DX*PAR(3))
    #            IF (IDERIV .GT. 0) THEN
    #               TERM(1) = (2.*P4FOD*DX**2/P1SQ-PROFIL)/PAR(1)
    #               TERM(2) = (2.*P4FOD*DY**2/P2SQ-PROFIL)/PAR(2)
    #               TERM(3) = - P4FOD*XY
    #C              TERM(4) = PROFIL*(1./(PAR(4)-1.)-ALOG(DENOM))
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            XSQ(IX) = X(IX)**2
    #            P1XSQ(IX) = XSQ(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            YSQ = Y**2
    #            P2YSQ = YSQ/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY=X(IX)*Y
    #               DENOM = 1. + ALPHA*(P1XSQ(IX) + P2YSQ + XY*PAR(3))
    #               FUNC = (PAR(4) - 1.) / (P1P2 * DENOM**PAR(4))
    #               P4FOD = PAR(4)*ALPHA*FUNC/DENOM
    #               WP4FOD = WT*P4FOD
    #               WF = WT*FUNC
    #               PROFIL = PROFIL + WF
    #               DHDXC = DHDXC + WP4FOD*(2.*X(IX)/P1SQ + Y*PAR(3))
    #               DHDYC = DHDYC + WP4FOD*(2.*Y/P2SQ + X(IX)*PAR(3))
    #               IF (IDERIV .GT. 0) THEN
    #                  TERM(1) = TERM(1) + 
    #     .                     (2.*WP4FOD*P1XSQ(IX)-WF)/PAR(1)
    #                  TERM(2) = TERM(2) + 
    #     .                     (2.*WP4FOD*P2YSQ-WF)/PAR(2)
    #                  TERM(3) = TERM(3) - WP4FOD*XY
    #C                 TERM(4) = TERM(4) + WF*(1./(PAR(4)-1.)-ALOG(DENOM))
    #               END IF
    #            END DO
    #         END DO
    #C
    #      ELSE IF (IPSTYP .EQ. 4) THEN
    #C
    #C MOFFAT function  BETA = 3.5
    #C
    #C                            BETA-1
    #C F = --------------------------------------------------------
    #C      Ax * Ay * [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]**BETA
    #C
    #C PAR(1) is the HWHM in x at y = 0: 
    #C
    #C             1/2 = 1/[1 + (PAR(1)/Ax)**2]**BETA
    #C so
    #C             2**(1/BETA) - 1 = (PAR(1)/Ax)**2
    #C
    #C             Ax**2 = PAR(1)**2/[2**(1/BETA) - 1]
    #C
    #C When BETA = 2.5, Ax**2 = 4.56593 * PAR(1)**2
    #C
    #C Hence, let us use
    #C
    #C                                  1
    #C F = ---------------------------------------------------------------
    #C     P(1)*P(2)*{1+0.2190137*[(X/P(1))**2+(Y/P(2))**2+(XY*P(3))]**2.5
    #C 
    #C neglecting a constant of proportionality.
    #C
    #         ALPHA = 0.2190137
    #         TALPHA = 0.4380273
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         P1P2 = PAR(1)*PAR(2)
    #         XY = DX*DY
    #C
    #         DENOM = 1. + ALPHA*(DX**2/P1SQ + DY**2/P2SQ + XY*PAR(3))
    #         IF (DENOM .GT. 1.E4) RETURN
    #         FUNC = 1. / (P1P2 * DENOM**PAR(4))
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = (PAR(4) - 1.) * FUNC
    #            P4FOD = PAR(4)*ALPHA*PROFIL/DENOM
    #            DHDXC = P4FOD*(2.*DX/P1SQ + DY*PAR(3))
    #            DHDYC = P4FOD*(2.*DY/P2SQ + DX*PAR(3))
    #            IF (IDERIV .GT. 0) THEN
    #               TERM(1) = (2.*P4FOD*DX**2/P1SQ-PROFIL)/PAR(1)
    #               TERM(2) = (2.*P4FOD*DY**2/P2SQ-PROFIL)/PAR(2)
    #               TERM(3) = - P4FOD*XY
    #C              TERM(4) = PROFIL*(1./(PAR(4)-1.)-ALOG(DENOM))
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            XSQ(IX) = X(IX)**2
    #            P1XSQ(IX) = XSQ(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            YSQ = Y**2
    #            P2YSQ = YSQ/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY=X(IX)*Y
    #               DENOM = 1. + ALPHA*(P1XSQ(IX) + P2YSQ + XY*PAR(3))
    #               FUNC = (PAR(4) - 1.) / (P1P2 * DENOM**PAR(4))
    #               P4FOD = PAR(4)*ALPHA*FUNC/DENOM
    #               WP4FOD = WT*P4FOD
    #               WF = WT*FUNC
    #               PROFIL = PROFIL + WF
    #               DHDXC = DHDXC + WP4FOD*(2.*X(IX)/P1SQ + Y*PAR(3))
    #               DHDYC = DHDYC + WP4FOD*(2.*Y/P2SQ + X(IX)*PAR(3))
    #               IF (IDERIV .GT. 0) THEN
    #                  TERM(1) = TERM(1) + 
    #     .                     (2.*WP4FOD*P1XSQ(IX)-WF)/PAR(1)
    #                  TERM(2) = TERM(2) + 
    #     .                     (2.*WP4FOD*P2YSQ-WF)/PAR(2)
    #                  TERM(3) = TERM(3) - WP4FOD*XY
    #C                 TERM(4) = TERM(4) + WF*(1./(PAR(4)-1.)-ALOG(DENOM))
    #               END IF
    #            END DO
    #         END DO
    #C
    #      ELSE IF (IPSTYP .EQ. 5) THEN
    #C
    #C LORENTZ function
    #C                      1
    #C F = --------------------------------------
    #C     [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]
    #C
    #C PAR(1) is the HWHM in x at y = 0.
    #C
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         P1P2 = PAR(1)*PAR(2)
    #         XY = DX*DY
    #C
    #         DENOM = 1. + DX**2/P1SQ + DY**2/P2SQ + XY*PAR(3)
    #         IF (DENOM .GT. 1.E10) RETURN
    #         FUNC = 1. / DENOM
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = FUNC
    #            WFSQ = FUNC**2
    #            DHDXC = WFSQ*(2.*DX/P1SQ + DY*PAR(3))
    #            DHDYC = WFSQ*(2.*DY/P2SQ + DX*PAR(3))
    #            IF (IDERIV .GT. 0) THEN
    #               TERM(1) = WFSQ*(2.*DX**2/P1SQ)/PAR(1)
    #               TERM(2) = WFSQ*(2.*DY**2/P2SQ)/PAR(2)
    #               TERM(3) = - WFSQ*XY
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            XSQ(IX) = X(IX)**2
    #            P1XSQ(IX) = XSQ(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            YSQ = Y**2
    #            P2YSQ = YSQ/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY=X(IX)*Y
    #               DENOM = 1. + P1XSQ(IX) + P2YSQ + XY*PAR(3)
    #               FUNC = 1. / DENOM
    #               WF = WT*FUNC
    #               WFSQ = WF*FUNC
    #               PROFIL = PROFIL + WF
    #               DHDXC = DHDXC + WFSQ*(2.*X(IX)/P1SQ + Y*PAR(3))
    #               DHDYC = DHDYC + WFSQ*(2.*Y/P2SQ + X(IX)*PAR(3))
    #               IF (IDERIV .GT. 0) THEN
    #                  TERM(1) = TERM(1) + WFSQ*(2.*P1XSQ(IX))/PAR(1)
    #                  TERM(2) = TERM(2) + WFSQ*(2.*P2YSQ)/PAR(2)
    #                  TERM(3) = TERM(3) - WFSQ*XY
    #               END IF
    #            END DO
    #         END DO
    #C
    #      ELSE IF (IPSTYP .EQ. 6) THEN
    #C
    #C Penny function --- Gaussian core plus Lorentzian wings.  The Lorentzian 
    #C is elongated along the x or y axis, the Gaussian may be tilted.
    #C
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         ONEMP3 = 1.-PAR(3)
    #         XY = DX*DY
    #C
    #         RSQ = DX**2/P1SQ + DY**2/P2SQ
    #         IF (RSQ .GT. 1.E10) RETURN
    #C
    #         F = 1./(1.+RSQ)
    #         RSQ = RSQ + XY*PAR(4)
    #         IF (RSQ .LT. 34.) THEN
    #            E = EXP(-0.6931472*RSQ)
    #            FUNC = PAR(3)*E + ONEMP3*F
    #         ELSE
    #            E = 0.
    #            FUNC = ONEMP3*F
    #         END IF
    #C
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = FUNC
    #            DFBY = ONEMP3*F**2
    #            DEBY = 0.6931472*PAR(3)*E
    #            DBYX0 = 2.*DX/P1SQ
    #            DBYY0 = 2.*DY/P2SQ
    #            DHDXC = DEBY*(DBYX0 + DY*PAR(4)) + DFBY*DBYX0
    #            DHDYC = DEBY*(DBYY0 + DX*PAR(4)) + DFBY*DBYY0
    #            IF (IDERIV .GT. 0) THEN
    #               DBYX0 = DBYX0*DX/PAR(1)
    #               DBYY0 = DBYY0*DY/PAR(2)
    #               DFBY = DFBY + DEBY
    #               TERM(1) = DFBY * DBYX0
    #               TERM(2) = DFBY * DBYY0
    #               TERM(3) = E - F
    #               TERM(4) = - DEBY * XY
    #     .              / (0.5 - ABS(PAR(4)))
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            P1XSQ(IX) = X(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            P2YSQ = Y/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY = X(IX)*Y
    #               RSQ = P1XSQ(IX)*X(IX) + P2YSQ*Y
    #               F = 1./(1.+RSQ)
    #               RSQ = RSQ + XY*PAR(4)
    #               IF (RSQ .LT. 34.) THEN
    #                  E = EXP(-0.6931472*RSQ)
    #                  FUNC = PAR(3)*E + ONEMP3*F
    #                  DEBY = 0.6931472*WT*PAR(3)*E
    #               ELSE
    #                  E = 0.
    #                  FUNC = ONEMP3*F
    #                  DEBY = 0.
    #               END IF
    #               PROFIL = PROFIL + WT*FUNC
    #               DFBY = WT*ONEMP3*F**2
    #               DBYX0 = 2.*P1XSQ(IX)
    #               DBYY0 = 2.*P2YSQ
    #               DHDXC = DHDXC + 
    #     .              DEBY*(DBYX0 + DY*PAR(4)) + DFBY*DBYX0
    #               DHDYC = DHDYC +
    #     .              DEBY*(DBYY0 + DX*PAR(4)) + DFBY*DBYY0
    #               IF (IDERIV .GT. 0) THEN
    #                  DBYX0 = DBYX0*DX/PAR(1)
    #                  DBYY0 = DBYY0*DY/PAR(2)
    #                  TERM(1) = TERM(1) + (DFBY+DEBY)*DBYX0
    #                  TERM(2) = TERM(2) + (DFBY+DEBY)*DBYY0
    #                  TERM(3) = TERM(3) + WT*(E-F)
    #                  TERM(4) = TERM(4) - DEBY * XY
    #               END IF
    #            END DO
    #         END DO
    #C
    #      ELSE IF (IPSTYP .EQ. 7) THEN
    #C
    #C Penny function --- Gaussian core plus Lorentzian wings.
    #C The Lorentzian and Gaussian may be tilted in different
    #C directions.
    #C
    #         P1SQ = PAR(1)**2
    #         P2SQ = PAR(2)**2
    #         ONEMP3 = 1.-PAR(3)
    #         XY = DX*DY
    #C
    #         RSQ = DX**2/P1SQ + DY**2/P2SQ
    #         DFBY = RSQ + PAR(5)*XY
    #         IF (DFBY .GT. 1.E10) RETURN
    #         F = 1./(1.+DFBY)
    #C
    #         DEBY = RSQ + PAR(4)*XY
    #         IF (DEBY .LT. 34.) THEN
    #            E = EXP(-0.6931472*DEBY)
    #         ELSE
    #            E = 0.
    #         END IF
    #C
    #         FUNC = PAR(3)*E + ONEMP3*F
    #         IF (FUNC .GE. 0.046) THEN
    #            NPT = 4
    #         ELSE IF (FUNC .GE. 0.0022) THEN
    #            NPT = 3
    #         ELSE IF (FUNC .GE. 0.0001) THEN
    #            NPT = 2
    #         ELSE IF (FUNC .GE. 1.E-10) THEN
    #            PROFIL = FUNC
    #            DFBY = ONEMP3*F**2
    #            DEBY = 0.6931472*PAR(3)*E
    #            DBYX0 = 2.*DX/P1SQ
    #            DBYY0 = 2.*DY/P2SQ
    #            DHDXC = DEBY*(DBYX0 + DY*PAR(4)) + 
    #     .              DFBY*(DBYX0 + DY*PAR(5))
    #            DHDYC = DEBY*(DBYY0 + DX*PAR(4)) + 
    #     .              DFBY*(DBYY0 + DX*PAR(5))
    #            IF (IDERIV .GT. 0) THEN
    #               DBYX0 = DBYX0*DX/PAR(1)
    #               DBYY0 = DBYY0*DY/PAR(2)
    #               TERM(5) = -DFBY * XY
    #               DFBY = DFBY + DEBY
    #               TERM(1) = DFBY * DBYX0
    #               TERM(2) = DFBY * DBYY0
    #               TERM(3) = E - F
    #               TERM(4) = - DEBY * XY
    #            END IF
    #            RETURN
    #         ELSE
    #            RETURN
    #         END IF
    #C
    #         DO IX=1,NPT
    #            X(IX) = DX+D(IX,NPT)
    #            P1XSQ(IX) = X(IX)/P1SQ
    #         END DO
    #C
    #         DO IY=1,NPT
    #            Y = DY+D(IY,NPT)
    #            P2YSQ = Y/P2SQ
    #            DO IX=1,NPT
    #               WT = W(IY,NPT)*W(IX,NPT)
    #               XY = X(IX)*Y
    #               RSQ = P1XSQ(IX)*X(IX) + P2YSQ*Y
    #               F = RSQ + PAR(5)*XY
    #               IF (F .LE. -1.) THEN
    #                  F = 0.
    #               ELSE
    #                  F = 1./(1.+F)
    #               END IF
    #               DEBY = RSQ + PAR(4)*XY
    #               IF (DEBY .LT. 34.) THEN
    #                  E = EXP(-0.6931472*DEBY)
    #                  FUNC = PAR(3)*E + ONEMP3*F
    #                  DEBY = 0.6931472*WT*PAR(3)*E
    #               ELSE
    #                  E = 0.
    #                  FUNC = ONEMP3*F
    #                  DEBY = 0.
    #               END IF
    #               PROFIL = PROFIL + WT*FUNC
    #               DFBY = WT*ONEMP3*F**2
    #               DBYX0 = 2.*P1XSQ(IX)
    #               DBYY0 = 2.*P2YSQ
    #               DHDXC = DHDXC + 
    #     .              DEBY*(DBYX0 + DY*PAR(4)) + 
    #     .              DFBY*(DBYX0 + DY*PAR(5))
    #               DHDYC = DHDYC +
    #     .              DEBY*(DBYY0 + DX*PAR(4)) + 
    #     .              DFBY*(DBYY0 + DX*PAR(5))
    #               IF (IDERIV .GT. 0) THEN
    #                  DBYX0 = DBYX0*DX/PAR(1)
    #                  DBYY0 = DBYY0*DY/PAR(2)
    #                  TERM(1) = TERM(1) + (DFBY+DEBY)*DBYX0
    #                  TERM(2) = TERM(2) + (DFBY+DEBY)*DBYY0
    #                  TERM(3) = TERM(3) + WT*(E-F)
    #                  TERM(4) = TERM(4) - DEBY * XY
    #                  TERM(5) = TERM(5) - DFBY * XY
    #               END IF
    #            END DO
    #         END DO
    #      ELSE
    #         CALL STUPID ('Invalid PSF type.')
    #         CALL BYEBYE
    #      END IF
    #      RETURN
    #      END


def psfnparam(ipstyp, fwhm, label, par, maxpar):
    #INTEGER FUNCTION  NPARAM  (IPSTYP, FWHM, LABEL, PAR, MAXPAR)
    #CHARACTER*8 LABEL
    #INTEGER MAXPAR
    #REAL PAR(MAXPAR)
    #PAR(1) = FWHM/2.
    #PAR(2) = PAR(1)

    par = np.zeros(maxpar,float)
    par[0] = fwhm/2
    par[1] = fwhm/2
    
    #IF (IPSTYP .EQ. 1) THEN
    #   NPARAM = 2
    #   LABEL = 'GAUSSIAN'
    if (ipstyp==1):
        nparam = 2
        label = 'GAUSSIAN'
    #ELSE IF (IPSTYP .EQ. 2) THEN
    #   NPARAM = 3
    #   PAR(3) = 0.
    #   PAR(4) = 1.5
    #   LABEL = 'MOFFAT15'
    elif (ipstyp==2):
        nparam = 3
        par[2] = 0.
        par[3] = 1.5
        label = 'MOFFAT15'
    #ELSE IF (IPSTYP .EQ. 3) THEN
    #   NPARAM = 3
    #   PAR(3) = 0.
    #   PAR(4) = 2.5
    #   LABEL = 'MOFFAT25'
    elif (ipstyp==3):
        nparam = 3
        par[2] = 0.
        par[3] = 2.5
        label = 'MOFFAT25'
    #ELSE IF (IPSTYP .EQ. 4) THEN
    #   NPARAM = 3
    #   PAR(3) = 0.
    #   PAR(4) = 3.5
    #   LABEL = 'MOFFAT35'        
    elif (ipstyp==4):
        nparam = 3
        par[2] = 0.
        par[3] = 3.5
        label = 'MOFFAT35'
    #ELSE IF (IPSTYP .EQ. 5) THEN
    #   NPARAM = 3
    #   PAR(3) = 0.
    #   LABEL = 'LORENTZ '
    elif (ipstyp==5):
        nparam = 3
        par[2] = 0.
        label = 'LORENTZ'
    #ELSE IF (IPSTYP .EQ. 6) THEN
    #   NPARAM = 4
    #   PAR(3) = 0.75
    #   PAR(4) = 0.0
    #   LABEL = 'PENNY1  '
    elif (ipstyp==5):
        nparam = 4
        par[2] = 0.75
        par[3] = 0.0
        label = 'PENNY1'
    #ELSE IF (IPSTYP .EQ. 7) THEN
    #   NPARAM = 5
    #   PAR(3) = 0.75
    #   PAR(4) = 0.0
    #   PAR(5) = 0.0
    #   LABEL = 'PENNY2  '
    elif (ipstyp==6):
        nparam = 5
        par[2] = 0.75
        par[3] = 0.0
        par[4] = 0.0        
        label = 'PENNY2'
    #ELSE
    #   CALL STUPID ('Invalid PSF type: '//CHAR(IPSTYP+48))
    #END IF
    #RETURN
    #END

    return (nparam,par)


def rdpsf(psffile):
    """ Load a DAOPHOT .psf file"""
    # Check if the file exists
    if os.path.exists(psffile) is False:
        raise ValueError(psffile+" NOT FOUND")

    # Check RDPSF in mathsubs.f

    #      INTEGER FUNCTION  RDPSF  (PSFFIL, IPSTYP, PAR, MAXPAR, NPAR,
    #     .     PSF, MAXPSF, MAXEXP, NPSF, NEXP, NFRAC, 
    #     .     PSFMAG, BRIGHT, XPSF, YPSF)
    #
    # Read in the point-spread function
    #
    #      IMPLICIT NONE
    #      INTEGER MAXPSF, MAXTYP, MAXPAR, MAXEXP
    #      PARAMETER (MAXTYP=7)
    #
    #      REAL PAR(MAXPAR), PSF(MAXPSF,MAXPSF,MAXEXP)
    #
    #      CHARACTER*30 PSFFIL
    #      CHARACTER*8 LABEL, CHECK
    #      REAL PSFMAG, BRIGHT, XPSF, YPSF
    #      INTEGER I, J, K, IPSTYP, NPSF, NPAR, NEXP, NFRAC, ISTAT
    #      INTEGER NTERM, NPARAM
    #
    #      CALL INFILE (3, PSFFIL, ISTAT)
    #      IF (ISTAT .NE. 0) THEN
    #         RDPSF = -1
    #         RETURN
    #      END IF

    # Read header line
    #      READ (3,302,IOSTAT=ISTAT) LABEL, NPSF, NPAR, NEXP, NFRAC, PSFMAG, 
    #     .     BRIGHT, XPSF, YPSF
    #  302 FORMAT (1X, A8, 4I5, F9.3, F15.3, 2F9.1)    
    lines = dln.readlines(psffile)
    nlines = len(lines)
    line1 = lines[0]
    fmt = '(1X, A8, 4I5, F9.3, F15.3, 2F9.1)'
    label,npsf,npar,nexp,nfrac,psfmag,bright,xpsf,ypsf = dln.fread(line1,fmt)
    label = label.strip()
    header = {'label':label, 'npsf':npsf, 'npar':npar, 'nexp':nexp, 'nfrac':nfrac,
              'psfmag':psfmag, 'bright':bright, 'xpsf':xpsf, 'ypsf':ypsf}
    
    # Checking something here, maybe that the parameters are okay
    #      DO IPSTYP=1,MAXTYP
    #         I = NPARAM(IPSTYP, 1., CHECK, PAR, MAXPAR)
    #         IF ((LABEL .EQ. CHECK) .AND. (I .EQ. NPAR)) GO TO 1100
    #      END DO
    #      CALL STUPID ('Inappropriate PSF: '//LABEL)


    # Read in the parameters
    # 1100 READ (3,301,IOSTAT=ISTAT) (PAR(I), I=1,NPAR)
    #  301 FORMAT (1X, 6E13.6)
    line1 = lines[1]
    par = np.zeros(npar,float)
    for i in range(npar):
        par[i] = float(line1[1+i*13:1+(i+1)*13])
    
    # Read in the data
    #      NTERM = NEXP+NFRAC
    #      IF (NTERM .GE. 1) THEN
    #         DO K=1,NTERM
    #            READ (3,311,IOSTAT=ISTAT) ((PSF(I,J,K), I=1,NPSF), J=1,NPSF)
    #  311       FORMAT (1X, 6E13.6)
    #         END DO
    #      END IF
    nterm = nexp + nfrac
    header['nterm'] = nterm
    if nterm>=1:
        # Put all of the data into one long array
        bigline = ''
        for i in range(nlines-2):
            bigline += lines[i+2][1:]
        nnum = int(len(bigline)/13)
        numbers = np.zeros(nnum,float)
        for i in range(nnum):
            numbers[i] = float(bigline[i*13:(i+1)*13])
        # Now put in 3D array
        psf = np.zeros((npsf,npsf,nterm),float)
        n2 = npsf*npsf
        for k in range(nterm):
            numbers1 = numbers[k*n2:(k+1)*n2]
            psf[:,:,k] = numbers1.reshape(npsf,npsf)

        # plt.imshow(psf[:,:,0],origin='lower',aspect='auto')   

    # header, par, psf
    return (header, par, psf)
        
    #return (ipstyp, par, maxpar, npar, psf, maxpsf, maxexp, npsf, nexp, nfrac, psfmag, bright, xpsf, ypsf)

    
    
class PSF:
    """ DAOPHOT PSF class."""
    
    def __init__(self,header,par,psf):
        # Initalize the psf object
        self.header = header
        self.par = par
        self.psf = psf

    def call(self,x,y,mag,full=False,deriv=False,origin=0):
        """ Create a PSF image."""

        # CHECK ADDSTAR.F

        # have option to return the derivative
        
    def __str__(self):
        pass

    @classmethod
    def read(self,filename):
        """ Read in a PSF from a .psf file."""
        # Check if the file exists
        if os.path.exists(filename) is False:
            raise ValueError(filename+" NOT FOUND")
        # Load the file
        header, par, psf = rdpsf(filename)
        # Initalize the psf object
        return PSF(header, par, psf)
        
    def write(self,outfile):
        """ Write the PSF to a file."""
        pass


