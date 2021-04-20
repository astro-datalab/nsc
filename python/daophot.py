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

    #      DATA D / 0.00000000,  0.0,        0.0       , 0.0       ,
    #     .        -0.28867513,  0.28867513, 0.0       , 0.0       ,
    #     .        -0.38729833,  0.00000000, 0.38729833, 0.0       ,
    #     .        -0.43056816, -0.16999052, 0.16999052, 0.43056816/
    #      DATA W / 1.00000000,  0.0       , 0.0       , 0.0       ,
    #     .         0.50000000,  0.50000000, 0.0       , 0.0       ,
    #     .         0.27777778,  0.44444444, 0.27777778, 0.0       ,
    #     .         0.17392742,  0.32607258, 0.32607258, 0.17392742/

PROFILE_DATAD = [[0.00000000,  0.0,        0.0       , 0.0       ],
                 [-0.28867513,  0.28867513, 0.0       , 0.0       ],
                 [-0.38729833,  0.00000000, 0.38729833, 0.0       ],
                 [-0.43056816, -0.16999052, 0.16999052, 0.43056816]]
PROFILE_DATAD = np.array(PROFILE_DATAD).T
PROFILE_DATAW = [[1.00000000,  0.0       , 0.0       , 0.0       ],
                 [0.50000000,  0.50000000, 0.0       , 0.0       ],
                 [0.27777778,  0.44444444, 0.27777778, 0.0       ],
                 [0.17392742,  0.32607258, 0.32607258, 0.17392742]]
PROFILE_DATAW = np.array(PROFILE_DATAW).T

def daoerf(xin,x0,beta):
    #      REAL FUNCTION  DAOERF  (XIN, XO, BETA, DFDXO, DFDBET)
    #C
    #C Numerically integrate a Gaussian function 
    #C
    #C          F = EXP {-0.5*[(x-XO)/SIGMA]**2 },
    #C
    #C from XIN-0.5 to XIN+0.5 using Gauss-Legendre integration.  BETA
    #C is the half-width at half-maximum, which is equal to 1.17741 * SIGMA.
    #C Thus,
    #C
    #C          F = EXP {-0.6931472*[(x-XO)/BETA]**2 }.
    #C
    #C Also: provide the first derivative of the integral with respect to 
    #C Xo and BETA.  Use Gauss-Legendre integration.
    #C
    #C-----------------------------------------------------------------------
    #C
    #      IMPLICIT NONE
    #      INTEGER MAXPT
    #      PARAMETER (MAXPT=4)
    #C
    #      REAL DX(MAXPT,MAXPT), WT(MAXPT,MAXPT)
    #C
    #      REAL EXP
    #C
    #      REAL X, XSQ
    #      REAL XIN, XO, BETA, DFDXO, DFDBET, BETASQ, DELTAX, F, WF
    #      INTEGER NPT, I
    #C
    #      DATA DX / 0.00000000,  0.0,        0.0       , 0.0       ,
    #     .         -0.28867513,  0.28867513, 0.0       , 0.0       ,
    #     .         -0.38729833,  0.00000000, 0.38729833, 0.0       ,
    #     .         -0.43056816, -0.16999052, 0.16999052, 0.43056816/
    #      DATA WT / 1.00000000,  0.0       , 0.0       , 0.0       ,
    #     .          0.50000000,  0.50000000, 0.0       , 0.0       ,
    #     .          0.27777778,  0.44444444, 0.27777778, 0.0       ,
    #     .          0.17392742,  0.32607258, 0.32607258, 0.17392742/
    datadx = [[0.00000000,  0.0,        0.0       , 0.0      ] ,
              [-0.28867513,  0.28867513, 0.0       , 0.0      ] ,
              [-0.38729833,  0.00000000, 0.38729833, 0.0      ] ,
              [-0.43056816, -0.16999052, 0.16999052, 0.43056816]]
    datadx = np.array(datadx).T
    datawt = [[1.00000000,  0.0       , 0.0       , 0.0       ],
              [0.50000000,  0.50000000, 0.0       , 0.0       ],
              [0.27777778,  0.44444444, 0.27777778, 0.0       ],
              [0.17392742,  0.32607258, 0.32607258, 0.17392742]]
    datawt = np.array(datawt).T

    #      DAOERF = 0.
    #      DFDXO = 0.
    #      DFDBET = 0.
    #      BETASQ=BETA**2
    #      DELTAX = XIN-XO
    #C

    daoerf = 0.0
    dfdx0 = 0.0
    dfdbet = 0.0
    betasq = beta**2
    deltax = xin-x0
    
    #      XSQ = DELTAX**2
    #      F = XSQ/BETASQ

    xsq = deltax**2
    f = xsq/betasq

    #      IF (F .GT. 34.) RETURN
    #      F = EXP(-0.6931472*F)
    #      IF (F .GE. 0.046) THEN
    #         NPT = 4
    #      ELSE IF (F .GE. 0.0022) THEN
    #         NPT = 3
    #      ELSE IF (F .GE. 0.0001) THEN
    #         NPT = 2
    #      ELSE IF (F .GE. 1.E-10) THEN
    #         DAOERF = F
    #         DFDXO = 1.3862944 * DELTAX * F / BETASQ
    #         DFDBET = 1.3862944 * XSQ * F / (BETASQ*BETA)
    #         RETURN
    #      ELSE
    #         RETURN
    #      END IF

    if (f > 34):
        return (daoerf,dfdx0,dfdbet)
    f = np.exp(-0.6931472*f)
    if (f >= 0.046):
        npt = 4
    elif (f >= 0.0022):
        npt = 3
    elif (f >= 0.0001):
        npt = 2
    elif (f >= 1e-10):
        daoerf = f
        dfdx0 = 1.3862944 * deltax * f / betasq
        dfdbet = 1.3862944 * xsq * f / (betasq*beta)
        return (daoerf,dfdx0,dfdbet)
    else:
        return (daoerf,dfdx0,dfdbet)
    
    #C
    #      DO I=1,NPT
    #         X = DELTAX + DX(I,NPT)
    #         XSQ = X**2
    #         F = EXP(-0.6931472*XSQ/BETASQ)
    #         WF = WT(I,NPT)*F
    #         DAOERF = DAOERF+WF
    #         DFDXO = DFDXO + X*WF
    #         DFDBET = DFDBET + XSQ*WF
    #      END DO
    #      DFDXO = 1.3862944*DFDXO/BETASQ
    #      DFDBET = 1.3862944*DFDBET/(BETASQ*BETA)

    for i in range(npt):
        #x = deltax + datadx[npt-1,i]
        x = deltax + datadx[i,npt-1]        
        xsq = x**2
        f = np.exp(-0.6931472*xsq/betasq)
        #wf = datawt[npt-1,i]*f
        wf = datawt[i,npt-1]*f        
        daoerf += wf
        dfdx0 += x*wf
        dfdbet += xsq*wf
    dfdx0 = 1.3862944*dfdx0/betasq
    dfdbet = 1.3862944*dfdbet/(betasq*beta)

    return (daoerf,dfdx0,dfdbet)

def erf(xin,x0,beta):
    """
    Pythonic version of DAOERF.
    xin = position of pixel
    x0 = center of Gaussian
    beta = half-width at half-maximum, 1.17741 * SIGMA

    Returns
    f - integrated Gaussian function for this pixel
    dfdx0 - derivative with respect to x0
    dfdbet - derivative with respect to beta

    """
    
    #C Numerically integrate a Gaussian function 
    #C
    #C          F = EXP {-0.5*[(x-XO)/SIGMA]**2 },
    #C
    #C from XIN-0.5 to XIN+0.5 using Gauss-Legendre integration.  BETA
    #C is the half-width at half-maximum, which is equal to 1.17741 * SIGMA.
    #C Thus,
    #C
    #C          F = EXP {-0.6931472*[(x-XO)/BETA]**2 }.
    #C
    #C Also: provide the first derivative of the integral with respect to 
    #C Xo and BETA.  Use Gauss-Legendre integration.

    # The Gaussian is FORCED to have an AMPLITUDE of 1.
    # A = ht*sigma*sqrt(2*pi)

    
    # Gaussian sigma
    sigma = beta / 1.17741
    factor = 0.6931472
    
    w = (xin-x0)/beta
    w0 = (xin-0.5-x0)/beta
    w1 = (xin+0.5-x0)/beta
    
    # scipy.special.erf() does:
    # 2/sqrt(pi)*integral(exp(-t**2), t=0..z).
    # goes from -1 to +1, area under Gaussian is 2
    # Gaussian area = height*sigma*sqrt(2*pi)    
    f0 = scipy.special.erf(np.sqrt(factor)*w0)
    f1 = scipy.special.erf(np.sqrt(factor)*w1)
    f = (f1-f0)*np.sqrt(2*np.pi)*sigma/2

    # Derivative wrt x0
    # the derivative removes the integral over x and we just evaluate at the two points
    g0 = np.exp(-factor*w0**2)
    g1 = np.exp(-factor*w1**2)
    dfdx = -(g1-g0)
    
    # Derivative wrt beta
    # phi = normal function
    # Phi = cumulative normal function
    # Integral of x^2*phi(x) = Phi(x) - x*phi(x)
    dfdbet = f/beta - (w1*g1) + (w0*g0)
    
    return f,dfdx,dfdbet


def profile_gaussian(dx,dy,par,ideriv=False,binned=True):
    """ GAUSSIAN PSF analytical profile."""

    # dx/dy can be scalars or 1D arrays
    
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
    n = dln.size(dx)
    p1p2 = par[0]*par[1]
    erfx,dhdxc,dhdsx = erf(dx,0.0,par[0])
    erfy,dhdyc,dhdsy = erf(dy,0.0,par[1])  
    profil = erfx*erfy/p1p2
    dhdxc *= erfy/p1p2
    dhdyc *= erfx/p1p2
    term = np.zeros((2,n),float)    
    if ideriv==True:
        term[0,:] = (dhdsx-erfx/par[0])*erfy/p1p2
        term[1,:] = (dhdsy-erfy/par[1])*erfx/p1p2

    return (profil,dhdxc,dhdyc,term)


def vecbin_moffat(alpha,par,x,y,wt,ideriv):
    # This performs the vectorized version of the binned MOFFAT15/25/35 functions
    # The only differences between the three is the alpha input value
    # x, y and wt should be (2,2,n), or (3,3,n) or (4,4,n)
    # where the first two dimensions are the oversampling of the X and Y dimensions
    p1sq = par[0]**2
    p2sq = par[1]**2
    p1p2 = par[0]*par[1]
    nx,ny,n = x.shape
    term = np.zeros((4,n),float)
    xsq = x**2
    p1xsq = xsq/p1sq
    ysq = y**2
    p2ysq = ysq/p2sq
    xy = x*y
    denom = 1.0 + alpha*(p1xsq + p2ysq + xy*par[2])
    func = (par[3]-1.0) / (p1p2 * denom**par[3])
    p4fod = par[3]*alpha*func/denom
    wp4fod = wt*p4fod
    wf = wt*func
    profil = np.sum(np.sum(wf,axis=0),axis=0)
    dhdxc = np.sum(np.sum(wp4fod*(2.0*x/p1sq + y*par[2]),axis=0),axis=0)
    dhdyc = np.sum(np.sum(wp4fod*(2.0*y/p2sq + x*par[2]),axis=0),axis=0)
    if ideriv==True:
        term[0,:] = np.sum(np.sum((2.0*wp4fod*p1xsq-wf)/par[0],axis=0),axis=0)
        term[1,:] = np.sum(np.sum((2.0*wp4fod*p2ysq-wf)/par[1],axis=0),axis=0)
        term[2,:] = np.sum(np.sum(-wp4fod/xy,axis=0),axis=0)
        term[3,:] = np.sum(np.sum(wf*(1.0/(par[3]-1.0)-np.log(denom)),axis=0),axis=0)
    return profil,dhdxc,dhdyc,term

def profile_moffat(dx,dy,par,beta,ideriv=False,binned=True):
    """ MOFFAT PSF analytical profile."""

    #C MOFFAT function   BETA can be  1.5, 2.5 or 3.5
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

    n = dln.size(dx)
    profil = np.zeros(n,float)
    dhdxc = np.zeros(n,float)
    dhdyc = np.zeros(n,float)
    term = np.zeros((4,n),float)

    if beta==1.5:
        alpha = 0.5874011
    elif beta==2.5:
        alpha = 0.3195079
    elif beta==3.5:
        alpha = 0.2190137
    else:
        raise ValueError('Beta can only be 1.5, 2.5, or 3.5')
        
    p1sq = par[0]**2
    p2sq = par[1]**2
    p1p2 = par[0]*par[1]

    xy = dx*dy
    denom = 1.0 + alpha*(dx**2/p1sq + dy**2/p2sq + xy*par[2])
    func = 1.0 / (p1p2*denom**par[3])
    
    # No binning
    if binned is False:
        profil = (par[3]-1.0)*func
        p4fod = par[3]*alpha*profil/denom
        dhdxc = p4fod*(2.0*dx/p1sq + dy*par[2])
        dhdyc = p4fod*(2.0*dy/p2sq + dx*par[2])        
        if (ideriv==True):
            term[0,:] = (2.0*p4fod*dx**2/p1sq-profil)/par[0]
            term[1,:] = (2.0*p4fod*dy**2/p2sq-profil)/par[1]
            term[2,:] = -p4fod*dy
            term[3,:] = profil*(1.0/(par[3]-1.0)-np.log(denom))
        return (profil,dhdxc,dhdyc,term)
            
    # Figure out the binned oversampling, 1-4 pixels
    npt = np.zeros(n,int)-1
    ind4, = np.where((denom<=5e6) & (func >= 0.046))
    nind4 = len(ind4)
    if nind4>0: npt[ind4] = 4
    ind3, = np.where((denom<=5e6) & (func >= 0.0022) & (func < 0.046))
    nind3 = len(ind3)
    if nind3>0: npt[ind3] = 3
    ind2, = np.where((denom<=5e6) & (func >= 0.0001) & (func < 0.0022))
    nind2 = len(ind2)
    if nind2>0: npt[ind2] = 2
    ind1, = np.where((denom<=5e6) & (func >= 1e-10)  & (func < 0.0001))
    nind1 = len(ind1)
    if nind1>0: npt[ind1] = 1
    ind0, = np.where(denom>5e6)
    nind0 = len(ind0)
    if nind0>0: npt[ind1] = 0    


    if nind1>0:
        profil[ind1] = (par[3]-1.0)*func[ind1]
        p4fod = par[3]*alpha*profil[ind1]/denom[ind1]
        dhdxc[ind1] = p4fod*(2.0*dx[ind1]/p1sq + dy[ind1]*par[2])
        dhdyc[ind1] = p4fod*(2.0*dy[ind1]/p2sq + dx[ind1]*par[2])        
        if (ideriv==True):
            term[0,ind1] = (2.0*p4fod*dx[ind1]**2/p1sq-profil[ind1])/par[0]
            term[1,ind1] = (2.0*p4fod*dy[ind1]**2/p2sq-profil[ind1])/par[1]
            term[2,ind1] = -p4fod*dy[ind1]
            term[3,ind1] = profil[ind1]*(1.0/(par[3]-1.0)-np.log(denom[ind1]))

    for pt in np.arange(3)+2:
        ind, = np.where(npt==pt)
        if len(ind)>0:
            datad2 = np.tile(PROFILE_DATAD[0:pt,pt-1],(pt,1))
            x = np.tile(dx[ind],(pt,pt,1)) + np.tile(datad2.T,(len(ind),1,1)).T
            y = np.tile(dy[ind],(pt,pt,1)) + np.tile(datad2,(len(ind),1,1)).T
            wt = np.tile(np.outer(PROFILE_DATAW[0:pt,pt-1],PROFILE_DATAW[0:pt,pt-1]),(len(ind),1,1)).T          
            profil1,dhdxc1,dhdyc1,term1 = vecbin_moffat(alpha,par,x,y,wt,ideriv)
            profil[ind] = profil1
            dhdxc[ind] = dhdxc1
            dhdyc[ind] = dhdyc1        
            term[:,ind] = term1

    return (profil,dhdxc,dhdyc,term)


def vecbin_lorentz(par,x,y,wt,ideriv):
    # This performs the vectorized version of the binned LORENTZ function
    # x, y and wt should be (2,2,n), or (3,3,n) or (4,4,n)
    # where the first two dimensions are the oversampling of the X and Y dimensions
    p1sq = par[0]**2
    p2sq = par[1]**2
    p1p2 = par[0]*par[1]
    nx,ny,n = x.shape
    term = np.zeros((3,n),float)
    xsq = x**2
    p1xsq = xsq/p1sq
    ysq = y**2
    p2ysq = ysq/p2sq
    xy = x*y
    denom = 1.0 + p1xsq + p2ysq + xy*par[2]
    func = 1.0 / denom
    wf = wt*func
    wfsq = wf*func
    profil = np.sum(np.sum(wf,axis=0),axis=0)
    dhdxc = np.sum(np.sum(wfsq*(2.0*x/p1sq + y*par[2]),axis=0),axis=0)
    dhdyc = np.sum(np.sum(wfsq*(2.0*y/p2sq + x*par[2]),axis=0),axis=0)
    if ideriv==True:
        term[0,:] = np.sum(np.sum(wfsq*(2.0*p1xsq)/par[0],axis=0),axis=0)
        term[1,:] = np.sum(np.sum(wfsq*(2.0*p2ysq)/par[1],axis=0),axis=0)
        term[2,:] = np.sum(np.sum(-wfsq/xy,axis=0),axis=0)
    return profil,dhdxc,dhdyc,term

def profile_lorentz(dx,dy,par,ideriv=False,binned=True):
    """ LORENTZ PSF analytical profile."""
    
    #      ELSE IF (IPSTYP .EQ. 5) THEN
    #C
    #C LORENTZ function
    #C                      1
    #C F = --------------------------------------
    #C     [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]
    #C
    #C PAR(1) is the HWHM in x at y = 0.
    #C

    n = dln.size(dx)
    profil = np.zeros(n,float)
    dhdxc = np.zeros(n,float)
    dhdyc = np.zeros(n,float)
    term = np.zeros((3,n),float)

    p1sq = par[0]**2
    p2sq = par[1]**2
    p1p2 = par[0]*par[1]

    xy = dx*dy
    denom = 1.0 + dx**2/p1sq + dy**2/p2sq + xy*par[2]
    func = 1.0 / denom
    
    # No binning
    if binned is False:
        profil = func
        wfsq = func**2
        dhdxc = wfsq*(2.0*dx/p1sq + dy*par[2])
        dhdyc = wfsq*(2.0*dy/p2sq + dx*par[2])
        if (ideriv==True):
            term[0,:] = wfsq*(2.0*dx**2/p1sq)/par[0]
            term[1,:] = wfsq*(2.0*dy**2/p2sq)/par[1]
            term[2,:] = -wfsq*xy
        return (profil,dhdxc,dhdyc,term)
            
    # Figure out the binned oversampling, 1-4 pixels
    npt = np.zeros(n,int)-1
    ind4, = np.where((denom<=1e10) & (func >= 0.046))
    nind4 = len(ind4)
    if nind4>0: npt[ind4] = 4
    ind3, = np.where((denom<=1e10) & (func >= 0.0022) & (func < 0.046))
    nind3 = len(ind3)
    if nind3>0: npt[ind3] = 3
    ind2, = np.where((denom<=1e10) & (func >= 0.0001) & (func < 0.0022))
    nind2 = len(ind2)
    if nind2>0: npt[ind2] = 2
    ind1, = np.where((denom<=1e10) & (func >= 1e-10)  & (func < 0.0001))
    nind1 = len(ind1)
    if nind1>0: npt[ind1] = 1
    ind0, = np.where(denom>1e10)
    nind0 = len(ind0)
    if nind0>0: npt[ind1] = 0    


    if nind1>0:
        profil[ind1] = func[ind1]
        wfsq = func[ind1]**2
        dhdxc = wfsq*(2.0*dx[ind1]/p1sq + dy[ind1]*par[2])
        dhdyc = wfsq*(2.0*dy[ind1]/p2sq + dx[ind1]*par[2])
        if (ideriv==True):
            term[0,:] = wfsq*(2.0*dx[ind1]**2/p1sq)/par[0]
            term[1,:] = wfsq*(2.0*dy[ind1]**2/p2sq)/par[1]
            term[2,:] = -wfsq*xy[ind1]

    for pt in np.arange(3)+2:
        ind, = np.where(npt==pt)
        if len(ind)>0:
            datad2 = np.tile(PROFILE_DATAD[0:pt,pt-1],(pt,1))
            x = np.tile(dx[ind],(pt,pt,1)) + np.tile(datad2.T,(len(ind),1,1)).T
            y = np.tile(dy[ind],(pt,pt,1)) + np.tile(datad2,(len(ind),1,1)).T
            wt = np.tile(np.outer(PROFILE_DATAW[0:pt,pt-1],PROFILE_DATAW[0:pt,pt-1]),(len(ind),1,1)).T          
            profil1,dhdxc1,dhdyc1,term1 = vecbin_lorentz(par,x,y,wt,ideriv)
            profil[ind] = profil1
            dhdxc[ind] = dhdxc1
            dhdyc[ind] = dhdyc1        
            term[:,ind] = term1

    return (profil,dhdxc,dhdyc,term)


def vecbin_penny1(par,dx,dy,x,y,wt,ideriv):
    # This performs the vectorized version of the binned PENNY1 function
    # dx, dy, x, y and wt should be (2,2,n), or (3,3,n) or (4,4,n)
    # where the first two dimensions are the oversampling of the X and Y dimensions
    p1sq = par[0]**2
    p2sq = par[1]**2
    onemp3 = 1.0-par[2]
    nx,ny,n = x.shape
    term = np.zeros((4,n),float)
    p1xsq = x/p1sq
    p2ysq = y/p2sq
    xy = x*y
    rsq = p1xsq*x + p2ysq*y
    f = 1.0/(1.0+rsq)
    rsq += xy*par[3]
    e = np.zeros((nx,ny,n),float)
    func = np.zeros((nx,ny,n),float)
    deby = np.zeros((nx,ny,n),float)        
    low = (rsq<34.0)
    if np.sum(low)>0:
        e[low] = np.exp(-0.6931472*rsq[low])
        func[low] = par[2]*e[low] + onemp3*f[low]
        deby[low] = 0.6931472*wt[low]*par[2]*e[low]
    if np.sum(~low)>0:
        #e[~low] = 0.0     # already zero
        func[~low] = onemp3*f[~low]
        #deby[~low] = 0.0  # already zero
    profil = np.sum(np.sum(wt*func,axis=0),axis=0)
    dfby = wt*onemp3*f**2
    dbyx0 = 2.0*p1xsq
    dbyy0 = 2.0*p2ysq
    dhdxc = np.sum(np.sum(deby*(dbyx0+dy*par[3])+dfby*dbyx0,axis=0),axis=0)
    dhdyc = np.sum(np.sum(deby*(dbyy0+dx*par[3])+dfby*dbyy0,axis=0),axis=0)
    if ideriv==True:
        dbyx0 = dbyx0*dx/par[0]
        dbyy0 = dbyy0*dy/par[1]
        term[0,:] = np.sum(np.sum((dfby+deby)*dbyx0,axis=0),axis=0)
        term[1,:] = np.sum(np.sum((dfby+deby)*dbyy0,axis=0),axis=0)
        term[2,:] = np.sum(np.sum(wt*(e-f),axis=0),axis=0)
        term[3,:] = np.sum(np.sum(-deby*xy,axis=0),axis=0)
    return profil,dhdxc,dhdyc,term

def profile_penny1(dx,dy,par,ideriv=False,binned=True):
    """ PENNY1 PSF analytical profile."""
       
    #      ELSE IF (IPSTYP .EQ. 6) THEN
    #C
    #C Penny function --- Gaussian core plus Lorentzian wings.  The Lorentzian 
    #C is elongated along the x or y axis, the Gaussian may be tilted.
    #C

    n = dln.size(dx)
    profil = np.zeros(n,float)
    dhdxc = np.zeros(n,float)
    dhdyc = np.zeros(n,float)
    term = np.zeros((4,n),float)

    p1sq = par[0]**2
    p2sq = par[1]**2
    onemp3 = 1.0-par[2]
    xy = dx*dy
    rsq = dx**2/p1sq + dy**2/p2sq

    f = 1.0/(1.0+rsq)
    rsq += xy*par[3]
    e = np.zeros(n,float)
    func = np.zeros(n,float)    
    low = (rsq<34.0)
    if np.sum(low)>0:
        e[low] = np.exp(-0.6931472*rsq[low])
        func[low] = par[2]*e[low] + onemp3*f[low]
    if np.sum(~low)>0:
        e[~low] = 0.0
        func[~low] = onemp3*f[~low]
    
    # No binning
    if binned is False:
        profil = func
        dfby = onemp3*f**2
        deby = 0.6931472*par[2]*e
        dbyx0 = 2.0*dx/p1sq
        dbyy0 = 2.0*dy/p2sq
        dhdxc = deby*(dbyx0 + dy*par[3]) + dfby*dbyx0
        dhdyc = deby*(dbyy0 + dx*par[3]) + dfby*dbyy0
        if (ideriv==True):
            dbyx0 = dbyx0*dx/par[0]
            dbyy0 = dbyy0*dy/par[1]
            dfby += deby
            term[0,:] = dfby * dbyx0
            term[1,:] = dfby * dbyy0
            term[2,:] = e - f
            term[3,:] = -deby * xy / (0.5 - np.abs(par[3]))
        return (profil,dhdxc,dhdyc,term)
            
    # Figure out the binned oversampling, 1-4 pixels
    npt = np.zeros(n,int)-1
    ind4, = np.where((rsq<=1e10) & (func >= 0.046))
    nind4 = len(ind4)
    if nind4>0: npt[ind4] = 4
    ind3, = np.where((rsq<=1e10) & (func >= 0.0022) & (func < 0.046))
    nind3 = len(ind3)
    if nind3>0: npt[ind3] = 3
    ind2, = np.where((rsq<=1e10) & (func >= 0.0001) & (func < 0.0022))
    nind2 = len(ind2)
    if nind2>0: npt[ind2] = 2
    ind1, = np.where((rsq<=1e10) & (func >= 1e-10)  & (func < 0.0001))
    nind1 = len(ind1)
    if nind1>0: npt[ind1] = 1
    ind0, = np.where(rsq>1e10)
    nind0 = len(ind0)
    if nind0>0: npt[ind1] = 0    


    if nind1>0:
        profil[ind1] = func[ind1]
        dfby = onemp3*f[ind1]**2
        deby = 0.6931472*par[2]*e[ind1]
        dbyx0 = 2.0*dx[ind1]/p1sq
        dbyy0 = 2.0*dy[ind1]/p2sq
        dhdxc = deby*(dbyx0 + dy[ind1]*par[3]) + dfby*dbyx0
        dhdyc = deby*(dbyy0 + dx[ind1]*par[3]) + dfby*dbyy0
        if (ideriv==True):
            dbyx0 = dbyx0*dx[ind1]/par[0]
            dbyy0 = dbyy0*dy[ind1]/par[1]
            dfby += deby
            term[0,:] = dfby * dbyx0
            term[1,:] = dfby * dbyy0
            term[2,:] = e[ind1] - f[ind1]
            term[3,:] = -deby * xy[ind1] / (0.5 - np.abs(par[3]))


    for pt in np.arange(3)+2:
        ind, = np.where(npt==pt)
        if len(ind)>0:
            datad2 = np.tile(PROFILE_DATAD[0:pt,pt-1],(pt,1))
            dx2 = np.tile(dx[ind],(pt,pt,1))
            dy2 = np.tile(dy[ind],(pt,pt,1))
            x = dx2 + np.tile(datad2.T,(len(ind),1,1)).T
            y = dy2 + np.tile(datad2,(len(ind),1,1)).T
            wt = np.tile(np.outer(PROFILE_DATAW[0:pt,pt-1],PROFILE_DATAW[0:pt,pt-1]),(len(ind),1,1)).T          
            profil1,dhdxc1,dhdyc1,term1 = vecbin_penny1(par,dx2,dy2,x,y,wt,ideriv)
            profil[ind] = profil1
            dhdxc[ind] = dhdxc1
            dhdyc[ind] = dhdyc1        
            term[:,ind] = term1

    return (profil,dhdxc,dhdyc,term)


def vecbin_penny2(par,dx,dy,x,y,wt,ideriv):
    # This performs the vectorized version of the binned PENNY2 function
    # dx, dy, x, y and wt should be (2,2,n), or (3,3,n) or (4,4,n)
    # where the first two dimensions are the oversampling of the X and Y dimensions
    p1sq = par[0]**2
    p2sq = par[1]**2
    onemp3 = 1.0-par[2]
    nx,ny,n = x.shape
    term = np.zeros((5,n),float)
    p1xsq = x/p1sq
    p2ysq = y/p2sq
    xy = x*y
    rsq = p1xsq*x + p2ysq*y
    f = rsq + par[4]*xy
    low = (f<=-1)
    if np.sum(low)>0:
        f[low] = 0.0
    if np.sum(~low)>0:
        f[~low] = 1.0/(1.0+f[~low])
    deby = rsq + par[3]*xy
    e = np.zeros((nx,ny,n),float)
    func = np.zeros((nx,ny,n),float)    
    low = (deby<34.0)
    if np.sum(low)>0:
        e[low] = np.exp(-0.6931472*deby[low])
        func[low] = par[2]*e[low] + onemp3*f[low]
        deby[low] = 0.6931472*wt[low]*par[2]*e[low]
    if np.sum(~low)>0:
        #e[~low] = 0.0            # already zero
        func[~low] = onemp3*f[~low]
        #deby[~low] = 0.0         $ already zero
    profil = np.sum(np.sum(wt*func,axis=0),axis=0)

    dfby = wt*onemp3*f**2
    dbyx0 = 2.0*p1xsq
    dbyy0 = 2.0*p2ysq
    dhdxc = np.sum(np.sum(deby*(dbyx0+dy*par[3])+dfby*(dbyx0+dy*par[4]),axis=0),axis=0)
    dhdyc = np.sum(np.sum(deby*(dbyy0+dx*par[3])+dfby*(dbyy0+dx*par[4]),axis=0),axis=0)
    if ideriv==True:
        dbyx0 = dbyx0*dx/par[0]
        dbyy0 = dbyy0*dy/par[1]
        term[0,:] = np.sum(np.sum((dfby+deby)*dbyx0,axis=0),axis=0)
        term[1,:] = np.sum(np.sum((dfby+deby)*dbyy0,axis=0),axis=0)
        term[2,:] = np.sum(np.sum(wt*(e-f),axis=0),axis=0)
        term[3,:] = np.sum(np.sum(-deby*xy,axis=0),axis=0)
        term[4,:] = np.sum(np.sum(-dfby*xy,axis=0),axis=0)
    return profil,dhdxc,dhdyc,term

def profile_penny2(dx,dy,par,ideriv=False,binned=True):
    """ PENNY2 PSF analytical profile."""

    #      ELSE IF (IPSTYP .EQ. 7) THEN
    #C
    #C Penny function --- Gaussian core plus Lorentzian wings.
    #C The Lorentzian and Gaussian may be tilted in different
    #C directions.
    #C

    n = dln.size(dx)
    profil = np.zeros(n,float)
    dhdxc = np.zeros(n,float)
    dhdyc = np.zeros(n,float)
    term = np.zeros((5,n),float)

    p1sq = par[0]**2
    p2sq = par[1]**2
    onemp3 = 1.0-par[2]
    xy = dx*dy
    rsq = dx**2/p1sq + dy**2/p2sq
    dfby = rsq + par[4]*xy
    f = 1.0/(1.0+dfby)
    deby = rsq + par[3]*xy

    e = np.zeros(n,float)
    low = (deby<34.0)
    if np.sum(low)>0:
        e[low] = np.exp(-0.6931472*deby[low])
    #if np.sum(~low)>0:  # already zero
    #    e[~low] = 0.0
    func = par[2]*e + onemp3*f

    
    # No binning
    if binned is False:
        profil = func
        dfby = onemp3*f**2
        deby = 0.6931472*par[2]*e
        dbyx0 = 2.0*dx/p1sq
        dbyy0 = 2.0*dy/p2sq
        dhdxc = deby*(dbyx0 + dy*par[3]) + dfby*(dbyx0 + dy*par[4])
        dhdyc = deby*(dbyy0 + dx*par[3]) + dfby*(dbyy0 + dx*par[4])
        if (ideriv==True):
            dbyx0 = dbyx0*dx/par[0]
            dbyy0 = dbyy0*dy/par[1]
            term[4,:] = -dfby * xy
            dfby += deby
            term[0,:] = dfby * dbyx0
            term[1,:] = dfby * dbyy0
            term[2,:] = e - f
            term[3,:] = -deby * xy 
        return (profil,dhdxc,dhdyc,term)
            
    # Figure out the binned oversampling, 1-4 pixels
    npt = np.zeros(n,int)-1
    ind4, = np.where((rsq<=1e10) & (func >= 0.046))
    nind4 = len(ind4)
    if nind4>0: npt[ind4] = 4
    ind3, = np.where((rsq<=1e10) & (func >= 0.0022) & (func < 0.046))
    nind3 = len(ind3)
    if nind3>0: npt[ind3] = 3
    ind2, = np.where((rsq<=1e10) & (func >= 0.0001) & (func < 0.0022))
    nind2 = len(ind2)
    if nind2>0: npt[ind2] = 2
    ind1, = np.where((rsq<=1e10) & (func >= 1e-10)  & (func < 0.0001))
    nind1 = len(ind1)
    if nind1>0: npt[ind1] = 1
    ind0, = np.where(rsq>1e10)
    nind0 = len(ind0)
    if nind0>0: npt[ind1] = 0    


    if nind1>0:
        profil[ind1] = func[ind1]
        dfby = onemp3*f[ind1]**2
        deby = 0.6931472*par[2]*e[ind1]
        dbyx0 = 2.0*dx[ind1]/p1sq
        dbyy0 = 2.0*dy[ind1]/p2sq
        dhdxc = deby*(dbyx0 + dy[ind1]*par[3]) + dfby*(dbyx0 + dy[ind1]*par[4])
        dhdyc = deby*(dbyy0 + dx[ind1]*par[3]) + dfby*(dbyy0 + dx[ind1]*par[4])
        if (ideriv==True):
            dbyx0 = dbyx0*dx[ind1]/par[0]
            dbyy0 = dbyy0*dy[ind1]/par[1]
            term[4,:] = -dfby * xy[ind1]
            dfby += deby
            term[0,:] = dfby * dbyx0
            term[1,:] = dfby * dbyy0
            term[2,:] = e[ind1] - f[ind1]
            term[3,:] = -deby * xy[ind1]

    for pt in np.arange(3)+2:
        ind, = np.where(npt==pt)
        if len(ind)>0:
            datad2 = np.tile(PROFILE_DATAD[0:pt,pt-1],(pt,1))
            dx2 = np.tile(dx[ind],(pt,pt,1))
            dy2 = np.tile(dy[ind],(pt,pt,1))
            x = dx2 + np.tile(datad2.T,(len(ind),1,1)).T
            y = dy2 + np.tile(datad2,(len(ind),1,1)).T
            wt = np.tile(np.outer(PROFILE_DATAW[0:pt,pt-1],PROFILE_DATAW[0:pt,pt-1]),(len(ind),1,1)).T          
            profil1,dhdxc1,dhdyc1,term1 = vecbin_penny2(par,dx2,dy2,x,y,wt,ideriv)
            profil[ind] = profil1
            dhdxc[ind] = dhdxc1
            dhdyc[ind] = dhdyc1        
            term[:,ind] = term1

    return (profil,dhdxc,dhdyc,term)


def profile(ipstyp,dx,dy,par,ideriv=False,binned=True):
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

    maxpar = len(par)
    term = np.zeros(maxpar,float)

    # Moffat
    if ipstyp>=2 and ipstyp<=4:
        beta = {2:1.5, 3:2.5, 4:3.5}[ipstyp]
        out = profile_moffat(dx,dy,par,beta,ideriv=ideriv,binned=binned)
        return out
    
    # Other functions
    profdict = {1:profile_gaussian, 2:profile_moffat15, 3:profile_moffat25,
                4:profile_moffat35, 5:profile_lorentz, 6:profile_penny1,
                7:profile_penny2}
    out = profdict[ipstyp](dx,dy,par,ideriv=ideriv,binned=binned)
    # returns (profil,dhdxc,dhdxy,term)
    return out


def bicubic(f, dx, dy):
    """
    Perform bicubic interpolation in a grid of values.
    And return the derivatives wrt x and y.
    """

    #C
    #      REAL FUNCTION  BICUBC  (F, NBOX, DX, DY, DFDX, DFDY)
    #C
    #C Perform a type of bicubic interpolation in a grid of values.
    #C For a point located DX, DY distant from the corner of the grid
    #C (defined to be 1,1), return both the interpolated value and
    #C its first derivatives with respect to x and y.
    #C
    #      IMPLICIT NONE
    #      INTEGER NBOX
    #      REAL F(NBOX,NBOX), TEMP(4), DFDXT(4)
    #C
    #      REAL DX, DY, DFDX, DFDY, C1, C2, C3, C4
    #      INTEGER JY
    #C
    #C By construction, the point at which we want to estimate the function
    #C will lie between the second and third columns, and between the second
    #C and third rows of F, at a distance of (DX,DY) from the (2,2) element
    #C of F.
    #C
    #      DO JY=1,4
    #         C1 = 0.5*(F(3,JY)-F(1,JY))
    #         C4 = F(3,JY) - F(2,JY) - C1
    #         C2 = 3.*C4 - 0.5*(F(4,JY)-F(2,JY)) + C1
    #         C3 = C4 - C2
    #         C4 = DX*C3
    #         TEMP(JY) = DX*(DX*(C4+C2)+C1)+F(2,JY)
    #         DFDXT(JY)= DX*(C4*3.+2.*C2)+C1
    #      END DO
    #      C1 = 0.5*(TEMP(3)-TEMP(1))
    #      C4 = TEMP(3) - TEMP(2) - C1
    #      C2 = 3.*C4 - 0.5*(TEMP(4)-TEMP(2)) + C1
    #      C3 = C4 - C2
    #      C4 = DY*C3
    #      BICUBC = DY*(DY*(C4+C2)+C1)+TEMP(2)
    #      DFDY = DY*(C4*3.+2.*C2)+C1
    #      C1 = 0.5*(DFDXT(3)-DFDXT(1))
    #      C4 = DFDXT(3) - DFDXT(2) - C1
    #      C2 = 3.*C4 - 0.5*(DFDXT(4)-DFDXT(2)) + C1
    #      C3 = C4 - C2
    #      DFDX = DY*(DY*(DY*C3+C2)+C1)+DFDXT(2)
    #      RETURN
    #      END

    temp = np.zeros(4,float)
    dfdxt = np.zeros(4,float)    
    
    for jy in range(4):
        c1 = 0.5*(f[2,jy]-f[0,jy])
        c4 = f[2,jy] - f[1,jy] - c1
        c2 = 3.0*c4 - 0.5*(f[3,jy]-f[1,jy]) + c1
        c3 = c4 - c2
        c4 = dx*c3
        temp[jy] = dx*(dx*(c4+c2)+c1)+f[1,jy]
        dfdxt[jy] = dx*(c4*3.0+2.0*c2)+c1
    c1 = 0.5*(temp[2]-temp[0])
    c4 = temp[2] - temp[1] - c1
    c2 = 3.0*c4 - 0.5*(temp[3]-temp[1]) + c1
    c3 = c4 - c2
    c4 = dy*c3
    bicubic = dy*(dy*(c4+c2)+c1)+temp[1]
    dfdy = dy*(c4*3.0+2.0*c2)+c1
    c1 = 0.5*(dfdxt[2]-dfdxt[0])
    c4 = dfdxt[2] - dfdxt[1] - c1
    c2 = 3.0*c4 - 0.5*(dfdxt[3]-dfdxt[1]) + c1
    c3 = c4 - c2
    dfdx = dy*(dy*(dy*c3+c2)+c1)+dfdxt[1]
    
    return bicubic, dfdx, dfdy


def usepsf(ipstyp,dx,dy,bright,par,psf,npsf,npar,nexp,nfrac,deltax,deltay,binned=True):
    """ Evaluate the PSF for a point."""
    #    C
    #      REAL FUNCTION  USEPSF  (IPSTYP, DX, DY, BRIGHT, PAR, PSF, 
    #     .     NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
    #C
    #C Evaluate the PSF for a point distant DX, DY from the center of a
    #C star located at relative frame coordinates DELTAX, DELTAY.
    #C
    #      IMPLICIT NONE
    #      INTEGER MAXPSF, MAXPAR, MAXEXP
    #      PARAMETER (MAXPSF=207, MAXPAR=6, MAXEXP=10)
    #C
    #      REAL PAR(*), PSF(MAXPSF,MAXPSF,*), JUNK(MAXEXP)
    #C
    #      REAL PROFIL, BICUBC
    #C
    #      REAL MIDDLE, BRIGHT, DVDXC, DVDYC, DELTAX, DELTAY, XX, YY
    #      REAL DX, DY, CORR, DFDX, DFDY
    #      INTEGER K, LX, LY, IPSTYP
    #      INTEGER NFRAC, NTERM, NPSF, NEXP, NPAR
    #C
    #      NTERM = NEXP + NFRAC
    #      USEPSF = BRIGHT*PROFIL(IPSTYP, DX, DY, PAR, DVDXC, DVDYC, 
    #     .     JUNK, 0)
    #      DVDXC = BRIGHT*DVDXC
    #      DVDYC = BRIGHT*DVDYC
    #      IF (NTERM .LT. 0) RETURN
    #      MIDDLE = (NPSF+1)/2

    maxpsf = 207
    maxpar = 6
    maxexp = 10
    
    n = dln.size(dx)
    #if n>1:
    #    upsf = np.zeros(n,float)
    #    dvdxc = np.zeros(n,float)
    #    dvdyc = np.zeros(n,float)
    #    for i in range(n):
    #        upsf1,dvdxc1,dvdyc1 = usepsf(ipstyp,dx[i],dy[i],bright,par,psf,npsf,npar,nexp,nfrac,deltax,deltay,binned=binned)
    #        upsf[i] = upsf1
    #        dvdxc[i] = dvdxc1
    #        dvdyc[i] = dvdyc1
    #    return upsf,dvdxc,dvdyc

    nterm = nexp + nfrac    
    upsf,dvdxc,dvdyc,junk = profile(ipstyp,dx,dy,par,binned=binned)
    upsf *= bright
    dvdxc *= bright
    dvdyc *= bright
    if nterm<0 or psf is None:
        return upsf,dvdxc,dvdyc

    middle = npsf//2

    import pdb; pdb.set_trace()
    
    #C
    #C The PSF look-up tables are centered at (MIDDLE, MIDDLE).
    #C
    #      IF (NEXP .GE. 0) THEN
    #         JUNK(1) = 1.
    #         IF (NEXP .GE. 2) THEN
    #            JUNK(2) = DELTAX
    #            JUNK(3) = DELTAY
    #            IF (NEXP .GE. 4) THEN
    #               JUNK(4) = 1.5*DELTAX**2-0.5
    #               JUNK(5) = DELTAX*DELTAY
    #               JUNK(6) = 1.5*DELTAY**2-0.5
    #               IF (NEXP .GE. 7) THEN
    #                  JUNK(7) = DELTAX*(5.*JUNK(4)-2.)/3.
    #                  JUNK(8) = JUNK(4)*DELTAY
    #                  JUNK(9) = DELTAX*JUNK(6)
    #                  JUNK(10) = DELTAY*(5.*JUNK(6)-2.)/3.
    #               END IF
    #            END IF
    #         END IF
    #      END IF
    #C

    junk = np.zeros(10,float)
    if (nexp >= 0):
        junk[0] = 1
        if (nexp >= 2):
            junk[1] = deltax
            junk[2] = deltay
            if (nexp >= 4):
                junk[3] = 1.5*deltax**2-0.5
                junk[4] = deltax*deltay
                junk[5] = 1.5*deltay**2-0.5
                if (nexp >= 7):
                    junk[6] = deltax*(5.0*junk[3]-2.0)/3.
                    junk[7] = junk[3]*deltay
                    junk[8] = deltax*junk[5]
                    junk[9] = deltay*(5.0*junk[5]-2.0)/3.      

                    
    #C     IF (NFRAC .GT. 0) THEN
    #C        J = NEXP+1
    #C        JUNK(J) = -2.*(DX - REAL(NINT(DX)))
    #C        J = J+1
    #C        JUNK(J) = -2.*(DY - REAL(NINT(DY)))
    #C        J = J+1
    #C        JUNK(J) = 1.5*JUNK(J-2)**2 - 0.5
    #C        J = J+1
    #C        JUNK(J) = JUNK(J-3)*JUNK(J-2)
    #C        J = J+1
    #C        JUNK(J) = 1.5*JUNK(J-3)**2 - 0.5
    #C     END IF
    #      XX = 2.*DX+MIDDLE
    #      LX = INT(XX)
    #      YY = 2.*DY+MIDDLE
    #      LY = INT(YY)

    xx = 2.0*dx+middle
    lx = int(xx)
    yy = 2.0*dy+middle
    ly = int(yy)

    #C
    #C This point in the stellar profile lies between columns LX and LX+1,
    #C and between rows LY and LY+1 in the look-up tables.
    #C
    #      DO K=1,NTERM
    #         CORR = BICUBC(PSF(LX-1,LY-1,K), MAXPSF, 
    #     .        XX-REAL(LX), YY-REAL(LY), DFDX, DFDY)
    #         USEPSF = USEPSF + JUNK(K)*CORR
    #         DVDXC = DVDXC-JUNK(K)*DFDX
    #         DVDYC = DVDYC-JUNK(K)*DFDY
    #      END DO
    #      RETURN
    #      END
    #C

    for k in range(nterm):
        corr, dfdx, dfdy = bicubic(psf[lx-1:,ly-1:,k], xx-float(lx), yy-float(ly))
        upsf += junk[k]*corr
        dvdxc -= junk[k]*dfdx
        dvdyc -= junk[k]*dfdy        
    
    return upsf,dvdxc,dvdyc


def psfnparam(ipstyp, fwhm):
    #INTEGER FUNCTION  NPARAM  (IPSTYP, FWHM, LABEL, PAR, MAXPAR)
    #CHARACTER*8 LABEL
    #INTEGER MAXPAR
    #REAL PAR(MAXPAR)
    #PAR(1) = FWHM/2.
    #PAR(2) = PAR(1)

    maxpar = 5
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
    elif (ipstyp==6):
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
    elif (ipstyp==7):
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

    return (nparam,par,label)


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
    maxtyp = 7
    #      REAL PAR(MAXPAR), PSF(MAXPSF,MAXPSF,MAXEXP)
    #
    #      CHARACTER*30 PSFFIL
    #      CHARACTER*8 LABEL, CHECK
    #      REAL PSFMAG, BRIGHT, XPSF, YPSF
    #      INTEGER I, J, K, IPSTYP, NPSF, NPAR, NEXP, NFRAC, ISTAT
    #      INTEGER NTERM, NPARAM


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
    ipstyp = {'GAUSSIAN':1, 'MOFFAT15':2, 'MOFFAT25':3, 'MOFFAT35':4,
              'LORENTZ':5, 'PENNY1':6, 'PENNY2':7}[label]
    header = {'label':label, 'ipstyp':ipstyp, 'npsf':npsf, 'npar':npar, 'nexp':nexp,
              'nfrac':nfrac, 'psfmag':psfmag, 'bright':bright, 'xpsf':xpsf, 'ypsf':ypsf,
              'binned':1}

    # Add some other values
    #----------------------
    # PSF radius
    # addstar.f line 74
    psfradius = (float(header['npsf']-1)/2.0 - 1.0)/2.
    header['psfradius'] = psfradius
    
    # psf.f line 477+478
    # XMID = REAL(NCOL-1)/2.  # same as XPSF
    # YMID = REAL(NROW-1)/2.  # same as YPSF
    ncol = int(header['xpsf']*2+1)
    nrow = int(header['ypsf']*2+1)
    header['ncol'] = ncol
    header['nrow'] = nrow    
    
    
    #      DO IPSTYP=1,MAXTYP
    #         I = NPARAM(IPSTYP, 1., CHECK, PAR, MAXPAR)
    #         IF ((LABEL .EQ. CHECK) .AND. (I .EQ. NPAR)) GO TO 1100
    #      END DO
    #      CALL STUPID ('Inappropriate PSF: '//LABEL)

    # Check that the number of input parameters is correct for this PSF type
    # and initial PSF with some extra parameters
    nparcheck,par,label = psfnparam(ipstyp, 1.0)
    if nparcheck != npar:
        raise ValueError('Inappropriate PSF: '+label)
    
    # Read in the parameters
    # 1100 READ (3,301,IOSTAT=ISTAT) (PAR(I), I=1,NPAR)
    #  301 FORMAT (1X, 6E13.6)
    line1 = lines[1]
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
            psf[:,:,k] = numbers1.reshape(npsf,npsf).T
    else:
        psf = None
        
    # header, par, psf
    return (header, par, psf)
        
    #return (ipstyp, par, maxpar, npar, psf, maxpsf, maxexp, npsf, nexp, nfrac, psfmag, bright, xpsf, ypsf)

def numinpvals(inpvals):
    """ Helper function to get the number of input values."""
    if (type(inpvals) is list) or (type(inpvals) is tuple):
        listoflists = False
        if (type(inpvals[0]) is list) or (type(inpvals[0]) is tuple):
            listoflists = True
        if listoflists:
            nstars = len(inpvals)
        else:
            nstars = 1
    else:
        nstars = len(inpvals)
    return nstars
        
def getinpvals(inpvals,i):
    """ Helper function get the next input values."""
    if type(inpvals) is list or type(inpvals) is tuple:
        listoflists = False
        if (type(inpvals[0]) is list) or (type(inpvals[0]) is tuple):
            listoflists = True
        if listoflists:
            x,y,mag = inpvals[i]
        else:
            x,y,mag = inpvals
    else:
        x = inpvals['x'][i]
        y = inpvals['y'][i]
        mag = inpvals['mag'][i]  
    return x,y,mag


        
class PSF:
    """ DAOPHOT PSF class."""
    
    def __init__(self,header,par,psf):
        # Initalize the psf object
        self.header = header
        self.par = par
        self.psf = psf

    def __call__(self,inpvals,xy=None,full=False,deriv=False,origin=0,binned=None):
        """ Create a PSF image."""

        nstars = numinpvals(inpvals)

        # Multiple stars
        #---------------
        if nstars>1:
            # Initialize the output image
            ncol = self.header['ncol']
            nrow = self.header['nrow']            
            image = np.zeros((ncol,nrow),float)
            # Loop over stars
            for i in range(nstars):
                x,y,mag = getinpvals(inpvals,i)
                # Get the PSF image for this star
                xy = self._getxyranges(x,y)  # get X/Y ranges
                upsf = self((x,y,mag),xy=xy,deriv=deriv,origin=origin,binned=binned)
                # Add to full image
                image[xy[0][0]:xy[0][1]+1,xy[1][0]:xy[1][1]+1] += upsf
            return image

        
        # Single star PSF
        #----------------
        x,y,mag = getinpvals(inpvals,0)
        
        # Full image
        if full is True:
            ncol = self.header['ncol']
            nrow = self.header['nrow']
            image = np.zeros((ncol,nrow),float)
                
        # Scale x/y values
        # addstar.f line 190-191
        deltax = (x-1.0)/self.header['xpsf'] - 1.0
        deltay = (y-1.0)/self.header['ypsf'] - 1.0        
        
        # X/Y pixel ranges
        psfradius = self.header['psfradius']        
        npix = 2*int(psfradius)-1
        if full is False:
            if xy is None:
                dx = np.arange(npix)-npix//2
                dx2 = np.repeat(dx,npix).reshape(npix,npix)
                dy = np.arange(npix)-npix//2
                dy2 = np.repeat(dy,npix).reshape(npix,npix).T
                nxpix = npix
                nypix = npix
            else:
                x0,x1 = xy[0]
                y0,y1 = xy[1]
                dx = np.arange(x0,x1+1).astype(float)-x
                nxpix = len(dx)
                dy = np.arange(y0,y1+1).astype(float)-y
                nypix = len(dy)
                dx2 = np.repeat(dx,nypix).reshape(nxpix,nypix)
                dy2 = np.repeat(dy,nxpix).reshape(nypix,nxpix).T                    
        else:
            xy = self._getxyranges(x,y)            
            x0,x1 = xy[0]
            dx = np.arange(x0,x1+1).astype(float)-x
            nxpix = len(dx)
            y0,y1 = xy[1]
            dy = np.arange(y0,y1+1).astype(float)-y
            nypix = len(dy)
            dx2 = np.repeat(dx,nypix).reshape(nxpix,nypix)
            dy2 = np.repeat(dy,nxpix).reshape(nypix,nxpix).T

        if binned is None:
            binned = self.header['binned']
        upsf,dvdxc,dvdyc = usepsf(self.header['ipstyp'],dx2.flatten(),dy2.flatten(),self.header['bright'],
                                  self.par,self.psf,self.header['npsf'],self.header['npar'],
                                  self.header['nexp'],self.header['nfrac'],deltax,deltay,
                                  binned=binned)
        upsf = upsf.reshape(nxpix,nypix)
        dvdxc = dvdxc.reshape(nxpix,nypix)
        dvdyc = dvdyc.reshape(nxpix,nypix)        

        # Impose the psf radius
        rad = np.sqrt(dx2**2+dy2**2)
        upsf[rad>psfradius] = 0.0
        
        # Scale it with the magnitude
        # from addstar.f line 196
        scale = 10.0**(0.4*(self.header['psfmag']-mag))
        upsf *= scale

        # Single PSF
        if full is False:
            # No derivatives
            if deriv is False:
                return upsf
            else:
                return upsf,dvdxc,dvdyc
        else:
            # Add to full image
            image[x0:x1+1,y0:y0+1] = upsf
            return image
            
        
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

    def _getxyranges(self,x,y):
        psfradius = self.header['psfradius']        
        npix = 2*int(psfradius)-1
        ncol = self.header['ncol']
        nrow = self.header['nrow']        
        x0 = int(np.maximum(np.round(x)-npix//2,0))
        x1 = int(np.minimum(np.round(x)+npix//2,ncol-1))
        y0 = int(np.maximum(np.round(y)-npix//2,0))
        y1 = int(np.minimum(np.round(y)+npix//2,nrow-1))
        xy = ((x0,x1),(y0,y1))
        return xy

