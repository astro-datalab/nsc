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
from scipy.interpolate import interp2d, RectBivariateSpline
import time
from dlnpyutils import utils as dln

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
    
    # Numerically integrate a Gaussian function 
    #
    #          F = EXP {-0.5*[(x-XO)/SIGMA]**2 },
    #
    # from XIN-0.5 to XIN+0.5 using Gauss-Legendre integration.  BETA
    # is the half-width at half-maximum, which is equal to 1.17741 * SIGMA.
    # Thus,
    #
    #          F = EXP {-0.6931472*[(x-XO)/BETA]**2 }.
    #
    # Also: provide the first derivative of the integral with respect to 
    # Xo and BETA.  Use Gauss-Legendre integration.

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
    
    # GAUSSIAN
    #
    #     F = ERFX * ERFY / (PAR(1) * PAR(2))
    #
    # PAR(1) is the HWHM in X; sigma(x) = 0.8493218 * HWHM
    # PAR(2) is the HWHM in Y; ditto

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
    
    # LORENTZ function
    #                      1
    # F = --------------------------------------
    #     [1 + (X/Ax)**2 + (Y/Ay)**2 + (XY*Axy)]
    #
    # PAR(1) is the HWHM in x at y = 0.

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
       
    # Penny function --- Gaussian core plus Lorentzian wings.  The Lorentzian 
    # is elongated along the x or y axis, the Gaussian may be tilted.

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

    # Penny function --- Gaussian core plus Lorentzian wings.
    # The Lorentzian and Gaussian may be tilted in different
    # directions.

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

    # Compute the value of an ANALYTIC prfile for a point DX,DY distant
    # from the centroid.  Return both the computed value and its
    # first derivatives with respect to x and y.  If IDERIV .NE. 0,
    # return also the first derivatives with respect to all the parameters
    # defining the profile.

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


def usepsf(ipstyp,dx,dy,bright,par,psf,npsf,npar,nexp,nfrac,deltax,deltay,
           ideriv=True,binned=True):
    """ Evaluate the PSF for a point."""

    # Evaluate the PSF for a point distant DX, DY from the center of a
    # star located at relative frame coordinates DELTAX, DELTAY.

    maxpsf = 207
    maxpar = 6
    maxexp = 10
    
    n = dln.size(dx)
    nterm = nexp + nfrac    
    upsf,dvdxc,dvdyc,junk = profile(ipstyp,dx,dy,par,ideriv=ideriv,binned=binned)
    upsf *= bright
    dvdxc *= bright
    dvdyc *= bright
    if nterm<0 or psf is None:
        return upsf,dvdxc,dvdyc

    middle = npsf//2

    # The PSF look-up tables are centered at (MIDDLE, MIDDLE).

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

    xx = 2.0*dx+middle
    lx = np.array(xx).astype(int)
    yy = 2.0*dy+middle
    ly = np.array(yy).astype(int)

    # This point in the stellar profile lies between columns LX and LX+1,
    # and between rows LY and LY+1 in the look-up tables.
    
    for k in range(nterm):
        f_psf = RectBivariateSpline(np.arange(npsf), np.arange(npsf), psf[:,:,k],kx=3,ky=3,s=0)
        corr = f_psf(xx,yy,grid=False)
        upsf += junk[k]*corr
        if ideriv is True:
            dfdx = f_psf(xx,yy,dx=1,grid=False)
            dfdy = f_psf(xx,yy,dy=1,grid=False)
            dvdxc -= junk[k]*dfdx
            dvdyc -= junk[k]*dfdy
        
    return upsf,dvdxc,dvdyc


def psfnparam(ipstyp, fwhm):
    """
    Figure out the label for the input PSF type and
    initialize the parameter array.
    """

    maxpar = 5
    par = np.zeros(maxpar,float)
    par[0] = fwhm/2
    par[1] = fwhm/2
    
    if (ipstyp==1):
        nparam = 2
        label = 'GAUSSIAN'
    elif (ipstyp==2):
        nparam = 3
        par[2] = 0.
        par[3] = 1.5
        label = 'MOFFAT15'
    elif (ipstyp==3):
        nparam = 3
        par[2] = 0.
        par[3] = 2.5
        label = 'MOFFAT25'
    elif (ipstyp==4):
        nparam = 3
        par[2] = 0.
        par[3] = 3.5
        label = 'MOFFAT35'
    elif (ipstyp==5):
        nparam = 3
        par[2] = 0.
        label = 'LORENTZ'
    elif (ipstyp==6):
        nparam = 4
        par[2] = 0.75
        par[3] = 0.0
        label = 'PENNY1'
    elif (ipstyp==7):
        nparam = 5
        par[2] = 0.75
        par[3] = 0.0
        par[4] = 0.0        
        label = 'PENNY2'
    else:
        raise ValueError(str(ipstyp)+' PSF type not supported')

    return (nparam,par,label)


def rdpsf(psffile):
    """ Load a DAOPHOT .psf file"""
    # Check if the file exists
    if os.path.exists(psffile) is False:
        raise ValueError(psffile+" NOT FOUND")

    # Read in the point-spread function
    maxtyp = 7

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

def numinpvals(inpvals):
    """ Helper function to get the number of input values to PSF()."""
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
    """ Helper function get the next input values to PSF()."""
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


# Fitting the PSF to selected stars
# psf.f
# getpsf() generates PSF from several stars
# fitana() fits the ANALYTIC profile to selected stars
#    works for a SINGLE ipstyp of PSF

# get subarray around each star and check it for invalid pixels

# After running fitana() it generates the look-up table
# of corrections.

# fitana() calls profil() and uses the derivatives dhdxc/dhdyc

# Can I use curve_fit() to do this?
