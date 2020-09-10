#!/usr/bin/env python

# Fit proper motion and parallax using ra/dec/mjd data

# Most of this code was taken from here:
# https://github.com/ctheissen/WISE_Parallaxes/blob/master/WISE_Parallax.py

import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
from astropy import units as u
from scipy.optimize import curve_fit, minimize
from astropy.time import Time
import astropy.coordinates as coords
from dlnpyutils import utils as dln, coords as dcoords

# Set some constants
d2a  = 3600.
d2ma = 3600000.
d2y  = 1/365.25 

def astrometryfunc(x, Delta1, Delta2, PMra, PMdec, pi):
    """ Compute proper motion and parallax model for a set of ra/dec/mjd values."""
    # x: input list of central RA and DEC positions and array of MJDs
    # Delta1: initial dRA position
    # Delta2: initial dDEC position
    # PMra: proper motion in RA (mas/yr)
    # PMdec: proper motion in DEC (mas/yr)
    # pi: parallax (mas)
    
    ra0, dec0, mjds = x
    n = len(mjds)
    years = (mjds - mjds[0])*d2y
    ras = np.zeros(n,np.float64)+ra0
    decs = np.zeros(n,np.float64)+dec0    
    
    bary = coords.get_body_barycentric('earth', Time(mjds, format='mjd'))    
    
    # Parallax factors   
    Fac1 = (bary.x * np.sin(ras*np.pi/180.) - bary.y * np.cos(ras*np.pi/180.) ) 
    Fac2 = bary.x * np.cos(ras*np.pi/180.) * np.sin(decs*np.pi/180.) + \
           bary.y * np.sin(ras*np.pi/180.) * np.sin(decs*np.pi/180.) - \
           bary.z * np.cos(decs*np.pi/180.)

    RAsend  = Delta1 + PMra  * years + pi * Fac1.value
    DECsend = Delta2 + PMdec * years + pi * Fac2.value

    return np.concatenate( [RAsend, DECsend]).flatten()


def fit(cat):
    """ Fit proper motion and parallax to ra/dec/mjd data in a table."""
    mjd = cat['mjd']
    ra = cat['ra']
    raerr = cat['raerr']
    dec = cat['dec']
    decerr = cat['decerr']

    # Compute relative positions
    cenra = np.mean(ra)
    cendec = np.mean(dec)
    lon,lat = dcoords.rotsphcen(ra,dec,cenra,cendec,gnomic=True)
    lon *= d2a
    lat *= d2a

    # Fit proper motion and parallax
    pars, cov = curve_fit(astrometryfunc, [ra, dec, mjd] , 
                            np.concatenate( [lon,lat] ).flatten(), 
                            sigma=np.concatenate( [ raerr, decerr ] ).flatten() )

    return pars,cov


def plotfit(cat,pars,cov,savefig=None):
    """ Plot a figure of the data and the proper motion/parallax fit."""
    
    plt.rcParams.update({'font.size': 12})

    # Compute relative positions
    cenra = np.mean(cat['ra'])
    cendec = np.mean(cat['dec'])
    lon,lat = dcoords.rotsphcen(cat['ra'],cat['dec'],cenra,cendec,gnomic=True)
    lon *= d2a
    lat *= d2a

    # Array of MJDs for model curve
    mjd = np.linspace(np.min(cat['mjd']),np.max(cat['mjd']),100)
    out = astrometryfunc([cenra,cendec,mjd],pars[0],pars[1],pars[2],pars[3],pars[4])
    ll = out[0:100]
    bb = out[100:]
    
    # Plot the model and data
    plt.plot(ll,bb)
    plt.errorbar(lon,lat,xerr=cat['raerr'],yerr=cat['decerr'],fmt='o',color='black',
                 markersize=5,ecolor='lightgray',elinewidth=2,linestyle='none',capsize=0)
    plt.xlabel('dRA (arcsec)')
    plt.ylabel('dDEC (arcsec)')
    xr = dln.minmax(np.concatenate((lon,ll)))
    xr = [xr[0]-0.05*dln.valrange(xr),xr[1]+0.05*dln.valrange(xr)]
    yr = dln.minmax(np.concatenate((lat,bb)))
    yr = [yr[0]-0.05*dln.valrange(yr),yr[1]+0.05*dln.valrange(yr)]
    plt.xlim(xr)
    plt.ylim(yr)
    perr = np.sqrt(np.diag(cov))
    plt.annotate(r'$\mu_\alpha$ = %5.3f $\pm$ %5.3f mas/yr' % (pars[2],perr[2]) + '\n' +
                 r'$\mu_\delta$ = %5.3f $\pm$ %5.3f mas/yr' % (pars[3],perr[3]) + '\n' + 
                 r'$\pi$ = %5.3f $\pm$ %5.3f mas' % (pars[4],perr[4]),
                 xy=(xr[0]+0.05*dln.valrange(xr),yr[1]-0.20*dln.valrange(yr)),ha='left')
    if savefig is not None:
        plt.savefig(savefig)
