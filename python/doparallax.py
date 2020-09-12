#!/usr/bin/env python

# Fit proper motion and parallax using ra/dec/mjd data

# Most of this code was taken from here:
# https://github.com/ctheissen/WISE_Parallaxes/blob/master/WISE_Parallax.py

import os, sys
import numpy as np
from astropy.table import Table, vstack, join
#import matplotlib.pyplot as plt
from astropy import units as u
from scipy.optimize import curve_fit, minimize
from astropy.time import Time
import astropy.coordinates as coords
from dlnpyutils import utils as dln, coords as dcoords
from argparse import ArgumentParser
import time
from dl import queryClient as qc
import psycopg2 as pq

# Set some constants
d2a  = 3600.
d2ma = 3600000.
d2y  = 1/365.25 

def astrometryfunc(x, Delta1, Delta2, PMra, PMdec, pi):
    """ Compute proper motion and parallax model for a set of ra/dec/mjd values."""
    # x: input list of central RA and DEC positions and array of MJDs
    # Delta1: initial dRA position
    # Delta2: initial dDEC position
    # PMra: proper motion in RA (arcsec/yr)
    # PMdec: proper motion in DEC (arcsec/yr)
    # pi: parallax (arcsec)
    
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
    plt.annotate(r'$\mu_\alpha$ = %5.3f $\pm$ %5.3f mas/yr' % (pars[2]*1e3,perr[2]*1e3) + '\n' +
                 r'$\mu_\delta$ = %5.3f $\pm$ %5.3f mas/yr' % (pars[3]*1e3,perr[3]*1e3) + '\n' + 
                 r'$\pi$ = %5.3f $\pm$ %5.3f mas' % (pars[4]*1e3,perr[4]*1e3),
                 xy=(xr[0]+0.05*dln.valrange(xr),yr[1]-0.20*dln.valrange(yr)),ha='left')
    if savefig is not None:
        plt.savefig(savefig)


# Main command-line program
if __name__ == "__main__":
    parser = ArgumentParser(description='Run Doppler fitting on spectra')
    parser.add_argument('healpix0', type=int, nargs=1, help='Starting healpix')
    parser.add_argument('healpix1', type=int, nargs=1, help='Ending healpix')
    args = parser.parse_args()

    t0 = time.time()
    pix0 = args.healpix0[0]
    pix1 = args.healpix1[0]

    connection = pq.connect(user="dlquery",host="db01.datalab.noao.edu",
                            password="",port = "5432",database = "tapdb")
    cur = connection.cursor()
    osql = '''select id,ra,dec,gmag,imag,ndet from nsc_dr2.hpm where 
              pix>=%d and pix<=%d and ndet>=20 and deltamjd > 1095 and
              (POWER(pmra/pmraerr,2) + POWER(pmdec/pmdecerr,2)) > 27.63''' % (pix0,pix1)    
    cur.execute(osql)
    data = cur.fetchall()
    # Convert to numpy structured array
    dtype = np.dtype([('id',np.str,50),('ra',np.float64),('dec',np.float64),
                      ('gmag',float),('imag',float),('ndet',int)])
    obj = np.zeros(len(data),dtype=dtype)
    obj[...] = data
    del(data)
    nobj = len(obj)
    if nobj==0:
        print('No objects found')
        sys.exit()
    print(str(nobj)+' total objects found')

    dt = np.dtype([('objectid',np.str,50),('nmeas',int),('chi2_motion',float),('deltamjd',float),
                   ('class_star',float),('gmag',float),('rmag',float),('imag',float),('zmag',float),
                   ('pars',np.float64,5),('perr',np.float64,5),('cov',np.float64,(5,5))])
    cat = np.zeros(nobj,dtype=dt)


    # HEALPix loop
    pix = np.arange(pix0,pix1+1)
    cnt = 0
    for p in pix:
        #print('Pix = '+str(p))
        osql1 = '''select id,ra,dec,pmra,pmraerr,pmdec,pmdecerr,gmag,rmag,imag,zmag,ndet,class_star,deltamjd
                   from nsc_dr2.hpm where pix=%d and ndet>=20 and deltamjd > 1095 and
                   (POWER(pmra/pmraerr,2) + POWER(pmdec/pmdecerr,2)) > 27.63''' % p
        cur.execute(osql1)
        data = cur.fetchall()
        # Convert to numpy structured array
        dtype = np.dtype([('id',np.str,50),('ra',np.float64),('dec',np.float64),
                          ('pmra',np.float64),('pmraerr',float),('pmdec',np.float64),('pmdecerr',float),
                          ('gmag',float),('rmag',float),('imag',float),('zmag',float),
                          ('ndet',int),('class_star',float),('deltamjd',np.float64)])
        obj1 = np.zeros(len(data),dtype=dtype)
        obj1[...] = data
        del(data)
        nobj1 = len(obj1)
        #print(str(nobj1)+' objects')

        if nobj1>0:
            msql = '''select meas.objectid,meas.ra,meas.raerr,meas.dec,meas.decerr,meas.mjd,meas.class_star
                      from nsc_dr2.meas as meas join nsc_dr2.hpm as obj on meas.objectid=obj.id where
                      obj.pix=%d and obj.ndet>=20 and obj.deltamjd > 1095 and 
                      (POWER(obj.pmra/obj.pmraerr,2) + POWER(obj.pmdec/obj.pmdecerr,2)) > 27.63''' % p
            cur.execute(msql)
            data = cur.fetchall()
            # Convert to numpy structured array
            dtype = np.dtype([('objectid',np.str,50),('ra',np.float64),('raerr',float),('dec',np.float64),
                              ('decerr',float),('mjd',np.float64),('class_star',float)])
            meas = np.zeros(len(data),dtype=dtype)
            meas[...] = data
            del(data)
            nmeas = len(meas)

            # Loop over objects
            for i in range(nobj1):
                ind, = np.where(meas['objectid']==obj1['id'][i])
                nind = len(ind)
                meas1 = meas[ind]
                pars, cov = fit(meas1)
                perr = np.sqrt(np.diag(cov))
                print(str(cnt)+' '+obj1['id'][i]+' '+str(nind)+' '+str(pars))
                cat['objectid'][cnt] = obj1['id'][i]
                cat['nmeas'][cnt] = nind
                cat['chi2_motion'][cnt] = (obj1['pmra'][i]/obj1['pmraerr'][i])**2 + (obj1['pmdec'][i]/obj1['pmdecerr'][i])**2
                cat['deltamjd'][cnt] = obj1['deltamjd'][i]
                cat['class_star'][cnt] = obj1['class_star'][i]
                cat['gmag'][cnt] = obj1['gmag'][i]
                cat['rmag'][cnt] = obj1['rmag'][i]
                cat['imag'][cnt] = obj1['imag'][i]
                cat['zmag'][cnt] = obj1['zmag'][i]
                cat['pars'][cnt] = pars
                cat['perr'][cnt] = perr
                cat['cov'][cnt] = cov
                cnt += 1
                                        
    cur.close()
    connection.close()

    # Write the output file
    outfile = '/net/dl2/dnidever/nsc/instcal/v3/parallax/plx_'+str(pix0)+'_'+str(pix1)+'.fits'
    print('Writing to '+outfile)
    Table(cat).write(outfile,overwrite=True)
    
    dt = time.time()-t0
    print('dt = '+str(dt)+' sec.')

