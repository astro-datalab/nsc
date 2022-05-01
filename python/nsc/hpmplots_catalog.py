#!/usr/bin/env python

# Script to make diagnostic plots for NSC DR2 high proper motion candidates

import os
from dl import queryClient as qc
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from scipy import stats
from dlnpyutils import utils as dln, coords
#from functions.dlnpyutils import coords
#import functions.parallax as parallax
from astropy.table import Table
from astropy.io import fits
import random
from time import perf_counter 
from argparse import ArgumentParser
import traceback

def make_4panelplot(obj,outdir):
    matplotlib.use('Agg')
    #params = {'tex.usetex': True}
    #plt.rcParams.update(params)
    #plt.rc(usetex = True)

    idv = obj['id']

    #gather all data used
    dataselect = "select mjd,ra,dec,mag_auto,raerr,decerr,filter,fwhm,exposure from nsc_dr2.meas where objectid='" + idv + "'"
    meas = qc.query(sql=dataselect,fmt='table',profile='db01')
    #obj = qc.query(sql="select id,pmra,pmdec from nsc_dr2.object where id='"+ idv +"'",fmt='table',profile='db01')
    #obj = qc.query(sql="select id,pmra,pmdec from nsc_dr2.object where id='"+ idv +"'",fmt='table')

    # Make cut on FWHM
    # maybe only use values for 0.5*fwhm_chip to 1.5*fwhm_chip
    sql = "select chip.* from nsc_dr2.chip as chip join nsc_dr2.meas as meas on chip.exposure=meas.exposure and chip.ccdnum=meas.ccdnum"
    sql += " where meas.objectid='"+idv+"'"
    chip = qc.query(sql=sql,fmt='table')
    ind3,ind4 = dln.match(chip['exposure'],meas['exposure'])
    si = np.argsort(ind4)   # sort by input meas catalog
    ind3 = ind3[si]
    ind4 = ind4[si]
    chip = chip[ind3]
    meas = meas[ind4]
    gdfwhm, = np.where((meas['fwhm'] > 0.2*chip['fwhm']) & (meas['fwhm'] < 2.0*chip['fwhm']))
    if len(gdfwhm)==0:
        print('All measurements have bad FWHM values')
        return
    if len(gdfwhm) < len(meas):
        print('Removing '+str(len(meas)-len(gdfwhm))+' measurements with bad FWHM values')
        meas = meas[gdfwhm]

    plt.subplots(2,2, figsize = (12,8))
    plt.subplots_adjust(hspace = .4, wspace = .4)
    plt.suptitle(idv)
    meas["mjd"] -= min(meas["mjd"])
    mjd = (meas["mjd"])
    cenra = np.mean(meas["ra"])
    cendec = np.mean(meas["dec"])
    mag = meas["mag_auto"]
    filters = meas["filter"]
    #colors = ["b", "g", "r", "c", "y"]
    colordict = {'u':'c', 'g':'b', 'r':'g', 'i':'y', 'z':'orange', 'Y':'r', 'VR':'purple'}
    ra = meas["ra"]
    raerr = meas['raerr']
    dra = ra - cenra
    dra *= 3600 * np.cos(np.deg2rad(cendec))
    dec = meas["dec"]
    decerr = meas['decerr']
    ddec = dec - cendec
    ddec *= 3600
    t = meas["mjd"] - np.mean(meas["mjd"])
    pmra = obj["pmra"]/1000/365.2425 # mas/yr -> as/day
    pmdec = obj["pmdec"]/1000/365.2425

    size = 15
        
    # filter some
    goodind = np.where(np.logical_and(abs(ddec) < 500, abs(dra) < 500))             
    ddec = ddec[goodind]
    dra = dra[goodind]
    raerr = raerr[goodind]
    decerr = decerr[goodind]
    mjd = mjd[goodind]
    filters = filters[goodind]
    mag = mag[goodind]
    meandra = np.mean(dra)
    meanddec = np.mean(ddec)

    # Unique filters
    # put them in this order [u,g,r,i,z,Y,VR]
    ufilter = []
    for filt in ['u','g','r','i','z','Y','VR']:
        if filt in filters:
            ufilter.append(filt)
    
    # ra dec plot
    plt.subplot(2,2,1)
    plt.errorbar(dra, ddec, raerr, decerr, fmt='none', ecolor='lightgray', elinewidth=1, capsize=0, alpha=0.5, zorder=0)
    plt.scatter(dra, ddec, c = mjd, s = size, zorder=1)
    diffra = max(dra) - min(dra)
    diffdec = max(ddec) - min(ddec)
    plt.xlim(min(dra) - diffra/4, max(dra) + diffra/4)
    plt.ylim(min(ddec) - diffdec/4, max(ddec) + diffdec/4)
    m, b = np.polyfit(dra, ddec, 1)
    plt.plot(dra, m*dra + b, c = "k")
    plt.colorbar(label=r'$\Delta$ MJD')
    plt.xlabel("$\Delta$ RA (arcsec)")
    plt.ylabel("$\Delta$ DEC (arcsec)")

    # ra mjd plot
    plt.subplot(2,2,2)
    plt.errorbar(mjd, dra, yerr=raerr, fmt='none', ecolor='lightgray', elinewidth=1, capsize=0, alpha=0.5, zorder=0)
    count = 0
    for fil in ufilter:
        filind = np.where(filters == fil)
        #plt.scatter(mjd[filind], dra[filind], c = colors[count], label = fil, s = size, zorder=1)
        plt.scatter(mjd[filind], dra[filind], c = colordict[fil], label = fil, s = size, zorder=1)
        count+=1
    plt.legend()
    diffmjd = max(mjd) - min(mjd)
    plt.xlim(min(mjd) - diffmjd/4,  max(mjd) + diffmjd/4)
    plt.ylim(min(dra) - diffra/4, max(dra) + diffra/4)
    m, b = np.polyfit(mjd, dra, 1)        
    plt.plot(mjd, mjd*pmra+b, c = "k")
    plt.xlabel(r'$\Delta$ MJD (days)')
    plt.ylabel(r'$\Delta$ RA (arcsec)')

    # dec mjd plot
    plt.subplot(2,2,4)
    plt.errorbar(mjd, ddec, yerr=decerr, fmt='none', ecolor='lightgray', elinewidth=1, capsize=0, alpha=0.5, zorder=0)
    count = 0
    for fil in ufilter:
        filind = np.where(filters == fil)
        plt.scatter(mjd[filind], ddec[filind], c = colordict[fil], label = fil, s = size, zorder=1)
        count+=1
    plt.legend()
    plt.xlim(min(mjd) - diffmjd/4,  max(mjd) + diffmjd/4)
    plt.ylim(min(ddec) - diffdec/4, max(ddec) + diffdec/4)
    m, b = np.polyfit(mjd, ddec, 1)
    plt.plot(mjd, mjd*pmdec+b, c = "k")
    plt.xlabel(r'$\Delta$ MJD (days)')
    plt.ylabel(r'$\Delta$ DEC (arcsec)')

    #magtime
    plt.subplot(2,2,3)
    count = 0
    for fil in ufilter:
        filind = np.where(filters == fil)
        plt.scatter(mjd[filind], mag[filind], c = colordict[fil], label = fil, s = size)
        count+=1
    plt.legend()
    diffmag = max(mag) - min(mag)
    plt.xlim(min(mjd) - diffmjd/4,  max(mjd) + diffmjd/4)
    plt.ylim(min(mag) - diffmag/4, max(mag) + diffmag/4)
    plt.xlabel(r'$\Delta$ MJD (days)')
    plt.ylabel("Magnitude")

    outfile = outdir+idv+'.png'
    #outfile = outdir+idv+'_4panel.png'
    print('Saving figure to '+outfile)
    plt.savefig(outfile,bbox_inches='tight')
    #plt.close(fig)

        #parallax
    #     plt.figure()
    #     pars, cov = parallax.fit(meas)
    #     parallax.plotfit(meas,pars, cov)
    #except:
    #    print("This field threw an error ", idv)


# Main command-line program
if __name__ == "__main__":
    parser = ArgumentParser(description='Make diagnostic HPM plots')
    parser.add_argument('catfile', type=str, nargs=1, help='Catalog file')
    parser.add_argument('--outdir', type=str, nargs=1, default='',help='Output directory')
    args = parser.parse_args()
    catfile = args.catfile[0]
    outdir = dln.first_el(args.outdir)
    if outdir=='':
        outdir = '/net/dl2/dnidever/nsc/instcal/v3/hpm2/4panel/'

    cat = fits.getdata(catfile,1)
    ncat = len(cat)

    for i in range(ncat):
        print(str(i+1)+' '+cat['id'][i])
        try:
            make_4panelplot(cat[i],outdir)
        except Exception as e:
            print('Failed on '+cat['id'][i]+' '+str(e))
            traceback.print_exc()
