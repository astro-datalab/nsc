#!/usr/bin/env python                                                                                                                                                              

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
from astropy.time import Time
import healpy as hp
from dlnpyutils import utils as dln, coords, bindata, db, job_daemon as jd
import subprocess
import time
from argparse import ArgumentParser
import socket
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
from sklearn.cluster import DBSCAN
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
import sqlite3
import gc
import psutil


def make_meas_missing(exposure):
    """ Make meas_missing.fits for exposures."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    version = 'v3'

    # Load the exposures table
    print('Loading exposure table')
    expstr = fits.getdata('/net/dl2/dnidever/nsc/instcal/'+version+'/lists/nsc_v3_exposure_table.fits.gz',1)

    # Make sure it's a list
    if type(exposure) is str: exposure=[exposure]

    # Match exposures to exposure catalog
    ind1,ind2 = dln.match(expstr['EXPOSURE'],exposure)
    nmatch = len(ind1)
    print(str(nmatch)+' exposures found in exposure table')

    # Loop over files
    for i in range(nmatch):
        exp1 = expstr['EXPOSURE'][ind1[i]]
        instcode = expstr['INSTRUMENT'][ind1[i]]
        dateobs = expstr['DATEOBS'][ind1[i]]
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        measfile = '/net/dl2/dnidever/nsc/instcal/v3/'+instcode+'/'+night+'/'+exp1+'/'+exp1+'_meas.fits.gz'
        if os.path.exists(measfile):
            meas = fits.getdata(measfile,1)
            objectid = np.char.array(meas['OBJECTID']).strip()
            miss, = np.where(objectid=='')
            if len(miss)>0:
                print(str(i+1)+' '+exp1+' '+str(len(miss))+' missing')
                missfile = '/net/dl2/dnidever/nsc/instcal/v3/'+instcode+'/'+night+'/'+exp1+'/'+exp1+'_meas_missing.fits'
                print('  Writing missing catalog to '+missfile)
                Table(meas[miss]).write(missfile,overwrite=True)
            else:
                print(str(i+1)+' '+exp1+' none missing')
        else:
            print(measfile+' NOT FOUND')

 
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('exposure', type=str, nargs=1, help='exposure')
    args = parser.parse_args()

    version = 'v3'

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    exposure = args.exposure[0]

    # Input is a list
    if exposure[0]=='@':
        listfile = exposure[1:]
        if os.path.exists(listfile): 
            exposure = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    make_meas_missing(exposure)
