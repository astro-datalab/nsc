#!/usr/bin/env python

# Update measurement catalogs using the broken up measid/objectid lists
# from nsc_instcal_combine_breakup_idstr.py

import os
import sys
import numpy as np
import time
from dlnpyutils import utils as dln, db
from astropy.table import Table
from astropy.io import fits
import sqlite3
import socket
from argparse import ArgumentParser
import logging
from glob import glob
import subprocess
import healpy as hp

def check_meas_objectid(measfiles):
    """ Check that the meas objectids are in the correct format."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    iddir = '/data0/dnidever/nsc/instcal/v3/idstr/'
    version = 'v3'

    if type(measfiles) is not list: measfiles=[measfiles]

    dt = np.dtype([('measfile',np.str,300),('exists',bool),('okay',bool),('ndot',int)])
    out = np.zeros(len(measfiles),dtype=dt)
    out['measfile'] = measfiles
    out['exists'] = False
    out['okay'] = False
    out['ndot'] = -1

    # Loop over files
    for i in range(len(measfiles)):
        print(str(i+1)+' '+measfiles[i])
        out['exists'][i] = os.path.exists(measfiles[i])
        if out['exists'][i]==False:
            print(measfiles[i]+' NOT FOUND')
            continue
        try:
            meas = fits.getdata(measfiles[i],1)
            num = np.max(meas['objectid'].count('.'))
            out['ndot'][i] = num
            if num==1:
                out['okay'][i] = True
                print('OKAY')
            else:
                out['okay'][i] = False
                print('PROBLEM')
        except:
            print('problem with '+measfiles[i])

    print('dt = %6.1f sec.' % (time.time()-t00))

    return out


if __name__ == "__main__":
    parser = ArgumentParser(description='Update measid in exposure.')
    parser.add_argument('measfiles', type=str, nargs=1, help='meas filenames')
    parser.add_argument('outfile', type=str, nargs=1, help='output file')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    measfiles = args.measfiles[0]
    outfile = args.outfile[0]

    # Input is a list
    if measfiles[0]=='@':
        listfile = measfiles[1:]
        if os.path.exists(listfile): 
            measfiles = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    # Check the measurement files
    out = check_meas_objectid(measfiles)

    # Write to output file
    print('Writing to '+outfile)
    if os.path.exists(outfile): os.remove(outfile)
    Table(out).write(outfile,overwrite=True)
