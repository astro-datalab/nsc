#!/usr/bin/env python

# Break up idstr file into separate measid/objectid lists per exposure on /data0

import os
import sys
import numpy as np
import time
from dlnpyutils import utils as dln, db
from astropy.io import fits
import sqlite3
import socket
from argparse import ArgumentParser


def breakup_idstr(dbfile):
    """ Break-up idstr file into separate measid/objectid lists per exposure on /data0."""

    t00 = time.time()

    outdir = '/data0/dnidever/nsc/instcal/v3/idstr/'

    # Load the exposures table
    expcat = fits.getdata('/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_v3_exposure_table.fits.gz',1)

    # Make sure it's a list
    if type(dbfile) is str: dbfile=[dbfile]

    print('Breaking up '+str(len(dbfile))+' database files')

    # Loop over files
    for i,dbfile1 in enumerate(dbfile):
        print(str(i+1)+' '+dbfile1)
        if os.path.exists(dbfile1):
            t0 = time.time()
            dbbase1 = os.path.basename(dbfile1)[0:-9]  # remove _idstr.db ending
            # Get existing index names for this database
            d = sqlite3.connect(dbfile1, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
            cur = d.cursor()
            cmd = 'select measid,exposure,objectid from idstr'
            t1 = time.time()
            data = cur.execute(cmd).fetchall()
            print('  '+str(len(data))+' rows read in %5.1f sec. ' % (time.time()-t1))
            # Break up data into lists
            measid,exposure,objectid = list(zip(*data))
            measid = np.array(measid)
            objectid = np.array(objectid)
            exposure = np.array(exposure)
            eindex = dln.create_index(exposure)
            # Match exposures to exposure catalog
            ind1,ind2 = dln.match(expcat['EXPOSURE'],eindex['value'])
            # Loop over exposures and write output files
            nexp = len(eindex['value'])
            print('  '+str(nexp)+' exposures')
            measid_maxlen = np.max(dln.strlen(measid))
            objectid_maxlen = np.max(dln.strlen(objectid))
            df = np.dtype([('measid',np.str,measid_maxlen+1),('objectid',np.str,objectid_maxlen+1)])
            # Loop over the exposures and write out the files
            for k in range(nexp):
                if nexp>100:
                    if k % 100 == 0: print('  '+str(k+1))
                ind = eindex['index'][eindex['lo'][k]:eindex['hi'][k]+1]
                cat = np.zeros(len(ind),dtype=df)
                cat['measid'] = measid[ind]
                cat['objectid'] = objectid[ind]
                instcode = expcat['INSTRUMENT'][ind1[k]]
                dateobs = expcat['DATEOBS'][ind1[k]]
                night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
                if os.path.exists(outdir+instcode+'/'+night+'/'+eindex['value'][k]) is False:
                    # Sometimes this crashes because another process is making the directory at the same time
                    try:
                        os.makedirs(outdir+instcode+'/'+night+'/'+eindex['value'][k])
                    except:
                        pass
                outfile = outdir+instcode+'/'+night+'/'+eindex['value'][k]+'/'+eindex['value'][k]+'__'+dbbase1+'.npy'
                np.save(outfile,cat)
            print('  dt = %6.1f sec. ' % (time.time()-t0))
        else:
            print('  '+dbfile1+' NOT FOUND')

    print('dt = %6.1f sec.' % (time.time()-t00))

if __name__ == "__main__":
    parser = ArgumentParser(description='Break up idstr into separate lists per exposure.')
    parser.add_argument('dbfile', type=str, nargs=1, help='Database filename')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    dbfile = args.dbfile[0]

    # Input is a list
    if dbfile[0]=='@':
        listfile = dbfile[1:]
        if os.path.exists(listfile): 
            dbfile = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    breakup_idstr(dbfile)
