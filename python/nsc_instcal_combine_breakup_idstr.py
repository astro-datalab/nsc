#!/usr/bin/env python

# Break up idstr file into separate measid/objectid lists per exposure on /data0

import os
import sys
import numpy as np
import time
from dlnpyutils import utils as dln, db
import sqlite3
import socket
from argparse import ArgumentParser


def breakup_idstr(dbfile):
    """ Break-up idstr file into separate measid/objectid lists per exposure on /data0."""

    t0 = time.time()

    # Make sure it's a list
    if type(dbfile) is str: dbfile=[dbfile]

    print('Breaking up '+str(len(dbfile))+' database files')

    # Loop over files
    for i,dbfile1 in enumerate(dbfile):
        print(str(i)+' '+dbfile1)
        if os.path.exists(dbfile1):
            dbbase1 = os.path.basename(dbfile1)[0:-9]  # remove _idstr.db ending
            # Get existing index names for this database
            d = sqlite3.connect(dbfile1, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
            cur = d.cursor()
            # Get number of rows
            cur.execute('select count(*) from idstr')
            nmeas = cur.fetchall()
            nmeas = nmeas[0][0]
            # Number of rows per exposure
            #  on larger db files this can take ~20 sec
            cur.execute('select exposure,count(measid) from idstr group by exposure')
            data = cur.fetchall()
            nexp = len(data)
            dt = np.dtype([('exposure',np.str,50),('nmeas',int)])
            expcat = np.zeros(nexp,dtype=dt)
            expcat[...] = data
            # Get minimum and maximum rowIDs
            minrowid = cur.execute('select min(rowid) from idstr').fetchall()[0][0]
            maxrowid = cur.execute('select max(rowid) from idstr').fetchall()[0][0]
            # Number of groups
            ngroups = int(np.ceil(nmeas/500000))
            rowidpergroup = int(np.ceil((maxrowid-minrowid+1)/ngroups))+1
            lastmaxrowid = minrowid-1
            for i in range(ngroups):
                #thisminrowid = i*rowidpergroup+minrowid
                #thismaxrowid = (i+1)*rowidpergroup+minrowid
                thisminrowid = lastmaxrowid+1
                thismaxrowid = thisminrowid+rowidpergroup
                if i<(ngroups-1):
                    data = cur.execute('select rowid,measid,exposure,objectid from idstr where rowid>='+str(thisminrowid)+' and rowid<'+str(thismaxrowid)).fetchall()
                else:
                    # last one doesn't not upper limit
                    data = cur.execute('select rowid,measid,exposure,objectid from idstr where rowid>='+str(thisminrowid)).fetchall()
                rowid,measid,exposure,objectid = list(zip(*data))
                eindex = dln.create_index(exposure)
                # resort by when they appear
                si = np.argsort(eindex['index'][eindex['lo']])
                eindex['value'] = eindex['value'][si]
                eindex['num'] = eindex['num'][si]
                eindex['lo'] = eindex['lo'][si]
                eindex['hi'] = eindex['hi'][si]

                # Remove last exposure since it might have gotten split across the query breaks
                #  not on last one
                if i<(ngroups-1):
                    old = eindex.copy()
                    eindex = {'index':old['index'],'value':old['value'][0:-1],
                              'num':old['num'][0:-1],'lo':old['lo'][0:-1],'hi':old['hi'][0:-1]}
                    ## DOUBLE-CHECK THAT THIS WORKS!!!
                    del old

                rowid = np.array(rowid)
                lastmaxrowid = np.max(rowid[eindex['index'][eindex['hi']]])
                # Loop over exposures and write output files
                nexp = len(eindex['value'])
                df = np.dtype([('measid',np.str,100),('objectid',np.str,100)])
                for j in range(nexp):
                    ind = eindex['index'][eindex['lo'][j]:eindex['hi'][j]+1]
                    cat = np.zeros(len(ind),dtype=df)
                    cat['measid'] = res['measid'][ind]
                    cat['objectid'] = res['objectid'][ind]
                    outfile = '/data0/dnidever/nsc/instcal/v3/tmp/writetest/'+dbbase1+'_'+eindex['value'][j]+'.npy'
                    np.save(outfile,cat)

    print('dt = '+str(time.time()-t0)+' sec.')

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
