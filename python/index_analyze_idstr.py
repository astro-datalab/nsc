#!/usr/bin/env python

# Index and analyze idstr database

import os
import sys
import numpy as np
import time
from dlnpyutils import utils as dln, db
import sqlite3
import socket
from argparse import ArgumentParser


def index_analyze_idstr(dbfile):
    """ Index and analyze idstr database."""

    # Make sure it's a list
    if type(dbfile) is str: dbfile=[dbfile]

    print('Indexing and analyzing '+str(len(dbfile))+' database files')

    # Loop over files
    for i,dbfile1 in enumerate(dbfile):
        print(str(i)+' '+dbfile1)
        if os.path.exists(dbfile1):
            # Get existing index names for this database
            d = sqlite3.connect(dbfile1, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
            cur = d.cursor()
            cur.execute('select name from sqlite_master where tbl_name="idstr" and type="index"')
            index_names = cur.fetchall()
            if len(index_names)>0:
                index_names = [i[0] for i in index_names]

            # Create objectid and exposure index, then analyze
            if 'idx_objectid_idstr' not in index_names:
                print('Creating OBJECTID index')
                db.createindex(dbfile1,'objectid',table='idstr',unique=False)
            if 'idx_exposure_idstr' not in index_names:
                print('Creating EXPOSURE index')
                db.createindex(dbfile1,'exposure',table='idstr',unique=False)
            print('Analyzing')
            db.analyzetable(dbfile1,'idstr')

if __name__ == "__main__":
    parser = ArgumentParser(description='Index and analyze idstr database.')
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

    index_analyze_idstr(dbfile)
