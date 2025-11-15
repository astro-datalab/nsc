#!/usr/bin/env python

import os
import sys
import numpy as np
from dlnpyutils import utils as dln
import time
import sqlite3

#def writecat2db(cat,dbfile,tablename):
#    """ Write a catalog to the database """
#    ncat = dln.size(cat)
#    sqlite3.register_adapter(np.int16, int)
#    sqlite3.register_adapter(np.int64, int)
#    sqlite3.register_adapter(np.float64, float)
#    sqlite3.register_adapter(np.float32, float)
#    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
#    #db = sqlite3.connect('test.db')
#    #db.text_factory = lambda x: str(x, 'latin1')
#    #db.row_factory = sqlite3.Row
#    c = db.cursor()
#    # Create the table
#    #   the primary key ROWID is automatically generated
#    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="meas"').fetchall()) < 1:
#        c.execute('''CREATE TABLE meas(measid TEXT, objlabel INTEGER, exposure TEXT, ccdnum INTEGER, filter TEXT, mjd REAL,
#                     ra REAL, raerr REAL, dec REAL, decerr REAL, mag_auto REAL, magerr_auto REAL, asemi REAL, asemierr REAL,
#                     bsemi REAL, bsemierr REAL, theta REAL, thetaerr REAL, fwhm REAL, flags INTEGER, class_star REAL)''')
#    data = list(zip(cat['measid'],np.zeros(ncat,int)-1,cat['exposure'],cat['ccdnum'],cat['filter'],cat['mjd'],cat['ra'],
#                    cat['raerr'],cat['dec'],cat['decerr'],cat['mag_auto'],cat['magerr_auto'],cat['asemi'],cat['asemierr']3,
#                    cat['bsemi'],cat['bsemierr'],cat['theta'],cat['thetaerr'],cat['fwhm'],cat['flags'],cat['class_star']))
#    c.executemany('''INSERT INTO meas(measid,objlabel,exposure,ccdnum,filter,mjd,ra,raerr,dec,decerr,mag_auto,magerr_auto,
#                     asemi,asemierr,bsemi,bsemierr,theta,thetaerr,fwhm,flags,class_star)
#                     VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', data)
#    db.commit()
#    db.close()


def writecat(cat,dbfile,table='meas'):
    """ Write a catalog to the database """
    ncat = dln.size(cat)
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"S":"TEXT", "i":"INTEGER", "f":"REAL"}

    # Get the column names
    cnames = cat.dtype.names
    cdict = dict(cat.dtype.fields)
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="'+table+'"').fetchall()) < 1:
        columns = cnames[0].lower()+' '+d2d[cdict[cnames[0]][0].kind]
        for n in cnames[1:]: columns+=', '+n.lower()+' '+d2d[cdict[n][0].kind]
        c.execute('CREATE TABLE '+table+'('+columns+')')
    # Insert statement
    columns = []
    for n in cnames: columns.append(n.lower())
    qmarks = np.repeat('?',dln.size(cnames))
    c.executemany('INSERT INTO '+table+'('+','.join(columns)+') VALUES('+','.join(qmarks)+')', list(cat))
    db.commit()
    db.close()


def getdbcoords(dbfile):
    """ Get the coordinates and ROWID from the database """
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    c.execute('''SELECT rowid,ra,dec FROM meas''')
    data = c.fetchall()
    db.close()

    # Convert to nump structured array
    dtype = np.dtype([('ROWID',int),('RA',np.float64),('DEC',np.float64)])
    cat = np.zeros(len(data),dtype=dtype)
    cat[...] = data
    del data

    return cat

def createindex(dbfile,col='measid',table='meas',unique=True,verbose=False):
    """ Index a column in the database """
    t0 = time.time()
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    index_name = 'idx_'+col+'_'+table
    # Check if the index exists first
    c.execute('select name from sqlite_master')
    d = c.fetchall()
    for nn in d:
        if nn[0]==index_name:
            print(index_name+' already exists')
            return
    # Create the index
    print('Indexing '+col)
    if unique:
        c.execute('CREATE UNIQUE INDEX '+index_name+' ON '+table+'('+col+')')
    else:
        c.execute('CREATE INDEX '+index_name+' ON '+table+'('+col+')')
    data = c.fetchall()
    db.close()
    if verbose:
        print('indexing done after '+str(time.time()-t0)+' sec')

def insertobjlabelsdb(rowid,labels,dbfile):
    """ Insert objectlabel values into the database """
    print('Inserting object labels')
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    data = list(zip(labels,rowid))
    c.executemany('''UPDATE meas SET objlabel=? WHERE rowid=?''', data) 
    db.commit() 
    db.close()
    print('inserting done after '+str(time.time()-t0)+' sec')

def updatecoldb(selcolname,selcoldata,updcolname,updcoldata,table,dbfile):
    """ Update column in database """
    print('Updating '+updcolname+' column in '+table+' table using '+selcolname)
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    data = list(zip(updcoldata,selcoldata))
    c.executemany('''UPDATE '''+table+''' SET '''+updcolname+'''=? WHERE '''+selcolname+'''=?''', data) 
    db.commit() 
    db.close()
    print('updating done after '+str(time.time()-t0)+' sec')    

def deleterowsdb(colname,coldata,table,dbfile):
    """ Delete rows from the database using rowid"""
    print('Deleting rows from '+table+' table using '+colname)
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    data = list(zip(coldata))
    c.executemany('''DELETE from '''+table+''' WHERE '''+colname+'''=?''', data) 
    db.commit() 
    db.close()
    print('deleting done after '+str(time.time()-t0)+' sec')

    
def writeidstr2db(cat,dbfile):
    """ Insert IDSTR database values """
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="idstr"').fetchall()) < 1:
        c.execute('''CREATE TABLE idstr(measid TEXT, exposure TEXT, objectid TEXT, objectindex INTEGER)''')
    data = list(zip(cat['measid'],cat['exposure'],cat['objectid'],cat['objectindex']))
    c.executemany('''INSERT INTO idstr(measid,exposure,objectid,objectindex)
                     VALUES(?,?,?,?)''', data)
    db.commit() 
    db.close()
    #print('inserting done after '+str(time.time()-t0)+' sec')

def querydb(dbfile,table='meas',cols='rowid,*',where=None):
    """ Query database table """
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()
    cmd = 'SELECT '+cols+' FROM '+table
    if where is not None: cmd += ' WHERE '+where
    cur.execute(cmd)
    data = cur.fetchall()
    db.close()

    # Return results
    return data

def executedb(dbfile,cmd):
    """ Execute a database command """
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()
    cur.execute(cmd)
    data = cur.fetchall()
    db.close()
    return data    

def getdatadb(dbfile,table='meas',cols='rowid,*',objlabel=None,rar=None,decr=None,verbose=False):
    """ Get measurements for an object(s) from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()
    cmd = 'SELECT '+cols+' FROM '+table
    # OBJLABEL constraints
    if objlabel is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        if dln.size(objlabel)==2:
            cmd += 'objlabel>='+str(objlabel[0])+' AND objlabel<='+str(objlabel[1])
        else:
            cmd += 'objlabel='+str(objlabel)
    # RA constraints
    if rar is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'ra>='+str(rar[0])+' AND ra<'+str(rar[1])
    # DEC constraints
    if decr is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'dec>='+str(decr[0])+' AND dec<'+str(decr[1])

    # Execute the select command
    #print('CMD = '+cmd)
    cur.execute(cmd)
    data = cur.fetchall()
    db.close()

    # No results
    if len(data)==0:
        return np.array([])

    # Convert to numpy structured array
    dtype_hicat = np.dtype([('ROWID',int),('MEASID',np.str,30),('OBJLABEL',int),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
                            ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
                            ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
                            ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])
    cat = np.zeros(len(data),dtype=dtype_hicat)
    cat[...] = data
    del data

    if verbose: print('got data in '+str(time.time()-t0)+' sec.')

    return cat

def getradecrangedb(dbfile):
    """ Get RA/DEC ranges from database """
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    c.execute('''SELECT MIN(ra),MAX(ra),MIN(dec),MAX(dec) FROM meas''')
    data = c.fetchall()
    db.close()

    return data[0]

def add_elements(cat,nnew=300000):
    """ Add more elements to a catalog"""
    ncat = len(cat)
    old = cat.copy()
    nnew = dln.gt(nnew,ncat)
    cat = np.zeros(ncat+nnew,dtype=old.dtype)
    cat[0:ncat] = old
    del old
    return cat    

