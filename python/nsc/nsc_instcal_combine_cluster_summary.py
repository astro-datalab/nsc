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
from dlnpyutils import utils as dln, coords, bindata
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

def writecat2db(cat,dbfile):
    """ Write a catalog to the database """
    ncat = dln.size(cat)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    #db = sqlite3.connect('test.db')
    #db.text_factory = lambda x: str(x, 'latin1')
    #db.row_factory = sqlite3.Row
    c = db.cursor()
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="meas"').fetchall()) < 1:
        c.execute('''CREATE TABLE meas(measid TEXT, objlabel INTEGER, exposure TEXT, ccdnum INTEGER, filter TEXT, mjd REAL,
                     ra REAL, raerr REAL, dec REAL, decerr REAL, mag_auto REAL, magerr_auto REAL, asemi REAL, asemierr REAL,
                     bsemi REAL, bsemierr REAL, theta REAL, thetaerr REAL, fwhm REAL, flags INTEGER, class_star REAL)''')
    data = list(zip(cat['measid'],np.zeros(ncat,int)-1,cat['exposure'],cat['ccdnum'],cat['filter'],cat['mjd'],cat['ra'],
                    cat['raerr'],cat['dec'],cat['decerr'],cat['mag_auto'],cat['magerr_auto'],cat['asemi'],cat['asemierr'],
                    cat['bsemi'],cat['bsemierr'],cat['theta'],cat['thetaerr'],cat['fwhm'],cat['flags'],cat['class_star']))
    c.executemany('''INSERT INTO meas(measid,objlabel,exposure,ccdnum,filter,mjd,ra,raerr,dec,decerr,mag_auto,magerr_auto,
                     asemi,asemierr,bsemi,bsemierr,theta,thetaerr,fwhm,flags,class_star)
                     VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', data)
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

def createindexdb(dbfile,col='measid',table='meas',unique=True):
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

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix pixel number')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    radeg = np.float64(180.00) / np.pi

    # Inputs
    pix = int(args.pix[0])
    version = args.version[0]

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        dir = "/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        dir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    outdir=dir+'combine/'
    subdir = str(int(pix)//1000)    # use the thousands to create subdirectory grouping
    
    # IDSTR schema
    dtype_idstr = np.dtype([('measid',np.str,200),('exposure',np.str,200),('objectid',np.str,200),('objectindex',int)])

    # OBJ schema
    dtype_obj = np.dtype([('objectid',np.str,100),('pix',int),('ra',np.float64),('dec',np.float64),('raerr',np.float32),('decerr',np.float32),
                          ('pmra',np.float32),('pmdec',np.float32),('pmraerr',np.float32),('pmdecerr',np.float32),('mjd',np.float64),
                          ('deltamjd',np.float32),('ndet',np.int16),('nphot',np.int16),
                          ('ndetu',np.int16),('nphotu',np.int16),('umag',np.float32),('urms',np.float32),('uerr',np.float32),
                             ('uasemi',np.float32),('ubsemi',np.float32),('utheta',np.float32),
                          ('ndetg',np.int16),('nphotg',np.int16),('gmag',np.float32),('grms',np.float32),('gerr',np.float32),
                             ('gasemi',np.float32),('gbsemi',np.float32),('gtheta',np.float32),
                          ('ndetr',np.int16),('nphotr',np.int16),('rmag',np.float32),('rrms',np.float32),('rerr',np.float32),
                             ('rasemi',np.float32),('rbsemi',np.float32),('rtheta',np.float32),
                          ('ndeti',np.int16),('nphoti',np.int16),('imag',np.float32),('irms',np.float32),('ierr',np.float32),
                             ('iasemi',np.float32),('ibsemi',np.float32),('itheta',np.float32),
                          ('ndetz',np.int16),('nphotz',np.int16),('zmag',np.float32),('zrms',np.float32),('zerr',np.float32),
                             ('zasemi',np.float32),('zbsemi',np.float32),('ztheta',np.float32),
                          ('ndety',np.int16),('nphoty',np.int16),('ymag',np.float32),('yrms',np.float32),('yerr',np.float32),
                             ('yasemi',np.float32),('ybsemi',np.float32),('ytheta',np.float32),
                          ('ndetvr',np.int16),('nphotvr',np.int16),('vrmag',np.float32),('vrrms',np.float32),('vrerr',np.float32),
                            ('vrasemi',np.float32),('vrbsemi',np.float32),('vrtheta',np.float32),
                          ('asemi',np.float32),('asemierr',np.float32),('bsemi',np.float32),('bsemierr',np.float32),
                          ('theta',np.float32),('thetaerr',np.float32),('fwhm',np.float32),('flags',np.int16),('class_star',np.float32),
                          ('ebv',np.float32),('rmsvar',np.float32),('madvar',np.float32),('iqrvar',np.float32),('etavar',np.float32),
                          ('jvar',np.float32),('kvar',np.float32),('chivar',np.float32),('romsvar',np.float32),
                          ('variable10sig',np.int16),('nsigvar',np.float32)])

    dbfile = tmproot+str(pix)+'_combine.db'

    # IDSTR database file
    dbfile_idstr = outdir+'/'+subdir+'/'+str(pix)+'_idstr.db'
    usedb = True

    # Load the object structured array
    obj = fits.getdata(outdir+'/'+subdir+'/'+str(pix)+'.fits.gz',2)
    nobj = len(obj)

    # Initialize the OBJ structured array
    #obj = np.zeros(nobj,dtype=dtype_obj)
    #obj['objectid'] = dln.strjoin( str(pix)+'.', ((np.arange(nobj)+1).astype(np.str)) )
    #obj['pix'] = pix
    ## all bad to start
    #for f in ['pmra','pmraerr','pmdec','pmdecerr','asemi','bsemi','theta','asemierr',
    #          'bsemierr','thetaerr','fwhm','class_star','rmsvar','madvar','iqrvar',
    #          'etavar','jvar','kvar','chivar','romsvar']: obj[f]=np.nan
    #for f in ['u','g','r','i','z','y','vr']:
    #    obj[f+'mag'] = 99.99
    #    obj[f+'err'] = 9.99
    #    obj[f+'rms'] = np.nan
    #    obj[f+'asemi'] = np.nan
    #    obj[f+'bsemi'] = np.nan
    #    obj[f+'theta'] = np.nan
    #obj['variable10sig'] = 0
    #obj['nsigvar'] = np.nan
    #idstr = np.zeros(ncat,dtype=dtype_idstr)

    ## Higher precision catalog
    #dtype_hicat = np.dtype([('MEASID',np.str,30),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
    #                        ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
    #                        ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
    #                        ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])

    ## Convert to nump structured array
    #dtype_hicatdb = np.dtype([('MEASID',np.str,30),('OBJLABEL',int),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
    #                          ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
    #                          ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
    #                          ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])

    dtype_sumstr = np.dtype([('objectid',np.str,100),('maxdist',float)])
    sumstr = np.zeros(nobj,dtype=dtype_sumstr)

    t1 = time.time()

    import pdb; pdb.set_trace()

    # Loop over the objects
    meascount = 0
    ngroup = -1
    grpcount = 0
    maxmeasload = 50000
    ngrpcat = 0
    ncat1 = 0
    fidmag = np.zeros(nobj,float)+np.nan  # fiducial magnitude
    for i,objid in enumerate(obj['objectid']):
        if (i % 1000)==0: print(i)

        # Get meas data for this object
        if usedb is False:
            oindx = np.arange(objstr['LO'][i],objstr['HI'][i]+1)  # this fails if start,stop are the same
            if objstr['NMEAS'][i]==1: oindx=np.atleast_1d(objstr['LO'][i])
            ncat1 = dln.size(oindx)
            cat1_orig = cat[oindx]
            # Upgrade precisions of catalog
            cat1 = np.zeros(ncat1,dtype=dtype_hicat)
            cat1[...] = cat1_orig   # stuff in the data
            #for n in dtype_hicat.names: cat1[n] = cat1_orig[n]
            del cat1_orig
        # Get from the database
        else:            
            # Get next group of object measurements
            if grpcount>=ngroup:
                # Use maxmeasload to figure out how many objects we can load
                if i==0:
                    ngroup = np.max(np.where(meascumcount[i:]<=maxmeasload)[0])+1
                else:
                    ngroup = np.max(np.where((meascumcount[i:]-meascumcount[i-1])<=maxmeasload)[0])+1
                ngroup = np.max([1,ngroup])   # need to load at least 1
                objid0 = objid
                objid1 = obj['objectid'][np.min([i+ngroup-1,nobj-1])]
                #lab1 = labelindex['value'][np.min([i+ngroup-1,nobj-1])]
                if ngrpcat>0: del grpcat
                if ncat1>0: del cat1
                grpcat = getdatadb(dbfile,objectid=[objid0,objid1])
                ngrpcat = dln.size(grpcat)
                grpindex = dln.create_index(grpcat['OBJLABEL'])                
                #ngroup = len(grpindex['value'])
                grpcount = 0
            # Get the measurement data for this object
            gindx = grpindex['index'][grpindex['lo'][grpcount]:grpindex['hi'][grpcount]+1]
            cat1 = np.atleast_1d(grpcat[gindx])
            ncat1 = len(cat1)
            grpcount += 1
            oindx = np.arange(ncat1)+meascount
            meascount += ncat1            

        import pdb; pdb.set_trace()

        # Compute maximum spherical distance of the measurements from the mean value        
        sumstr['objectid'][i] = obj['objectid'][i]
        sumstr['maxdist'][i] = np.max(dist)
        #obj['ndet'][i] = ncat1


    # Write the output file
    outfile = outdir+'/'+subdir+'/'+str(pix)+'_summary.fits'
    print('Writing combined catalog to '+outfile)
    Table(sumstr).write(outfile)               # first, summary table

    dt = time.time()-t0
    print('dt = '+str(dt)+' sec.')
