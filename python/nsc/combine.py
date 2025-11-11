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
from . import utils

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

    
def writeidtab2db(cat,dbfile):
    """ Insert IDTAB database values """
    t0 = time.time()
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="idtab"').fetchall()) < 1:
        c.execute('''CREATE TABLE idtab(measid TEXT, exposure TEXT, objectid TEXT, objectindex INTEGER)''')
    data = list(zip(cat['measid'],cat['exposure'],cat['objectid'],cat['objectindex']))
    c.executemany('''INSERT INTO idtab(measid,exposure,objectid,objectindex)
                     VALUES(?,?,?,?)''', data)
    db.commit() 
    db.close()
    #print('inserting done after '+str(time.time()-t0)+' sec')

def readidtabdb(dbfile):
    """ Get data from IDTAB database"""
    data = querydb(dbfile,table='idtab',cols='*')
    # Put in catalog
    dtype_idtab = np.dtype([('measid',str,200),('exposure',str,200),('objectid',str,200),('objectindex',int)])
    cat = np.zeros(len(data),dtype=dtype_idtab)
    cat[...] = data
    del data    
    return cat

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
    dtype_hicat = np.dtype([('ROWID',int),('MEASID',str,30),('OBJLABEL',int),('EXPOSURE',str,40),('CCDNUM',int),('FILTER',str,3),
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

def add_elements(meas,nnew=300000):
    """ Add more elements to a catalog"""
    nmeas = len(meas)
    old = meas.copy()
    nnew = dln.gt(nnew,nmeas)
    meas = np.zeros(nmeas+nnew,dtype=old.dtype)
    meas[0:nmeas] = old
    del old
    return meas    

def seqcluster(meas,dcr=0.5,iter=False,inpobj=None,trim=False):
    """ Sequential clustering of measurements in exposures.  This was the old method."""

    nmeas = len(meas)
    labels = np.zeros(nmeas)-1

    # Iterate
    if iter is not False:
        done = False
        niter = 1
        maxiter = 10
        lastlabels = np.zeros(nmeas)-1
        while (done is False):
            # First time
            if niter==1:
                inpobj = None
            # Second and up
            else:
                del labels1, obj1
                inpobj = obj2
            # Cluster
            labels1,obj1 = seqcluster(meas,dcr=dcr,iter=False,inpobj=inpobj)
            print('Iter='+str(niter)+' '+str(int(np.max(labels1)))+' clusters')
            # Calculate average ra/dec
            obj2 = propermotion(meas,labels1)
            #print(labels1-lastlabels)
            # Are we done?
            if (niter==maxiter) | (np.sum(labels1-lastlabels)==0): done=True
            lastlabels = labels1
            niter += 1
        return labels1, obj2

    # Create exposures index
    index = dln.create_index(meas['exposure'])
    nexp = len(index['value'])

    # Create object catalog
    dtype_obj = np.dtype([('label',int),('ra',np.float64),('dec',np.float64),('ndet',int)])    
    # Is there an input object catalog that we are starting with?
    if inpobj is not None:
        obj = inpobj
        cnt = len(obj)
    else:
        obj = np.zeros(np.min([500000,nmeas]),dtype=dtype_obj)
        cnt = 0
    nobj = len(obj)

    # Loop over exposures
    for i in range(nexp):
        #print(str(i)+' '+index['value'][i])
        indx = index['index'][index['lo'][i]:index['hi'][i]+1]
        meas1 = meas[indx]
        nmeas1 = len(meas1)
        if dln.size(dcr)>1:
            dcr1 = dcr[indx]
        else:
            dcr1 = dcr
        
        # First exposure
        if cnt==0:
            ind1 = np.arange(nmeas1)
            obj['label'][ind1] = ind1
            obj['ra'][ind1] = meas1['ra']
            obj['dec'][ind1] = meas1['dec']
            obj['ndet'][ind1] = 1
            labels[indx] = ind1
            cnt += nmeas1

        # Second and up
        else:
            #  Match new sources to the objects
            #ind1,ind2,dist = coords.xmatch(obj[0:cnt]['ra'],obj[0:cnt]['dec'],meas1['RA'],meas1['DEC'],dcr,unique=True)
            ind2,ind1,dist = coords.xmatch(meas1['ra'],meas1['dec'],obj[0:cnt]['ra'],obj[0:cnt]['dec'],dcr1,unique=True)            
            nmatch = dln.size(ind1)
            #  Some matches, add data to existing record for these sources
            if nmatch>0:
                obj['ndet'][ind1] += 1
                labels[indx[ind2]] = ind1
                if nmatch<nmeas1:
                    indx0 = indx.copy()
                    indx = np.delete(indx,ind2)
                    meas1 = np.delete(meas1,ind2)
                    nmeas1 = dln.size(meas1)
                else:
                    meas1 = np.array([])
                    nmeas1 = 0

            # Some left, add records for these sources
            if nmeas1>0:
                # Add new elements
                if (cnt+nmeas1)>nobj:
                    obj = add_elements(obj)
                    nobj = len(obj)
                ind1 = np.arange(nmeas1)+cnt
                obj['label'][ind1] = ind1
                obj['ra'][ind1] = meas1['ra']
                obj['dec'][ind1] = meas1['dec']
                obj['ndet'][ind1] = 1
                labels[indx] = ind1

                cnt += nmeas1
    # Trim off the excess elements
    obj = obj[0:cnt]
    # Trim off any objects that do not have any detections
    #  could happen if an object catalog was input
    if trim is True:
        bd, nbd = dln.where(obj['ndet']<1)
        if nbd>0: obj = np.delete(obj,bd)
    
    # Maybe iterate
    # -measure mean ra/dec for each object and go through the process again

    return labels, obj

def seqpms(obj):
    """
    Calculate mean coordinates and slopes.
    """

    ndet = obj['ndet']
    sumra = obj['sumra']
    sumdec = obj['sumdec']
    sumt = obj['sumt']
    sumt2 = obj['sumt2']
    sumtra = obj['sumtra']
    sumtdec = obj['sumtdec']
    
    # Calculate mean coordinates
    mnra = sumra/ndet
    mndec = sumdec/ndet
    mnt = sumt/ndet
    # Calculate proper motions
    slpra = np.zeros(nmatch,float)
    slpdec = np.zeros(nmatch,float)
    twodet = (ndet > 1)
    denom = (sumt2[twodet]/ndet[twodet]-mnt[twodet]**2)
    slpra[twodet] = (sumtra[twodet]/ndet[twodet]-mnra[twodet]*mnt[twodet]) / denom
    slpdec[twodet] = (sumtdec[twodet]/ndet[twodet]-mndec[twodet]*mnt[twodet]) / denom
    # weighted mean and slope
    # mnx = np.sum(wt*x)/totwt
    # mny = np.sum(wt*y)/totwt
    # mnxerr = np.sqrt( np.sum( ((x-mnx)**2)*wt)*n / ((n-1)*np.sum(wt))) / np.sqrt(n)
    #        ~ stddev/sqrt(n)
    #        = np.sqrt( np.sum(wt*(x**2-2*x*mnx+mnx**2))*n / ((n-1)*np.sum(wt))) / np.sqrt(n)
    #        = np.sqrt( np.sum(wt*x**2)-2*mnx*np.sum(wt*x)+np.sum(wt*mnx**2) ) *
    #                 sqrt(n/((n-1)*np.sum(wt))) / np.sqrt(n)
    #        = np.sqrt( np.sum(wt*x**2)-2*mnx*np.sum(wt*x)+np.sum(wt)*mnx**2 ) *
    #                 sqrt(1/((n-1)*np.sum(wt))
    # wtx = (np.sum(wt*x*y)/totwt-mnx*mny)/(np.sum(wt*x**2)/totwt-mnx**2)
    # wtxerr = 1.0/np.sqrt( np.sum(wt*x**2)-mnx**2 * np.sum(wt))
    wmnt = obj['sumwt']/obj['sumw']
    wmnra = obj['sumwra']/obj['sumw']
    wmndec = obj['sumwdec']/obj['sumw']
    denom = (obj['sumwt2']/obj['sumw']-wmnt**2)
    wslpra = (obj['sumwtra']/obj['sumw']-wmnt*wmnra)/denom
    wslpdec = (obj['sumwtdec']/obj['sumw']-wmnt*wmndec)/denom
    wslpraerr = 1.0/np.sqrt( obj['sumwra2']-wmnra**2 * obj['sumw'])
    wslpraerr = 1.0/np.sqrt( obj['sumwdec2']-wmndec**2 * obj['sumw'])
    wmnraerr = (np.sqrt((obj['sumwra2']-2*wmnra*obj['sumwra']+obj['sumw']*wmnra**2) / 
                        (np.sqrt(obj['ndet']-1)*obj['sumw'])))
    wmndecerr = (np.sqrt((obj['sumwdec2']-2*wmndec*obj['sumwdec']+obj['sumw']*wmndec**2) / 
                        (np.sqrt(obj['ndet']-1)*obj['sumw'])))

    return (mnt,mnra,mndec,slpra,slpdec,
            wmnt,wmnra,wmndec,wslpra,wslpdec,
            wmnraerr,wmndecerr,wslpraerr,wslpdecerr)

def seqclusterpm(meas,dcr=0.5,doiter=False,inpobj=None,calcpm=True,trim=False,minmeaspm=5):
    """
    Sequential clustering of measurements in exposures with proper motion.
    If you are rerunning with a previous object table and do NOT want to
    to recalculate the proper motions (as we go), then set calcpm=False
    """

    nmeas = len(meas)
    labels = np.zeros(nmeas)-1   # object label (also its index) for all the measurements

    # Create exposures index
    index = dln.create_index(meas['exposure'])
    nexp = len(index['value'])

    # Create object catalog
    dtype_obj = np.dtype([('label',int),('ra',np.float64),('dec',np.float64),('ndet',int),
                          ('sumt',float),('sumt2',float),('sumra',float),('sumdec',float),
                          ('sumtra',float),('sumtdec',float),('slpra',float),('slpdec',float),('mnt',float),
                          ('sumw',float),('sumwt',float),('sumwt2',float),
                          ('sumwra',float),('sumwdec',float),('sumwtra',float),('sumwtdec',float),
                          ('sumwra2',float),('sumwdec2',float)])

    # weighted slope
    # mnx = np.sum(wt*x)/totwt
    # mny = np.sum(wt*y)/totwt
    # wtx = (np.sum(wt*x*y)/totwt-mnx*mny)/(np.sum(wt*x**2)/totwt-mnx**2)

    # Is there an input object catalog that we are starting with?
    if inpobj is not None:
        obj = inpobj.copy()
        cnt = len(obj)
    else:
        obj = np.zeros(np.min([500000,nmeas]),dtype=dtype_obj)
        cnt = 0
    nobj = len(obj)

    # Loop over exposures
    for i in range(nexp):
        #print(str(i)+' '+index['value'][i])
        indx = index['index'][index['lo'][i]:index['hi'][i]+1]
        meas1 = meas[indx]
        nmeas1 = len(meas1)
        measmjd = meas1['mjd'][0]
        if dln.size(dcr)>1:
            dcr1 = dcr[indx]
        else:
            dcr1 = dcr
        
        # First exposure
        if cnt==0:
            ind1 = np.arange(nmeas1)
            obj['label'][ind1] = ind1
            obj['ra'][ind1] = meas1['ra']
            obj['dec'][ind1] = meas1['dec']
            obj['ndet'][ind1] = 1
            obj['sumt'][ind1] = meas1['mjd']
            obj['sumt2'][ind1] = meas1['mjd']
            obj['sumra'][ind1] = meas1['ra']
            obj['sumdec'][ind1] = meas1['dec']
            obj['sumra'][ind1] = meas1['ra']
            obj['sumdec'][ind1] = meas1['dec']
            obj['sumtra'][ind1] = meas1['mjd']*meas1['ra']
            obj['sumtdec'][ind1] = meas1['mjd']*meas1['dec']
            obj['slpra'][ind1] = 0.0
            obj['slpdec'][ind1] = 0.0
            obj['mnt'][ind1] = measmjd
            # weighted values
            wt = 1.0/meas1['raerr']**2  # ra/dec uncertainty are the same
            obj['sumw'][ind1] = wt
            obj['sumwt'][ind1] = wt*meas1['mjd']
            obj['sumwt2'][ind1] = wt*meas1['mjd']**2
            obj['sumwra'][ind1] = wt*meas1['ra']
            obj['sumwdec'][ind1] = wt*meas1['dec']
            obj['sumwra2'][ind1] = wt*meas1['ra']**2
            obj['sumwdec2'][ind1] = wt*meas1['dec']**2
            obj['sumwtra'][ind1] = wt*meas1['mjd']*meas1['ra']
            obj['sumwtdec'][ind1] = wt*meas1['mjd']*meas1['dec']

            labels[indx] = ind1
            cnt += nmeas1

        # Second and up, or object catalog input
        else:
            # Predict current coordinates with linear fit
            # this is helpful if there is a big temporal gap between
            # two exposures
            predra = obj['ra'][:cnt] + obj['slpra'][:cnt]*(measmjd-obj['mnt'][:cnt])
            preddec = obj['dec'][:cnt] + obj['slpdec'][:cnt]*(measmjd-obj['mnt'][:cnt])

            ind2,ind1,dist = coords.xmatch(meas1['ra'],meas1['dec'],predra,
                                           preddec,dcr1,unique=True)
            #ind2,ind1,dist = coords.xmatch(meas1['RA'],meas1['DEC'],obj[:cnt]['ra'],
            #                               obj[:cnt]['dec'],dcr1,unique=True)
            nmatch = dln.size(ind1)
            #  Some matches, add data to existing records for these measurements
            if nmatch>0:
                obj['ndet'][ind1] += 1
                obj['sumt'][ind1] += meas1['mjd'][ind2]
                obj['sumt2'][ind1] += meas1['mjd'][ind2]
                obj['sumra'][ind1] += meas1['ra'][ind2]
                obj['sumdec'][ind1] += meas1['dec'][ind2]
                obj['sumtra'][ind1] += meas1['mjd'][ind2]*meas1['ra'][ind2]
                obj['sumtdec'][ind1] += meas1['mjd'][ind2]*meas1['dec'][ind2]
                # weighted values
                wt = 1.0/meas1['RAERR'][ind2]**2  # ra/dec uncertainty are the same
                obj['sumw'][ind1] += wt
                obj['sumwt'][ind1] += wt*meas1['mjd'][ind2]
                obj['sumwt2'][ind1] += wt*meas1['mjd'][ind2]**2
                obj['sumwra'][ind1] += wt*meas1['ra'][ind2]
                obj['sumwdec'][ind1] += wt*meas1['dec'][ind2]
                obj['sumwra2'][ind1] += wt*meas1['ra'][ind2]**2
                obj['sumwdec2'][ind1] += wt*meas1['dec'][ind2]**2
                obj['sumwtra'][ind1] = wt*meas1['mjd'][ind2]*meas1['ra'][ind2]
                obj['sumwtdec'][ind1] = wt*meas1['mjd'][ind2]*meas1['dec'][ind2]
                # Calculate mean coordinates and proper motions
                if calcpm:
                    mnt,mnra,mndec,slpra,slpdec = seqpms(obj[ind1])
                    pm = np.sqrt(slpra**2+slpdec**2)
                    # Limit pm mesurements to good ones
                    gd, = np.where((obj['ndet'][ind1]>minmeaspm) &
                                   (pm/pmerr > 3))
                    if len(gd)>0:
                        obj['slpra'][ind1][gd] = slpra[gd]
                        obj['slpdec'][ind1][gd] = slpdec[gd]
                    obj['ra'][ind1] = mnra
                    obj['dec'][ind1] = mndec
                    obj['mnt'][ind1] = mnt
                labels[indx[ind2]] = ind1
                if nmatch<nmeas1:
                    indx0 = indx.copy()
                    indx = np.delete(indx,ind2)
                    meas1 = np.delete(meas1,ind2)
                    nmeas1 = dln.size(meas1)
                else:
                    meas1 = np.array([])
                    nmeas1 = 0

            # Some left, add records for these sources
            if nmeas1>0:
                # Add new elements
                if (cnt+nmeas1)>nobj:
                    obj = add_elements(obj)
                    nobj = len(obj)
                ind1 = np.arange(nmeas1)+cnt
                obj['label'][ind1] = ind1
                obj['ra'][ind1] = meas1['ra']
                obj['dec'][ind1] = meas1['dec']
                obj['sumt'][ind1] = meas1['mjd']
                obj['sumt2'][ind1] = meas1['mjd']
                obj['sumra'][ind1] = meas1['ra']
                obj['sumdec'][ind1] = meas1['dec']
                obj['sumtra'][ind1] = meas1['mjd']*meas1['ra']
                obj['sumtdec'][ind1] = meas1['mjd']*meas1['dec']
                obj['ndet'][ind1] = 1
                obj['slpra'][ind1] = 0.0
                obj['slpdec'][ind1] = 0.0
                obj['mnt'][ind1] = measmjd
                # weighted values
                wt = 1.0/meas1['raerr']**2  # ra/dec uncertainty are the same
                obj['sumw'][ind1] = wt
                obj['sumwt'][ind1] = wt*meas1['mjd']
                obj['sumwt2'][ind1] = wt*meas1['mjd']**2
                obj['sumwra'][ind1] = wt*meas1['ra']
                obj['sumwdec'][ind1] = wt*meas1['dec']
                obj['sumwra2'][ind1] = wt*meas1['ra']**2
                obj['sumwdec2'][ind1] = wt*meas1['dec']**2
                obj['sumwtra'][ind1] = wt*meas1['mjd']*meas1['ra']
                obj['sumwtdec'][ind1] = wt*meas1['mjd']*meas1['dec']

                labels[indx] = ind1

                cnt += nmeas1
    # Trim off the excess elements
    obj = obj[:cnt]
    # Trim off any objects that do not have any detections
    #   could happen if an object catalog was input
    if trim is True:
        bd, nbd = dln.where(obj['ndet']<1)
        if nbd>0: obj = np.delete(obj,bd)
    
    # Maybe iterate
    # -measure mean ra/dec for each object and go through the process again

    return labels, obj

def meancoords(meas,labels):
    """ Measure mean RA/DEC."""

    # Make object index
    index = dln.create_index(labels)
    nobj = len(index['value'])
    radeg = np.float64(180.00) / np.pi
    
    dtype_obj = np.dtype([('label',int),('ndet',int),('ra',np.float64),('dec',np.float64),('raerr',np.float32),
                          ('decerr',np.float32),('asemi',np.float32),('bsemi',np.float32),
                          ('theta',np.float32),('fwhm',np.float32)])
    obj = np.zeros(nobj,dtype=dtype_obj)

    # Loop over the objects
    for i in range(nobj):
        indx = index['index'][index['lo'][i]:index['hi'][i]+1]
        nmeas1 = dln.size(indx)
        obj['label'][i] = index['value'][i]
        obj['ndet'][i] = nmeas1
        
        # Computing quantities
        # Mean RA/DEC, RAERR/DECERR
        if nmeas1>1:
            wt_ra = 1.0/meas['raerr'][indx]**2
            wt_dec = 1.0/meas['decerr'][indx]**2
            obj['ra'][i] = np.sum(meas['ra'][indx]*wt_ra)/np.sum(wt_ra)
            obj['raerr'][i] = np.sqrt(1.0/np.sum(wt_ra))
            obj['dec'][i] = np.sum(meas['dec'][indx]*wt_dec)/np.sum(wt_dec)
            obj['decerr'][i] = np.sqrt(1.0/np.sum(wt_dec))
        else:
            obj['ra'][i] = meas['ra'][indx]
            obj['dec'][i] = meas['dec'][indx]
            obj['raerr'][i] = meas['raerr'][indx]
            obj['decerr'][i] = meas['decerr'][indx]

        # Compute median FWHM
        if nmeas1>1:
            obj['asemi'][i] = np.median(meas['asemi'][indx])
            obj['bsemi'][i] = np.median(meas['bsemi'][indx])
            obj['theta'][i] = np.median(meas['theta'][indx])
            obj['fwhm'][i] = np.median(meas['fwhm'][indx])
        else:
            obj['asemi'][i] = meas['asemi'][indx]
            obj['bsemi'][i] = meas['bsemi'][indx]
            obj['theta'][i] = meas['theta'][indx]
            obj['fwhm'][i] = meas['fwhm'][indx]

    return obj
            
    
def propermotion(meas,labels):
    """ Measure proper motions."""
    
    # Make object index
    index = dln.create_index(labels)
    nobj = len(index['value'])
    radeg = np.float64(180.00) / np.pi

    obj = meancoords(meas,labels)
    dtype_pm = np.dtype([('pmra',np.float32),('pmdec',np.float32),('pmraerr',np.float32),
                         ('pmdecerr',np.float32),('mjd',np.float64)])
    obj = dln.addcatcols(obj,dtype_pm)

    # Loop over the objects
    for i in range(nobj):
        indx = index['index'][index['lo'][i]:index['hi'][i]+1]
        nmeas1 = dln.size(indx)

        # Mean proper motion and errors
        if nmeas1>1:
            raerr = np.array(meas['raerr'][indx]*1e3,np.float64)    # milli arcsec
            ra = np.array(meas['ra'][indx],np.float64)
            ra -= np.mean(ra)
            ra *= 3600*1e3 * np.cos(obj['dec'][i]/radeg)     # convert to true angle, milli arcsec
            t = meas['mjd'][indx].copy()
            t -= np.mean(t)
            t /= 365.2425                          # convert to year
            # Calculate robust slope
            pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
            obj['pmra'][i] = pmra                 # mas/yr
            obj['pmraerr'][i] = pmraerr           # mas/yr

            decerr = np.array(meas['DECERR'][indx]*1e3,np.float64)   # milli arcsec
            dec = np.array(meas['DEC'][indx],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            # Calculate robust slope
            pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
            obj['pmdec'][i] = pmdec               # mas/yr
            obj['pmdecerr'][i] = pmdecerr         # mas/yr

    return obj

def moments(meas,labels):
    # Measure XX, YY, XY comments of multiple measurements of an object:

    # Make object index
    index = dln.create_index(labels)
    nobj = len(index['value'])
    radeg = np.float64(180.00) / np.pi

    obj = meancoords(meas,labels)
    dtype_mom = np.dtype([('x2',np.float32),('y2',np.float32),('xy',np.float32),
                          ('asemi',np.float32),('bsemi',np.float32),('theta',np.float32)])
    obj = dln.addcatcols(obj,dtype_mom)

    # Loop over the objects
    for i in range(nobj):
        indx = index['index'][index['lo'][i]:index['hi'][i]+1]
        nmeas1 = dln.size(indx)
        
        # Measure moments
        if nmeas1>1:
            # See sextractor.pdf pg. 30
            x2 = np.sum( ((meas['ra'][indx]-obj['ra'][i])*np.cos(np.deg2rad(obj['dec'][i])))**2 ) / (nmeas1-1) * 3600**2
            y2 = np.sum( (meas['DEC'][indx]-obj['dec'][i])**2 ) / (nmeas1-1) * 3600**2
            xy = np.sum( (meas['ra'][indx]-obj['ra'][i])*np.cos(np.deg2rad(obj['dec'][i])) * (meas['DEC'][indx]-obj['dec'][i]) ) / (nmeas1-1) * 3600**2
            obj['x2'][i] = x2
            obj['y2'][i] = y2
            obj['xy'][i] = xy
            # See sextractor.pdf pg. 31
            obj['asemi'][i] = np.sqrt( 0.5*(x2+y2) + np.sqrt(((x2-y2)*0.5)**2 + xy**2) )
            obj['bsemi'][i] = np.sqrt( 0.5*(x2+y2) - np.sqrt(((x2-y2)*0.5)**2 + xy**2) )
            if (x2==y2):
                obj['theta'][i] = 0.0
            else:
                obj['theta'][i] = np.rad2deg(np.arctan(2*xy/(x2-y2))*0.5)
        else:
            obj['x2'][i] = obj['raerr'][i]**2
            obj['y2'][i] = obj['decerr'][i]**2
            obj['xy'][i] = 0.0
            obj['asemi'][i] = obj['x2'][i]
            obj['bsemi'][i] = obj['y2'][i]
            obj['theta'][i] = 0.0

    return obj

def ellipsecoords(pars,npoints=100):
    """ Create coordinates of an ellipse."""
    # [x,y,asemi,bsemi,theta]
    # copied from ellipsecoords.pro
    xc = pars[0]
    yc = pars[1]
    asemi = pars[2]
    bsemi = pars[3]
    pos_ang = pars[4]
    phi = 2*np.pi*(np.arange(npoints,dtype=float)/(npoints-1))   # Divide circle into Npoints
    ang = np.deg2rad(pos_ang)                             # Position angle in radians
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    x =  asemi*np.cos(phi)                              # Parameterized equation of ellipse
    y =  bsemi*np.sin(phi)
    xprime = xc + x*cosang - y*sinang               # Rotate to desired position angle
    yprime = yc + x*sinang + y*cosang
    return xprime, yprime

def checkboundaryoverlap(metafiles,buffdict,verbose=False):
    """ Check a list of fits files against a buffer and return metadata of overlapping exposures."""

    # New meta-data format
    dtype_meta = np.dtype([('file',str,500),('base',str,200),('instrument',str,3),('expnum',int),('ra',np.float64),
                           ('dec',np.float64),('dateobs',str,100),('mjd',np.float64),('filter',str,50),
                           ('exptime',float),('airmass',float),('nsources',int),('fwhm',float),
                           ('nchips',int),('badchip31',bool),('rarms',float),('decrms',float),
                           ('ebv',float),('gaianmatch',int),('zpterm',float),('zptermerr',float),
                           ('zptermsig',float),('refmatch',int)])

    allmeta = None
    for m,mfile in enumerate(np.atleast_1d(metafiles)):
        noverlap = 0
        if os.path.exists(mfile) is False and os.path.exists(mfile+'.gz') is False:
            if verbose: print(mfile+' NOT FOUND')
            continue
        if os.path.exists(mfile):
            meta = fits.getdata(mfile,1)
        else:
            meta = fits.getdata(mfile+'.gz',1)
        if verbose: print(str(m+1)+' Loading '+mfile)
        t = Time(meta['dateobs'], format='isot', scale='utc')
        meta['mjd'] = t.mjd                    # recompute because some MJD are bad
        chmeta = fits.getdata(mfile,2)      # chip-level meta-data structure

        # Convert META to new format
        newmeta = np.zeros(1,dtype=dtype_meta)
        # Copy over the meta information
        for c in newmeta.dtype.names:
            if c in meta.dtype.names: newmeta[c]=meta[c]

        # Get the name
        fdir = os.path.dirname(mfile)
        fbase, ext = os.path.splitext(os.path.basename(mfile))
        fbase = fbase[:-5]   # remove _meta at end
        
        # Loop over the chip files
        for j in range(len(chmeta)):
            # Check that this overlaps the healpix region
            inside = True
            vra = chmeta['vra'][j]
            vdec = chmeta['vdec'][j]
            vlon, vlat = coords.rotsphcen(vra,vdec,buffdict['cenra'],buffdict['cendec'],gnomic=True)
            if coords.doPolygonsOverlap(buffdict['lon'],buffdict['lat'],vlon,vlat) is False:
                if verbose: print('This chip does NOT overlap the HEALPix region+buffer')
                inside = False
            if inside is True:
                #chfile1 = chmeta['FILENAME'][j]
                #if os.path.exists(chfile1) is True: chfiles.append(chfile1)
                noverlap += 1

        if verbose: print('  FILTER='+meta['filter'][0]+'  EXPTIME='+str(meta['exptime'][0])+' sec  '+str(noverlap)+' chips overlap')
        if noverlap>0:
            if allmeta is None:
                allmeta = newmeta
            else:
                allmeta = np.hstack((allmeta,newmeta))
    
    if allmeta is None:
        nallmeta = 0
    else:
        nallmeta = len(allmeta)
    if verbose: print(str(nallmeta)+' exposures overlap')

    return allmeta


def find_obj_parent(obj):
    """ Find objects that have other objects "inside" them. """

    # Use crossmatch

    X1 = np.vstack((obj['ra'],obj['dec'])).T
    X2 = np.vstack((obj['ra'],obj['dec'])).T
    
    X1 = X1 * (np.pi / 180.)
    X2 = X2 * (np.pi / 180.)
    max_distance = (np.max(obj['fwhm']) / 3600) * (np.pi / 180.)

    # Convert 2D RA/DEC to 3D cartesian coordinates
    Y1 = np.transpose(np.vstack([np.cos(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 0]) * np.cos(X1[:, 1]),
                                 np.sin(X1[:, 1])]))
    Y2 = np.transpose(np.vstack([np.cos(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 1])]))

    # law of cosines to compute 3D distance
    max_y = np.sqrt(2 - 2 * np.cos(max_distance))
    dist, ind = coords.crossmatch(Y1, Y2, max_y, k=2)

    # convert distances back to angles using the law of tangents
    not_inf = ~np.isinf(dist)
    x = 0.5 * dist[not_inf]
    dist[not_inf] = (180. / np.pi * 2 * np.arctan2(x,
                                  np.sqrt(np.maximum(0, 1 - x ** 2))))
    dist[not_inf] *= 3600.0      # in arcsec

    # Add "parent" column if necessary
    if 'parent' not in obj.dtype.names:
        obj = dln.addcatcols(obj,np.dtype([('parent',bool)]))
    
    # Check if there are any objects within FWHM
    #  the closest object will be itself, so check the second one
    bd,nbd = dln.where( dist[:,1] <= np.minimum(0.5*obj['fwhm'],obj['asemi']))

    # Check that they are inside their ellipse footprint
    obj['parent'] = False    # all false to start
    if nbd>0:
        for i in range(nbd):
            ind1 = bd[i]
            ind2 = ind[bd[i],1]
            lon1,lat1 = (0.0, 0.0)
            cenra = obj['ra'][ind1]
            cendec = obj['dec'][ind1]
            lon2,lat2 = coords.rotsphcen(obj['ra'][ind2],obj['dec'][ind2],cenra,cendec,gnomic=True)
            pars = [lon1*3600,lon1*3600,obj['asemi'][ind1],obj['bsemi'][ind1],obj['theta'][ind1]]
            ll,bb = ellipsecoords(pars,npoints=10)
            obj['parent'][ind1] = coords.isPointInPolygon(ll,bb,np.atleast_1d(lon2*3600),np.atleast_1d(lat2*3600))
            #obj['parent'][ind1] = coords.doPolygonsOverlap(ll,bb,np.atleast_1d(lon2*3600),np.atleast_1d(lat2*3600))
    
    return obj


def hybridcluster(meas):
    """ use both DBSCAN and sequential clustering to cluster the data"""

    # Hybrid clustering algorithm
    # 1) Find "object" centers by using DBSCAN with a smallish eps (~0.2-0.3") and maybe minclusters of 2-3
    # 2) Do sequential clustering using the object centers on the leftover measurements.

    # Empty catalog input
    if len(meas)==0:
        return np.array([]), np.array([])

    # Only one exposure, don't cluster
    expindex = dln.create_index(meas['exposure'])
    nexp = len(expindex['value'])
    if nexp==1:
        print('Only one exposure. Do not need to cluster')
        labels = np.arange(len(meas))
        obj = np.zeros(len(meas),dtype=np.dtype([('label',int),('ndet',int),('ra',np.float64),('dec',np.float64),
                                                ('raerr',np.float32),('decerr',np.float32),('asemi',np.float32),
                                                ('bsemi',np.float32),('theta',np.float32),('fwhm',np.float32)]))
        obj['label'] = labels
        obj['ndet'] = 1
        for n in ['ra','dec','raerr','decerr','asemi','bsemi','theta','fwhm']: obj[n]=meas[n.upper()]
        return labels, obj

    
    # Step 1: Find object centers using DBSCAN with a small eps
    t0 = time.time()
    # DBSCAN does not deal with cos(dec), convert to a different projection
    cenra = np.mean(meas['ra'])
    cendec = np.mean(meas['dec'])
    # Deal with RA=0 wrap
    if (np.max(meas['ra'])-np.min(meas['ra']))>100:
        rr = meas['ra']
        bb,nbb = dln.where(rr>180)
        if nbb>0: rr[bb]-=360
        cenra = np.mean(rr)
        if cenra<0: cenra+=360
    lon,lat = coords.rotsphcen(meas['ra'],meas['dec'],cenra,cendec,gnomic=True)
    X1 = np.column_stack((lon,lat))
    err = np.sqrt(meas['raerr']**2+meas['decerr']**2)
    eps = np.maximum(3*np.median(err),0.3)
    print('DBSCAN eps=%4.2f' % eps)
    # Minimum number of measurements needed to define a cluster/object
    minsamples = 3
    if nexp<3: minsamples=nexp
    dbs1 = DBSCAN(eps=eps/3600, min_samples=minsamples).fit(X1)
    gdb,ngdb,bdb,nbdb = dln.where(dbs1.labels_ >= 0,comp=True)
    # No clusters, lower minsamples
    while (ngdb==0):
        minsamples -= 1
        print('No clusters. Lowering min_samples to '+str(minsamples))
        dbs1 = DBSCAN(eps=eps/3600, min_samples=minsamples).fit(X1)
        gdb,ngdb,bdb,nbdb = dln.where(dbs1.labels_ >= 0,comp=True)
    print('DBSCAN after %5.2f sec. ' % (time.time()-t0))

    # Get mean coordinates for each object
    #   only use the measurements that were clustered
    obj1 = meancoords(meas[gdb],dbs1.labels_[gdb])
    inpobj = obj1
    print(str(ngdb)+' measurements clustered into '+str(len(obj1))+' objects. '+str(nbdb)+' remaining.')
    
    # Step 2: sequential clustering with original list of objects with the outliers
    #  this allows more distance measurements with larger errors to be clustered as well
    #  the RA/DEC uncertainties can be very small, set a lower threshold of EPS
    if (nbdb>0):
        print('Sequential Clustering the remaining measurements')
        dcr = np.maximum(3*err[bdb],eps)
        measrem = meas[bdb]
        labels2, obj2 = seqcluster(measrem,dcr=dcr,inpobj=inpobj)
        # Add these new labels to the original list
        #  offset the numbers so they don't overlap
        labels = dbs1.labels_
        labels[bdb] = labels2+np.max(labels)+1
        obj = meancoords(meas,labels)    # Get mean coordinates again
    else:
        obj = obj1
        labels = dbs1.labels_

    print(str(len(obj))+' final objects')
    
    return labels, obj

def calibmeas(meas1,chmeta1,meta,version='v4',verbose=False):
    """ Apply astrometric and photometric calibrations to raw measurement data """
    
    nmeas1 = len(meas1)
    mdt = [('measid',(str,100)),('objectid',(str,100)),('exposure',(str,50)),
           ('ccdnum',np.int8),('filter',(str,50)),('mjd',float),('x',np.float32),('y',np.float32),
           ('ra',float),('raerr',np.float32),('dec',float),('decerr',np.float32),('mag_auto',np.float32),
           ('magerr_auto',np.float32),('mag_aper1',np.float32),('magerr_aper1',np.float32),('mag_aper2',np.float32),
           ('magerr_aper2',np.float32),('mag_aper4',np.float32),('magerr_aper4',np.float32),
           ('mag_aper6',np.float32),('magerr_aper6',np.float32),('mag_aper8',np.float32),
           ('magerr_aper8',np.float32),('mag_iso',np.float32),('magerr_iso',np.float32),
           ('kron_radius',np.float32),('background',np.float32),('threshold',np.float32),('isoarea_image',np.float32),
           ('isoarea_world',np.float32),('x2_world',np.float32),('y2_world',np.float32),('xy_world',np.float32),
           ('asemi',np.float32),('asemierr',np.float32),('bsemi',np.float32),('bsemierr',np.float32),
           ('theta',np.float32),('thetaerr',np.float32),('ellipticity',np.float32),
           ('errx2_world',np.float32),('erry2_world',np.float32),('errxy_world',np.float32),
           ('fwhm',np.float32),('flags',np.int16),('imaflags_iso',np.int32),('nimaflags_iso',np.int32),
           ('class_star',np.float32),
           ('ndet_iter',int),('repeat',int),('xpsf',float),('ypsf',float),('magpsf',float),
           ('errpsf',float),('skypsf',float),('iter',int),('chi',float),('sharp',float),
           ('rapsf',float),('decpsf',float),('ebv',np.float32),('haspsf',bool),('snr',np.float32),
           ('badflag',int)]
    final = np.zeros(nmeas1,dtype=np.dtype(mdt))
    final = Table(final)
    # Add units
    for c in ['ra','dec','rapsf','decpsf','theta','thetaerr']:
        final[c].unit = 'degree'
    for c in ['raerr','decerr','fwhm','asemi','asemierr','bsemi','bsemierr','kron_radius']:
        final[c].unit = 'arcsec'
    for c in ['mag_auto','magerr_auto','mag_aper1','magerr_aper1','mag_aper2','magerr_aper2',
              'mag_aper4','magerr_aper4','mag_aper6','magerr_aper6','mag_aper8','magerr_aper8',
              'mag_iso','magerr_iso','magpsf','errpsf','ebv']:
        final[c].unit = 'magnitude'
    for c in ['x','y','xpsf','ypsf']:
        final[c].unit = 'pixel'
    for c in ['background','skypsf']:
        final[c].unit = 'ADU'
    final['mjd'].unit = 'day'
    # Copy over all columns in common
    for n in final.colnames:
        if n in meas1.colnames:
            final[n][:] = meas1[n]
    # Special handling of some columns
    final['exposure'][:] = meta['base'][0]
    final['x'][:] = meas1['x_image']
    final['y'][:] = meas1['y_image']
    # Calibrate the photometry
    #  only calibrate good photometric measurements
    zpoff = 2.5*np.log10(meta['exptime'][0]) + chmeta1['zpterm'][0]
    for c in ['mag_auto','mag_aper1','mag_aper2','mag_aper4','mag_aper6','mag_aper8','magpsf']:
        if c[:8]=='mag_aper':
            indx = {'1':0,'2':1,'4':2,'6':3,'8':4}[c[-1]]
            mag = meas1['mag_aper'][:,indx]
        else:
            mag = meas1[c]
        gdmag = ((mag < 50) & np.isfinite(mag))
        final[c][:] = np.nan   # bad by default
        final[c][gdmag] = mag[gdmag]+zpoff
    final['magerr_auto'][:] = meas1['magerr_auto']
    final['magerr_aper1'][:] = meas1['magerr_aper'][:,0] 
    final['magerr_aper2'][:] = meas1['magerr_aper'][:,1] 
    final['magerr_aper4'][:] = meas1['magerr_aper'][:,2]
    final['magerr_aper6'][:] = meas1['magerr_aper'][:,3] 
    final['magerr_aper8'][:] = meas1['magerr_aper'][:,4]
    final['asemi'][:] = meas1['a_world'] * 3600.         # convert to arcsec 
    final['asemierr'][:] = meas1['erra_world'] * 3600.   # convert to arcsec 
    final['bsemi'][:] = meas1['b_world'] * 3600.         # convert to arcsec 
    final['bsemierr'][:] = meas1['errb_world'] * 3600.   # convert to arcsec 
    final['theta'][:] = 90-meas1['theta_world']          # make CCW E of N 
    final['thetaerr'][:] = meas1['errtheta_world'] 
    final['fwhm'][:] = meas1['fwhm_world'] * 3600.       # convert to arcsec 
    final['skypsf'][:] = meas1['sky']
    final['haspsf'][:] = meas1['magpsf'] < 50
    final['snr'][:] = 1.087/meas1['magerr_auto']

    # measid, objectid, ccdnum, filter, mjd
    
    # Set some catalog values
    filt = meta['filter'][0]
    mjd = meta['mjd'][0]
    instrument = meta['instrument'][0]
    expnum = meta['expnum'][0]
    ccdnum = chmeta1['ccdnum'][0]
    final['ccdnum'] = ccdnum
    final['filter'] = filt
    final['mjd'] = mjd 
    final['measid'] = instrument+'.'+str(expnum)+'.'+np.char.array(final['ccdnum'].astype(str))+'.'+np.char.array(meas1['number'].astype(str))
    
    if version < 'v4':
        racol = 'alpha_j2000'
        deccol = 'delta_j2000'
        xcol = 'x_image'
        ycol = 'y_image'
    else:
        racol = 'rapsf'
        deccol = 'decpsf'
        xcol = 'xpsf'
        ycol = 'ypsf'
        
    meas1['ccdnum'] = chmeta1['ccdnum'][0]
        
    # Apply QA cuts 
    #---------------- 
    badflag = np.zeros(len(meas1),int)

    # Mask bad chip data 
    # Half of chip 31 for MJD>56660 
    #  c4d_131123_025436_ooi_r_v2 with MJD=56619 has problems too 
    #  if the background b/w the left and right half is large then BAD
    lft31, = np.where((meas1[xcol] < 1024) & (meas1['ccdnum'] == 31))
    rt31, = np.where((meas1[xcol] >= 1024) & (meas1['ccdnum'] == 31))
    if len(lft31) > 10 and len(rt31) > 10: 
        lftback = np.median(meas1['background'][lft31]) 
        rtback = np.median(meas1['background'][rt31]) 
        mnback = 0.5*(lftback+rtback) 
        sigback = dln.mad(meas1['background']) 
        if np.abs(lftback-rtback) > (np.sqrt(mnback)>sigback): 
            jump31 = True
        else: 
            jump31 = False
        #logger.info('  Big jump in CCDNUM 31 background levels' 
    else:
        jump31 = False 
    if meta['mjd'] > 56600 or jump31 == 1: 
        badchip31 = True 
        # Remove bad measurements 
        # X: 1-1024 okay 
        # X: 1025-2049 bad 
        # use 1000 as the boundary since sometimes there's a sharp drop 
        # at the boundary that causes problem sources with SExtractor 
        bdind, = np.where((meas1[xcol] > 1000) & (meas1['ccdnum'] == 31))
        ngdind = len(meas1)-len(bdind)
        if len(bdind) > 0:  # some bad ones found
            if verbose: print('  '+str(len(bdind))+' bad chip 31 measurements')
            badflag[bdind] += 1
        else:
            badchip31 = False  # chip 31 
                         
    # Mask sourcds with bad quality mask flag (IMAFLAGS_ISO) 
    bdmeas, = np.where(meas1['imaflags_iso'] > 0)
    if len(bdmeas) > 0: 
        if verbose: print('  '+str(len(bdmeas))+' sources with bad CP flags.')
        badflag[bdmeas] += 2 
            
    # Mask sources with bad SE FLAGS 
    #   this masks problematic truncatd sources near chip edges 
    bdseflags, = np.where( ((meas1['flags'] & 8) > 0) |   # object truncated 
                           ((meas1['flags'] & 16) > 0))   # aperture truncate 
    if len(bdseflags) > 0: 
        if verbose: print('  '+str(len(bdseflags))+' truncated sources')
        badflag[bdseflags] += 4
            
    # Mask low-S/N sources 
    #  snr = 1.087/err 
    snrcut = 5.0 
    bdsnr, = np.where(1.087/meas1['magerr_auto'] < snrcut) 
    if len(bdsnr) > 0: 
        if verbose: print('  '+str(len(bdsnr))+' sources with S/N<'+str(snrcut))
        badflag[bdsnr] += 8
    
    final['badflag'][:] = badflag
    # Other columns we are keeping in case we need to recreate the original measurement file
    # already copied over in for loop above

    # Calibrate coordinates
    racoef = chmeta1['racoef'][0]
    deccoef = chmeta1['deccoef'][0]
    rarms = chmeta1['rarms'][0]
    decrms = chmeta1['decrms'][0]
    lon,lat = coords.rotsphcen(meas1['alpha_j2000'],meas1['delta_j2000'],
                               chmeta1['cenra'][0],chmeta1['cendec'][0],gnomic=True)
    lon2 = lon + utils.func_poly2d(lon,lat,*racoef)
    lat2 = lat + utils.func_poly2d(lon,lat,*deccoef) 
    ra2,dec2 = coords.rotsphcen(lon2,lat2,chmeta1['cenra'][0],chmeta1['cendec'][0],reverse=True,gnomic=True)
    final['ra'] = ra2 
    final['dec'] = dec2
    if version >= 'v4':
        gdmeas, = np.where(np.isfinite(meas1['rapsf']) & np.isfinite(meas1['decpsf']))
        lon,lat = coords.rotsphcen(meas1['rapsf'][gdmeas],meas1['decpsf'][gdmeas],
                                   chmeta1['cenra'][0],chmeta1['cendec'][0],gnomic=True)
        lon2 = lon + utils.func_poly2d(lon,lat,*racoef)
        lat2 = lat + utils.func_poly2d(lon,lat,*deccoef) 
        ra2,dec2 = coords.rotsphcen(lon2,lat2,chmeta1['cenra'][0],chmeta1['cendec'][0],reverse=True,gnomic=True)
        final['rapsf'][gdmeas] = ra2 
        final['decpsf'][gdmeas] = dec2
        
    # Add to astrometric errors 
    final['raerr'] = np.sqrt(final['raerr']**2 + rarms**2) 
    final['decerr'] = np.sqrt(final['decerr']**2 + decrms**2) 

    # Get reddening
    coo = SkyCoord(ra=final['ra'],dec=final['dec'],unit='deg')
    sfd = SFDQuery()
    ebv = sfd(coo)                
    final['ebv'] = ebv 
    
    return final
    

def loadmeas(metafile=None,buffdict=None,dbfile=None,verbose=False):

    t0 = time.time()

    if metafile is None:
        print('Need metafile')
        return np.array([]), np.array([])

    # New meta-data format
    dtype_meta = np.dtype([('file',str,500),('base',str,200),('expnum',int),('ra',np.float64),
                           ('dec',np.float64),('dateobs',str,100),('mjd',np.float64),('filter',str,50),
                           ('exptime',float),('airmass',float),('nsources',int),('fwhm',float),
                           ('nchips',int),('badchip31',bool),('rarms',float),('decrms',float),
                           ('ebv',float),('gaianmatch',int),('zpterm',float),('zptermerr',float),
                           ('zptermsig',float),('refmatch',int)])

    # All columns in MEAS catalogs (32)
    dtype_meas = np.dtype([('measid',str,30),('objectid',str,50),('exposure',str,40),('ccdnum',np.int8),('filter',str,3),
                           ('mjd',float),('ra',float),('raerr',np.float32),('dec',float),('decerr',np.float32),
                           ('mag_auto',np.float32),('magerr_auto',np.float32),('asemi',np.float32),('asemierr',np.float32),
                           ('bsemi',np.float32),('bsemierr',np.float32),('theta',np.float32),('thetaerr',np.float32),
                           ('fwhm',np.float32),('flags',np.int16),('class_star',np.float32),('magpsf',np.float32),
                           ('errpsf',np.float32),('skypsf',np.float32),('chi',np.float32),('sharp',np.float32),
                           ('rapsf',np.float64),('decpsf',np.float64),('haspsf',bool),('snr',np.float32),('badflag',bool)])
    
    #  Loop over exposures
    meas = None
    nmeas = 0
    allmeta = None
    meascount = 0
    metafile = np.atleast_1d(metafile)
    for m,mfile in enumerate(metafile):
        expmeascount = 0
        if os.path.exists(mfile)==False and os.path.exists(mfile+'.gz')==False:
            print(mfile+' NOT FOUND')
            continue
        if os.path.exists(mfile):
            mhdu = fits.open(mfile)
        else:
            mhdu = fits.open(mfile+'.gz')
        print(str(m+1)+' Loading '+mfile)
        meta = mhdu[1].data
        t = Time(meta['dateobs'], format='isot', scale='utc')
        meta['mjd'] = t.mjd                    # recompute because some MJD are bad

        chmeta = []
        for k in range(2,len(mhdu)):
            chmeta.append(np.array(mhdu[k].data))
        chmeta = Table(np.array(chmeta))
        #chmeta = fits.getdata(mfile,2)      # chip-level meta-data table
        print('  FILTER='+meta['filter'][0]+'  EXPTIME='+str(meta['exptime'][0])+' sec')

        v = psutil.virtual_memory()
        process = psutil.Process(os.getpid())
        print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

        # Convert META to new format
        newmeta = np.zeros(1,dtype=dtype_meta)
        # Copy over the meta information
        for c in newmeta.dtype.names:
            if c in meta.dtype.names: newmeta[c]=meta[c]

        # Get the name
        fdir = os.path.dirname(mfile)
        fbase, ext = os.path.splitext(os.path.basename(mfile))
        fbase = fbase[:-5]   # remove _meta at end

        measfile = os.path.join(fdir,fbase+'_meas.fits')
        if os.path.exists(measfile)==False:
            measfile += '.gz'
        if os.path.exists(measfile)==False:
            print(measfile,'NOT FOUND')
            import pdb; pdb.set_trace()
        hdu = fits.open(measfile)
        expmeta = meta[0]

        version = mfile.split('/')[-6]
        
        # Loop over the chip files
        for j in range(len(chmeta)):
            # Check that this chip was astrometrically calibrated
            #   and falls in to HEALPix region
            # Also check for issues with my astrometric corrections
            astokay = True
            if (chmeta['ngaiamatch'][j] == 0) | (np.max(np.abs(chmeta['racoef'][j]))>1) | (np.max(np.abs(chmeta['deccoef'][j]))>1):
                if verbose: print('This chip was not astrometrically calibrated or has astrometric issues')
                astokay = False
                
            # Check that this overlaps the healpix region
            inside = True
            if buffdict is not None:
                vra = chmeta['vra'][j].flatten()
                vdec = chmeta['vdec'][j].flatten()
                vlon, vlat = coords.rotsphcen(vra,vdec,buffdict['cenra'],buffdict['cendec'],gnomic=True)
                if coords.doPolygonsOverlap(buffdict['lon'],buffdict['lat'],vlon,vlat) is False:
                    if verbose: print('This chip does NOT overlap the HEALPix region+buffer')
                    inside = False

            # Check if the chip-level file exists
            #chfile = fdir+'/'+fbase+'_'+str(chmeta['ccdnum'][j])+'_meas.fits'
            #chfile_exists = os.path.exists(chfile)
            #if chfile_exists is False:
            #    print(chfile+' NOT FOUND')

            # Load this one
            #if (chfile_exists is True) and (inside is True) and (astokay is True):
            if (inside is True) and (astokay is True):
                # Load the chip-level catalog
                #meas1 = fits.getdata(chfile,1)
                meas1 = Table(hdu[j+1].data)
                for c in meas1.colnames: meas1[c].name = c.lower()
                nmeas1 = len(meas1)
                chmeta1 = chmeta[j]
                #print('  chip '+str(chmeta[j]['ccdnum'])+'  '+str(nmeas1)+' sources')
                
                # Calibrate the data, if necessary
                if 'measid' not in meas1.colnames:
                    print('calibrating chip '+str(chmeta1['ccdnum'][0]))
                    meas1 = calibmeas(meas1,chmeta1,meta,version=version,verbose=False)
                    
                # Fix negative FWHM values
                #  use A_WORLD and B_WORLD which are never negative
                bd,nbd = dln.where(meas1['fwhm']<0.1)
                if nbd>0:
                    meas1['fwhm'][bd] = np.sqrt(meas1['asemi'][bd]**2+meas1['bsemi'][bd]**2)*2.35
                # Fix RAERR=DECERR=0
                bd,nbd = dln.where(meas1['raerr']<0.0001)
                if nbd>0:
                    snr = 1.087/meas1['magerr_auto'][bd]
                    coorderr = 0.664*meas1['fwhm'][bd]/snr
                    meas1['raerr'][bd] = coorderr
                    meas1['decerr'][bd] = coorderr

                # Make sure it's in the right format
                #if len(meas1.dtype.fields) != 32:    # v3
                if len(meas1.dtype.fields) != 65:     # v4
                    if verbose: print('  This catalog does not have the right format. Skipping')
                    del meas1
                    nmeas1 = 0

                # Only include sources inside Boundary+Buffer zone
                #  -use ROI_CUT
                #  -reproject to tangent plane first so we don't have to deal
                #     with RA=0 wrapping or pole issues
                if buffdict is not None:
                    lon, lat = coords.rotsphcen(meas1['ra'],meas1['dec'],buffdict['cenra'],buffdict['cendec'],gnomic=True)
                    ind_out, ind_in = dln.roi_cut(buffdict['lon'],buffdict['lat'],lon,lat)
                    nmatch = dln.size(ind_in)
                    # Only want source inside this pixel
                    if nmatch>0:
                        meas1 = meas1[ind_in]
                        nmeas1 = len(meas1)
                    else:
                        meas1 = None
                        nmeas1 = 0
                    #if verbose: print('  '+str(nmatch)+' sources are inside this pixel')
                    
                # Combine the catalogs
                if nmeas1 > 0:
                    # Keep it all in memory
                    if dbfile is None:
                        if meas is None:
                            #dtype_meas = meas1.dtype
                            #nmeas_init = np.sum(chmeta['nsources'])*dln.size(metafile)
                            nmeas_init = np.maximum(100000,nmeas1)
                            meas = np.zeros(nmeas_init,dtype=dtype_meas)
                            meascount = 0
                        # Add more elements if necessary
                        if (meascount+nmeas1)>nmeas:
                            meas = add_elements(meas,np.maximum(100000,nmeas1))
                            nmeas = len(meas)

                        # Add it to the main MEAS measalog
                        for c in dtype_meas.names: meas[c][meascount:meascount+nmeas1] = meas1[c]
                    # Use the database
                    else:
                        writecat2db(meas1,dbfile)

                    if verbose: print('  chip '+str(chmeta['ccdnum'][j])+'  '+str(nmeas1)+' measurements')

                    meascount += nmeas1
                    expmeascount += nmeas1


        # Add metadata to ALLMETA, only if some measurements overlap
        if expmeascount>0:
            if allmeta is None:
                allmeta = newmeta
            else:
                allmeta = np.hstack((allmeta,newmeta))
        # Total measurements for this exposure
        print('  '+str(expmeascount)+' measurements')
        print(str(meascount)+' measurements total so far')
        
    #print('all exposures loaded. trimming now')
    if (meas is not None) & (meascount<nmeas): meas=meas[0:meascount]   # delete excess elements
    if meas is None: meas=np.array([])         # empty meas
    if allmeta is None: allmeta=np.array([])

    print('loading measurements done after '+str(time.time()-t0))

    return meas, meascount, allmeta

def clusterdata(meas,nmeas,dbfile=None):
    """ Perform spatial clustering """

    t00 = time.time()
    print('Spatial clustering')    
    # Divide into subregions
    if (nmeas>1000000) & (dbfile is not None):
        print('Dividing clustering problem into subregions')
        # Index RA and DEC
        createindexdb(dbfile,'ra',unique=False)
        createindexdb(dbfile,'dec',unique=False)
        db.analyzetable(dbfile,'meas')
        # Subdivide
        nsub = int(np.ceil(nmeas/100000))
        print(str(nsub)+' sub regions')
        nx = int(np.ceil(np.sqrt(nsub)))  # divide RA and DEC intro nx regions
        # Get RA/DEC ranges from the database
        ranges = getradecrangedb(dbfile)  # [min(ra),max(ra),min(dec),max(dec)]
        xr = [ranges[0]-0.001, ranges[1]+0.001]  # extend slightly 
        print('RA: '+str(xr[0])+' '+str(xr[1]))
        dx = (xr[1]-xr[0])/nx
        if (xr[1]-xr[0])>180:   # across RA=0
            dx = (xr[0]-(xr[1]-360))/nx
        yr = [ranges[2]-0.001, ranges[3]+0.001]  # extend slightly 
        mndec = np.mean(yr)
        print('DEC: '+str(yr[0])+' '+str(yr[1]))
        dy = (yr[1]-yr[0])/nx
        buff = 10./3600.0  # buffer in arc seconds
        rabuff = buff/np.cos(np.deg2rad(mndec))  # correct for cos(dec)
        objtab = np.zeros(100000,dtype=np.dtype([('OBJLABEL',int),('RA',float),('DEC',float),('NMEAS',int)]))
        nobjtab = len(objtab)
        # Loop over sub regions
        lastobjlabel = -1
        objcount = 0
        # RA loop
        for r in range(nx):
            r0 = xr[0]+r*dx
            r1 = xr[0]+(r+1)*dx
            # DEC loop
            for d in range(nx):
                d0 = yr[0]+d*dy
                d1 = yr[0]+(d+1)*dy
                print(str(r+1)+' '+str(d+1))
                print('RA: '+str(r0)+' '+str(r1)+'  DEC: '+str(d0)+' '+str(d1))
                meas1 = getdatadb(dbfile,rar=[r0-rabuff,r1+rabuff],decr=[d0-buff,d1+buff],verbose=True)
                nmeas1 = len(meas1)
                if nmeas1>0:
                    gmeas1,ngmeas1 = dln.where(meas1['OBJLABEL']==-1)  # only want ones that haven't been taken yet
                    if ngmeas1>0:
                        meas1 = meas1[gmeas1]
                    nmeas1 = len(meas1)
                    print(str(nmeas1)+' measurements with no labels')

                v = psutil.virtual_memory()
                process = psutil.Process(os.getpid())
                print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

                # Some measurements to work with
                if nmeas1>0:
                    # Cluster
                    t0 = time.time()
                    # Cluster labels are integers and in ascending order, but there are gaps
                    objlabels1, initobj1 = hybridcluster(meas1)
                    objlabels1 += lastobjlabel+1                 # add offset to labels
                    labelindex1 = dln.create_index(objlabels1)   # create inex
                    nobj1 = len(labelindex1['value'])
                    print(str(nmeas1)+' measurements for '+str(nobj1)+' objects')
                    # Compute mean positions
                    obj1 = np.zeros(nobj1,dtype=np.dtype([('objlabel',int),('ra',float),('dec',float),('nmeas',int)]))
                    obj1['objlabel'] = labelindex1['value']
                    obj1['nmeas'] = labelindex1['num']
                    for k in range(nobj1):
                        indx = labelindex1['index'][labelindex1['lo'][k]:labelindex1['hi'][k]+1]
                        wt_ra = 1.0/meas1['raerr'][indx]**2
                        wt_dec = 1.0/meas1['decerr'][indx]**2
                        obj1['ra'][k] = np.sum(meas1['ra'][indx]*wt_ra)/np.sum(wt_ra)
                        obj1['dec'][k] = np.sum(meas1['dec'][indx]*wt_dec)/np.sum(wt_dec)
                    # Only keep objects (and measurements) inside the box region
                    #  keep objects on LOWER boundary in RA/DEC
                    gdobj, ngdobj = dln.where((obj1['ra']>=r0) & (obj1['ra']<r1) & (obj1['dec']>=d0) & (obj1['dec']<d1))
                    print(str(ngdobj)+' objects are inside the boundary')
                    # Some objects in the region
                    if ngdobj>0:
                        obj1 = obj1[gdobj]
                        nobj1 = ngdobj
                        # Arrays of rowid and objlabels to add
                        add_rowid1 = np.zeros(np.sum(labelindex1['num'][gdobj]),int)
                        add_objlabels1 = np.zeros(np.sum(labelindex1['num'][gdobj]),int)
                        cnt1 = 0
                        for k in range(ngdobj):
                            indx = labelindex1['index'][labelindex1['lo'][gdobj[k]]:labelindex1['hi'][gdobj[k]]+1]
                            nmeas1 = labelindex1['num'][gdobj[k]]
                            add_rowid1[cnt1:cnt1+nmeas1] = meas1['rowid'][indx]
                            add_objlabels1[cnt1:cnt1+nmeas1] = labelindex1['value'][gdobj[k]]
                            cnt1 += nmeas1

                        # Add the object labels into the database
                        #  much faster if in rowid order
                        si = np.argsort(add_rowid1)
                        insertobjlabelsdb(add_rowid1[si],add_objlabels1[si],dbfile)

                        # Add OBJ1 to OBJTAB
                        if (objcount+nobj1>nobjtab):    # add new elements
                            print('Adding more elements to OBJTAB')
                            t1 = time.time()
                            objtab = add_elements(objtab,np.max([nobj1,100000]))
                            nobjtab = len(objtab)
                            print('more elements added in '+str(time.time()-t1)+' sec.')
                        objtab[objcount:objcount+nobj1] = obj1
                        objcount += nobj1

                        # Keep track of last label
                        lastobjlabel = np.max(obj1['objlabel'])

                        #import pdb; pdb.set_trace()

        # Trim extra elements
        if nobjtab>objcount:
            objtab = objtab[0:objcount]

    # No subdividing
    else:
        # Get MEASID, RA, DEC from database
        if dbfile is not None:
            #meas = getdbcoords(dbfile)
            meas = getdatadb(dbfile,verbose=True)
        objlabels, initobj = hybridcluster(meas)
        labelindex = dln.create_index(objlabels)   # create index
        nobj = len(labelindex['value'])
        print(str(nmeas)+' measurements for '+str(nobj)+' objects')
        # Make structure
        objtab = np.zeros(nobj,dtype=np.dtype([('objlabel',int),('nmeas',int),('lo',int),('hi',int)]))
        objtab['objlabel'] = labelindex['value']
        objtab['nmeas'] = labelindex['num']
        nobjtab = len(objtab)
        # Insert object label into database
        if dbfile is not None:
            insertobjlabelsdb(meas['rowid'],objlabels,dbfile)
        # Resort MEAS, and use index LO/HI
        meas = meas[labelindex['index']]
        objtab['lo'] = labelindex['lo']
        objtab['hi'] = labelindex['hi']

    print(str(len(objtab))+' final objects')

    # Index objlabel in database
    if dbfile is not None:
        createindexdb(dbfile,'objlabel',unique=False)

    print('clustering done after '+str(time.time()-t00)+' sec.')

    return objtab, meas


def breakup_idtab(dbfile):
    """ Break-up idtab file into separate measid/objectid lists per exposure on /data0."""

    t00 = time.time()

    outdir = '/data0/dnidever/nsc/instcal/v3/idtab/'

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
            dbbase1 = os.path.basename(dbfile1)[0:-9]  # remove _idtab.db ending
            # Get existing index names for this database
            d = sqlite3.connect(dbfile1, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
            cur = d.cursor()
            cmd = 'select measid,exposure,objectid from idtab'
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
            df = np.dtype([('measid',str,measid_maxlen+1),('objectid',str,objectid_maxlen+1)])
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




# Combine data for one NSC healpix region
def combine(pix,version,nside=128,redo=False,verbose=False,multilevel=True,outdir=None,nmulti=None):

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    radeg = np.float64(180.00) / np.pi

    tmpdir = '/tmp/'  # default
    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        dir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
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

    # Only nside>=128 supported right now
    if nside<128:
        print('Only nside=>128 supported')
        sys.exit()

    #print('*** KLUDGE: Forcing output to /scratch1 ***')
    #outdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/combine/'
    #outdir = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/combine/'
    outdir = '/home/group/davidnidever/nsc/instcal/v4/combine/'
    if os.path.exists(outdir) is False: os.mkdir(outdir)

    # nside>128
    if nside > 128:
        # Get parent nside=128 pixel
        pra,pdec = hp.pix2ang(nside,pix,lonlat=True)
        parentpix = hp.ang2pix(128,pra,pdec,lonlat=True)
        print('The nside=128 parent pixel is '+str(parentpix))
        # Output filenames
        outbase = str(parentpix)+'_n'+str(int(nside))+'_'+str(pix)
        subdir = str(int(parentpix)//1000)    # use the thousands to create subdirectory grouping
        if os.path.exists(outdir+'/'+subdir) is False: os.mkdir(outdir+'/'+subdir)
        outfile = outdir+'/'+subdir+'/'+outbase+'.fits'

    # nside=128
    else:
        # Output filenames
        outbase = str(pix)
        subdir = str(int(pix)//1000)    # use the thousands to create subdirectory grouping
        if os.path.exists(outdir+'/'+subdir) is False: os.mkdir(outdir+'/'+subdir)
        outfile = outdir+'/'+subdir+'/'+str(pix)+'.fits'

    # Check if output file already exists
    if (os.path.exists(outfile) or os.path.exists(outfile+'.gz')) & (not redo):
        print(outfile+' EXISTS already and REDO not set')
        sys.exit()

    print("Combining InstCal catalogs for Healpix pixel = "+str(pix))

    # Use the healpix list, nside=128
    #listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_instcal_combine_healpix_list.db'
    #listfile = '/home1/09970/dnidever/scratch1/nsc/instcal/'+version+'/lists/nsc_instcal_combine_healpix_list.db'
    #listfile = '/corral/projects/NOIRLab/nsc/instcal/'+version+'/lists/nsc_instcal_combine_healpix_list.db'
    listfile = '/home/group/davidnidever/nsc/instcal/'+version+'/lists/nsc_instcal_combine_healpix_list.db'
    if os.path.exists(listfile) is False:
        print(listfile+" NOT FOUND")
        sys.exit()


    # nside>128
    if nside > 128:
        # Find our pixel
        hlist = db.query(listfile,'hlist',where='PIX='+str(parentpix))
        nlist = len(hlist)
        if nlist == 0:
            print("No entries for Healpix pixel '"+str(parentpix)+"' in the list")
            sys.exit()
        hlist = Table(hlist)
        # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
        #  so we can deal with the edge cases
        neipix = hp.get_all_neighbours(128,parentpix)
        for neip in neipix:
            hlist1 = db.query(listfile,'hlist',where='PIX='+str(neip))
            nhlist1 = len(hlist1)
            if nhlist1>0:
                hlist1 = Table(hlist1)
                hlist = vstack([hlist,hlist1])

    # nside=128
    else:
        parentpix = pix
        # Find our pixel
        hlist = db.query(listfile,'hlist',where='PIX='+str(pix))
        nlist = len(hlist)
        if nlist == 0:
            print("No entries for Healpix pixel '"+str(pix)+"' in the list")
            sys.exit()
        hlist = Table(hlist)
        # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
        #  so we can deal with the edge cases
        neipix = hp.get_all_neighbours(nside,pix)
        for neip in neipix:
            hlist1 = db.query(listfile,'hlist',where='PIX='+str(neip))
            nhlist1 = len(hlist1)
            if nhlist1>0:
                hlist1 = Table(hlist1)
                hlist = vstack([hlist,hlist1])

    # Fix filenames for tempest
    hlist['measfile'] = [f.replace('/home1/09970/dnidever/scratch1/','/home/group/davidnidever/')+'.gz' for f in hlist['measfile']]

    
    # Rename to be consistent with the FITS file
    hlist['measfile'].name = 'FILE'
    hlist['base'].name = 'BASE'
    hlist['pix'].name = 'PIX'

    # Use entire exposure files
    # Get unique values
    u, ui = np.unique(hlist['FILE'],return_index=True)
    hlist = hlist[ui]
    nhlist = len(hlist)
    print(str(nhlist)+' exposures that overlap this pixel and neighbors')

    # Get the boundary coordinates
    #   healpy.boundaries but not sure how to do it in IDL
    #   pix2vec_ring/nest can optionally return vertices but only 4
    #     maybe subsample myself between the vectors
    # Expand the boundary to include a "buffer" zone
    #  to deal with edge cases
    vecbound = hp.boundaries(nside,pix,step=100)
    rabound, decbound = hp.vec2ang(np.transpose(vecbound),lonlat=True)

    # Expand the boundary by the buffer size
    cenra, cendec = hp.pix2ang(nside,pix,lonlat=True)
    # reproject onto tangent plane
    lonbound, latbound = coords.rotsphcen(rabound,decbound,cenra,cendec,gnomic=True)
    # expand by a fraction, it's not an extact boundary but good enough
    buffsize = 10.0/3600. # in deg
    radbound = np.sqrt(lonbound**2+latbound**2)
    frac = 1.0 + 1.5*np.max(buffsize/radbound)
    lonbuff = lonbound*frac
    latbuff = latbound*frac
    rabuff, decbuff = coords.rotsphcen(lonbuff,latbuff,cenra,cendec,gnomic=True,reverse=True)
    if (np.max(rabuff)-np.min(rabuff))>100:  # deal with RA=0 wraparound
        bd,nbd = dln.where(rabuff>180)
        if nbd>0:rabuff[bd] -=360.0
    buffdict = {'cenra':cenra,'cendec':cendec,'rar':dln.minmax(rabuff),'decr':dln.minmax(decbuff),
                'ra':rabuff,'dec':decbuff,'lon':lonbuff,'lat':latbuff,
                'lr':dln.minmax(lonbuff),'br':dln.minmax(latbuff)}

    # IDTAB schema
    dtype_idtab = np.dtype([('measid',str,200),('exposure',str,200),('objectid',str,200),('objectindex',int)])

    # OBJ schema
    dtype_obj = np.dtype([('objectid',str,100),('pix',int),('ra',np.float64),('dec',np.float64),('raerr',np.float32),
                          ('decerr',np.float32),('pmra',np.float32),('pmdec',np.float32),('pmraerr',np.float32),
                          ('pmdecerr',np.float32),('mjd',np.float64),('deltamjd',np.float32),('ndet',np.int16),
                          ('nphot',np.int16),('ndetu',np.int16),('nphotu',np.int16),('umag',np.float32),
                          ('urms',np.float32),('uerr',np.float32),('uasemi',np.float32),('ubsemi',np.float32),
                          ('utheta',np.float32),('ndetg',np.int16),('nphotg',np.int16),('gmag',np.float32),
                          ('grms',np.float32),('gerr',np.float32),('gasemi',np.float32),('gbsemi',np.float32),
                          ('gtheta',np.float32),('ndetr',np.int16),('nphotr',np.int16),('rmag',np.float32),
                          ('rrms',np.float32),('rerr',np.float32),('rasemi',np.float32),('rbsemi',np.float32),
                          ('rtheta',np.float32),('ndeti',np.int16),('nphoti',np.int16),('imag',np.float32),
                          ('irms',np.float32),('ierr',np.float32),('iasemi',np.float32),('ibsemi',np.float32),
                          ('itheta',np.float32),('ndetz',np.int16),('nphotz',np.int16),('zmag',np.float32),
                          ('zrms',np.float32),('zerr',np.float32),('zasemi',np.float32),('zbsemi',np.float32),
                          ('ztheta',np.float32),('ndety',np.int16),('nphoty',np.int16),('ymag',np.float32),
                          ('yrms',np.float32),('yerr',np.float32),('yasemi',np.float32),('ybsemi',np.float32),
                          ('ytheta',np.float32),('ndetvr',np.int16),('nphotvr',np.int16),('vrmag',np.float32),
                          ('vrrms',np.float32),('vrerr',np.float32),('vrasemi',np.float32),('vrbsemi',np.float32),
                          ('vrtheta',np.float32),('asemi',np.float32),('asemierr',np.float32),('bsemi',np.float32),
                          ('bsemierr',np.float32),('theta',np.float32),('thetaerr',np.float32),('fwhm',np.float32),
                          ('flags',np.int16),('class_star',np.float32),('ebv',np.float32),('rmsvar',np.float32),
                          ('madvar',np.float32),('iqrvar',np.float32),('etavar',np.float32),
                          ('jvar',np.float32),('kvar',np.float32),('chivar',np.float32),('romsvar',np.float32),
                          ('variable10sig',np.int16),('nsigvar',np.float32),('overlap',bool)])

    # Estimate number of measurements in pixel
    #metafiles = [m.replace('_cat','_meta').strip() for m in hlist['FILE']]
    metafiles = [m.replace('_meas','_meta').strip() for m in hlist['FILE']]
    metafiles = [m.replace('.fits.gz','.fits') for m in metafiles]
    metatab = checkboundaryoverlap(metafiles,buffdict,verbose=True)  # verbose=False
    nmeasperarea = np.zeros(dln.size(metatab),int)
    areadict = {'c4d':3.0, 'k4m':0.3, 'ksb':1.0}  # total area
    for j in range(dln.size(metatab)):
        nmeasperarea[j] = metatab['nsources'][j]/areadict[metatab['instrument'][j]]
    pixarea = hp.nside2pixarea(nside,degrees=True)
    nmeasperpix = nmeasperarea * pixarea
    totmeasest = np.sum(nmeasperpix)


    # Break into smaller healpix regions
    if (multilevel is True) & (nside == 128):
        nsub = int(np.ceil(totmeasest/500000))
        bestval,bestind = dln.closest([1,4,16,64],nsub)
        hinside = [128,256,512,1024][bestind]
        # Break into multiple smaller healpix
        if hinside>128:
            print('')
            print('----- Breaking into smaller HEALPix using nside='+str(hinside)+' ------')
            vecbound = hp.boundaries(nside,pix)
            allpix = hp.query_polygon(hinside,np.transpose(vecbound))
            print('Pix = '+','.join(allpix.astype(str)))
            outfiles = []
            # Check if any healpix need to be run/rerun
            dopix = []
            for i in range(len(allpix)):
                pix1 = allpix[i]
                # check the output file
                outbase1 = str(parentpix)+'_n'+str(int(hinside))+'_'+str(pix1)
                subdir1 = str(int(parentpix)//1000)    # use the thousands to create subdirectory grouping
                outfile1 = outdir+'/'+subdir1+'/'+outbase1+'.fits.gz'
                outfiles.append(outfile1)
                if (os.path.exists(outfile1) is False) | redo:
                    dopix.append(pix1)
            print(str(len(dopix))+' nside='+str(hinside)+' healpix to run')

            # Some healpix to run
            if len(dopix)>0:
                # Single process, just use subprocess
                if nmulti==1:
                    for i in range(len(dopix)):
                        pix1 = dopix[i]
                        print('')
                        print('########### '+str(i+1)+' '+str(pix1)+' ###########')
                        print('')
                        # check the output file
                        outbase1 = str(parentpix)+'_n'+str(int(hinside))+'_'+str(pix1)
                        subdir1 = str(int(parentpix)//1000)    # use the thousands to create subdirectory grouping
                        outfile1 = outdir+'/'+subdir1+'/'+outbase1+'.fits.gz'
                        if redo is True:
                            retcode = subprocess.call(['python',os.path.abspath(__file__),str(pix1),version,'--nside',str(hinside),'-r'],shell=False)
                        else:
                            retcode = subprocess.call(['python',os.path.abspath(__file__),str(pix1),version,'--nside',str(hinside)],shell=False)
                # Multiple parallel processes, Running job daemon
                else:
                    cmd = []
                    for i in range(len(dopix)):
                        cmd1 = os.path.abspath(__file__)+' '+str(dopix[i])+' '+version+' --nside '+str(hinside)
                        if redo: cmd1 = cmd1+' -r'
                        cmd.append(cmd1)
                    dirs = np.zeros(len(dopix),(str,200))
                    dirs[:] = tmpdir
                    jobs = jd.job_daemon(cmd,dirs,hyperthread=True,prefix='nsccmb',nmulti=nmulti)

            # Load and concatenate all of the files
            print('Combining all of the object catalogs')
            allmeta = None
            allobj = None
            nobjects = []
            totobjects = 0
            for i in range(len(allpix)):
                pix1 = allpix[i]
                outfile1 = outfiles[i]
                if os.path.exists(outfile1) is False:
                    print(outfile1+' NOT FOUND')
                    sys.exit()
                # meta columns different: nobjects   there'll be repeats
                meta1 = fits.getdata(outfile1,1)
                if allmeta is None:
                    allmeta = meta1
                else:
                    allmeta = np.hstack((allmeta,meta1))
                hd1 = fits.getheader(outfile1,2)
                print(str(i+1)+' '+outfile1+' '+str(hd1['naxis2']))
                obj1 = fits.getdata(outfile1,2)
                nobj1 = len(obj1)

                # Update the objectIDs
                dbfile_idtab1 = outfile1.replace('.fits.gz','_idtab.db')
                objectid_orig = obj1['objectid']
                objectid_new = dln.strjoin( str(parentpix)+'.', ((np.arange(nobj1)+1+totobjects).astype(str)) )
                #updatecoldb(selcolname,selcoldata,updcolname,updcoldata,table,dbfile):
                updatecoldb('objectid',objectid_orig,'objectid',objectid_new,'idtab',dbfile_idtab1)
                # Update objectIDs in catalog
                obj1['objectid'] = objectid_new

                # Update objectIDs in high resolution HEALPix output file
                print('Updating objectIDs in '+outfile1)
                outfile1fits = outfile1.replace('.fits.gz','.fits')
                if os.path.exists(outfile1fits): os.remove(outfile1fits)
                Table(meta1).write(outfile1fits)               # first, summary table
                #  append other fits binary tables
                hdulist = fits.open(outfile1fits)
                hdu = fits.table_to_hdu(Table(obj1))        # second, catalog
                hdulist.append(hdu)
                hdulist.writeto(outfile1fits,overwrite=True)
                hdulist.close()
                if os.path.exists(outfile1): os.remove(outfile1)
                ret = subprocess.call(['gzip',outfile1fits])    # compress final catalog

                if allobj is None:
                    allobj = obj1.copy()
                else:
                    allobj = np.hstack((allobj,obj1.copy()))
                nobjects.append(nobj1)
                totobjects += nobj1

            # Deal with duplicate metas
            metaindex = dln.create_index(allmeta['base'])
            for i in range(len(metaindex['value'])):
                indx = metaindex['index'][metaindex['lo'][i]:metaindex['hi'][i]+1]
                meta1 = allmeta[indx[0]].copy()
                if len(indx)>1:
                    meta1['nobjects'] = np.sum(allmeta['nobjects'][indx])
                if i==0:
                    sumstr = meta1
                else:
                    sumstr = np.hstack((sumstr,meta1))
            sumstr = Table(sumstr)

            # Write the output file
            print('Writing combined catalog to '+outfile)
            if os.path.exists(outfile): os.remove(outfile)
            sumstr.write(outfile)               # first, summary table
            #  append other fits binary tables
            hdulist = fits.open(outfile)
            hdu = fits.table_to_hdu(Table(allobj))        # second, catalog
            hdulist.append(hdu)
            hdulist.writeto(outfile,overwrite=True)
            hdulist.close()
            if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
            ret = subprocess.call(['gzip',outfile])    # compress final catalog

            dt = time.time()-t0
            print('dt = '+str(dt)+' sec.')

            print('Breaking-up IDTAB information')
            dbfiles_idtab = []
            for i in range(len(allpix)):
                outfile1 = outfiles[i]
                dbfile_idtab1 = outfile1.replace('.fits.gz','_idtab.db')
                dbfiles_idtab.append(dbfile_idtab1)
            breakup_idtab(dbfiles_idtab)            

            sys.exit()


    # Decide whether to load everything into RAM or use temporary database
    usedb = False
    if totmeasest>500000: usedb=True
    dbfile = None
    if usedb:
        dbfile = tmproot+outbase+'_combine.db'
        print('Using temporary database file = '+dbfile)
        if os.path.exists(dbfile): os.remove(dbfile)
    else:
        print('Keeping all measurement data in memory')

    #import pdb; pdb.set_trace()

    # IDTAB database file
    dbfile_idtab = outdir+'/'+subdir+'/'+outbase+'_idtab.db'
    if os.path.exists(dbfile_idtab): os.remove(dbfile_idtab)

    # Load the measurement catalog
    #  this will contain excess rows at the end, if all in RAM
    #  if using database, MEAS is empty
    meas, meascount, allmeta = loadmeas(metafiles,buffdict,dbfile=dbfile)
    nmeas = meascount
    print(str(nmeas))


    #import pdb; pdb.set_trace()
    
    # No measurements
    if nmeas==0:
        print('No measurements for this healpix')
        if (dbfile is not None):
            if os.path.exists(dbfile): os.remove(dbfile)
        if os.path.exists(dbfile_idtab): os.remove(dbfile_idtab)
        print('Writing blank output file to '+outfile)
        fits.PrimaryHDU().writeto(outfile)
        if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
        ret = subprocess.call(['gzip',outfile])    # compress final catalog
        sys.exit()

    # Spatially cluster the measurements with DBSCAN
    #   this might also resort MEAS
    objtab, meas = clusterdata(meas,nmeas,dbfile=dbfile)
    nobj = dln.size(objtab)
    meascumcount = np.cumsum(objtab['nmeas'])
    print(str(nobj)+' unique objects clustered')

    # Initialize the OBJ structured array
    obj = np.zeros(nobj,dtype=dtype_obj)
    # if nside>128 then we need unique IDs, so use PIX and *not* PARENTPIX
    #  add nside as well to make it truly unique
    if nside>128:
        obj['objectid'] = dln.strjoin( str(nside)+'.'+str(pix)+'.', ((np.arange(nobj)+1).astype(str)) )
    else:
        obj['objectid'] = dln.strjoin( str(pix)+'.', ((np.arange(nobj)+1).astype(str)) )
    obj['pix'] = parentpix    # use PARENTPIX
    # all bad to start
    for f in ['pmra','pmraerr','pmdec','pmdecerr','asemi','bsemi','theta','asemierr',
              'bsemierr','thetaerr','fwhm','class_star','rmsvar','madvar','iqrvar',
              'etavar','jvar','kvar','chivar','romsvar']: obj[f]=np.nan
    for f in ['u','g','r','i','z','y','vr']:
        obj[f+'mag'] = 99.99
        obj[f+'err'] = 9.99
        obj[f+'rms'] = np.nan
        obj[f+'asemi'] = np.nan
        obj[f+'bsemi'] = np.nan
        obj[f+'theta'] = np.nan
    obj['variable10sig'] = 0
    obj['nsigvar'] = np.nan
    #idtab = np.zeros(ncat,dtype=dtype_idtab)

    # Initialize temporary IDTAB structure
    idtab = np.zeros(100000,dtype=dtype_idtab)
    nidtab = dln.size(idtab)

    # Higher precision catalog
    dtype_himeas = np.dtype([('measid',str,30),('objectid',str,50),('exposure',str,40),('ccdnum',int),('filter',str,3),
                             ('mjd',float),('ra',float),('raerr',float),('dec',float),('decerr',float),
                             ('mag_auto',float),('magerr_auto',float),('asemi',float),('asemierr',float),
                             ('bsemi',float),('bsemierr',float),('theta',float),('thetaerr',float),
                             ('fwhm',float),('flags',int),('class_star',float),('magpsf',float),
                             ('errpsf',float),('skypsf',float),('chi',float),('sharp',float),
                             ('rapsf',float),('decpsf',float),('haspsf',bool),('snr',float),('badflag',bool)])

    # Convert to numpy structured array
    dtype_himeasdb = np.dtype([('measid',str,30),('objlabel',int),('exposure',str,40),('ccdnum',int),('filter',str,3),
                               ('mjd',float),('ra',float),('raerr',float),('dec',float),('decerr',float),
                               ('mag_auto',float),('magerr_auto',float),('asemi',float),('asemierr',float),
                               ('basemi',float),('bsemierr',float),('theta',float),('thetaerr',float),
                               ('fwhm',float),('flags',int),('class_star',float),('magpsf',float),
                               ('errpsf',float),('skypsf',float),('chi',float),('sharp',float),
                               ('rapsf',float),('decpsf',float),('haspsf',bool),('snr',float),('badflag',bool)])


    #dtype_meas = np.dtype([('measid',str,30),('objectid',str,50),('exposure',str,40),('ccdnum',np.int8),('filter',str,3),
    #                       ('mjd',float),('ra',float),('raerr',np.float32),('dec',float),('decerr',np.float32),
    #                       ('mag_auto',np.float32),('magerr_auto',np.float32),('asemi',np.float32),('asemierr',np.float32),
    #                       ('bsemi',np.float32),('bsemierr',np.float32),('theta',np.float32),('thetaerr',np.float32),
    #                       ('fwhm',np.float32),('flags',np.int16),('class_star',np.float32),('magpsf',np.float32),
    #                       ('errpsf',np.float32),('skypsf',np.float32),('chi',np.float32),('sharp',np.float32),
    #                       ('rapsf',np.float64),('decpsf',np.float64),('haspsf',bool),('snr',np.float32),('badflag',bool)])



    t1 = time.time()

    # Loop over the objects
    meascount = 0
    ngroup = -1
    grpcount = 0
    maxmeasload = 50000
    ngrpmeas = 0
    nmeas1 = 0
    idtab_count = 0
    idtab_grpcount = 0
    fidmag = np.zeros(nobj,float)+np.nan  # fiducial magnitude
    for i,lab in enumerate(objtab['objlabel']):
        if (i % 1000)==0: print(i)

        if (i % 1000)==0:
            v = psutil.virtual_memory()
            process = psutil.Process(os.getpid())
            print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

        # Get meas data for this object
        if usedb is False:
            oindx = np.arange(objtab['lo'][i],objtab['hi'][i]+1)  # this fails if start,stop are the same
            if objtab['nmeas'][i]==1: oindx=np.atleast_1d(objtab['lo'][i])
            nmeas1 = dln.size(oindx)
            meas1_orig = meas[oindx]
            # Upgrade precisions of catalog
            meas1 = np.zeros(nmeas1,dtype=dtype_himeas)
            meas1[...] = meas1_orig   # stuff in the data
            #for n in dtype_himeas.names: meas1[n] = meas1_orig[n]
            del meas1_orig
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
                lab0 = lab
                lab1 = objtab['objlabel'][np.min([i+ngroup-1,nobj-1])]
                #lab1 = labelindex['value'][np.min([i+ngroup-1,nobj-1])]
                if ngrpmeas>0: del grpmeas
                if nmeas1>0: del meas1
                grpmeas = getdatadb(dbfile,objlabel=[lab0,lab1])
                ngrpmeas = dln.size(grpmeas)
                grpindex = dln.create_index(grpmeas['objlabel'])
                #ngroup = len(grpindex['value'])
                grpcount = 0
            # Get the measurement data for this object
            gindx = grpindex['index'][grpindex['lo'][grpcount]:grpindex['hi'][grpcount]+1]
            meas1 = np.atleast_1d(grpmeas[gindx])
            nmeas1 = len(meas1)
            grpcount += 1
            oindx = np.arange(nmeas1)+meascount
            meascount += nmeas1            

        obj['ndet'][i] = nmeas1


        # Add IDTAB information to IDTAB structure/database
        #  update in groups to database so it takes less time
        if idtab_count+nmeas1 > nidtab:
            print('  Adding more elements to temporary IDTAB structure')
            idtab = add_elements(idtab,50000)  # add more elements if necessary
        # Add information to temporary IDTAB structure for this object
        idtab['measid'][idtab_count:idtab_count+nmeas1] = meas1['measid']
        idtab['exposure'][idtab_count:idtab_count+nmeas1] = meas1['exposure']
        idtab['objectid'][idtab_count:idtab_count+nmeas1] = obj['objectid'][i]
        idtab['objectindex'][idtab_count:idtab_count+nmeas1] = i
        idtab_count += nmeas1
        idtab_grpcount += 1
        # Write to database and reinitialize the temporary IDTAB structure
        if (idtab_grpcount>5000) | (idtab_count>30000) |  (i==(nobj-1)):
            print('  Writing data to IDTAB database')
            writeidtab2db(idtab[0:idtab_count],dbfile_idtab)
            idtab = np.zeros(100000,dtype=dtype_idtab)
            nidtab = dln.size(idtab)
            idtab_count = 0
            idtab_grpcount = 0

        # Computing quantities
        # Mean RA/DEC, RAERR/DECERR
        if nmeas1>1:
            wt_ra = 1.0/meas1['raerr']**2
            wt_dec = 1.0/meas1['decerr']**2
            obj['ra'][i] = np.sum(meas1['ra']*wt_ra)/np.sum(wt_ra)
            obj['raerr'][i] = np.sqrt(1.0/np.sum(wt_ra))
            obj['dec'][i] = np.sum(meas1['dec']*wt_dec)/np.sum(wt_dec)
            obj['decerr'][i] = np.sqrt(1.0/np.sum(wt_dec))
            obj['mjd'][i] = np.mean(meas1['mjd'])
            obj['deltamjd'][i] = np.max(meas1['mjd'])-np.min(meas1['mjd'])
        else:
            obj['ra'][i] = meas1['ra']
            obj['dec'][i] = meas1['dec']
            obj['raerr'][i] = meas1['raerr']
            obj['decerr'][i] = meas1['decerr']
            obj['mjd'][i] = meas1['mjd']
            obj['deltamjd'][i] = 0

        # Check for negative RA values
        if obj['ra'][i] < 0:
            obj['ra'][i] += 360

        # Mean proper motion and errors
        if nmeas1>1:
            raerr = np.array(meas1['raerr']*1e3,np.float64)    # milli arcsec
            ra = np.array(meas1['ra'],np.float64)
            ra -= np.mean(ra)
            ra *= 3600*1e3 * np.cos(obj['dec'][i]/radeg)     # convert to true angle, milli arcsec
            t = meas1['mjd'].copy()
            t -= np.mean(t)
            t /= 365.2425                          # convert to year
            # Calculate robust slope
            pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
            obj['pmra'][i] = pmra                 # mas/yr
            obj['pmraerr'][i] = pmraerr           # mas/yr

            decerr = np.array(meas1['decerr']*1e3,np.float64)   # milli arcsec
            dec = np.array(meas1['dec'],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            # Calculate robust slope
            pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
            obj['pmdec'][i] = pmdec               # mas/yr
            obj['pmdecerr'][i] = pmdecerr         # mas/yr

        # Mean magnitudes
        # Convert totalwt and totalfluxwt to MAG and ERR
        #  and average the morphology parameters PER FILTER
        filtindex = dln.create_index(meas1['filter'].astype(str))
        nfilters = len(filtindex['value'])
        resid = np.zeros(nmeas1)+np.nan     # residual mag
        relresid = np.zeros(nmeas1)+np.nan  # residual mag relative to the uncertainty
        for f in range(nfilters):
            filt = filtindex['value'][f].lower()
            findx = filtindex['index'][filtindex['lo'][f]:filtindex['hi'][f]+1]
            obj['ndet'+filt][i] = filtindex['num'][f]
            gph,ngph = dln.where(meas1['mag_auto'][findx]<50)
            obj['nphot'+filt][i] = ngph
            if ngph==1:
                obj[filt+'mag'][i] = meas1['mag_auto'][findx[gph]]
                obj[filt+'err'][i] = meas1['magerr_auto'][findx[gph]]
            if ngph>1:
                newmag, newerr = dln.wtmean(meas1['mag_auto'][findx[gph]], meas1['magerr_auto'][findx[gph]],
                                            magnitude=True,reweight=True,error=True)
                obj[filt+'mag'][i] = newmag
                obj[filt+'err'][i] = newerr
                # Calculate RMS
                obj[filt+'rms'][i] = np.sqrt(np.mean((meas1['mag_auto'][findx[gph]]-newmag)**2))
                # Residual mag
                resid[findx[gph]] = meas1['mag_auto'][findx[gph]]-newmag
                # Residual mag relative to the uncertainty
                #  set a lower threshold of 0.02 in the uncertainty
                relresid[findx[gph]] = np.sqrt(ngph/(ngph-1)) * (meas1['mag_auto'][findx[gph]]-newmag)/np.maximum(meas1['magerr_auto'][findx[gph]],0.02)

            # Calculate mean morphology parameters
            obj[filt+'asemi'][i] = np.mean(meas1['asemi'][findx])
            obj[filt+'bsemi'][i] = np.mean(meas1['bsemi'][findx])
            obj[filt+'theta'][i] = np.mean(meas1['theta'][findx])

        # Calculate variability indices
        gdresid = np.isfinite(resid)
        ngdresid = np.sum(gdresid)
        if ngdresid>0:
            resid2 = resid[gdresid]
            sumresidsq = np.sum(resid2**2)
            tsi = np.argsort(meas1['mjd'][gdresid])
            resid2tsi = resid2[tsi]
            quartiles = np.percentile(resid2,[25,50,75])
            # RMS
            rms = np.sqrt(sumresidsq/ngdresid)
            # MAD
            madvar = 1.4826*np.median(np.abs(resid2-quartiles[1]))
            # IQR
            iqrvar = 0.741289*(quartiles[2]-quartiles[0])
            # 1/eta
            etavar = sumresidsq / np.sum((resid2tsi[1:]-resid2tsi[0:-1])**2)
            obj['rmsvar'][i] = rms
            obj['madvar'][i] = madvar
            obj['iqrvar'][i] = iqrvar
            obj['etavar'][i] = etavar

        # Calculate variability indices wrt to uncertainties
        gdrelresid = np.isfinite(relresid)
        ngdrelresid = np.sum(gdrelresid)
        if ngdrelresid>0:
            relresid2 = relresid[gdrelresid]
            pk = relresid2**2-1
            jvar = np.sum( np.sign(pk)*np.sqrt(np.abs(pk)) )/ngdrelresid
            #avgrelvar = np.mean(np.abs(relresid2))    # average of absolute relative residuals
            chivar = np.sqrt(np.sum(relresid2**2))/ngdrelresid
            kdenom = np.sqrt(np.sum(relresid2**2)/ngdrelresid)
            if kdenom!=0:
                kvar = (np.sum(np.abs(relresid2))/ngdrelresid) / kdenom
            else:
                kvar = np.nan
            # RoMS
            romsvar = np.sum(np.abs(relresid2))/(ngdrelresid-1)
            obj['jvar'][i] = jvar
            obj['kvar'][i] = kvar
            #obj['avgrelvar'][i] = avgrelvar
            obj['chivar'][i] = chivar
            obj['romsvar'][i] = romsvar
            #if chivar>50: import pdb; pdb.set_trace()

        # Make NPHOT from NPHOTX
        obj['nphot'][i] = obj['nphotu'][i]+obj['nphotg'][i]+obj['nphotr'][i]+obj['nphoti'][i]+obj['nphotz'][i]+obj['nphoty'][i]+obj['nphotvr'][i]

        # Fiducial magnitude, used to select variables below
        #  order of priority: r,g,i,z,Y,VR,u
        if obj['nphot'][i]>0:
            magarr = np.zeros(7,float)
            for ii,nn in enumerate(['rmag','gmag','imag','zmag','ymag','vrmag','umag']): magarr[ii]=obj[nn][i]
            gfid,ngfid = dln.where(magarr<50)
            if ngfid>0: fidmag[i]=magarr[gfid[0]]

        # Mean morphology parameters
        obj['asemi'][i] = np.mean(meas1['asemi'])
        obj['bsemi'][i] = np.mean(meas1['bsemi'])
        obj['theta'][i] = np.mean(meas1['theta'])
        obj['asemierr'][i] = np.sqrt(np.sum(meas1['asemierr']**2)) / nmeas1
        obj['bsemierr'][i] = np.sqrt(np.sum(meas1['bsemierr']**2)) / nmeas1
        obj['thetaerr'][i] = np.sqrt(np.sum(meas1['thetaerr']**2)) / nmeas1
        obj['fwhm'][i] = np.mean(meas1['fwhm'])
        obj['class_star'][i] = np.mean(meas1['class_star'])
        obj['flags'][i] = np.bitwise_or.reduce(meas1['flags'])  # OR combine


    v = psutil.virtual_memory()
    process = psutil.Process(os.getpid())
    print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

    # Created OBJECTID index in IDTAB database
    createindexdb(dbfile_idtab,'objectid',table='idtab',unique=False)
    createindexdb(dbfile_idtab,'exposure',table='idtab',unique=False)
    db.analyzetable(dbfile_idtab,'idtab')


    # Select Variables
    #  1) Construct fiducial magnitude (done in loop above)
    #  2) Construct median VAR and sigma VAR versus magnitude
    #  3) Find objects that Nsigma above the median VAR line
    si = np.argsort(fidmag)   # NaNs are at end
    varcol = 'madvar'
    gdvar,ngdvar,bdvar,nbdvar = dln.where(np.isfinite(obj[varcol]) & np.isfinite(fidmag),comp=True)
    if ngdvar>0:
        nbins = np.ceil((np.max(fidmag[gdvar])-np.min(fidmag[gdvar]))/0.25)
        nbins = int(np.max([2,nbins]))
        fidmagmed, bin_edges1, binnumber1 = bindata.binned_statistic(fidmag[gdvar],fidmag[gdvar],statistic='nanmedian',bins=nbins)
        numhist, _, _ = bindata.binned_statistic(fidmag[gdvar],fidmag[gdvar],statistic='count',bins=nbins)
        # Fix NaNs in fidmagmed
        bdfidmagmed,nbdfidmagmed = dln.where(np.isfinite(fidmagmed)==False)
        if nbdfidmagmed>0:
            fidmagmed_bins = 0.5*(bin_edges1[0:-1]+bin_edges1[1:])
            fidmagmed[bdfidmagmed] = fidmagmed_bins[bdfidmagmed]
        # Median metric
        varmed, bin_edges2, binnumber2 = bindata.binned_statistic(fidmag[gdvar],obj[varcol][gdvar],statistic='nanmedian',bins=nbins)
        # Smooth, it handles NaNs well
        smlen = 5
        smvarmed = dln.gsmooth(varmed,smlen)
        bdsmvarmed,nbdsmvarmed = dln.where(np.isfinite(smvarmed)==False)
        if nbdsmvarmed>0:
            smvarmed[bdsmvarmed] = np.nanmedian(smvarmed)
        # Interpolate to all the objects
        gv,ngv,bv,nbv = dln.where(np.isfinite(smvarmed),comp=True)
        fvarmed = interp1d(fidmagmed[gv],smvarmed[gv],kind='linear',bounds_error=False,
                           fill_value=(smvarmed[0],smvarmed[-1]),assume_sorted=True)
        objvarmed = np.zeros(nobj,float)
        objvarmed[gdvar] = fvarmed(fidmag[gdvar])
        objvarmed[gdvar] = np.maximum(np.min(smvarmed[gv]),objvarmed[gdvar])   # lower limit
        if nbdvar>0: objvarmed[bdvar]=smvarmed[gv[-1]]   # objects with bad fidmag, set to last value
        # Scatter in metric around median
        #  calculate MAD ourselves so that it's around our computed median metric line
        varsig, bin_edges3, binnumber3 = bindata.binned_statistic(fidmag[gdvar],np.abs(obj[varcol][gdvar]-objvarmed[gdvar]),
                                                                  statistic='nanmedian',bins=nbins)
        varsig *= 1.4826   # scale MAD to stddev
        # Fix values for bins with few points
        bdhist,nbdhist,gdhist,ngdhist = dln.where(numhist<3,comp=True)
        if nbdhist>0:
            if ngdhist>0:
                varsig[bdhist] = np.nanmedian(varsig[gdhist])
            else:
                varsig[:] = 0.02
            
        # Smooth
        smvarsig = dln.gsmooth(varsig,smlen)
        # Interpolate to all the objects
        gv,ngv,bv,nbv = dln.where(np.isfinite(smvarsig),comp=True)
        fvarsig = interp1d(fidmagmed[gv],smvarsig[gv],kind='linear',bounds_error=False,
                           fill_value=(smvarsig[gv[0]],smvarsig[gv[-1]]),assume_sorted=True)
        objvarsig = np.zeros(nobj,float)
        objvarsig[gdvar] = fvarsig(fidmag[gdvar])
        objvarsig[gdvar] = np.maximum(np.min(smvarsig[gv]),objvarsig[gdvar])   # lower limit
        if nbdvar>0: objvarsig[bdvar]=smvarsig[gv[-1]]   # objects with bad fidmag, set to last value
        # Detect positive outliers
        nsigvarthresh = 10.0
        nsigvar = (obj[varcol]-objvarmed)/objvarsig
        obj['nsigvar'][gdvar] = nsigvar[gdvar]
        isvar,nisvar = dln.where(nsigvar[gdvar]>nsigvarthresh)
        print(str(nisvar)+' variables detected')
        if nisvar>0:
            obj['variable10sig'][gdvar[isvar]] = 1

    # Add E(B-V)
    print('Getting E(B-V)')
    sfd = SFDQuery()
    c = SkyCoord(obj['ra'],obj['dec'],frame='icrs',unit='deg')
    #c = SkyCoord('05h00m00.00000s','+30d00m00.0000s', frame='icrs') 
    ebv = sfd(c)
    obj['ebv'] = ebv

    
    # FIGURE OUT IF THERE ARE OBJECTS **INSIDE** OTHER OBJECTS!!
    #   could be a deblending problem, extended galaxy that was shredded, or asteroids going through
    obj = find_obj_parent(obj)
    bd,nbd = dln.where(obj['parent']==True)
    print(str(nbd)+' objects have other objects inside their footprint')

    
    # ONLY INCLUDE OBJECTS WITH AVERAGE RA/DEC
    # WITHIN THE BOUNDARY OF THE HEALPIX PIXEL!!!
    ipring = hp.pixelfunc.ang2pix(nside,obj['ra'],obj['dec'],lonlat=True)
    ind1,nmatch = dln.where(ipring == pix)
    if nmatch==0:
        print('None of the final objects fall inside the pixel')
        if (dbfile is not None):
            if os.path.exists(dbfile): os.remove(dbfile)
        if os.path.exists(dbfile_idtab): os.remove(dbfile_idtab)
        print('Writing blank output file to '+outfile)
        fits.PrimaryHDU().writeto(outfile)
        if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
        ret = subprocess.call(['gzip',outfile])    # compress final catalog
        sys.exit()
    # Get trimmed objects and indices
    objtokeep = np.zeros(nobj,bool)         # boolean to keep or trim objects
    objtokeep[ind1] = True
    if nmatch<nobj:                         # some to trim
        trimind = np.arange(nobj)
        trimind = np.delete(trimind,ind1)
        trimobj = obj[trimind]          # trimmed objects
    newobjindex = np.zeros(nobj,int)-1    # new indices
    newobjindex[ind1] = np.arange(nmatch)
    # Keep the objects inside the Healpix
    obj = obj[ind1]
    print(str(nmatch)+' final objects fall inside the pixel')

    #import pdb; pdb.set_trace()

    # Remove trimmed objects from IDTAB database
    if nmatch<nobj:
        # Delete measurements for the objects that we are trimming
        deleterowsdb('objectid',trimobj['objectid'],'idtab',dbfile_idtab)
        # Update OBJECTINDEX for the objects that we are keeping
        updatecoldb('objectid',obj['objectid'],'objectindex',np.arange(nmatch),'idtab',dbfile_idtab)

    v = psutil.virtual_memory()
    process = psutil.Process(os.getpid())
    print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

    # Get unique exposures in IDTAB database
    uexposure = executedb(dbfile_idtab,'SELECT DISTINCT exposure from idtab')
    # this returns a list of tuples, unpack
    uexposure = [i[0] for i in uexposure]
    # create sumstr for these using allmeta
    nuexposure = len(uexposure)
    ind1,ind2 = dln.match(allmeta['base'],uexposure)
    nmatch = len(ind1)
    sumstr = Table(allmeta[ind1])
    col_nobj = Column(name='nobjects', dtype=int, length=len(sumstr))
    col_healpix = Column(name='healpix', dtype=int, length=len(sumstr))
    sumstr.add_columns([col_nobj, col_healpix])
    sumstr['nobjects'] = 0
    sumstr['healpix'] = parentpix   # use PARENTPIX
    # get number of objects per exposure
    data = executedb(dbfile_idtab,'SELECT exposure, count(DISTINCT objectid) from idtab GROUP BY exposure')
    out = np.zeros(len(data),dtype=np.dtype([('exposure',str,40),('nobjects',int)]))
    out[...] = data
    ind1,ind2 = dln.match(sumstr['base'],out['exposure'])
    sumstr['nobjects'][ind1] = out['nobjects'][ind2]


    v = psutil.virtual_memory()
    process = psutil.Process(os.getpid())
    print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

    # save the measurement data to a file
    #outmeasfile = outdir+'/'+subdir+'/'+str(pix)+'_meas.fits'
    #if os.path.exists(outmeasfile): os.remove(outmeasfile)
    #Table(cat).write(outmeasfile)

    # Write the output file
    print('Writing combined catalog to '+outfile)
    if os.path.exists(outfile): os.remove(outfile)
    sumstr.write(outfile)               # first, summary table
    #  append other fits binary tables
    hdulist = fits.open(outfile)
    hdu = fits.table_to_hdu(Table(obj))        # second, catalog
    hdulist.append(hdu)
    # The IDTAB table is now in a stand-alone sqlite3 database called PIX_idtab.db
    #hdu = fits.table_to_hdu(Table(idtab))      # third, ID table
    #hdulist.append(hdu)    
    hdulist.writeto(outfile,overwrite=True)
    hdulist.close()
    if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
    ret = subprocess.call(['gzip',outfile])    # compress final catalog

    dt = time.time()-t0
    print('dt = '+str(dt)+' sec.')
    print('dt = ',str(time.time()-t1)+' sec. after loading the catalogs')

    if dbfile is not None:
        print('Deleting temporary database file '+dbfile)
        os.remove(dbfile)

    # Delete all arrays before we quit
    del sumstr
    del obj
    # garbage collection
    gc.collect()

    # Breaking up idtab information
    if nside==128:
        print('Breaking-up IDTAB information')
        breakup_idtab(dbfile_idtab)
