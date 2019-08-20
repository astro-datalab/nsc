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

def createindexdb(dbfile,col='measid',unique=True):
    """ Index a column in the database """
    t0 = time.time()
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    index_name = 'idx_'+col+'_meas'
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
        c.execute('CREATE UNIQUE INDEX '+index_name+' ON meas('+col+')')
    else:
        c.execute('CREATE INDEX '+index_name+' ON meas('+col+')')
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

    # Convert to nump structured array
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

def loadmeas(metafile=None,buffdict=None,dbfile=None,verbose=False):

    t0 = time.time()

    if metafile is None:
        print('Need metafile')
        return np.array([]), np.array([])

    # New meta-data format
    dtype_meta = np.dtype([('file',np.str,500),('base',np.str,200),('expnum',int),('ra',np.float64),
                           ('dec',np.float64),('dateobs',np.str,100),('mjd',np.float64),('filter',np.str,50),
                           ('exptime',float),('airmass',float),('nsources',int),('fwhm',float),
                           ('nchips',int),('badchip31',bool),('rarms',float),('decrms',float),
                           ('ebv',float),('gaianmatch',int),('zpterm',float),('zptermerr',float),
                           ('zptermsig',float),('refmatch',int)])

    # All columns in MEAS catalogs (32)
    #dtype_cat = np.dtype([('MEASID',np.str,200),('OBJECTID',np.str,200),('EXPOSURE',np.str,200),('CCDNUM',int),('FILTER',np.str,10),
    #                      ('MJD',float),('X',float),('Y',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
    #                      ('MAG_AUTO',float),('MAGERR_AUTO',float),('MAG_APER1',float),('MAGERR_APER1',float),('MAG_APER2',float),
    #                      ('MAGERR_APER2',float),('MAG_APER4',float),('MAGERR_APER4',float),('MAG_APER8',float),('MAGERR_APER8',float),
    #                      ('KRON_RADIUS',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),('THETA',float),
    #                      ('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])
    # All the columns that we need (20)
    #dtype_cat = np.dtype([('MEASID',np.str,30),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
    #                      ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
    #                      ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
    #                      ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])
    dtype_cat = np.dtype([('MEASID',np.str,30),('EXPOSURE',np.str,40),('CCDNUM',np.int8),('FILTER',np.str,3),
                          ('MJD',float),('RA',float),('RAERR',np.float16),('DEC',float),('DECERR',np.float16),
                          ('MAG_AUTO',np.float16),('MAGERR_AUTO',np.float16),('ASEMI',np.float16),('ASEMIERR',np.float16),
                          ('BSEMI',np.float16),('BSEMIERR',np.float16),('THETA',np.float16),('THETAERR',np.float16),
                          ('FWHM',np.float16),('FLAGS',np.int16),('CLASS_STAR',np.float16)])

    #  Loop over exposures
    cat = None
    ncat = 0
    allmeta = None
    catcount = 0
    metafile = np.atleast_1d(metafile)
    for m,mfile in enumerate(metafile):
        expcatcount = 0
        if os.path.exists(mfile) is False:
            print(mfile+' NOT FOUND')
            continue
        meta = fits.getdata(mfile,1)
        print(str(m+1)+' Loading '+mfile)
        t = Time(meta['dateobs'], format='isot', scale='utc')
        meta['mjd'] = t.mjd                    # recompute because some MJD are bad
        chmeta = fits.getdata(mfile,2)      # chip-level meta-data structure
        print('  FILTER='+meta['filter'][0]+'  EXPTIME='+str(meta['exptime'][0])+' sec')

        # Convert META to new format
        newmeta = np.zeros(1,dtype=dtype_meta)
        # Copy over the meta information
        for n in newmeta.dtype.names:
            if n.upper() in meta.dtype.names: newmeta[n]=meta[n]

        # Get the name
        fdir = os.path.dirname(mfile)
        fbase, ext = os.path.splitext(os.path.basename(mfile))
        fbase = fbase[:-5]   # remove _meta at end
        
        # Loop over the chip files
        for j in range(len(chmeta)):
            # Check that this chip was astrometrically calibrated
            #   and falls in to HEALPix region
            if chmeta['ngaiamatch'][j] == 0:
                if verbose: print('This chip was not astrometrically calibrate')

            # Check that this overlaps the healpix region
            inside = True
            if buffdict is not None:
                vra = chmeta['vra'][j]
                vdec = chmeta['vdec'][j]
                if (np.max(vra)-np.min(vra)) > 100:    # deal with RA=0 wrapround
                    bd,nbd = dln.where(vra>180)
                    if nbd>0: vra[bd] -= 360
                if coords.doPolygonsOverlap(buffdict['ra'],buffdict['dec'],vra,vdec) is False:
                    if verbose: print('This chip does NOT overlap the HEALPix region+buffer')
                    inside = False

            # Check if the chip-level file exists
            chfile = fdir+'/'+fbase+'_'+str(chmeta['ccdnum'][j])+'_meas.fits'
            if os.path.exists(chfile) is False:
                print(chfile+' NOT FOUND')

            # Load this one
            if (os.path.exists(chfile) is True) and (inside is True) and (chmeta['ngaiamatch'][j]>1):
                # Load the chip-level catalog
                cat1 = fits.getdata(chfile,1)
                ncat1 = len(cat1)
                #print('  chip '+str(chmeta[j]['ccdnum'])+'  '+str(ncat1)+' sources')

                # Make sure it's in the right format
                if len(cat1.dtype.fields) != 32:
                    if verbose: print('  This catalog does not have the right format. Skipping')
                    del cat1
                    ncat1 = 0

                # Only include sources inside Boundary+Buffer zone
                #  -use ROI_CUT
                #  -reproject to tangent plane first so we don't have to deal
                #     with RA=0 wrapping or pol issues
                if buffdict is not None:
                    lon, lat = coords.rotsphcen(cat1['ra'],cat1['dec'],buffdict['cenra'],buffdict['cendec'],gnomic=True)
                    ind_out, ind_in = dln.roi_cut(buffdict['lon'],buffdict['lat'],lon,lat)
                    nmatch = dln.size(ind_in)
                    # Only want source inside this pixel
                    if nmatch>0:
                        cat1 = cat1[ind_in]
                        ncat1 = len(cat1)
                    else:
                        cat1 = None
                        ncat1 = 0
                    #if verbose: print('  '+str(nmatch)+' sources are inside this pixel')

                # Combine the catalogs
                if ncat1 > 0:
                    # Keep it all in memory
                    if dbfile is None:
                        if cat is None:
                            #dtype_cat = cat1.dtype
                            #ncat_init = np.sum(chmeta['nsources'])*dln.size(metafile)
                            ncat_init = np.maximum(100000,ncat1)
                            cat = np.zeros(ncat_init,dtype=dtype_cat)
                            catcount = 0
                        # Add more elements if necessary
                        if (catcount+ncat1)>ncat:
                            cat = add_elements(cat,np.maximum(100000,ncat1))
                            ncat = len(cat)

                        # Add it to the main CAT catalog
                        for n in dtype_cat.names: cat[n][catcount:catcount+ncat1] = cat1[n.upper()]
                    # Use the database
                    else:
                        writecat2db(cat1,dbfile)

                    if verbose: print('  chip '+str(chmeta['ccdnum'][j])+'  '+str(ncat1)+' measurements')

                    catcount += ncat1
                    expcatcount += ncat1

        # Add metadata to ALLMETA, only if some measurements overlap
        if expcatcount>0:
            if allmeta is None:
                allmeta = newmeta
            else:
                allmeta = np.hstack((allmeta,newmeta))
        # Total measurements for this exposure
        print('  '+str(expcatcount)+' measurements')
        print(str(catcount)+' measurements total so far')

    #print('all exposures loaded. trimming now')
    if (cat is not None) & (catcount<ncat): cat=cat[0:catcount]   # delete excess elements
    if cat is None: cat=np.array([])         # empty cat
    if allmeta is None: allmeta=np.array([])

    print('loading measurements done after '+str(time.time()-t0))

    return cat, catcount, allmeta

def clusterdata(cat,ncat,dbfile=None):
    """ Perform spatial clustering """

    t00 = time.time()
    print('Spatial clustering with DBSCAN')    
    # Divide into subregions
    if (ncat>1000000) & (dbfile is not None):
        print('Dividing clustering problem into subregions')
        # Index RA and DEC
        createindexdb(dbfile,'ra',unique=False)
        createindexdb(dbfile,'dec',unique=False)
        # Subdivide
        nsub = int(np.ceil(ncat/100000))
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
        objstr = np.zeros(100000,dtype=np.dtype([('OBJLABEL',int),('RA',float),('DEC',float),('NMEAS',int)]))
        nobjstr = len(objstr)
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
                #import pdb; pdb.set_trace()
                cat1 = getdatadb(dbfile,rar=[r0-rabuff,r1+rabuff],decr=[d0-buff,d1+buff],verbose=True)
                ncat1 = len(cat1)
                print(str(ncat1)+' measurements')

                # Some measurements to work with
                if ncat1>0:
                    # Run DBSCAN
                    #import pdb; pdb.set_trace()
                    t0 = time.time()
                    X1 = np.column_stack((np.array(cat1['RA']),np.array(cat1['DEC'])))
                    dbs1 = DBSCAN(eps=0.5/3600, min_samples=1).fit(X1)
                    print('DBSCAN after '+str(time.time()-t0)+' sec.')
                    # Cluster labels are integers and in ascending order, but there are gaps
                    objlabels1 = dbs1.labels_
                    objlabels1 += lastobjlabel+1                 # add offset to labels
                    labelindex1 = dln.create_index(objlabels1)   # create inex
                    nobj1 = len(labelindex1['value'])
                    print(str(ncat1)+' measurements for '+str(nobj1)+' objects')
                    # Compute mean positions
                    obj1 = np.zeros(nobj1,dtype=np.dtype([('OBJLABEL',int),('RA',float),('DEC',float),('NMEAS',int)]))
                    obj1['OBJLABEL'] = labelindex1['value']
                    obj1['NMEAS'] = labelindex1['num']
                    for k in range(nobj1):
                        indx = labelindex1['index'][labelindex1['lo'][k]:labelindex1['hi'][k]+1]
                        wt_ra = 1.0/cat1['RAERR'][indx]**2
                        wt_dec = 1.0/cat1['DECERR'][indx]**2
                        obj1['RA'][k] = np.sum(cat1['RA'][indx]*wt_ra)/np.sum(wt_ra)
                        obj1['DEC'][k] = np.sum(cat1['DEC'][indx]*wt_dec)/np.sum(wt_dec)
                    # Only keep objects (and measurements) inside the box region
                    #  keep objects on LOWER boundary in RA/DEC
                    gdobj, ngdobj = dln.where((obj1['RA']>=r0) & (obj1['RA']<r1) & (obj1['DEC']>=d0) & (obj1['DEC']<d1))
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
                            add_rowid1[cnt1:cnt1+nmeas1] = cat1['ROWID'][indx]
                            add_objlabels1[cnt1:cnt1+nmeas1] = labelindex1['value'][gdobj[k]]
                            cnt1 += nmeas1

                        # Add the object labels into the database
                        #  much faster if in rowid order
                        si = np.argsort(add_rowid1)
                        insertobjlabelsdb(add_rowid1[si],add_objlabels1[si],dbfile)

                        # Add OBJ1 to OBJSTR
                        if (objcount+nobj1>nobjstr):    # add new elements
                            print('Adding more elements to OBSTR')
                            t1 = time.time()
                            objstr = add_elements(objstr,np.max([nobj1,100000]))
                            nobjstr = len(objstr)
                            print('more elements added in '+str(time.time()-t1)+' sec.')
                        objstr[objcount:objcount+nobj1] = obj1
                        objcount += nobj1

                        # Keep track of last label
                        lastobjlabel = np.max(obj1['OBJLABEL'])

                        #import pdb; pdb.set_trace()

        # Trim extra elements
        if nobjstr>objcount:
            objstr = objstr[0:objcount]

    # No subdividing
    else:
        # Get MEASID, RA, DEC from database
        if dbfile is not None:
            cat = getdbcoords(dbfile)
        # Spatially cluster the measurements with DBSCAN
        # coordinates of measurement
        X = np.column_stack((np.array(cat['RA']),np.array(cat['DEC'])))
        # Compute DBSCAN on all measurements
        dbs = DBSCAN(eps=0.5/3600, min_samples=1).fit(X)
        # Cluster labels are integers and in ascending order, but there are gaps
        objlabels = dbs.labels_
        labelindex = dln.create_index(objlabels)   # create inex
        nobj = len(labelindex['value'])
        print(str(ncat)+' measurements for '+str(nobj)+' objects')
        # Make structure
        objstr = np.zeros(nobj,dtype=np.dtype([('OBJLABEL',int),('NMEAS',int),('LO',int),('HI',int)]))
        objstr['OBJLABEL'] = labelindex['value']
        objstr['NMEAS'] = labelindex['num']
        nobjstr = len(objstr)
        # Insert object label into database
        if dbfile is not None:
            insertobjlabelsdb(cat['ROWID'],objlabels,dbfile)
        # Resort CAT, and use index LO/HI
        cat = cat[labelindex['index']]
        objstr['LO'] = labelindex['lo']
        objstr['HI'] = labelindex['hi']

    print(str(len(objstr))+' final objects')

    # Index objlabel in database
    if dbfile is not None:
        createindexdb(dbfile,'objlabel',unique=False)

    print('clustering done after '+str(time.time()-t00)+' sec.')

    return objstr, cat


# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix pixel number')
    parser.add_argument('version', type=str, nargs=1, help='Version number')
    parser.add_argument('--nside', type=int, default=128, help='HEALPix Nside')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    parser.add_argument('-v','--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--outdir', type=str, default='', help='Output directory')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    radeg = np.float64(180.00) / np.pi

    # Inputs
    pix = int(args.pix[0])
    version = args.version[0]
    verbose = args.verbose
    nside = args.nside
    redo = args.redo
    outdir = args.outdir

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

    # Check if output file already exists
    if outdir == '': outdir=dir+'combine/'
    subdir = str(int(pix)//1000)    # use the thousands to create subdirectory grouping
    outfile = outdir+'/'+subdir+'/'+str(pix)+'.fits'
    if (os.path.exists(outfile) or os.path.exists(outfile+'.gz')) & (not redo):
        print(outfile+' EXISTS already and REDO not set')
        sys.exit()

    print("Combining InstCal SExtractor catalogs for Healpix pixel = "+str(pix))

    # Load the list
    listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_instcal_combine_healpix_list.fits.gz'
    if os.path.exists(listfile) is False:
        print(listfile+" NOT FOUND")
        sys.exit()
    healstr = Table(fits.getdata(listfile,1))
    index = Table(fits.getdata(listfile,2))
    # Find our pixel
    ind,nind = dln.where(index['PIX'] == pix)
    if nind == 0:
        print("No entries for Healpix pixel '"+str(pix)+"' in the list")
        sys.exit()
    ind = ind[0]
    hlist = healstr[index['LO'][ind]:index['HI'][ind]+1]
    nlist = len(hlist)
    # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
    #  so we can deal with the edge cases
    neipix = hp.get_all_neighbours(nside,pix)
    for neip in neipix:
        ind1,nind1 = dln.where(index['PIX'] == neip)
        if nind1>0:
            ind1 = ind1[0]
            hlist1 = healstr[index[ind1]['LO']:index[ind1]['HI']+1]
            hlist = vstack([hlist,hlist1])

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
    buffdict = {'cenra':cenra,'cendec':cendec,'rar':dln.minmax(rabuff),'decr':dln.minmax(decbuff),'ra':rabuff,'dec':decbuff,\
                'lon':lonbuff,'lat':latbuff,'lr':dln.minmax(lonbuff),'br':dln.minmax(latbuff)}

    
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

    # Decide whether to load everything into RAM or use temporary database
    metafiles = [m.replace('_cat','_meta').strip() for m in hlist['FILE']]
    nmeasperchip = np.zeros(dln.size(metafiles),int)
    for i,m in enumerate(metafiles):
        expstr = fits.getdata(m,1)
        nmeasperchip[i] = expstr['NMEAS']/expstr['NCHIPS']
    totmeasest = np.sum(nmeasperchip)
    usedb = False
    if totmeasest>500000: usedb=True
    dbfile = None
    if usedb:
        dbfile = tmproot+str(pix)+'_combine.db'
        print('Using temporary database file = '+dbfile)
        if os.path.exists(dbfile): os.remove(dbfile)

    # Load the measurement catalog
    #  this will contain excess rows at the end, if all in RAM
    #  if using database, CAT is empty
    cat, catcount, allmeta = loadmeas(metafiles,buffdict,dbfile=dbfile)
    ncat = catcount
    print(str(ncat))

    # Spatially cluster the measurements with DBSCAN
    #   this might also resort CAT
    objstr, cat = clusterdata(cat,ncat,dbfile=dbfile)
    nobj = dln.size(objstr)
    meascumcount = np.cumsum(objstr['NMEAS'])
    print(str(nobj)+' unique objects clustered within 0.5 arcsec')

    # Initialize the OBJ structured array
    obj = np.zeros(nobj,dtype=dtype_obj)
    obj['objectid'] = dln.strjoin( str(pix)+'.', ((np.arange(nobj)+1).astype(np.str)) )
    obj['pix'] = pix
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
    idstr = np.zeros(ncat,dtype=dtype_idstr)

    # Higher precision catalog
    dtype_hicat = np.dtype([('MEASID',np.str,30),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
                            ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
                            ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
                            ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])

    # Convert to nump structured array
    dtype_hicatdb = np.dtype([('MEASID',np.str,30),('OBJLABEL',int),('EXPOSURE',np.str,40),('CCDNUM',int),('FILTER',np.str,3),
                              ('MJD',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
                              ('MAG_AUTO',float),('MAGERR_AUTO',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),
                              ('THETA',float),('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])

    t1 = time.time()

    # Loop over the objects
    meascount = 0
    ngroup = -1
    grpcount = 0
    maxmeasload = 50000
    fidmag = np.zeros(nobj,float)+np.nan  # fiducial magnitude
    for i,lab in enumerate(objstr['OBJLABEL']):
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
                # Use maxmeaslead to figure out how many objects we can load
                if i==0:
                    ngroup = np.max(np.where(meascumcount[i:]<=maxmeasload)[0])+1
                else:
                    ngroup = np.max(np.where((meascumcount[i:]-meascumcount[i-1])<=maxmeasload)[0])+1
                ngroup = np.max([1,ngroup])   # need to load at least 1
                lab0 = lab
                lab1 = objstr['OBJLABEL'][np.min([i+ngroup-1,nobj-1])]
                #lab1 = labelindex['value'][np.min([i+ngroup-1,nobj-1])]
                grpcat = getdatadb(dbfile,objlabel=[lab0,lab1])
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

        obj['ndet'][i] = ncat1

        # Add in IDSTR information
        idstr['measid'][oindx] = cat1['MEASID']
        idstr['exposure'][oindx] = cat1['EXPOSURE']
        idstr['objectid'][oindx] = obj['objectid'][i]
        idstr['objectindex'][oindx] = i

        # Computing quantities
        # Mean RA/DEC, RAERR/DECERR
        if ncat1>1:
            wt_ra = 1.0/cat1['RAERR']**2
            wt_dec = 1.0/cat1['DECERR']**2
            obj['ra'][i] = np.sum(cat1['RA']*wt_ra)/np.sum(wt_ra)
            obj['raerr'][i] = np.sqrt(1.0/np.sum(wt_ra))
            obj['dec'][i] = np.sum(cat1['DEC']*wt_dec)/np.sum(wt_dec)
            obj['decerr'][i] = np.sqrt(1.0/np.sum(wt_dec))
            obj['mjd'][i] = np.mean(cat1['MJD'])
            obj['deltamjd'][i] = np.max(cat1['MJD'])-np.min(cat1['MJD'])
        else:
            obj['ra'][i] = cat1['RA']
            obj['dec'][i] = cat1['DEC']
            obj['raerr'][i] = cat1['RAERR']
            obj['decerr'][i] = cat1['DECERR']
            obj['mjd'][i] = cat1['MJD']
            obj['deltamjd'][i] = 0

        # Mean proper motion and errors
        if ncat1>1:
            raerr = np.array(cat1['RAERR']*1e3,np.float64)    # milli arcsec
            ra = np.array(cat1['RA'],np.float64)
            ra -= np.mean(ra)
            ra *= 3600*1e3 * np.cos(obj['dec'][i]/radeg)     # convert to true angle, milli arcsec
            t = cat1['MJD'].copy()
            t -= np.mean(t)
            t /= 365.2425                          # convert to year
            # Calculate robust slope
            pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
            obj['pmra'][i] = pmra                 # mas/yr
            obj['pmraerr'][i] = pmraerr           # mas/yr

            decerr = np.array(cat1['DECERR']*1e3,np.float64)   # milli arcsec
            dec = np.array(cat1['DEC'],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            # Calculate robust slope
            pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
            obj['pmdec'][i] = pmdec               # mas/yr
            obj['pmdecerr'][i] = pmdecerr         # mas/yr

        # Mean magnitudes
        # Convert totalwt and totalfluxwt to MAG and ERR
        #  and average the morphology parameters PER FILTER
        filtindex = dln.create_index(cat1['FILTER'].astype(np.str))
        nfilters = len(filtindex['value'])
        resid = np.zeros(ncat1)+np.nan     # residual mag
        relresid = np.zeros(ncat1)+np.nan  # residual mag relative to the uncertainty
        for f in range(nfilters):
            filt = filtindex['value'][f].lower()
            findx = filtindex['index'][filtindex['lo'][f]:filtindex['hi'][f]+1]
            obj['ndet'+filt][i] = filtindex['num'][f]
            gph,ngph = dln.where(cat1['MAG_AUTO'][findx]<50)
            obj['nphot'+filt][i] = ngph
            if ngph==1:
                obj[filt+'mag'][i] = cat1['MAG_AUTO'][findx[gph]]
                obj[filt+'err'][i] = cat1['MAGERR_AUTO'][findx[gph]]
            if ngph>1:
                newmag, newerr = dln.wtmean(cat1['MAG_AUTO'][findx[gph]], cat1['MAGERR_AUTO'][findx[gph]],magnitude=True,reweight=True,error=True)
                obj[filt+'mag'][i] = newmag
                obj[filt+'err'][i] = newerr
                # Calculate RMS
                obj[filt+'rms'][i] = np.sqrt(np.mean((cat1['MAG_AUTO'][findx[gph]]-newmag)**2))
                # Residual mag
                resid[findx[gph]] = cat1['MAG_AUTO'][findx[gph]]-newmag
                # Residual mag relative to the uncertainty
                #  set a lower threshold of 0.02 in the uncertainty
                relresid[findx[gph]] = np.sqrt(ngph/(ngph-1)) * (cat1['MAG_AUTO'][findx[gph]]-newmag)/np.maximum(cat1['MAGERR_AUTO'][findx[gph]],0.02)

            # Calculate mean morphology parameters
            obj[filt+'asemi'][i] = np.mean(cat1['ASEMI'][findx])
            obj[filt+'bsemi'][i] = np.mean(cat1['BSEMI'][findx])
            obj[filt+'theta'][i] = np.mean(cat1['THETA'][findx])

        # Calculate variability indices
        gdresid = np.isfinite(resid)
        ngdresid = np.sum(gdresid)
        if ngdresid>0:
            resid2 = resid[gdresid]
            sumresidsq = np.sum(resid2**2)
            tsi = np.argsort(cat1['MJD'][gdresid])
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
        obj['asemi'][i] = np.mean(cat1['ASEMI'])
        obj['bsemi'][i] = np.mean(cat1['BSEMI'])
        obj['theta'][i] = np.mean(cat1['THETA'])
        obj['asemierr'][i] = np.sqrt(np.sum(cat1['ASEMIERR']**2)) / ncat1
        obj['bsemierr'][i] = np.sqrt(np.sum(cat1['BSEMIERR']**2)) / ncat1
        obj['thetaerr'][i] = np.sqrt(np.sum(cat1['THETAERR']**2)) / ncat1
        obj['fwhm'][i] = np.mean(cat1['FWHM'])
        obj['class_star'][i] = np.mean(cat1['CLASS_STAR'])
        obj['flags'][i] = np.bitwise_or.reduce(cat1['FLAGS'])  # OR combine

    # Select Variables
    #  1) Construct fiducial magnitude (done in loop above)
    #  2) Construct median VAR and sigma VAR versus magnitude
    #  3) Find objects that Nsigma above the median VAR line
    si = np.argsort(fidmag)   # NaNs are at end
    varcol = 'madvar'
    gdvar,ngdvar,bdvar,nbdvar = dln.where(np.isfinite(obj[varcol]) & np.isfinite(fidmag),comp=True)
    nbins = np.ceil((np.max(fidmag[gdvar])-np.min(fidmag[gdvar]))/0.25)
    nbins = int(np.max([2,nbins]))
    if ngdvar>0:
        fidmagmed, bin_edges1, binnumber1 = bindata.binned_statistic(fidmag[gdvar],fidmag[gdvar],statistic='nanmedian',bins=nbins)
        # Median metric
        varmed, bin_edges2, binnumber2 = bindata.binned_statistic(fidmag[gdvar],obj[varcol][gdvar],statistic='nanmedian',bins=nbins)
        # Smooth, it handles NaNs well
        smlen = 5
        smvarmed = dln.gsmooth(varmed,smlen)
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


    # ONLY INCLUDE OBJECTS WITH AVERAGE RA/DEC
    # WITHIN THE BOUNDARY OF THE HEALPIX PIXEL!!!
    ipring = hp.pixelfunc.ang2pix(nside,obj['ra'],obj['dec'],lonlat=True)
    ind1,nmatch = dln.where(ipring == pix)
    if nmatch==0:
        print('None of the final objects fall inside the pixel')
        sys.exit()
    # Get trimmed objects and indices
    objtokeep = np.zeros(nobj,bool)         # boolean to keep or trim objects
    objtokeep[ind1] = True
    if nmatch<nobj:
        trimind = np.arange(nobj)
        trimind = np.delete(trimind,ind1)
        #trimind = dln.remove_indices(trimind,ind1)
        trimobj = obj[trimind]          # trimmed objects
    newobjindex = np.zeros(nobj,int)-1    # new indices
    newobjindex[ind1] = np.arange(nmatch)
    # Keep the objects inside the Healpix
    obj = obj[ind1]
    print(str(nmatch)+' final objects fall inside the pixel')

    # Remove trimmed objects from IDSTR
    totrim,ntotrim = dln.where(~objtokeep[idstr['objectindex']])  #using old index
    if ntotrim>0:
        # Trim objects
        idstr = np.delete(idstr,totrim)
        #idstr = dln.remove_indices(idstr,totrim)
        # Update IDSTR.objectindex
        old_idstr_objectindex = idstr['objectindex']
        idstr['objectindex'] = newobjindex[old_idstr_objectindex]

    # Create final summary structure from ALLMETA
    #  get exposures that are in IDSTR
    #  sometimes EXPNUM numbers have the leading 0s removed
    #  and sometimes not, so turn to LONG to match
    dum, uiexposure = np.unique(idstr['exposure'],return_index=True)
    uexposure = idstr['exposure'][uiexposure]
    nuexposure = len(uexposure)
    ind1,ind2 = dln.match(allmeta['base'],uexposure)
    nmatch = len(ind1)
    sumstr = Table(allmeta[ind1])
    col_nobj = Column(name='nobjects', dtype=np.int, length=len(sumstr))
    col_healpix = Column(name='healpix', dtype=np.int, length=len(sumstr))
    sumstr.add_columns([col_nobj, col_healpix])
    sumstr['nobjects'] = 0
    sumstr['healpix'] = pix
    # get number of objects per exposure
    exposure = idstr['exposure']
    siexp = np.argsort(exposure)
    exposure = exposure[siexp]
    if nuexposure>1:
        brklo,nbrk = dln.where(exposure != np.roll(exposure,1))
        brkhi = np.hstack((brklo[1:nbrk],len(exposure)))
        numobjexp = brkhi-brklo+1
    else:
        numobjexp=len(exposure)
    ind1,ind2 = dln.match(sumstr['base'],uexposure)
    nmatch = len(ind1)
    sumstr['nobjects'][ind1] = numobjexp

    # Write the output file
    print('Writing combined catalog to '+outfile)
    if os.path.exists(outdir) is False: os.mkdir(outdir)
    if os.path.exists(outdir+'/'+subdir) is False: os.mkdir(outdir+'/'+subdir)
    if os.path.exists(outfile): os.remove(outfile)
    sumstr.write(outfile)               # first, summary table
    #  append other fits binary tables
    hdulist = fits.open(outfile)
    hdu = fits.table_to_hdu(Table(obj))        # second, catalog
    hdulist.append(hdu)
    hdu = fits.table_to_hdu(Table(idstr))      # third, ID table
    hdulist.append(hdu)    
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
    del idstr
