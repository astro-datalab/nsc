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
from glob import glob

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
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix pixel number')
    args = parser.parse_args()

    parentpix = args.pix[0]
    version = 'v3'
    nside = 128
    t0 = time.time()

    # Output filename
    outdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/combine/'
    outbase = str(parentpix)
    subdir = str(int(parentpix)//1000)    # use the thousands to create subdirectory grouping    
    if os.path.exists(outdir+'/'+subdir) is False: os.mkdir(outdir+'/'+subdir)
    outfile = outdir+'/'+subdir+'/'+str(parentpix)+'.fits'

    # Get higher-resolution object filesnames
    outfiles = glob('/net/dl2/dnidever/nsc/instcal/v3/combine/'+str(int(parentpix)//1000)+'/'+str(parentpix)+'_n*.fits.gz')
    outfiles = np.array(outfiles)
    si = np.argsort(outfiles)
    outfiles = outfiles[si]
    base = [os.path.basename(f)[:-8] for f in outfiles]
    hinside = int(base[0].split('_')[1][1:])
    allpix = [b.split('_')[2] for b in base]


    # Load and concatenate all of the files                                                                                                                                
    print('Combining all of the object catalogs for '+parentpix)
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
        dbfile_idstr1 = outfile1.replace('.fits.gz','_idstr.db')
        objectid_orig = obj1['objectid']
        objectid_new = dln.strjoin( str(parentpix)+'.', ((np.arange(nobj1)+1+totobjects).astype(np.str)) )
        # only update if they are different
        diff, = np.where(objectid_orig != objectid_new)
        if len(diff)>0:
            print('Updating objectid in '+dbfile_idstr1)
            updatecoldb('objectid',objectid_orig,'objectid',objectid_new,'idstr',dbfile_idstr1)
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

    print('Breaking-up IDSTR information')
    dbfiles_idstr = []
    for i in range(len(allpix)):
        outfile1 = outfiles[i]
        dbfile_idstr1 = outfile1.replace('.fits.gz','_idstr.db')
        dbfiles_idstr.append(dbfile_idstr1)
    breakup_idstr(dbfiles_idstr)

    sys.exit()
