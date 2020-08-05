#!/usr/bin/env python

# Update measurement catalogs using the broken up measid/objectid lists
# from nsc_instcal_combine_breakup_idstr.py

import os
import sys
import numpy as np
import shutil
import time
from dlnpyutils import utils as dln, coords, db
from astropy.table import Table
from astropy.io import fits
import sqlite3
import socket
from argparse import ArgumentParser
import logging
from glob import glob
import subprocess
import healpy as hp
#import tempfile
import psycopg2 as pq
import psutil

def get_meas(pix,nside=128):
    """ Get the measurements for a particular healpix."""
    # objid, ra, raerr, dec, decerr, mjd

    t0 = time.time()

    connection = pq.connect(user="dlquery",host="db01.datalab.noao.edu",
                            password="",port = "5432",database = "tapdb")
    cur = connection.cursor()

    if nside==128:
        cmd = """SELECT m.objectid,m.ra,m.raerr,m.dec,m.decerr,m.mjd from nsc_dr2.meas as m join
             nsc_dr2.object as obj on m.objectid=obj.objectid where obj.pix={0};""".format(pix)
    if nside==256:
        cmd = """SELECT m.objectid,m.ra,m.raerr,m.dec,m.decerr,m.mjd from nsc_dr2.meas as m join
             nsc_dr2.object as obj on m.objectid=obj.objectid where obj.ring256={0};""".format(pix)

    if (nside==128) | (nside==256):
        cur.execute(cmd)
        data = cur.fetchall()
        # Convert to numpy structured array
        dtype = np.dtype([('objectid',np.str,50),('ra',np.float64),('raerr',float),
                          ('dec',np.float64),('decerr',float),('mjd',np.float64)])
        meas = np.zeros(len(data),dtype=dtype)
        meas[...] = data
        del(data)

    # nside>256
    else:
        ra,dec = hp.pix2ang(nside,pix,lonlat=True)
        radius = hp.nside2resol(nside,arcmin=True)/60.*1.5
        # First get object ra/dec and figure out the radius
        #cmd = """SELECT objectid,ra,dec from nsc_dr2.object where
        #         q3c_radial_query(ra,dec,{0},{1},{2});""".format(ra,dec,radius)
        #cur.execute(cmd)
        #data = cur.fetchall()
        ## Convert to numpy structured array
        #dtype = np.dtype([('objectid',np.str,50),('ra',np.float64),('dec',np.float64)])
        #obj = np.zeros(len(data),dtype=dtype)
        #obj[...] = data
        #del(data)

        # https://github.com/segasai/q3c
        # The polygonal query, i.e. the query of the objects which lie inside the region bounded by the polygon on the sphere.
        # To query the objects in the polygon ((0,0),(2,0),(2,1),(0,1)) ) (this is the spherical polygon with following vertices:
        #  (ra=0, dec=0) ; (ra=2, dec=0); (ra=2, dec=1); (ra=0, dec=1)):
        # my_db# SELECT * FROM mytable WHERE q3c_poly_query(ra, dec, '{0, 0, 2, 0, 2, 1, 0, 1}');
        
        vecbound = hp.boundaries(nside,pix)
        rabound, decbound = hp.vec2ang(np.transpose(vecbound),lonlat=True)
        # Expand the boundary by the buffer size
        cenra, cendec = hp.pix2ang(nside,pix,lonlat=True)
        # reproject onto tangent plane
        lonbound, latbound = coords.rotsphcen(rabound,decbound,cenra,cendec,gnomic=True)
        # expand by a fraction, it's not an exact boundary but good enough
        frac = 1.05
        lonbuff = lonbound*frac
        latbuff = latbound*frac
        rabuff, decbuff = coords.rotsphcen(lonbuff,latbuff,cenra,cendec,gnomic=True,reverse=True)
        if (np.max(rabuff)-np.min(rabuff))>100:  # deal with RA=0 wraparound
            bd,nbd = dln.where(rabuff>180)
            if nbd>0:rabuff[bd] -=360.0
        bnd = (rabuff[0],decbuff[0], rabuff[1],decbuff[1], rabuff[2],decbuff[2], rabuff[3],decbuff[3])

        cmd = "SELECT m.objectid,m.ra,m.raerr,m.dec,m.decerr,m.mjd from nsc_dr2.meas as m join "
        cmd += "nsc_dr2.object as obj on m.objectid=obj.objectid where "
        cmd += "q3c_poly_query(obj.ra,obj.dec,'{%10.6f,%10.6f, %10.6f,%10.6f, %10.6f,%10.6f, %10.6f,%10.6f}'::double precision[]);" % (bnd)

        print(bnd)

        #cmd = """SELECT m.objectid,m.ra,m.raerr,m.dec,m.decerr,m.mjd from nsc_dr2.meas as m join
        #     nsc_dr2.object as obj on m.objectid=obj.objectid where q3c_radial_query(obj.ra,obj.dec,{0},{1},{2});""".format(ra,dec,radius)
        #print("""RA={0} DEC={1} RADIUS={2}""".format(ra,dec,radius))
        
        cur.execute(cmd)
        data = cur.fetchall()
        # Convert to numpy structured array
        dtype = np.dtype([('objectid',np.str,50),('ra',np.float64),('raerr',float),
                          ('dec',np.float64),('decerr',float),('mjd',np.float64)])
        meas = np.zeros(len(data),dtype=dtype)
        meas[...] = data
        del(data)

    cur.close()
    connection.close()

    dt = time.time()-t0
    print('Retrieved '+str(len(meas))+' measurements in '+str(dt)+' seconds')

    return meas


def fix_pms(pix):
    """ Correct the proper motions in the healpix object catalog."""

    t00 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    version = 'v3'
    nside = 128
    radeg = np.float64(180.00) / np.pi

    hdir = '/net/dl2/dnidever/nsc/instcal/'+version+'/combine/'+str(int(pix)//1000)+'/'
    objfile = hdir+str(pix)+'.fits.gz'
    outfile = hdir+str(pix)+'_pmcorr.fits'
    
    print('Correcting proper motions for '+str(pix))

    # Check that the object file exists
    if os.path.exists(objfile) is False:
        print(objfile+' NOT FOUND')
        return

    # Check fixed file  
    if os.path.exists(outfile+'.gz') == True:
        print(str(pix)+' already fixed')
        return

    # Load the object file
    #meta = fits.getdata(objfile,1)
    #obj = fits.getdata(objfile,2)
    meta = Table.read(objfile,1)
    obj = Table.read(objfile,2)
    nobj = len(obj)
    print(str(nobj)+' objects with '+str(np.sum(obj['ndet']))+' measurements')
    print('KLUDGE!!! MAKING COPY OF OBJ!!!')
    orig = obj.copy()

    v = psutil.virtual_memory()
    process = psutil.Process(os.getpid())
    print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

    # Break up into subregions
    totmeas = np.sum(obj['ndet'])
    nsub,bestind = dln.closest([1,4,16,64],int(np.ceil(totmeas/500000)))
    hinside = [128,256,512,1024][bestind]
    vecbound = hp.boundaries(nside,int(pix))
    allpix = hp.query_polygon(hinside,np.transpose(vecbound))
    allra,alldec = hp.pix2ang(hinside,allpix,lonlat=True)
    print(str(nsub)+' sub regions')

    # Get the objects within this subpixel
    objpix = hp.ang2pix(hinside,obj['ra'],obj['dec'],lonlat=True)

    ndet = np.zeros(nobj,int)
    #allpmra_old = np.zeros(nobj,float)
    #allpmdec_old = np.zeros(nobj,float)
    #allpmra_linefit = np.zeros(nobj,float)

    # Loop over subpixels
    for i in range(nsub):
        pix1 = allpix[i]
        print(str(i+1)+' '+str(pix1))

        # Get the measurements
        meas = get_meas(pix1,nside=hinside)
        nmeas = len(meas)

        v = psutil.virtual_memory()
        process = psutil.Process(os.getpid())
        print('%6.1f Percent of memory used. %6.1f GB available.  Process is using %6.2f GB of memory.' % (v.percent,v.available/1e9,process.memory_info()[0]/1e9))

        # Get the objects within this subpixel
        objind, = np.where(objpix==pix1)
        obj1 = obj[objind]
        nobj1 = len(obj1)
        print('  '+str(nobj1)+' objects in this subregion')

        idindex = dln.create_index(meas['objectid'])
        ## Not all matched
        #if len(idindex['value']) != nobj:
        #    print('Number of unique OBJECTIDs in object and meas catalogs do not match')
        #    return
        ind1,ind2 = dln.match(obj1['objectid'],idindex['value'])
        # Not all matched
        if len(ind1) != nobj1:
            print(str(len(obj1))+' objects in this sub healpix but only measurements for '+str(len(ind1)))
            #print('Some objects are missing measurements')
            #return
        # sort by object index
        si = np.argsort(ind1)
        ind1 = ind1[si]
        ind2 = ind2[si]

        # Loop over
        ndet1 = np.zeros(nobj1,int)
        #allpmra_old1 = np.zeros(nobj1,float)
        #allpmdec_old1 = np.zeros(nobj1,float)
        #allpmra_linefit1 = np.zeros(nobj1,float)
        for j in range(len(ind1)):
            if (j % 1000)==0: print('  '+str(j))
            k = ind1[j]  # object index
            # Calculate the proper motions
            mind = idindex['index'][idindex['lo'][ind2[j]]:idindex['hi'][ind2[j]]+1]
            cat1 = meas[mind]
            ncat1 = len(cat1)
            ndet1[k] = ncat1
            if ncat1>1:
                raerr = np.array(cat1['raerr']*1e3,np.float64)    # milli arcsec
                ra = np.array(cat1['ra'],np.float64)
                ra -= np.mean(ra)
                ra *= 3600*1e3 * np.cos(obj1['dec'][k]/radeg)     # convert to true angle, milli arcsec
                t = cat1['mjd'].copy()
                t -= np.mean(t)
                t /= 365.2425                          # convert to year
                # Calculate robust slope
                try:
                    pmra, pmraerr = dln.robust_slope(t,ra,raerr,reweight=True)
                    #pmra_old, pmraerr_old = dln.robust_slope_old(t,ra,raerr,reweight=True)
                    #pmra_linefit = dln.poly_fit(t,ra,2,robust=True,sigma=raerr,initpar=pmra)
                except:
                    print('problem')
                    import pdb; pdb.set_trace()
                obj1['pmra'][k] = pmra                 # mas/yr
                obj1['pmraerr'][k] = pmraerr           # mas/yr
                #allpmra_old1[k] = pmra_old
                #allpmra_linefit1[k] = pmra_linefit

                decerr = np.array(cat1['decerr']*1e3,np.float64)   # milli arcsec
                dec = np.array(cat1['dec'],np.float64)
                dec -= np.mean(dec)
                dec *= 3600*1e3                         # convert to milli arcsec
                # Calculate robust slope
                try:
                    pmdec, pmdecerr = dln.robust_slope(t,dec,decerr,reweight=True)
                    #pmdec_old, pmdecerr_old = dln.robust_slope_old(t,dec,decerr,reweight=True)
                except:
                    print('problem')
                    import pdb; pdb.set_trace()
                obj1['pmdec'][k] = pmdec               # mas/yr
                obj1['pmdecerr'][k] = pmdecerr         # mas/yr
                #allpmdec_old1[k] = pmdec_old

        # Stuff subregion object back into big one
        obj[objind] = obj1
        ndet[objind] = ndet1
        #allpmra_old[objind] = allpmra_old1
        #allpmdec_old[objind] = allpmdec_old1
        #allpmra_linefit[objind] = allpmra_linefit1

        #import pdb; pdb.set_trace()

    #np.save(hdir+str(pix)+'_pmraold.npy',allpmra_old)
    #np.save(hdir+str(pix)+'_pmdecold.npy',allpmdec_old)
    #np.save(hdir+str(pix)+'_pmralinefit.npy',allpmra_linefit)

    #import pdb; pdb.set_trace()


    # Save the new version of obj
    # Write the output file
    print('Writing combined catalog to '+outfile)
    if os.path.exists(outfile): os.remove(outfile)
    #Table(meta).write(outfile)               # first, summary table
    meta.write(outfile)               # first, summary table
    #  append other fits binary tables
    hdulist = fits.open(outfile)
    #hdu = fits.table_to_hdu(Table(obj))        # second, catalog
    hdu = fits.table_to_hdu(obj)        # second, catalog
    hdulist.append(hdu)
    hdulist.writeto(outfile,overwrite=True)
    hdulist.close()
    if os.path.exists(outfile+'.gz'): os.remove(outfile+'.gz')
    ret = subprocess.call(['gzip',outfile])    # compress final catalog

    print('dt = %6.1f sec.' % (time.time()-t00))


if __name__ == "__main__":
    parser = ArgumentParser(description='Fix pms in healpix object catalogs.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix')
    args = parser.parse_args()

    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    pix = args.pix[0]

    # Input is a list
    if pix[0]=='@':
        listfile = pix[1:]
        if os.path.exists(listfile): 
            pix = dln.readlines(listfile)
        else:
            print(listfile+' NOT FOUND')
            sys.exit()

    # Fix the pms in healpix object catalogs
    if type(pix) is not list: pix=[pix]
    npix = len(pix)
    print('Correcting PMs for '+str(npix)+' HEALPix')
    for i in range(len(pix)):
        fix_pms(pix[i])


