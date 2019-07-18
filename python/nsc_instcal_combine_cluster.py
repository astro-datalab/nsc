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
from dlnpyutils import utils, coords
import subprocess
import time
from argparse import ArgumentParser
import socket
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
from sklearn.cluster import DBSCAN
from scipy.optimize import least_squares

def add_elements(cat,nnew=300000):
    """ Add more elements to a catalog"""
    ncat = len(cat)
    old = cat.copy()
    nnew = utils.gt(nnew,ncat)
    cat = np.zeros(ncat+nnew,dtype=old.dtype)
    cat[0:ncat] = old
    del(old)
    return cat
    
def add_cat(obj,totobj,idstr,idcnt,ind1,cat,meta):
    """ Add object information from a new meas catalog of matched objects"""

    ncat = len(cat)
    f = meta['filter'].lower().strip()[0]
    # Copy to final structure
    obj['ra'][ind1] = cat['RA']
    obj['dec'][ind1] = cat['DEC']
    obj['raerr'][ind1] += 1.0/cat['RAERR']**2                           # total(ra_wt)
    obj['decerr'][ind1] += 1.0/cat['DECERR']**2                         # total(dec_wt)
    obj['pmra'][ind1] += (1.0/cat['RAERR']**2) * meta['mjd']*cat['RA']     # total(wt*mjd*ra)
    obj['pmdec'][ind1] += (1.0/cat['DECERR']**2) * meta['mjd']*cat['DEC']  # total(wt*mjd*dec)
    obj['mjd'][ind1] += meta['mjd']                                 # total(mjd)
    obj['ndet'][ind1] += 1
    # Detection and morphology parameters for this FILTER
    obj['ndet'+f][ind1] += 1
    obj[f+'asemi'][ind1] += cat['ASEMI']
    obj[f+'bsemi'][ind1] += cat['BSEMI']
    obj[f+'theta'][ind1] += cat['THETA']
    # Good photometry for this FILTER
    gdmag, = np.where(cat['MAG_AUTO']<50)
    if len(gdmag)>0:
      obj[f+'mag'][ind1[gdmag]] += 2.5118864**cat['MAG_AUTO'][gdmag] * (1.0/cat['MAGERR_AUTO'][gdmag]**2)
      obj[f+'err'][ind1[gdmag]] += 1.0/cat['MAGERR_AUTO'][gdmag]**2
      obj['nphot'+f][ind1[gdmag]] += 1
    obj['asemi'][ind1] += cat['ASEMI']
    obj['asemierr'][ind1] += cat['ASEMIERR']**2
    obj['bsemi'][ind1] += cat['BSEMI']
    obj['bsemierr'][ind1] += cat['BSEMIERR']**2
    obj['theta'][ind1] += cat['THETA']
    obj['thetaerr'][ind1] += cat['THETAERR']**2
    obj['fwhm'][ind1] += cat['FWHM']  # in arcsec
    obj['flags'][ind1] += cat['FLAGS']
    obj['class_star'][ind1] += cat['CLASS_STAR']
    totobj['ra'][ind1] += cat['RA'] * (1.0/cat['RAERR']**2)             # total(ra*wt)
    totobj['dec'][ind1] += cat['DEC'] * (1.0/cat['DECERR']**2)          # total(dec*wt)
    totobj['ramjd'][ind1] += (1.0/cat['RAERR']**2) * meta['mjd']        # total(wt_ra*mjd)
    totobj['decmjd'][ind1] += (1.0/cat['DECERR']**2) * meta['mjd']      # total(wt_dec*mjd)
    totobj['ramjd2'][ind1] += (1.0/cat['RAERR']**2) * meta['mjd']**2    # total(wt_ra*mjd**2)
    totobj['decmjd2'][ind1] += (1.0/cat['DECERR']**2) * meta['mjd']**2  # total(wt_dec*mjd**2)
    totobj['minmjd'][ind1] = np.minimum( meta['mjd'][0], totobj['minmjd'][ind1] )
    totobj['maxmjd'][ind1] = np.maximum( meta['mjd'][0], totobj['maxmjd'][ind1] )
    if len(gdmag)>0:
        totobj[f+'tot'][ind1[gdmag]] += cat['MAG_AUTO'][gdmag]       # sum(mag)
        totobj[f+'mag2'][ind1[gdmag]] += np.float64(cat['MAG_AUTO'][gdmag])**2   # sum(mag**2), need dbl to precent underflow

    # Add new elements to IDSTR
    if idcnt+ncat > len(idstr):
        idstr = add_elements(idstr)
        nidstr = len(idstr)

    # Add to IDSTR
    idstr['measid'][idcnt:idcnt+ncat] = cat['MEASID']
    idstr['exposure'][idcnt:idcnt+ncat] = meta['base']
    idstr['expnum'][idcnt:idcnt+ncat] = meta['expnum']
    idstr['objectid'][idcnt:idcnt+ncat] = obj[ind1]['objectid']
    idstr['objectindex'][idcnt:idcnt+ncat] = ind1
    idcnt += ncat

    return obj,totobj,idstr,idcnt
    

def loadmeas(metafile=None,buffdict=None,verbose=False):

    if metafile is None:
        print('Need metafile')
        return np.array([]), np.array([])

    # New meta-data format
    dtype_meta = np.dtype([('file',np.str,500),('base',np.str,200),('expnum',int),('ra',np.float64),
                           ('dec',np.float64),('dateobs',np.str,100),('mjd',np.float64),('filter',np.str,50),
                           ('exptime',float),('airmass',float),('nsources',int),('fwhm',float),
                           ('nchips',int),('badchip31',bool),('rarms',float),('decrms',float),
                           ('ebv',float),('gaianmatch',int),('zpterm',float),('zptermerr',float),
                           ('zptermsig',float),('refmatch',int),('nmeas',int),('catindex',int)])

    #  Loop over exposures
    cat = None
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
        # Add index in CAT
        newmeta['catindex'] = catcount

        # Get the name
        fdir = os.path.dirname(mfile)
        fbase, ext = os.path.splitext(os.path.basename(mfile))
        fbase = fbase[:-5]   # remove _meta at end
        
        # Loop over the chip files
        for j in range(len(chmeta)):
            # Check that this chip was astrometrically calibrated
            #   and falls in to HEALPix region
            if chmeta[j]['ngaiamatch'] == 0:
                if verbose: print('This chip was not astrometrically calibrate')

            # Check that this overlaps the healpix region
            inside = True
            if buffdict is not None:
                vra = chmeta[j]['vra']
                vdec = chmeta[j]['vdec']
                if (np.max(vra)-np.min(vra)) > 100:    # deal with RA=0 wrapround
                    bd, = np.where(vra>180)
                    if len(bd)>0: vra[bd] -= 360
                if coords.doPolygonsOverlap(buffdict['ra'],buffdict['dec'],vra,vdec) is False:
                    if verbose: print('This chip does NOT overlap the HEALPix region+buffer')
                    inside = False

            # Check if the chip-level file exists
            chfile = fdir+'/'+fbase+'_'+str(chmeta[j]['ccdnum'])+'_meas.fits'
            if os.path.exists(chfile) is False:
                print(chfile+' NOT FOUND')

            # Load this one
            if (os.path.exists(chfile) is True) and (inside is True) and (chmeta[j]['ngaiamatch']>1):
                # Load the chip-level catalog
                cat1 = fits.getdata(chfile,1)
                ncat1 = len(cat1)
                print('  '+str(ncat1)+' sources')

                # Make sure it's in the right format
                if len(cat1.dtype.fields) != 32:
                    if verbose: print('  This catalog does not have the right format. Skipping')
                    del(cat1)
                    ncat1 = 0

                # Only include sources inside Boundary+Buffer zone
                #  -use ROI_CUT
                #  -reproject to tangent plane first so we don't have to deal
                #     with RA=0 wrapping or pol issues
                if buffdict is not None:
                    lon, lat = coords.rotsphcen(cat1['ra'],cat1['dec'],buffdict['cenra'],buffdict['cendec'],gnomic=True)
                    ind0, ind1 = utils.roi_cut(buffdict['lon'],buffdict['lat'],lon,lat)
                    nmatch = len(ind1)
                    # Only want source inside this pixel
                    if nmatch>0:
                        cat1 = cat1[ind1]
                    ncat1 = len(cat1)
                    if verbose: print('  '+str(nmatch)+' sources are inside this pixel')

                # Combine the catalogs
                if ncat1 > 0:
                    if cat is None:
                        dtype_cat = cat1.dtype
                        cat = np.zeros(np.sum(chmeta['nsources']),dtype=dtype_cat)
                        catcount = 0
                    cat[catcount:catcount+ncat1] = cat1
                    catcount += ncat1
                    expcatcount += ncat1

        # Add total number of measurements for this exposure
        newmeta['nmeas'] = expcatcount
        # Add metadata to ALLMETA, only if some measurements overlap
        if expcatcount>0:
            if allmeta is None:
                allmeta = newmeta
            else:
                allmeta = np.hstack((allmeta,newmeta))
            

    if cat is not None: cat=cat[0:catcount]  # trim excess
    if cat is None: cat=np.array([])         # empty cat
    if allmeta is None: allmeta=np.array([])

    return cat, allmeta
    

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix pixel number')
    parser.add_argument('--version', type=str, default='v3', help='Version number')
    parser.add_argument('--nside', type=int, default=128, help='HEALPix Nside')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    parser.add_argument('--outdir', type=str, default='', help='Output directory')
    #parser.add_argument('--filesexist', type=float, default=0.2, help='Time to wait between checking the status of running jobs')
    #parser.add_argument('--pixfiles', type=str, default=False, help='IDL program')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]
    radeg = np.float64(180.00) / np.pi

    # Inputs
    pix = int(args.pix[0])
    version = args.version
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
    if (os.path.exists(outfile) or os.path.exists(outfile+'.gz')) & ~redo:
        print(outfile+' EXISTS already and REDO not set')
        sys.exit()

    print("Combining InstCal SExtractor catalogs for Healpix pixel = "+str(pix))

    # Load the list
    listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
    if os.path.exists(listfile) is False:
        print(listfile+" NOT FOUND")
        sys.exit()
    healstr = Table(fits.getdata(listfile,1))
    index = Table(fits.getdata(listfile,2))
    # Find our pixel
    ind, = np.where(index['PIX'] == pix)
    nind = len(ind)
    if nind == 0:
        print("No entries for Healpix pixel '"+str(pix)+"' in the list")
        sys.exit()
    ind = ind[0]
    hlist = healstr[index[ind]['LO']:index[ind]['HI']+1]
    nlist = len(hlist)
    # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
    #  so we can deal with the edge cases
    neipix = hp.get_all_neighbours(nside,pix)
    for neip in neipix:
        ind1, = np.where(index['PIX'] == neip)
        nind1 = len(ind1)
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
        bd, = np.where(rabuff>180)
        if len(bd)>0:rabuff[bd] -=360.0
    buffdict = {'cenra':cenra,'cendec':cendec,'rar':utils.minmax(rabuff),'decr':utils.minmax(decbuff),'ra':rabuff,'dec':decbuff,\
                'lon':lonbuff,'lat':latbuff,'lr':utils.minmax(lonbuff),'br':utils.minmax(latbuff)}

    
    # Initialize the ID structure
    # this will contain the MeasID, Exposure name, ObjectID
    dtype_idstr = np.dtype([('measid',np.str,200),('exposure',np.str,200),('expnum',np.str,200),('objectid',np.str,200),('objectindex',int)])
    idstr = np.zeros(1000000,dtype=dtype_idstr)
    nidstr = len(idstr)
    idcnt = 0

    # Initialize the object structure
    dtype_obj = np.dtype([('objectid',np.str,100),('pix',int),('ra',np.float64),('dec',np.float64),('raerr',float),('decerr',float),
                          ('pmra',float),('pmdec',float),('pmraerr',float),('pmdecerr',float),('mjd',np.float64),
                          ('deltamjd',float),('ndet',int),('nphot',int),
                          ('ndetu',int),('nphotu',int),('umag',float),('urms',float),('uerr',float),('uasemi',float),('ubsemi',float),('utheta',float),
                          ('ndetg',int),('nphotg',int),('gmag',float),('grms',float),('gerr',float),('gasemi',float),('gbsemi',float),('gtheta',float),
                          ('ndetr',int),('nphotr',int),('rmag',float),('rrms',float),('rerr',float),('rasemi',float),('rbsemi',float),('rtheta',float),
                          ('ndeti',int),('nphoti',int),('imag',float),('irms',float),('ierr',float),('iasemi',float),('ibsemi',float),('itheta',float),
                          ('ndetz',int),('nphotz',int),('zmag',float),('zrms',float),('zerr',float),('zasemi',float),('zbsemi',float),('ztheta',float),
                          ('ndety',int),('nphoty',int),('ymag',float),('yrms',float),('yerr',float),('yasemi',float),('ybsemi',float),('ytheta',float),
                          ('ndetvr',int),('nphotvr',int),('vrmag',float),('vrrms',float),('vrerr',float),('vrasemi',float),('vrbsemi',float),('vrtheta',float),
                          ('asemi',float),('asemierr',float),('bsemi',float),('bsemierr',float),('theta',float),('thetaerr',float),
                          ('fwhm',float),('flags',int),('class_star',float),('ebv',float)])
    tags = dtype_obj.names
    obj = np.zeros(500000,dtype=dtype_obj)
    obj['pix'] = pix
    nobj = len(obj)
    dtype_totobj = np.dtype([('ra',np.float64),('dec',np.float64),('ramjd',np.float64),('decmjd',np.float64),('ramjd2',np.float64),
                             ('decmjd2',np.float64),('minmjd',np.float64),('maxmjd',np.float64),('umag2',np.float64),('gmag2',np.float64),
                             ('rmag2',np.float64),('imag2',np.float64),('zmag2',np.float64),('ymag2',np.float64),('vrmag2',np.float64),
                             ('utot',np.float64),('gtot',np.float64),('rtot',np.float64),('itot',np.float64),('ztot',np.float64),
                             ('ytot',np.float64),('vrtot',np.float64)])
    totags = dtype_totobj.names
    totobj = np.zeros(nobj,dtype=dtype_totobj)
    totobj['minmjd'] = 999999.0
    totobj['maxmjd'] = -999999.0    
    cnt = 0


    # Load the measurement catalog
    metafiles = [m.replace('_cat','_meta').strip() for m in hlist['FILE']]
    cat, allmeta = loadmeas(metafiles,buffdict)
    ncat = utils.size(cat)
    # currently ALLMETA includes exposures with no measurements in this healpix!!!

    # coordinates of measurement
    X = np.column_stack((np.array(cat['RA']),np.array(cat['DEC'])))
    # Compute DBSCAN on all measurements
    db = DBSCAN(eps=0.5/3600, min_samples=1).fit(X)
    labelindex = utils.create_index(db.labels_)
    nobj = len(labelindex['value'])
    print(str(nobj)+' unique objects clustered within 0.5 arcsec')

    # Initialize the OBJ structured arra
    obj = np.zeros(nobj,dtype=dtype_obj)
    obj['objectid'] = utils.strjoin( str(pix)+'.', ((np.arange(nobj)+1).astype(np.str)) )
    obj['pix'] = pix
    # all bad to start
    for f in ['pmra','pmraerr','pmdec','pmdecerr','asemi','bsemi','theta','asemierr','bsemierr','thetaerr','fwhm','class_star']: obj[f]=np.nan
    for f in ['u','g','r','i','z','y','vr']:
        obj[f+'mag'] = 99.99
        obj[f+'err'] = 9.99
        obj[f+'rms'] = np.nan
        obj[f+'asemi'] = np.nan
        obj[f+'bsemi'] = np.nan
        obj[f+'theta'] = np.nan

    # Loop over the objects
    for i,lab in enumerate(labelindex['value']):
        if (i % 100)==0: print(i)
        oindx = labelindex['index'][labelindex['lo'][i]:labelindex['hi'][i]+1]
        cat1 = cat[oindx]
        ncat1 = len(cat1)
        obj['ndet'][i] = ncat1

        # Computing quantities
        # Mean RA/DEC, RAERR/DECERR
        if ncat1>1:
            wt_ra = 1.0/cat1['RAERR']**2
            wt_dec = 1.0/cat1['DECERR']**2
            obj['ra'][i] = np.sum(cat1['RA']*wt_ra)/np.sum(wt_ra)
            obj['raerr'][i] = np.sum(1.0/np.sum(wt_ra))
            obj['dec'][i] = np.sum(cat1['DEC']*wt_dec)/np.sum(wt_dec)
            obj['raerr'][i] = np.sum(1.0/np.sum(wt_dec))
            obj['mjd'][i] = np.mean(cat1['MJD'])
            obj['deltamjd'][i] = np.max(cat1['MJD'])-np.min(cat1['MJD'])
        else:
            obj['ra'][i] = cat1['RA']
            obj['dec'][i] = cat1['DEC']
            obj['raerr'][i] = cat1['RAERR']
            obj['decerr'][i] = cat1['DECERR']
            obj['mjd'][i] = cat1['MJD']
            obj['deltamjd'][i] = 0

        #if ncat1>4:
        #    import pdb; pdb.set_trace()            

        # Mean proper motion and errors
        # fit robust line to RA values vs. time
        if ncat1>1:
            #raerr = cat1['RAERR'] / (3600*np.cos(obj['dec'][i]/radeg))   # convert to arcsec (true angle) to RA deg
            raerr = np.array(cat1['RAERR']*1e3,np.float64)    # milli arcsec
            ra = np.array(cat1['RA'],np.float64)
            ra -= np.mean(ra)
            ra *= 3600*1e3 * np.cos(obj['dec'][i]/radeg)     # convert to true angle, milli arcsec
            t = cat1['MJD']
            t -= np.mean(t)
            ra_coef, ra_coeferr = utils.poly_fit(t,ra,1,sigma=raerr,robust=True)
            pmra = ra_coef[0] * 365.2425           # mas/year
            pmraerr = ra_coeferr[0] * 365.2425     # mas/yr
            #ra_coef, ra_coeferr = utils.poly_fit(cat1['MJD'],cat1['RA'],1,sigma=raerr,robust=True)
            #pmra = ra_coef[0]                      # deg[ra]/day
            #pmra *= (3600*1e3)*365.2425            # mas/year
            #pmra *= np.cos(obj['dec'][i]/radeg)    # mas/year, true angle
            #pmraerr = ra_coeferr[0]                # deg[ra]/day
            #pmraerr *= (3600*1e3)*365.2425         # mas/year
            #pmraerr *= np.cos(obj['dec'][i]/radeg) # mas/year, true angle
            obj['pmra'][i] = pmra
            obj['pmraerr'][i] = pmraerr

            decerr = np.array(cat1['DECERR']*1e3,np.float64)   # milli arcsec
            dec = np.array(cat1['DEC'],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            dec_coef, dec_coeferr = utils.poly_fit(t,dec,1,sigma=decerr,robust=True)
            pmdec = dec_coef[0] * 365.2425          # mas/year
            pmdecerr = dec_coeferr[0] * 365.2425    # mas/year
            #decerr = cat1['DECERR'] / 3600.0       # convert from arcsec to deg
            #dec_coef, dec_coeferr = utils.poly_fit(cat1['MJD'],cat1['DEC'],1,sigma=decerr,robust=True)
            #pmdec = dec_coef[0]                    # deg[dec]/day
            #pmdec *= (3600*1e3)*365.2425           # mas/year
            #pmdecerr = dec_coeferr[0]             # deg[dec]/day
            #pmdecerr *= (3600*1e3)*365.2425        # mas/year
            obj['pmdec'][i] = pmdec
            obj['pmdecerr'][i] = pmdecerr

        # Mean magnitudes
        # Convert totalwt and totalfluxwt to MAG and ERR
        #  and average the morphology parameters PER FILTER
        filtindex = utils.create_index(cat1['FILTER'].astype(np.str))
        nfilters = len(filtindex['value'])
        for f in range(nfilters):
            filt = filtindex['value'][f].lower()
            findx = filtindex['index'][filtindex['lo'][f]:filtindex['hi'][f]+1]
            obj['ndet'+filt][i] = filtindex['num'][f]
            gph, = np.where(cat['MAG_AUTO'][findx]<50)
            ngph = len(gph)
            obj['nphot'+filt][i] = ngph
            if ngph>0:
                wt = 1.0/cat1['MAGERR_AUTO'][findx[gph]]**2
                newflux = np.sum( 2.5118864**cat1['MAG_AUTO'][findx[gph]] * wt) / np.sum(wt)
                newmag = 2.50*np.log10(newflux)
                newerr = np.sqrt(1.0/np.sum(wt))
                obj[filt+'mag'][i] = newmag
                obj[filt+'err'][i] = newerr
                # Calculate RMS
                if ngph>1: obj[filt+'rms'][i] = np.sqrt(np.mean((cat1['MAG_AUTO'][findx[gph]]-newmag)**2))

            #import pdb; pdb.set_trace()

            # Calculate mean morphology parameters
            obj[filt+'asemi'][i] = np.mean(cat1['ASEMI'][findx])
            obj[filt+'bsemi'][i] = np.mean(cat1['BSEMI'][findx])
            obj[filt+'theta'][i] = np.mean(cat1['THETA'][findx])
            
            #import pdb; pdb.set_trace()

        #import pdb; pdb.set_trace()

        # Make NPHOT from NPHOTX
        obj['nphot'][i] = obj['nphotu'][i]+obj['nphotg'][i]+obj['nphotr'][i]+obj['nphoti'][i]+obj['nphotz'][i]+obj['nphoty'][i]+obj['nphotvr'][i]

        # Mean morphology parameters
        obj['asemi'][i] = np.mean(cat1['ASEMI'])
        obj['bsemi'][i] = np.mean(cat1['BSEMI'])
        obj['theta'][i] = np.mean(cat1['THETA'])
        obj['asemierr'][i] = np.sqrt(np.sum(cat1['ASEMIERR']**2)) / ncat1
        obj['bsemierr'][i] = np.sqrt(np.sum(cat1['BSEMIERR']**2)) / ncat1
        obj['thetaerr'][i] = np.sqrt(np.sum(cat1['THETAERR']**2)) / ncat1
        obj['fwhm'][i] = np.mean(cat1['FWHM'])
        obj['class_star'][i] = np.mean(cat['CLASS_STAR'])

        # Stuff information into IDSTR


    # Add E(B-V)
    print('Getting E(B-V)')
    sfd = SFDQuery()
    c = SkyCoord(obj['ra'],obj['dec'],frame='icrs',unit='deg')
    #c = SkyCoord('05h00m00.00000s','+30d00m00.0000s', frame='icrs') 
    ebv = sfd(c)
    obj['ebv'] = ebv

    import pdb; pdb.set_trace()

    # ONLY INCLUDE OBJECTS WITH AVERAGE RA/DEC
    # WITHIN THE BOUNDARY OF THE HEALPIX PIXEL!!!
    ipring = hp.pixelfunc.ang2pix(nside,obj['ra'],obj['dec'],lonlat=True)
    ind1, = np.where(ipring == pix)
    nmatch = len(ind1)
    if nmatch==0:
        print('None of the final objects fall inside the pixel')
        sys.exit()
    # Get trimmed objects and indices
    objtokeep = np.zeros(nobj,bool)         # boolean to keep or trim objects
    objtokeep[ind1] = True
    if nmatch<nobj:
        trimind = np.arange(nobj)
        trimind = np.delete(trimind,ind1)
        #trimind = utils.remove_indices(trimind,ind1)
        trimobj = obj[trimind]          # trimmed objects
    newobjindex = np.zeros(nobj,int)-1    # new indices
    newobjindex[ind1] = np.arange(nmatch)
    # Keep the objects inside the Healpix
    obj = obj[ind1]
    print(str(nmatch)+' final objects fall inside the pixel')

    # Remove trimmed objects from IDSTR
    totrim, = np.where(objtokeep[idstr['objectindex']]==0)  #using old index
    if len(totrim)>0:
        # Trim objects
        idstr = np.delete(idstr,totrim)
        #idstr = utils.remove_indices(idstr,totrim)
        # Update IDSTR.objectindex
        old_idstr_objectindex = idstr['objectindex']
        idstr['objectindex'] = newobjindex[old_idstr_objectindex]

    # Create final summary structure from ALLMETA
    #  get exposures that are in IDSTR
    #  sometimes EXPNUM numbers have the leading 0s removed
    #  and sometimes not, so turn to LONG to match
    dum, uiexpnum = np.unique(idstr['expnum'].astype(int),return_index=True)
    uexpnum = idstr[uiexpnum]['expnum'].astype(int)
    nuexpnum = len(uexpnum)
    ind1,ind2 = utils.match(allmeta['expnum'].astype(int),uexpnum)
    nmatch = len(ind1)
    sumstr = Table(allmeta[ind1])
    col_nobj = Column(name='nobjects', dtype=np.int, length=len(sumstr))
    col_healpix = Column(name='healpix', dtype=np.int, length=len(sumstr))
    sumstr.add_columns([col_nobj, col_healpix])
    sumstr['nobjects'] = 0
    sumstr['healpix'] = pix
    # get number of objects per exposure
    expnum = idstr['expnum'].astype(int)
    siexp = np.argsort(expnum)
    expnum = expnum[siexp]
    if nuexpnum>1:
        brklo, = np.where(expnum != np.roll(expnum,1))
        nbrk = len(brklo)
        brkhi = np.hstack((brklo[1:nbrk],len(expnum)))
        numobjexp = brkhi-brklo+1
    else:
        numobjexp=len(expnum)
    ind1,ind2 = utils.match(sumstr['expnum'].astype(int),uexpnum)
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
