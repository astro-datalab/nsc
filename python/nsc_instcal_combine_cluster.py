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
from scipy.interpolate import interp1d

def add_elements(cat,nnew=300000):
    """ Add more elements to a catalog"""
    ncat = len(cat)
    old = cat.copy()
    nnew = utils.gt(nnew,ncat)
    cat = np.zeros(ncat+nnew,dtype=old.dtype)
    cat[0:ncat] = old
    del(old)
    return cat    

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
                           ('zptermsig',float),('refmatch',int)])

    dtype_cat = np.dtype([('MEASID',np.str,200),('OBJECTID',np.str,200),('EXPOSURE',np.str,200),('CCDNUM',int),('FILTER',np.str,10),
                          ('MJD',float),('X',float),('Y',float),('RA',float),('RAERR',float),('DEC',float),('DECERR',float),
                          ('MAG_AUTO',float),('MAGERR_AUTO',float),('MAG_APER1',float),('MAGERR_APER1',float),('MAG_APER2',float),
                          ('MAGERR_APER2',float),('MAG_APER4',float),('MAGERR_APER4',float),('MAG_APER8',float),('MAGERR_APER8',float),
                          ('KRON_RADIUS',float),('ASEMI',float),('ASEMIERR',float),('BSEMI',float),('BSEMIERR',float),('THETA',float),
                          ('THETAERR',float),('FWHM',float),('FLAGS',int),('CLASS_STAR',float)])

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
                print('  chip '+str(chmeta[j]['ccdnum'])+'  '+str(ncat1)+' sources')

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
                        #dtype_cat = cat1.dtype
                        #ncat_init = np.sum(chmeta['nsources'])*utils.size(metafile)
                        ncat_init = np.maximum(100000,ncat1)
                        cat = np.zeros(ncat_init,dtype=dtype_cat)
                        catcount = 0
                    # Add more elements
                    if (catcount+ncat1)>ncat:
                        cat = add_elements(cat,np.maximum(100000,ncat1))
                        ncat = len(cat)
                    cat[catcount:catcount+ncat1] = cat1
                    catcount += ncat1
                    expcatcount += ncat1

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

    
    # IDSTR schema
    dtype_idstr = np.dtype([('measid',np.str,200),('exposure',np.str,200),('objectid',np.str,200),('objectindex',int)])

    # OBJ schema
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
                          ('fwhm',float),('flags',int),('class_star',float),('ebv',float),('jvar',float),('kvar',float),('avgvar',float),
                          ('chivar',float),('phrms',float)])

    # Load the measurement catalog
    metafiles = [m.replace('_cat','_meta').strip() for m in hlist['FILE']]
    cat, allmeta = loadmeas(metafiles,buffdict)
    #import pdb; pdb.set_trace()
    # KLUDGE
    #cat = np.array(fits.getdata('60025_cat.fits',1))
    #allmeta = np.array(fits.getdata('60025_allmeta.fits',1))
    #cat = np.array(fits.getdata('148487_cat.fits',1))
    #allmeta = np.array(fits.getdata('148487_allmeta.fits',1))
    ncat = utils.size(cat)

    t1 = time.time()



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
    for f in ['pmra','pmraerr','pmdec','pmdecerr','asemi','bsemi','theta','asemierr',
              'bsemierr','thetaerr','fwhm','class_star','jvar','kvar','avgvar','chivar']: obj[f]=np.nan
    for f in ['u','g','r','i','z','y','vr']:
        obj[f+'mag'] = 99.99
        obj[f+'err'] = 9.99
        obj[f+'rms'] = np.nan
        obj[f+'asemi'] = np.nan
        obj[f+'bsemi'] = np.nan
        obj[f+'theta'] = np.nan
    idstr = np.zeros(ncat,dtype=dtype_idstr)

    # Loop over the objects
    for i,lab in enumerate(labelindex['value']):
        if (i % 1000)==0: print(i)
        oindx = labelindex['index'][labelindex['lo'][i]:labelindex['hi'][i]+1]
        cat1 = cat[oindx]
        ncat1 = len(cat1)
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
            t = cat1['MJD']
            t -= np.mean(t)
            t /= 365.2425                          # convert to year
            # Calculate robust slope
            pmra, pmraerr = utils.robust_slope(t,ra,raerr,reweight=True)
            obj['pmra'][i] = pmra                 # mas/yr
            obj['pmraerr'][i] = pmraerr           # mas/yr

            decerr = np.array(cat1['DECERR']*1e3,np.float64)   # milli arcsec
            dec = np.array(cat1['DEC'],np.float64)
            dec -= np.mean(dec)
            dec *= 3600*1e3                         # convert to milli arcsec
            # Calculate robust slope
            pmdec, pmdecerr = utils.robust_slope(t,dec,decerr,reweight=True)
            obj['pmdec'][i] = pmdec               # mas/yr
            obj['pmdecerr'][i] = pmdecerr         # mas/yr

        # Mean magnitudes
        # Convert totalwt and totalfluxwt to MAG and ERR
        #  and average the morphology parameters PER FILTER
        filtindex = utils.create_index(cat1['FILTER'].astype(np.str))
        nfilters = len(filtindex['value'])
        resid = np.zeros(ncat1)+np.nan     # residual mag
        relresid = np.zeros(ncat1)+np.nan  # residual mag relative to the uncertainty
        for f in range(nfilters):
            filt = filtindex['value'][f].lower()
            findx = filtindex['index'][filtindex['lo'][f]:filtindex['hi'][f]+1]
            obj['ndet'+filt][i] = filtindex['num'][f]
            gph, = np.where(cat['MAG_AUTO'][findx]<50)
            ngph = len(gph)
            obj['nphot'+filt][i] = ngph
            if ngph==1:
                obj[filt+'mag'][i] = cat1['MAG_AUTO'][findx[gph]]
                obj[filt+'err'][i] = cat1['MAGERR_AUTO'][findx[gph]]
            if ngph>1:
                newmag, newerr = utils.wtmean(cat1['MAG_AUTO'][findx[gph]], cat1['MAGERR_AUTO'][findx[gph]],magnitude=True,reweight=True,error=True)
                obj[filt+'mag'][i] = newmag
                obj[filt+'err'][i] = newerr
                # Calculate RMS
                obj[filt+'rms'][i] = np.sqrt(np.mean((cat1['MAG_AUTO'][findx[gph]]-newmag)**2))
                # Residual mag
                resid[findx[gph]] = cat1['MAG_AUTO'][findx[gph]]-newmag
                # Residual mag relative to the uncertainty
                #  set a lower threshold of 0.02 in the uncertainty
                relresid[findx[gph]] = np.sqrt(ngph/(ngph-1)) * (cat1['MAG_AUTO'][findx[gph]]-newmag)/np.sqrt(cat1['MAGERR_AUTO'][findx[gph]]**2+0.02**2)


            # Calculate mean morphology parameters
            obj[filt+'asemi'][i] = np.mean(cat1['ASEMI'][findx])
            obj[filt+'bsemi'][i] = np.mean(cat1['BSEMI'][findx])
            obj[filt+'theta'][i] = np.mean(cat1['THETA'][findx])

        # Calculate photometric RMS across all bands
        gdresid = np.isfinite(resid)
        ngdresid = np.sum(gdresid)
        if ngdresid>0:
            resid2 = resid[gdresid]
            rms = np.sqrt(np.mean(resid2**2))
            obj['phrms'][i] = rms
                       

        # Calculate variability indices
        gdrelresid = np.isfinite(relresid)
        ngdrelresid = np.sum(gdrelresid)
        if ngdrelresid>0:
            relresid2 = relresid[gdrelresid]
            pk = relresid2**2-1
            jvar = np.sum( np.sign(pk)*np.sqrt(np.abs(pk)) )/ngdrelresid
            avgvar = np.mean(relresid2)    # average of relative residuals
            chivar = np.sqrt(np.sum(relresid2**2))/ngdrelresid
            kdenom = np.sqrt(np.sum(relresid2**2)/ngdrelresid)
            if kdenom!=0:
                kvar = (np.sum(np.abs(relresid2))/ngdrelresid) / kdenom
            else:
                kvar = 0.0
            obj['jvar'][i] = jvar
            obj['kvar'][i] = kvar
            obj['avgvar'][i] = avgvar
            obj['chivar'][i] = chivar
            #if chivar>50: import pdb; pdb.set_trace()

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
        obj['class_star'][i] = np.mean(cat1['CLASS_STAR'])
        obj['flags'][i] = np.bitwise_or.reduce(cat1['FLAGS'])  # OR combine


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
    totrim, = np.where(~objtokeep[idstr['objectindex']])  #using old index
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
    dum, uiexposure = np.unique(idstr['exposure'],return_index=True)
    uexposure = idstr['exposure'][uiexposure]
    nuexposure = len(uexposure)
    ind1,ind2 = utils.match(allmeta['base'],uexposure)
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
        brklo, = np.where(exposure != np.roll(exposure,1))
        nbrk = len(brklo)
        brkhi = np.hstack((brklo[1:nbrk],len(exposure)))
        numobjexp = brkhi-brklo+1
    else:
        numobjexp=len(exposure)
    ind1,ind2 = utils.match(sumstr['base'],uexposure)
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
