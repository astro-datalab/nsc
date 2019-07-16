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
from dlnpyutils import utils
import subprocess
time

def add_elements(cat,nnew=300000L):
    """ Add more elements to a catalog"""
    ncat = len(cat)
    old = cat.copy()
    nnew = utils.gt(nnew,ncat)
    cat = np.zeros(ncat+nnew,dtype=old.dtype)
    cat[0:ncat] = old
    del(old)
    return cat
    
def add_cat(obj,totobj,idstr,cnt,idcnt,ind1,cat,meta):
    """ Add object information from a new meas catalog of matched objects"""

    f = meta['filter'].lower()
    # Copy to final structure
    obj[ind1]['ra'] += cat['ra']
    obj[ind1]['dec'] += cat['dec']
    obj[ind1]['raerr'] += 1.0/cat['raerr']**2                           # total(ra_wt)
    obj[ind1]['decerr'] += 1.0/cat['decerr']**2                         # total(dec_wt)
    obj[ind1]['pmra'] += (1.0/cat['raerr']**2) * meta.mjd*cat['ra']     # total(wt*mjd*ra)
    obj[ind1]['pmdec'] += (1.0/cat['decerr']**2) * meta.mjd*cat['dec']  # total(wt*mjd*dec)
    obj[ind1]['mjd'] += meta['mjd']                                 # total(mjd)
    obj[ind1]['ndet'] += 1
    # Detection and morphology parameters for this FILTER
    obj[obj]['ndet'+f] += 1
    obj[obj][f+'asemi'] += cat['asemi']
    obj[obj][f+'bsemi'] += cat['bsemi']
    obj[obj][f+'theta'] += cat['theta']
    # Good photometry for this FILTER
    gdmag, = np.where(cat['mag_auto']<50)
    if len(gdmag)>0:
      obj[ind1[gdmag]][f+'mag'] += 2.5118864d**cat[gdmag]['mag_auto'] * (1.0/cat[gdmag]['magerr_auto']**2)
      obj[ind1[gdmag]][f+'err'] += 1.0/cat[gdmag]['magerr_auto']**2
      obj[ind1[gdmag]]['nphot'+f] += 1
    obj[ind1]['asemi'] += cat['asemi']
    obj[ind1]['asemierr'] += cat['asemierr']**2
    obj[ind1]['bsemi'] += cat['bsemi']
    obj[ind1]['bsemierr'] += cat['bsemierr']**2
    obj[ind1]['theta'] += cat['theta']
    obj[ind1]['thetaerr'] += cat['thetaerr']**2
    obj[ind1]['fwhm'] += cat['fwhm']  # in arcsec
    obj[ind1]['flags'] += cat['flags']
    obj[ind1]['class_star'] += cat['class_star']
    totobj[ind1]['ra'] += cat['ra'] * (1.0/cat['raerr']**2)             # total(ra*wt)
    totobj[ind1]['dec'] += cat['dec'] * (1.0/cat['decerr']**2)          # total(dec*wt)
    totobj[ind1]['ramjd'] += (1.0/cat['raerr']**2) * meta['mjd']        # total(wt_ra*mjd)
    totobj[ind1]['decmjd'] += (1.0/cat['decerr']**2) * meta['mjd']      # total(wt_dec*mjd)
    totobj[ind1]['ramjd2'] += (1.0/cat['raerr']**2) * meta['mjd']**2    # total(wt_ra*mjd**2)
    totobj[ind1]['decmjd2'] += (1.0/cat['decerr']**2) * meta['mjd*']**2  # total(wt_dec*mjd**2)
    totobj[ind1]['minmjd'] = utils.lt( meta['mjd'], totobj[ind1]['minmjd'] )
    totobj[ind1]['maxmjd'] = utils.gt( meta['mjd'], totobj[ind1]['maxmjd'] )
    if len(gdmag)>0:
        totobj[ind1[gdmag]][f+'tot'] += cat[gdmag]['mag_auto']       # sum(mag)
        totobj[ind1[gdmag]][f+'mag2'] += np.float64(cat[gdmag]['mag_auto'])**2   # sum(mag**2), need dbl to precent underflow
    cnt += ncat

    # Add new elements to IDSTR
    if idcnt+ncat gt nidstr:
        idstr = add_elements(idstr)
        nidstr = len(idstr)

    # Add to IDSTR
    idstr[idcnt:idcnt+ncat]['measid'] = cat['measid']
    idstr[idcnt:idcnt+ncat]['exposure'] = meta['base']
    idstr[idcnt:idcnt+ncat]['expnum'] = meta['expnum']
    idstr[idcnt:idcnt+ncat]['objectid'] = newexp['objectid']
    idstr[idcnt:idcnt+ncat]['objectindex'] = ind1
    idcnt += ncat

    return obj,totobj,idstr,cnt,idcnt
    

def loadmeas(metafile,buffdict=None):

    if os.path.exists(metafile) is False:
        print(metafile+' NOT FOUND')
        return np.array([])
    meta = fits.getdata(metafile,1)
    chmeta = fits.getdata(metafile,2)    

    fdir = os.path.dirname(metafile)
    fbase, ext = os.path.splitext(os.path.basename(metafile))
    fbase = fbase[:-5]   # remove _meta at end

    # Loop over the chip files
    cat = None
    for j in range(len(chmeta)):
        # Check that this chip was astrometrically calibrated
        #   and falls in to HEALPix region
        if chmeta[j]['ngaiamatch'] == 0:
            print('This chip was not astrometrically calibrate')

        # Check that this overlaps the healpix region
        inside = True
        if buffdict is not None:
            vra = chmeta[j]['vra']
            vdec = chmeta[j]['vdec']
            if (np.max(vra)-np.min(vra)) > 100:    # deal with RA=0 wrapround
                bd, = np.where(vra>180)
                if len(bd)>0: vra[bd] -= 360
            if dopolygonsoverlap(buffdict['ra'],buffdict['dec'],vra,vdec):
                print('This chip does NOT overlap the HEALPix region+buffer')
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
                print('  This catalog does not have the right format. Skipping')
                del(cat1)
                ncat1 = 0

            # Only include sources inside Boundary+Buffer zone
            #  -use ROI_CUT
            #  -reproject to tangent plane first so we don't have to deal
            #     with RA=0 wrapping or pol issues
            if buffdict is not None:
                lon, lat = rotsphcen(cat1['ra'],cat1['dec'],buffdict['cenra'],buffdict['cendec'],gnomic=True)
                ind0, ind1 = utils.roi_cut(buffdict['lon'],buffdict['lat'],lon,lat)
                nmatch = len(ind1)
                # Only want source inside this pixel
                if nmatch>0:
                    cat1 = cat1[ind1]
                ncat1 = len(cat1)
                print('  '+str(nmatch)+' sources are inside this pixel')

            # Combine the catalogs
            if ncat1 > 0:
                if cat is None:
                    dtype_cat = cat1.dtype
                    cat = np.zeros(np.sum(chmeta['nsources']),dtype=dtype_cat)
                    catcount = 0L
                cat[catcount:catcount+ncat1] = cat1
                catcount += ncat1

              #BOMB1:
    if cat is not None: cat=cat[0:catcount]  # trim excess
    if cat is None: cat=np.array([])         # empty cat

    return cat
    

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC data for one healpix region.')
    parser.add_argument('pix', type=str, nargs=1, help='HEALPix pixel number')
    parser.add_argument('--version', type=str, default='v3', help='Version number')
    parser.add_argument('--nside', type=int, default=128, help='HEALPix Nside')
    parser.add_argument('--redo', type=str, default='No', help='Redo this HEALPIX')
    parser.add_argument('--outdir', type=str, default='', help='Output directory')
    #parser.add_argument('--filesexist', type=float, default=0.2, help='Time to wait between checking the status of running jobs')
    #parser.add_argument('--pixfiles', type=str, default=False, help='IDL program')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    pix = args.pix
    version = args.version
    nside = args.nside
    redo = args.redo
    if ((redo=='True') or (redo=='TRUE') or (redo=='Yes') or (redo=='YES') or (redo=='Y') or (redo=='y')):
        redo = True
    else:
        redo = False
    outdir = args.outdir

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        dir = "/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        dir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"?"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    # Not enough inputs
    n = len(sys.argv)
    if n < 1:
        print "Syntax - nsc_instcal_combine pix --version v# --nside ### --redo YES/NO --outdir OUTDIR"
        sys.exit()

    # Check if output file already exists
    if outdir == '': outdir=dir+'combine/'
    subdir = str(int(pix)/1000)    # use the thousands to create subdirectory grouping
    outfile = outdir+'/'+subdir+'/'+str(pix)+'.fits'
    if (os.path.exists(outfile) or os.path.exists(outfile+'.gz')) and ~redo:
        print(outfile+' EXISTS already and REDO not set')
        sys.exit()

    print("Combining InstCal SExtractor catalogs for Healpix pixel = "+str(pix,2))

    # Load the list
    listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
    if os.path.exists(list) is False:
        print(listfile+" NOT FOUND")
        sys.exist()
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
    lonbound, latbound = utils.rotsphcen(rabound,decbound,cenra,cendec,gnomic=True)
    # expand by a fraction, it's not an extact boundary but good enough
    buffsize = 10.0/3600. # in deg
    radbound = sqrt(lonbound**2+latbound**2)
    frac = 1.0 + 1.5*np.max(buffsize/radbound)
    lonbuff = lonbound*frac
    latbuff = latbound*frac
    rabuff, decbuff = utils.rotsphcen(lonbuff,latbuff,cenra,cendec,gnomic=True,reverse=True)
    if (np.max(rabuff)-np.min(rabuff))>100:  # deal with RA=0 wraparound
        bd, = np.where(rabuff>180)
        if len(bd)>0:rabuff[bd] -=360.0
    buffdict = {'cenra':cenra,'cendec':cendec,'rar':utils.minmax(rabuff),'decr':utils.minmax(decbuff),'ra':rabuff,'dec':decbuff,\
                'lon':lonbuff,'lat':latbuff,'lr':utils.minmax(lonbuff),'br':utils.minmax(latbuff)}

    
    # Initialize the ID structure
    # this will contain the MeasID, Exposure name, ObjectID
    dtype_idstr = np.dtype([('measid',np.str),('exposure',np.str),('expnum',np.str),('objectid',np.str),('objectindex',long)])
    idstr = np.zeros(1000000,dtype=dtype_idstr)
    nidstr = len(idstr)
    idcnt = 0L

    # Initialize the object structure
    dtype_obj = np.dtype([('objectid',np.str),('pix',long),('ra',np.float64),('dec',np.float64),('raerr',float),('decerr',float),
                          ('pmra',float),('pmdec',float),('pmraerr',float),('pmdecerr',float),('mjd',np.float64),
                          ('deltamjd',float),('ndet',long),('nphot',long),
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
    cnt = 0L

    # Loop over the exposures
    allmeta = []
    for i in range(nlist):
        print(str(i+1)+' Loading '+list[i]['FILE'])

        # Load meta data file first
        metafile = repstr(list[i]['FILE'],'_cat','_meta')
        if os.path.exists(metafile) is False:
            print(metafile+' NOT FOUND')
            #goto,BOMB
        meta = fits.getdata(metafile,1)
        t = Time(times, format='isot', scale='utc')
        meta['mjd'] = t.mjd                    # recompute because some MJD are bad
        chmeta = fits.getdata(metafile,2)      # chip-level meta-data structure
        print('  FILTER='+meta['filter']+'  EXPTIME='+str(meta['exptime'])+' sec')

        # Loop over the chip files
        cat = None
        catcount = 0L
        for j in range(len(chmeta)):
              # Check that this chip was astrometrically calibrated
              if chmeta[j]['ngaiamatch'] == 0:
                  print('This chip was not astrometrically calibrate')
                  #goto,BOMB1

              # Check that this overlaps the healpix region
              vra = chmeta[j]['vra']
              vdec = chmeta[j]['vdec']
              if (np.max(vra)-np.min(vra)) > 100:    # deal with RA=0 wrapround
                  bd, = np.where(vra>180)
                  if len(bd)>0: vra[bd] -= 360
              if dopolygonsoverlap(buffdict['ra'],buffdict['dec'],vra,vdec):
                  print('This chip does NOT overlap the HEALPix region+buffer')
                  #goto,BOMB1

              # Load the chip-level catalog
              chfile = os.path.dirname(list[i]['file'])+'/'+list[i]['base']+'_'+str(chmeta[j]['ccdnum'])+'_meas.fits'
              if os.path.exists(chfile) is False:
                  print(chfile+' NOT FOUND')
                  #goto,BOMB1
              cat1 = fits.getdata(chfile,1)
              ncat1 = len(cat1)
              print('  '+str(ncat1)+' sources')

              # Make sure it's in the right format
              if len(cat1.dtype.fields) != 32:
                  print('  This catalog does not have the right format. Skipping')
                  #goto,BOMB1

              # Only include sources inside Boundary+Buffer zone
              #  -use ROI_CUT
              #  -reproject to tangent plane first so we don't have to deal
              #     with RA=0 wrapping or pol issues
              lon, lat = rotsphcen(cat1['ra'],cat1['dec'],buffdict['cenra'],buffdict['cendec'],gnomic=True)
              ind0, ind1 = utils.roi_cut(buffdict['lon'],buffdict['lat'],lon,lat)
              nmatch = len(ind1)
              # Only want source inside this pixel
              if nmatch==0:
                  print('  No sources inside this pixel')
                  #goto,BOMB
              print('  '+str(nmatch)+' sources are inside this pixel')
              cat1 = cat1[ind1]
              ncat1 = nmatch

              # Combine the catalogs
              if cat is None:
                  dtype_cat = cat1.dtype
                  cat = np.zeros(np.sum(chmeta['nsources']),dtype=dtype_cat)
              cat[catcount:catcount+ncat1] = cat1
              catcount += ncat1

              #BOMB1:
        if cat is not None: cat=cat[0:catcount]  # trim excess
        ncat = utils.size(cat)
        if ncat==0:
              print('This exposure does NOT cover the HEALPix')
              #goto,BOMB

        # Add metadata to ALLMETA
        #  Make sure it's in the right format
        dtype_meta = np.dtype([('file',np.str),('base',np.str),('expnum',long),('ra',np.float64),
                               ('dec',np.float64),('dateobs',np.str),('mjd',np.float64),('filter',np.str),
                               ('exptime',float),('airmass',float),('nsources',long),('fwhm',float),
                               ('nchips',long),('badchip31',bool),('rarms',float),('decrms',float),
                               ('ebv',float),('gaianmatch',long),('zpterm',float),('zptermerr',float),
                               ('refmatch',long)])
        newmeta = np.zero(1,dtype=dtype_meta)
        # Copy over the meta information
        for n in newmeta.names:
            if n in meta.names: newmeta[n]=meta[n]
        allmeta.append(newmeta)

        # Combine the data
        #-----------------
        # First catalog
        if cnt==0:
            ind1 = np.arange(len(cat))
            obj,totobj,idstr,cnt,idcnt = add_cat(obj,totobj,idstr,cnt,idcnt,ind1,cat,meta)
            obj[ind1]['objectid'] = str(pix)+'.'+str(np.arange(ncat)+1)


        # Second and up
        else:
            #  Match new sources to the objects
            ind1,ind2,dist = utils.xmatch(obj[0:cnt]['ra'],obj[0:cnt]['dec'],cat['ra'],cat['dec'],0.5)
            nmatch = len(ind1)
            #  Some matches, add data to existing record for these sources
            if nmatch>0:
                obj,totobj,idstr,cnt,idcnt = add_cat(obj,totobj,idstr,cnt,idcnt,ind1,cat[ind2],meta)
                if nmatch<ncat:
                    cat = np.delete(cat,ind2)
                    ncat = len(cat)
                else:
                    cat = np.array([])
                    ncat = 0

            # Some left, add records for these sources
            if ncat>0:
                print('  '+str(ncat)+' sources left to add')
                # Add new elements
                if (cnt+ncat)>nobj:
                    obj = add_elements(obj)
                    nobj = len(obj)
                ind1 = np.arange(ncat)
                obj,totobj,idstr,cnt,idcnt = add_cat(obj,totobj,idstr,cnt,idcnt,ind1,cat,meta)
                obj[ind1]['objectid'] = str(pix)+'.'+str(cnt+np.arange(ncat)+1)

    # No sources
    if cnt==0:
        print('No sources in this pixel')
        sys.exit()
    # Trim off the excess elements
    obj = obj[0:cnt]
    totobj = totobj[0:cnt]
    nobj = len(obj)
    print(str(nobj)+' final objects')
    idstr = idstr[0:idcnt]

    # Make NPHOT from NPHOTX
    obj['nphot'] = obj['nphotu']+obj['nphotg']+obj['nphotr']+obj['nphoti']+obj['nphotz']+obj['nphoty']+obj['nphotvr']

    # Convert total(mjd*ra) to true proper motion values
    #  the slope of RA vs. MJD is
    #  pmra=(total(wt*mjd*ra)/total(wt)-<mjd>*<ra>)/(total(wt*mjd^2)/total(wt)-<mjd>^2)
    #  we are creating the totals cumulatively as we go
    totobj['ra'] /= obj['raerr']        # wt mean RA (totalrawt/totalwt)
    totobj['dec'] /= obj['decerr']      # wt mean DEC (totaldecwt/totalwt)
    obj['mjd'] /= obj['ndet']           # mean MJD
    totobj['ramjd'] /= obj['raerr']     # wt_ra mean MJD
    totobj['decmjd'] /= obj['decerr']   # wt_dec mean MJD
    pmra = (obj['pmra']/obj['raerr']-totobj['ramjd']*totobj['ra'])/(totobj['ramjd2']/obj['raerr']-totobj['ramjd']**2)   # deg[ra]/day
    pmra *= (3600*1e3)*365.2425     # mas/year
    pmra *= np.cos(obj['dec']/radeg)      # mas/year, true angle
    pmdec = (obj['pmdec']/obj['decerr']-totobj['decmjd']*totobj['dec'])/(totobj['decmjd2']/obj['decerr']-totobj['decmjd']**2)  # deg/day
    pmdec *= (3600*1e3)*365.2425    # mas/year
    # Proper motion errors
    # pmerr = 1/sqrt( sum(wt*mjd^2) - <mjd>^2 * sum(wt) )
    #   if wt=1/err^2 with err in degrees, but we are using arcsec
    #   Need to divide by 3600 for PMDECERR and 3600*cos(dec) for PMRAERR
    pmraerr = 1.0/np.sqrt( totobj['ramjd2'] - totobj['ramjd']**2 * obj['raerr'] )
    pmraerr /= (3600*np.cos(totobj['dec']/radeg))    # correction for raerr in arcsec
    pmraerr *= (3600*1e3)*365.2425     # mas/year
    pmraerr *= np.cos(obj['dec']/radeg)      # mas/year, true angle
    pmdecerr = 1.0/np.sqrt( totobj['decmjd2'] - totobj['decmjd']**2 * obj['decerr'] )
    pmdecerr /= 3600                   # correction for decerr in arcsec
    pmdecerr *= (3600*1e3)*365.2425    # mas/year
    gdet, = np.where(obj['ndet']>1)
    ngdet = len(gdet)
    if ngdet>0:
        obj[gdet]['pmra'] = pmra[gdet]
        obj[gdet]['pmdec'] = pmdec[gdet]
        obj[gdet]['pmraerr'] = pmraerr[gdet]
        obj[gdet]['pmdecerr'] = pmdecerr[gdet]
    bdet, = np.where((obj['ndet']<2) | ~np.isfinite(pmra))
    nbdet = len(bdet)
    # sometimes it happens that the denominator is 0.0 
    #  when there are few closely spaced points
    #  nothing we can do, just mark as bad
    if nbdet>0:
        obj[bdet]['pmra'] = 999999.0
        obj[bdet]['pmdec'] = 999999.0
        obj[bdet]['pmraerr'] = 999999.0
        obj[bdet]['pmdecerr'] = 999999.0
    obj['deltamjd'] = totobj['maxmjd']-totobj['minmjd']
    # Average coordinates
    obj['ra'] = totobj['ra']   # now stuff in the average coordinates
    obj['dec'] = totobj['dec']
    obj['raerr'] = np.sqrt(1.0/obj['raerr'])    # err in wt mean RA, arcsec
    obj['decerr'] = np.sqrt(1.0/obj['decerr'])  # err in wt mean DEC, arcsec


    # Convert totalwt and totalfluxwt to MAG and ERR
    #  and average the morphology parameters PER FILTER
    filters = ['u','g','r','i','z','y','vr']
    for f in filters:
        newflux = obj[f+'mag'] / obj[f+'err']
        newmag = 2.50*np.log10(newflux)
        newerr = np.sqrt(1.0/obj[f+'err'])
        obj[f+'mag'] = newmag
        obj[f+'err'] = newerr
        bdmag, = np.where(~np.isfinite(newmag))
        nbdmag = len(bdmag)
        if nbdmag>0:
            obj[bdmag][f+'mag'] = 99.99
            obj[bdmag][f+'err'] = 9.99

        # Calculate RMS scatter
        #  RMS^2 * N = sum(mag^2) - 2*<mag>*sum(mag) + N*<mag>^2
        #   where <mag> is a weighted average
        #  RMS = sqrt( sum(mag^2)/N - 2*<mag>*sum(mag)/N + <mag>^2 )
        #  sum(mag^2) is in the MAG2 column and sum(mag) is in TOT
        rms = np.zeros(nobj,float)
        gdrms, = np.where(obj['nphot'+f]>1)
        ngdrms = len(gdrms)
        bdrms, = np.where(obj['nphot'+f]<=1)
        nbdrms = len(bdrms)
        if ngdrms>0:
           rms[gdrms] = np.sqrt( totobj[gdrms][f+'mag2']/obj[gdrms]['nphot'+f] - $
                                 2*obj[gdrms][f+'mag']*totobj[gdrms][f+'tot']/obj[gdrms]['nphot'+f] + obj[gdrms][f+'mag'])**2 )
        if nbdrms>0: rms[bdrms] = 999999.
        obj[f+'rms'] = rms

        # Average the morphology parameters PER FILTER
        gdet, = np.where(obj['ndet'+f]>0)
        ngdet = len(gdet)
        bdet, = np.where(obj['ndet'+f]==0)
        nbdet = len(bdet)        
        if ngdet>0:
            obj[gdet][f+'asemi'] /= obj[gdet]['ndet'+f]
            obj[gdet][f+'bsemi'] /= obj[gdet]['ndet'+f]
            obj[gdet][f+'theta'] /= obj[gdet]['ndet'+f]
        if nbdet>0:
            obj[bdet][f+'asemi'] = 99.99
            obj[bdet][f+'bsemi'] = 99.99
            obj[bdet][f+'theta'] = 99.99

    # Average the morphology parameters, Need a separate counter for that maybe?
    mtags = ['asemi','bsemi','theta','fwhm','class_star']
    gdet, = np.where(obj['ndet']>0)
    ngdet = len(gdet)
    bdet, = np.where(obj['ndet']==0)
    nbdet = len(bdet)    
    for m in mtags:
        # Divide by the number of detections
        if ngdet>0: obj[gdet][m] /= obj[gdet]['ndet']
        if nbdet>0: obj[bdet][m] = 99.99   # no good detections

    # Get the average error
    metags = ['asemierr','bsemierr','thetaerr']
    for m in metags:
        # Just take the sqrt to complete the addition in quadrature
        if ngdet>0: obj[gdet][m] = np.sqrt(obj[gdet][m])
        if nbdet>0: obj[bdet][m] = 99.99

    # Add E(B-V)
    print('Getting E(B-V)')
    #glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
    #obj['ebv'] = DUST_GETVAL(glon,glat,/noloop,/interp)

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
    newobjindex = np.zeros(nobj,long)-1    # new indices
    newobjindex[ind1] = np.arange(nmatch)
    # Keep the objects inside the Healpix
    obj = obj[ind1]
    print(str(nmatch)+' final objects fall inside the pixel')

    # Remove trimmed objects from IDSTR
    totrim, = np.where(objtokeep[idstr.objectindex] eq 0,ntotrim)  #using old index
    if ntotrim>0:
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
    dum, uiexpnum = np.unique(np.long(idstr['expnum']),return_index=True)
    uexpnum = np.long(idstr[uiexpnum]['expnum'])
    nuexpnum = len(uexpnum)
    ind1,ind2 = utils.match(np.long(allmeta.expnum),uexpnum)
    nmatch = len(ind1)
    sumstr = Table(allmeta[ind1])
    col_nobj = Column(name='nobjects', dtype=np.long)
    col_healpix = Column(name='healpix', dtype=np.long)
    sumstr.add_columns([col_nobj col_healpix])
    sumstr['nobjects'] =
    sumstr['healpix'] = pix
    # get number of objects per exposure
    expnum = np.long(idstr.expnum)
    siexp = np.sort(expnum)
    expnum = expnum[siexp]
    if nuexpnum>1:
        brklo, = np.where(expnum != np.roll(expnum,1))
        nbrk = len(brklo)
        brkhi = np.hstack((brklo[1:nbrk-1],len(expnum)))
        numobjexp = brkhi-brklo+1
    else:
        numobjexp=len(expnum)
    ind1,ind2 = utils.match(np.long(sumstr.expnum),uexpnum)
    nmatch = len(ind1)
    sumstr[ind1]['nobjects'] = numobjexp

    # Write the output file
    print('Writing combined catalog to '+outfile)
    if os.path.exists(outdir) is False: os.mkdir(outdir)
    if os.path.exists(outdir+'/'+subdir) is False: os.mkdir(outdir+'/'+subdir)
    if os.path.exists(outfile): os.delete(outfile)
    sumstr.write(outfile)               # first, summary table
    #  append other fits binary tables
    hdulist = fits.open(outfile)
    hdu = fits.table_to_hdu(obj)        # second, catalog
    hdulist.append(hdu)
    hdu = fits.table_to_hdu(idstr)      # third, ID table
    hdulist.append(hdu)    
    hdulist.writeto(outfile,overwrite=True)
    hdulist.close()
    #MWRFITS,sumstr,outfile,/create      # first, summary table
    #MWRFITS,obj,outfile,/silent         # second, catalog
    #MWRFITS,idstr,outfile,/silent       # third, ID table
    if os.path.exists(outfile+'.gz'): os.delete(outfile)    

    ret = subprocess.call(['gzip',outfile])    # compress final catalog

    dt = time.time()-t0
    print('dt = '+str(dt)+' sec.')
