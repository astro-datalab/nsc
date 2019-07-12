#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack
from astropy.time import Time
import healpy as hp
import utils

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
    buffsize = 10.0/3600. ; in deg
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
        print('  FILTER='+meta['filter']+'  EXPTIME='+str(meta['exptime'])+' sec'

        # Loop over the chip files
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
                  bd, = np.where(vra gt 180)
                  if len(bd)>0: vra[bd] -= 360
              if dopolygonsoverlap(buffdict['ra'],buffdict['dec'],vra,vdec):
                  #print,'This chip does NOT overlap the HEALPix region+buffer'
                  #goto,BOMB1

              # Load the chip-level catalog
              chfile = file_dirname(list[i].file)+'/'+list[i].base+'_'+strtrim(chmeta[j].ccdnum,2)+'_meas.fits'
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
#GOT TO HERE
              if running_gdl() eq 0 then begin
                  ROI_CUT,buffer.lon,buffer.lat,lon,lat,ind0,ind1,fac=1000,/silent
                  nmatch = n_elements(ind1)
              endif else begin
                  ; first use WHERE with X/Y limits
                  in1 = where(lon ge buffer.lr[0] and lon le buffer.lr[1] and $
                              lat ge buffer.br[0] and lat le buffer.br[1],nin1)
                  if nin1 gt 0 then begin
                      inmask = INSIDE(lon[in1],lat[in1],buffer.lon,buffer.lat)
                      in2 = where(inmask eq 1,nin2)
                      if nin2 gt 0 then ind1=in1[in2] else undefine,ind1
                      nmatch = nin2
                  endif else begin
                      undefine,ind1
                      nmatch = 0
              # Only want source inside this pixel
              if nmatch eq 0:
                  print,'  No sources inside this pixel'
                  goto,BOMB
              print,'  ',strtrim(nmatch,2),' sources are inside this pixel'
              cat1 = cat1[ind1]
              ncat1 = nmatch

              # Combine the catalogs
              if n_elements(cat) eq 0:
                  cat_schema = cat1[0]
                  struct_assign,{dum:''},cat_schema
                  cat = replicate(cat_schema,total(chmeta.nsources))
              cat[catcount:catcount+ncat1-1] = cat1
              catcount += ncat1

              BOMB1:
        if n_elements(cat) gt 0 then cat = cat[0:catcount-1]  ; trim excess
        ncat = n_elements(cat)
        if ncat eq 0:
              print,'This exposure does NOT cover the HEALPix'
              goto,BOMB

        # Add metadata to ALLMETA
        #  Make sure it's in the right format
        newmeta = {file:'',base:'',expnum:0L,ra:0.0d0,dec:0.0d0,dateobs:'',mjd:0.0d0,filter:'',exptime:0.0,$
             airmass:0.0,nsources:0L,fwhm:0.0,nchips:0L,badchip31:0B,rarms:0.0,decrms:0.0,ebv:0.0,gaianmatch:0L,$
             zpterm:0.0,zptermerr:0.0,zptermsig:0.0,nrefmatch:0L}
        STRUCT_ASSIGN,meta,newmeta
        PUSH,allmeta,newmeta
