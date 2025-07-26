import os
import numpy as np
from glob import glob
from dlnpyutils import utils as dln,coords
from astropy.table import Table,vstack
from astropy.io import fits
import healpy as hp
from . import utils

def coverage(pix,version='v4',clobber=False):
    """
    Make the coverage map for a single nside=128 NSC healpix
    """

    # Combine all of the data
    dldir,mssdir,localdir = utils.nscrootdirs()
    dir = '/net/dl2/dnidever/nsc/instcal/'+version+'/'
    nside = 128
    nside2 = 4096
    radeg = 180.0d0 / !dpi

    # Does the coverage map already exist
    covfile = dir+'combine/coverage/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'_coverage.fits'
    if os.path.exists(covfile) and clobber==False:
        print(covfile,' EXISTS and clobber==False')
        return

    # Get healpix boundary coordinates
    PIX2VEC_RING,nside,pix,vec,vertex
    vertex = transpose(reform(vertex))  # [1,3,4] -> [4,3]
    VEC2ANG,vec,hcendec,hcenra,/astro
    VEC2ANG,vertex,hdec,hra,/astro

    # Rotate to tangent plane
    hlon,hlat = coords.rotsphcen(hra,hdec,hcenra,hcendec,gnomic=True)
    mmhlon = [np.min(hlon),np.max(hlon)]
    mmhlat = [np.min(hlat),np.max(hlat)]

    # Get the pixel numbers for nside=4096 healpix that are within
    # this larger pixel
    QUERY_POLYGON,nside2,vertex,listpix,nlistpix

    # Check that they belong inside this healpix
    PIX2ANG_RING,nside2,listpix,listtheta,listphi
    ANG2PIX_RING,nside,listtheta,listphi,listpix1
    bdpix, = np.where(listpix1 != pix)

    step = 100
    v = hp.pix2vec(nside,pix)
    radius = hp.nside2resol(nside)
    pix2 = hp.query_disc(nside2,v,radius=2*radius)
    theta,phi = hp.pix2ang(nside2,pix2)
    pix1 = hp.ang2pix(nside,theta,phi)
    gd = (pix1 == pix)
    pix = pix2[gd]
    nlistpix = len(listpix)


    # Initialize the coverage structure
    print('Creating coverage structure for pixel ',pix)
    dt = [('pix',int),('pix128',int),('ra',float),('dec',float),('nobj',int),('coverage',float),
          ('nexp',int),('ucoverage',float),('unexp',int),('udepth',float),
          ('gcoverage',float),('gnexp',int),('gdepth',float),('rcoverage',float),('rnexp',int),('rdepth',float),
          ('icoverage',float),('inexp',int),('idepth',float),('zcoverage',float),('znexp',int),('zdepth',float),
          ('ycoverage',float),('ynexp',int),('ydepth',float),('vrcoverage',float),('vrnexp',int),('vrdepth',float)]
    covtab = np.zeros(nlistpix,dtype=np.dtype(dt))
    covtab['pix'] = listpix
    PIX2ANG_RING,nside2,covtab.pix,theta,phi
    covtab['ra'] = phi*radeg
    covtab['dec'] = 90-theta*radeg
    covtags = covtab.dtype.names

    # Does the combined object file exist?
    objfile = dir+'combine/'+str(int(pix)//1000)+'/'+str(pix)+'.fits.gz'
    if os.path.exists(objfile)==False:
        print(objfile,' NOT FOUND')
        goto,SAVEFILE

    # Check that it's not blank
    head0 = fits.getheader(objfile,1,errmsg=errmsg)
    if errmsg != '':
        print(objfile,' IS BLANK')
        goto,SAVEFILE

    # Load the list of exposures
    exptab = Table.read(objfile,1)
    nexptab = len(exptab)
    exptab['file'] = exptab['file'].strip()
    exptab['base'] = exptab['base'].strip()
    exptab['dateobs'] = exptab['dateobs'].strip()
    exptab['filter'] = exptab['filter'].strip()
    exptab['success'] = False
    exptab['hoverlap'] = False
    exptab['chipindx'] = -1

    # Load the object table
    obj = Table.read(objfile,2)
    #ophi = np.deg2rad(obj['ra'])
    #otheta = np.deg2rad((90-obj['dec']))
    obj_pix4096 = hp.ang2pix(nside2,obj['ra'],obj['dec'],nest=False,lonlat=True)
    #ANG2PIX_RING,nside2,otheta,ophi,obj_pix4096

    # Load the chip summary structure for each exposure
    print('Loading the chip summary information')
    print(' Number                Exposure          Nchips  Noverlap')
    nchtab = 0
    cnt = 0
    for i in range(nexptab):
        # Construct metadata filename
        base = exptab['base'][i]
        instrument = 'c4d'
        if base.find('k4m') > -1:
            instrument = 'k4m'
        if base.find('ksb') > -1:
            instrument = 'ksb'
        dateobs = exptab['dateobs'][i]
        night = dateobs[:4]+dateobs[5:7]+dateobs[8:10]
        #night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
        metafile = dir+instrument+'/'+night+'/'+base+'/'+base+'_meta.fits'

        chtab1 = None
        if os.path.exists(metafile):
            # Load the chip summary structure
            chtab1 = Table.read(metafile,2)
            nchtab1 = len(chtab1)
            chtab1['expdir'] = chtab1['expdir'].strip()
            chtab1['filename'] = chtab1['filename'].strip()
            chtab1['hoverlap'] = False      # does it overlap healpix
            chtab1['filter'] = exptab['filter'][0].strip()
            exptab['success'][i] = True

            # Check if the chips overlap the healpix
            for j in range(nchtab1):
                lon,lat = coords.rotsphcen(chtab1['vra'][j],chtab1['vdec'][j],hcenra,hcendec,gnomic=True)
                chtab1['hoverlap'][j] = DOPOLYGONSOVERLAP(hlon,hlat,lon,lat)      
            exptab['hoverlap'][i] = np.max(chtab1.hoverlap)

            # Only keep overlapping chips
            gdchtab1, = np.where(chtab1['hoverlap']==True)
            # ,ngdchtab1,comp=bdchtab1,ncomp=nbdchtab1)
            print('{:5d} {:35s} {:6d} {:6d}'.format(i+1,base,nchtab1,ngdchtab1))
            if ngdchtab1 == 0:
                print('No chips overlap this healpix')
                exptab['hoverlap'][i] = 0
                exptab['nchips'][i] = 0
                exptab['chipindx'][i] = -1
                goto,BOMB
            # Remove non-overlapping chips
            if nbdchtab1 > 0:
                chtab1 = np.delete(chtab1,bdchtab1)
                nchtab1 = n_elements(chtab1)
            exptab['chipindx'][i] = cnt
            exptab['nchips'][i] = nchtab1

            # Start CHTAB structure
            if nchtab == 0:
                chtab = np.zeros(10000,dtype=np.dtype(dt))
                nchtab = len(chtab)

            # Add new elements
            if cnt+nchtab1 > nchtab:
                oldchtab = chtab.copy()
                chtab = np.zeros(nchtab+10000,dtype=np.dtype(dt))
                chtab[:nchtab] = oldchtab
                nchtab = len(chtab)
                del oldchtab

            # Stuff into the CHTAB structure
            newchtab1 = np.zeros(nchtab1,dtype=np.dtype(dt))
            for c in chtab1.dtype.names:
                newchtab1[c] = chtab1[c]
            #struct_assign,chtab1,newchtab1
            chtab[cnt:cnt+nchtab1] = newchtab1  
            cnt += nchtab1

        # Metadata file not found
        else:
            print(metafile,' NOT FOUND')


    # Trim off the extra elements
    chtab = chtab[:cnt]
    nchtab = len(chtab)

    # Trim off any exposures with no overlapping chips
    exptab0 = exptab.copy()
    bdexp, = np.where(exptab['hoverlap']==False)
    if len(bdexp) > 0:
        exptab = np.delete(exptab,bdexp)

    filters = ['u','g','r','i','z','Y','VR']
    nfilters = len(filters)

    # Loop over the small pixels and figure out the coverage and depth
    # and number of exposures
    for i in range(nlistpix):
        # Get the objects for this pixel
        _,pind1,pind2 = np.intersect1d(listpix[i],obj_pix4096,return_indices=True)
        nobjmatch = len(pind1)
        #MATCH,listpix[i],obj_pix4096,pind1,pind2,/sort,count=nobjmatch
        covtab['nobj'][i] = nobjmatch

        # Get healpix boundary coordinates
        vec1 = hp.pix2vec(nside2,listpix[i],nest=False)
        hcenra1,hcendec1 = hp.vec2ang(vec1,lonlat=True)
        hra1,hdec1 = hp.vec2ang(vertex1,lonlat=True)
        #PIX2VEC_RING,nside2,listpix[i],vec1,vertex1
        #vertex1 = transpose(reform(vertex1))  # [1,3,4] -> [4,3]
        #VEC2ANG,vec1,hcendec1,hcenra1,/astro
        #VEC2ANG,vertex1,hdec1,hra1,/astro

        # Rotate to tangent plane
        hlon1,hlat1 = coords.rotsphcen(hra1,hdec1,hcenra1,hcendec1,gnomic=True)

        dx = 1e-3  # gives ~20x20 pixel
        rlon = [np.min(hlon1),np.max(hlon1)]
        rlat = [np.min(hlat1),np.max(hlat1)]
        lon0 = rlon[0]
        lat0 = rlat[0]
        hx = (hlon1-lon0)/dx
        hy = (hlat1-lat0)/dx
        nx = int(np.ceil((rlon[1]-lon0)/dx))+1
        ny = int(np.ceil((rlat[1]-lat0)/dx))+1

        # Mask image for which pixels are in the healpix region
        mask = np.zeros((ny,nx),int)
        inm = polyfillv(hx,hy,nx,ny)
        mask[inm] = 1
        maskpix = np.sum(mask)

        # Loop over each filter
        for f in range(nfilters):
            filtind, = np.where(chstr['filter'] == filters[f])
            nfiltind = len(filtind)

            # Get columns for this filter
            covcol = filters[f]+'coverage'
            numcol = filters[f]+'nexp'
            depcol = filters[f]+'depth'

            # Coverage map for this filter
            numim = np.zeros((ny,nx),int)
            depthim = np.zeros((ny,nx),float)

            # Now loop over each chip
            alloverlap = 0
            for c in range(nfiltind):
                j = filtind[c]
                lon1,lat1 = coords.rotsphcen(chstr['vra'][j],chstr['vdec'][j],hcenra1,hcendec1,gnomic=True)
                # Check if they overlap
                overlap = DOPOLYGONSOVERLAP(hlon1,hlat1,lon1,lat1)
                if overlap == 1:
                    # Get the chip overlap image
                    # Transform to pixel-based tangent plane system
                    vx = (lon1-lon0)/dx
                    vy = (lat1-lat0)/dx
                    # Pixels in the chip
                    cin = POLYFILLV(vx,vy,nx,ny)
                    cmask = np.zeros((ny,nx),int)
                    if len(cin) > 1 or cin[0] != -1:
                        cmask[cin] = 1  # some good overlap
                    cmask = cmask * mask   # pixels in healpix region

                    numim[:,:] += cmask                        # number of chips that overlap
                    depthim[:,:] += cmask*chstr['depth95'][j]  # sum of depth of all pixels

            # Calculate coverage
            gdpix, = np.where(numim > 0)
            ngdpix = len(gdpix)
            overlapfrac = float(ngdpix) / maskpix
            # Calculate average depth image
            mndepthim = depthim/(numim > 1)
            if ngdpix > 0:
                depth = np.median(mndepthim[gdpix])
            else:
                depth = -9999.0
            # Number of chips
            nchipoverlap = np.max(numim)

            # Stuff into the coverage structure
            covtab[covcol][i] = overlapfrac
            covtab[depcol][i] = depth
            covtab[numcol][i] = nchipoverlap

            # Add to total coverage and exposure for this pixel
            covtab['coverage'][i] >= overlapfrac
            covtab['nexp'][i] += nchipoverlap
            
    # Save the coverage map
    print('Writing coverage information to ',covfile)
    if os.makedirs(os.path.dirname(covfile),exist_ok=True)
    covtab.write(covfile,overwrite=True)
