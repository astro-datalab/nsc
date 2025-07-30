import os
import numpy as np
from dlnpyutils import utils as dln
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from glob import glob
import healpy as hp

def get_sorted_indices(input_list):
    """
    Returns a list of indices that would sort the input_list in ascending order.

    Args:
        input_list: The list to be sorted by index.

    Returns:
        A list of integers representing the sorted indices.
    """
    # Create a list of (value, original_index) tuples
    indexed_list = []
    for i, value in enumerate(input_list):
        indexed_list.append((value, i))

    # Sort the list of tuples based on the value (first element of the tuple)
    indexed_list.sort()

    # Extract the original indices from the sorted tuples
    sorted_indices = [index for value, index in indexed_list]

    return sorted_indices

def makelist():

    tab = Table.read('decam_instcal_list_exptime10sec_20240727.fits.gz')
    # 614776

    lines = dln.readlines('/corral/projects/NOIRLab/nsc/instcal/v4/c4d/allmeas_expdir.txt')
    # 632166

    left = dln.readlines('/corral/projects/NOIRLab/nsc/instcal/v4/c4d/allmeas_expdir_left_072825.txt')
    # 46558
    lbase = [os.path.basename(f) for f in left]
    _,ind1,ind2 = np.intersect1d(tab['base'],lbase,return_indices=True)
    # 38611

    # Remove duplicates from "left" list
    left = np.unique(left)
    lbase = [os.path.basename(f) for f in left]
    # 46498
    # Remove non-standard filters from "left" list

    filt = []
    for b in lbase:
        if b[:3]=='c4d':
            filt.append(b.split('_')[-2])
        else:
            filt.append(b)
    filt = np.array(filt)
    #u   1100
    #g   8003
    #r  10000
    #i   6060
    #z   9587
    #Y   2546
    #VR  2076
    ll = [len(f) for f in filt]
    ll = np.array(ll)

    gd, = np.where(ll<=2)
    gdleft = left[gd]
    gdleft = np.unique(gdleft)
    # 39312
    gdlbase = [os.path.basename(f) for f in gdleft]
    gdlbase = np.array(gdlbase)

    # we have multiple versions of the same exposure, get the best one
    gdlbase2 = []
    for b in gdlbase:
        if b[:3]=='c4d':
            gdlbase2.append('_'.join(b.split('_')[:-1]))
        else:
            gdlbase2.append(b)
    gdlbase2 = np.array(gdlbase2)

    #In [81]: len(np.unique(gdlbase2))
    #Out[81]: 38001

    index = dln.create_index(gdlbase2)

    for i in range(len(index['value'])):
        ind = index['index'][index['lo'][i]:index['hi'][i]+1]
        nind = len(ind)
        if nind > 1:
            bases = gdlbase[ind]

    dd = []
    plver = []
    for i in range(len(tab)): 
        if tab['base'][i][:3]=='c4d': 
            dd.append(tab['base'][i].split('_')[-1]) 
            plver.append(tab['plver'][i])
    dd = np.array(dd)
    plver = np.array(plver)
    _,ui = np.unique(dd,return_index=True)



    uindex = dln.create_index(dd)
    uplver = len(uindex['value'])*[None]
    udd = uindex['value']
    for i in range(len(uindex['value'])):
        ind = uindex['index'][uindex['lo'][i]:uindex['hi'][i]+1]
        nindex = len(ind)
        plver1 = [plver[ii] for ii in ind]
        plver1.sort()
        mplver1 = plver1[-1]
        uplver[i] = mplver1
        print(i+1,uindex['value'][i],mplver1)

    si = get_sorted_indices(uplver)
    for i in range(len(uplver)):
        print(i+1,udd[si[i]],uplver[si[i]])

    # plver reduction versions for each CP extension type
    verplver = {item[0]:item[1] for item in zip(udd,uplver)}
    verplver['ls10.2'] = 'V5.2.4LS'

    # get the best reduction version for each exposure
    gdldd = []
    gdlubase = []
    gdlplver = []
    for i in range(len(gdlbase)):
        gdlubase.append('_'.join(gdlbase[i].split('_')[:-1]))
        gdldd.append(gdlbase[i].split('_')[-1])
        gdlplver.append(verplver[gdlbase[i].split('_')[-1]])
    gdlubase = np.array(gdlubase)
    gdldd = np.array(gdldd)
    gdlplver = np.array(gdlplver)
    gdlubase_index = dln.create_index(gdlubase)
    gdfinalbase = []
    gdfinalleft = []
    for i in range(len(gdlubase_index['value'])):
        ind = gdlubase_index['index'][gdlubase_index['lo'][i]:gdlubase_index['hi'][i]+1]
        nind = len(ind)
        if nind>1:
            gdldd1 = gdldd[ind]
            gdlplver1 = gdlplver[ind]
            si1 = get_sorted_indices(gdlplver1)
            useind1 = ind[si1[-1]]
            gdfinalbase.append(gdlbase[useind1])
            gdfinalleft.append(gdleft[useind1])
        else:
            gdfinalbase.append(gdlbase[ind[0]])
            gdfinalleft.append(gdleft[ind[0]])


    _,ind1,ind2 = np.intersect1d(tab['base'],gdfinalbase,return_indices=True)
    # 34750

    # Get information for all of our unique exposures
    tab['ubase'] = tab['base'].copy()
    for i in range(len(tab)):
        if tab['base'][i][:3].strip()=='c4d':
            tab['ubase'][i] = '_'.join(tab['base'][i].strip().split('_')[:-1])

    _,ind1,ind2 = np.intersect1d(tab['base'],gdfinalbase,return_indices=True)
    dt = [('base',str,50),('expdir',str,300),('ra',float),('dec',float),
          ('filter',str,10),('exptime',float),('npix',int)]
    out = np.zeros(len(gdfinalleft),dtype=np.dtype(dt))
    out['base'] = gdfinalbase
    out['expdir'] = gdfinalleft
    out['ra'][ind2] = tab['ra'][ind1]
    out['dec'][ind2] = tab['dec'][ind1]
    out['exptime'][ind2] = tab['exposure'][ind1]
    for i in range(len(gdfinalleft)):
        out['filter'][i] = out['base'][i].split('_')[-2]
        if out['exptime'][i]==0 or out['ra'][i]<0:
            ubase1 = '_'.join(out['base'][i].strip().split('_')[:-1])
            ind, = np.where(tab['ubase']==ubase1)
            if len(ind)>0:
                out['ra'][i] = tab['ra'][ind[0]]
                out['dec'][i] = tab['dec'][ind[0]]
                out['exptime'][i] = tab['exposure'][ind[0]]
        if out['exptime'][i]==0 or np.isfinite(out['ra'][i])==False or out['ra'][i]<0:
            headfile = os.path.join(out['expdir'][i],out['base'][i]+'_header.fits')
            if os.path.exists(headfile):
                head = fits.getheader(headfile,1)
                out['ra'][i] = head['crval1']
                out['dec'][i] = head['crval2']
                out['exptime'][i] = head['exptime']
                if out['exptime'][i]==0.0:
                    out['exptime'][i] = head['EXPREQ']
            if os.path.exists(headfile)==False:
                # use instcal files on tacc                                                        
                headfile = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/images/'
                headfile += '/'.join(out['expdir'][i].split('/')[-4:])+'.fits.fz'
                if os.path.exists(headfile):
                    head = fits.getheader(headfile,0)
                    ra = head.get('ra')
                    dec = head.get('dec')
                    if ra is None:
                        ra = head.get('telra')
                        dec = head.get('teldec')
                    if int(ra.split(':')[0])>24:
                        ra = '{:02d}:'.format(int(ra.split(':')[0])-24)+':'.join(ra.split(':')[1:])
                    coo = SkyCoord(ra,dec,unit=('hour','deg'),frame='icrs')
                    out['ra'][i] = coo.ra.degree
                    out['dec'][i] = coo.dec.degree
                    out['exptime'][i] = head['exptime']
                    if out['exptime'][i]==0.0:
                        out['exptime'][i] = head['EXPREQ']
            if os.path.exists(headfile)==False:
                # sometimes there's a different version, i.e. _d2 instead of _ls11                 
                base2 = '_'.join(out['base'][i].split('_')[:-1])
                night = '20'+base2[4:10]
                headfile = '/corral/projects/NOIRLab/nsc/instcal/v4/header'
                headfile = glob(os.path.join(headfile,'c4d',night,base2+'_*.hdr'))
                if len(headfile)>0:
                    headlines = dln.readlines(headfile[0])
                    lo = dln.grep(headlines,'^SIMPLE  =',index=True)
                    begind = dln.grep(headlines,'^XTENSION',index=True)
                    begind = lo+begind
                    endind = dln.grep(headlines,'^END',index=True)
                    hlines = headlines[begind[0]:endind[0]+1]
                    head = fits.Header.fromstring('\n'.join(hlines),sep='\n')
                    ra = head['ra']
                    if int(ra.split(':')[0])>24:
                        ra = '{:02d}:'.format(int(ra.split(':')[0])-24)+':'.join(ra.split(':')[1:])
                    coo = SkyCoord(ra,head['dec'],unit=('hour','deg'),frame='icrs')
                    out['ra'][i] = coo.ra.degree
                    out['dec'][i] = coo.dec.degree
                    out['exptime'][i] = head['exptime']
                    if out['exptime'][i]==0.0:
                        out['exptime'][i] = head['EXPREQ']

        print(i+1,out['base'][i])

    out['npix'] = hp.ang2pix(128,out['ra'],out['dec'],nest=False,lonlat=True)

    # there are some short exposures
    gd, = np.where(out['exptime']>=10)
    # 36334 of 38001
    Table(out).write('healpix_list_exptime10sec_left_073025.fits',overwrite=True)

    import pdb; pdb.set_trace()
