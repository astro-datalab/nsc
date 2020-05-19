#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
import healpy as hp
from dlnpyutils import utils as dln, coords, bindata, db
from astropy.coordinates import SkyCoord
from scipy.optimize import least_squares
from scipy.interpolate import interp1d


def varmetric(inpmeas):
    """ Compute photometric variability metrics."""
    # meas is a catalog of measurements for a single object

    nmeas = len(inpmeas)
    
    # Need the catalog to be a numpy array
    if isinstance(inpmeas,np.ndarray):
        meas = inpmeas
    else:
        meas = np.array(inpmeas)
        
    filtcol = 'FILTER'
    if filtcol not in meas.dtype.names:
        filtcol = 'filter'
    if filtcol not in meas.dtype.names:
        raise ValueError('No filter column')
    magcol = 'MAG_AUTO'
    if magcol not in meas.dtype.names:
        magcol = 'mag_auto'
    if magcol not in meas.dtype.names:
        raise ValueError('No mag_auto column')
    errcol = 'MAGERR_AUTO'
    if errcol not in meas.dtype.names:
        errcol = 'magerr_auto'
    if errcol not in meas.dtype.names:
        raise ValueError('No magerr_auto column')
    mjdcol = 'MJD'
    if mjdcol not in meas.dtype.names:
        mjdcol = 'mjd'
    if mjdcol not in meas.dtype.names:
        raise ValueError('No mjd column')
    

    # OBJ schema
    dtype_obj = np.dtype([('deltamjd',np.float32),('ndet',np.int16),('nphot',np.int16),
                          ('ndetu',np.int16),('nphotu',np.int16),('umag',np.float32),('urms',np.float32),('uerr',np.float32),
                          ('ndetg',np.int16),('nphotg',np.int16),('gmag',np.float32),('grms',np.float32),('gerr',np.float32),
                          ('ndetr',np.int16),('nphotr',np.int16),('rmag',np.float32),('rrms',np.float32),('rerr',np.float32),
                          ('ndeti',np.int16),('nphoti',np.int16),('imag',np.float32),('irms',np.float32),('ierr',np.float32),
                          ('ndetz',np.int16),('nphotz',np.int16),('zmag',np.float32),('zrms',np.float32),('zerr',np.float32),
                          ('ndety',np.int16),('nphoty',np.int16),('ymag',np.float32),('yrms',np.float32),('yerr',np.float32),
                          ('ndetvr',np.int16),('nphotvr',np.int16),('vrmag',np.float32),('vrrms',np.float32),('vrerr',np.float32),
                          ('rmsvar',np.float32),('madvar',np.float32),('iqrvar',np.float32),('etavar',np.float32),
                          ('jvar',np.float32),('kvar',np.float32),('chivar',np.float32),('romsvar',np.float32),
                          ('variable10sig',np.int16),('nsigvar',np.float32)])

    
    obj = np.zeros(1,dtype=dtype_obj)


    # Initialize the OBJ structured array
    obj = np.zeros(1,dtype=dtype_obj)
    # all bad to start
    for f in ['rmsvar','madvar','iqrvar','etavar','jvar','kvar','chivar','romsvar']: obj[f]=np.nan
    for f in ['u','g','r','i','z','y','vr']:
        obj[f+'mag'] = 99.99
        obj[f+'err'] = 9.99
        obj[f+'rms'] = np.nan

    obj['ndet'] = nmeas
    obj['deltamjd'] = np.max(meas[mjdcol])-np.min(meas[mjdcol])
    
    # Mean magnitudes
    # Convert totalwt and totalfluxwt to MAG and ERR
    #  and average the morphology parameters PER FILTER
    filtindex = dln.create_index(meas[filtcol].astype(np.str))
    nfilters = len(filtindex['value'])
    resid = np.zeros(nmeas)+np.nan     # residual mag
    relresid = np.zeros(nmeas)+np.nan  # residual mag relative to the uncertainty
    for f in range(nfilters):
        filt = filtindex['value'][f].lower()
        findx = filtindex['index'][filtindex['lo'][f]:filtindex['hi'][f]+1]
        obj['ndet'+filt] = filtindex['num'][f]
        gph,ngph = dln.where(meas[magcol][findx]<50)
        obj['nphot'+filt] = ngph
        if ngph==1:
            obj[filt+'mag'] = meas[magcol][findx[gph]]
            obj[filt+'err'] = meas[errcol][findx[gph]]
        if ngph>1:
            newmag, newerr = dln.wtmean(meas[magcol][findx[gph]], meas[errcol][findx[gph]],magnitude=True,reweight=True,error=True)
            obj[filt+'mag'] = newmag
            obj[filt+'err'] = newerr
            # Calculate RMS
            obj[filt+'rms'] = np.sqrt(np.mean((meas[magcol][findx[gph]]-newmag)**2))
            # Residual mag
            resid[findx[gph]] = meas[magcol][findx[gph]]-newmag
            # Residual mag relative to the uncertainty
            #  set a lower threshold of 0.02 in the uncertainty
            relresid[findx[gph]] = np.sqrt(ngph/(ngph-1)) * (meas[magcol][findx[gph]]-newmag)/np.maximum(meas[errcol][findx[gph]],0.02)

    # Calculate variability indices
    gdresid = np.isfinite(resid)
    ngdresid = np.sum(gdresid)
    if ngdresid>0:
        resid2 = resid[gdresid]
        sumresidsq = np.sum(resid2**2)
        tsi = np.argsort(meas[mjdcol][gdresid])
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
        obj['rmsvar'] = rms
        obj['madvar'] = madvar
        obj['iqrvar'] = iqrvar
        obj['etavar'] = etavar

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
        obj['jvar'] = jvar
        obj['kvar'] = kvar
        #obj['avgrelvar'] = avgrelvar
        obj['chivar'] = chivar
        obj['romsvar'] = romsvar
        #if chivar>50: import pdb; pdb.set_trace()

    # Make NPHOT from NPHOTX
    obj['nphot'] = obj['nphotu']+obj['nphotg']+obj['nphotr']+obj['nphoti']+obj['nphotz']+obj['nphoty']+obj['nphotvr']

    # Fiducial magnitude, used to select variables below
    #  order of priority: r,g,i,z,Y,VR,u
    if obj['nphot']>0:
        magarr = np.zeros(7,float)
        for ii,nn in enumerate(['rmag','gmag','imag','zmag','ymag','vrmag','umag']): magarr[ii]=obj[nn]
        gfid,ngfid = dln.where(magarr<50)
        if ngfid>0: fidmag=magarr[gfid[0]]

        
    return obj


def selectvariables(obj):
    """ Select variables using photometric variability indices."""

    nobj = len(obj)
    fidmag = np.zeros(nobj,float)+np.nan  # fiducial magnitude
    # Loop over the objects
    for i in range(nobj):
        # Fiducial magnitude, used to select variables below
        #  order of priority: r,g,i,z,Y,VR,u
        if obj['nphot'][i]>0:
            magarr = np.zeros(7,float)
            for ii,nn in enumerate(['rmag','gmag','imag','zmag','ymag','vrmag','umag']): magarr[ii]=obj[nn][i]
            gfid,ngfid = dln.where(magarr<50)
            if ngfid>0: fidmag[i]=magarr[gfid[0]]

    
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


    return obj
