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
#import subprocess
import time
from argparse import ArgumentParser
import socket
#from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
#from sklearn.cluster import DBSCAN
#from scipy.optimize import least_squares
#from scipy.interpolate import interp1d

# Combine data for one NSC healpix region
if __name__ == "__main__":
    parser = ArgumentParser(description='Combine NSC Instcal Catalogs.')
    parser.add_argument('--version', type=str, default='v3', help='Version number')
    parser.add_argument('--makelist', action='store_true', help='Make healpix list')
    parser.add_argument('-r','--redo', action='store_true', help='Redo this HEALPIX')
    parser.add_argument('--nmulti', type=int, default=20, help='Number of jobs to run')
    parser.add_argument('--nocuts', action='store_true', help='Do not apply any quality cuts')
    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = args.version
    redo = args.redo
    makelist = args.makelist
    nmulti = args.nmulti
    nocuts = args.nocuts
    nside = 128
    radeg = 180 / np.pi

    # on thing/hulk use
    if (host == "thing") or (host == "hulk"):
        basedir = "/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/mss1/"
        localdir = "/d0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"
    # on gp09 use
    if (host == "gp09") or (host == "gp08") or (host == "gp07") or (host == "gp06") or (host == "gp05"):
        basedir = "/net/dl1/users/dnidever/nsc/instcal/"+version+"/"
        mssdir = "/net/mss1/"
        localdir = "/data0/"
        tmproot = localdir+"dnidever/nsc/instcal/"+version+"/tmp/"

    t0 = time.time()

    # Combine all of the data
    if ~os.path.exists(basedir+'combine'): os.mkdir(basedir+'combine/')
    if ~os.path.exists(basedir+'combine/logs/'): os.mkdir(basedir+'combine/logs/')
    if ~os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/'): os.mkdir(localdir+'dnidever/nsc/instcal/'+version+'/')
    plotsdir = basedir+'plots/'
    if ~os.path.exists(plotsdir): os.mkdir(plotsdir)


    # Log file
    #------------------
    # format is nsc_combine_main.DATETIME.log
    jd = systime(/julian)
    caldat,jd,month,day,year,hour,minute,second
    smonth = strtrim(month,2)
    if month lt 10 then smonth = '0'+smonth
    sday = strtrim(day,2)
    if day lt 10 then sday = '0'+sday
    syear = strmid(strtrim(year,2),2,2)
    shour = strtrim(hour,2)
    if hour lt 10 then shour='0'+shour
    sminute = strtrim(minute,2)
    if minute lt 10 then sminute='0'+sminute
    ssecond = strtrim(round(second),2)
    if second lt 10 then ssecond='0'+ssecond
    logfile = dir+'combine/logs/nsc_instcal_combine_main.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
    #JOURNAL,logfile

    print("Combining NOAO InstCal catalogs")

    #goto,STARTRUNNING

    # Restore the calibration summary file
    temp = fits.getdata(basedir+'lists/nsc_instcal_calibrate.fits',1)
    schema = temp[0]
    struct_assign,{dum:''},schema
    schema = create_struct(schema,'chipindx',-1LL,'NGOODCHIPWCS',0)
    calstr = replicate(schema,n_elements(temp))
    struct_assign,temp,calstr,/nozero
    # Add WCSCAL and TELSTAT information
    add_tag,calstr,'wcscal','',calstr
    add_tag,calstr,'telstat','',calstr
    coords = fits.getdata(basedir+'lists/allcoords.fits',1)
    fluxfile = calstr['file']
    g, = np.where(strmid(fluxfile,0,4) eq '/net',ng)
    if ng>0: fluxfile[g]=strmid(fluxfile[g],4)
    MATCH,fluxfile,coords.file,ind1,ind2,/sort
    calstr['wcscal'][ind1] = coords['wcscal'][ind2]    # Failed (3153), Poor (14), Successful (308190)
    calstr['telstat'][ind1] = coords['telstat'][ind2]  # NAN (68188), Not (1222), Track (241826), UNKNOWN (116), Unknown (5)
    # the 2054 failed exposures did not match b/c no fluxfile info
    # Only want exposures with successful SE processing
    gd, = np.where(calstr['success']==1)
    ncalstr = len(gd)
    calstr = calstr[gd]
    si = np.argsort(calstr['expdir'])
    calstr = calstr[si]
    chstr = fits.getdata(basedir+'lists/nsc_instcal_calibrate.fits',2)
    nchstr = len(chstr)
    # Get indices for CHSTR
    siexp = np.argsort(chstr['expdir'])
    chstr = chstr[siexp]
    expdir = chstr['expdir']
    brklo, = np.where(expdir != np.roll(expdir,1))
    nbrk = len(brklo)
    brkhi = [brklo[1:nbrk-1]-1,len(expdir)-1]
    nchexp = brkhi-brklo+1
    if ncalstr==len(brklo):
        Exception('number of exposures in CALSTR and CHSTR do not match')
    calstr['chipindx'] = brklo
    calstr['nchips'] = nchexp
    # Getting number of good chip WCS for each exposures
    for i in range(len(calstr): calstr[i].ngoodchipwcs = total(chstr[brklo[i]:brkhi[i]].ngaiamatch gt 0)
    # Fixing absolute paths of flux filename
    file = calstr['file']
    g1, = np.where(stregex(file,'/net/mss1/',/boolean)==1)
    ng1 = len(g1)
    if ng1>0: file[g1] = strmid(file[g1],10)
    g2, = np.where(stregex(file,'/mss1/',/boolean)==1)
    ng2 = len(g2)
    if ng2>0: file[g2] = strmid(file[g2],6)
    # Fixing very negative RAs
    print('FIXING NEGATIVE RAs in CALSTR and CHSTR')
    #bdra, = np.where(chstr.cenra lt -180,nbdra)
    bdra, = np.where(chstr['cenra']<0)
    nbdra = len(bdra)
    dum,uibd = np.unique(chstr['expdir'][bdra],return_indices=True)
    ind1,ind2 = utils.match(calstr['expdir'],chstr['expdir'][bdra[uibd]])
    nmatch = len(ind1)
    for i in range(nmatch):
        ind3,ind4 = utils.match(chstr['expdir'][bdra],calstr['expdir'][ind1[i]])
        # Fix CALSTR RA
        chra = chstr[bdra[ind3]].cenra
        bd1, = np.where(chra lt -180,nbd1)
        if nbd1>0: chra[bd1]+=360
        cenra = np.mean(utils.minmax(chra))
        if cenra<0: cenra+=360
        calstr['ra'][ind1[i]] = cenra
        # Fix CHSTR CENRA
        bd2, = np.where(chra<0)
        nbd2 = len(bd2)
        if nbd2>0: chra[bd2]+=360
        chstr['cenra']][bdra[ind3]] = chra
        # Fix CHSTR VRA
        vra = chstr['vra'][bdra[ind3]]
        bd3, = np.where(vra<0)
        nbd3 = len(bd3)
        if nbd3>0: vra[bd3]+=360
        chstr['vra'][bdra[ind3]] = vra

    # Fix instrument in STR and CHSTR
    print('FIXING INSTRUMENT IN STR AND CHSTR')
    type = ['c4d','k4m','ksb']
    for i=0,len(type)-1:
        gd, = np.where(stregex(calstr.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
        if ngd gt 0 then calstr[gd].instrument=type[i]
        gd, = np.where(stregex(chstr.expdir,'/'+type[i]+'/',/boolean) eq 1,ngd)
        if ngd gt 0 then chstr[gd].instrument=type[i]

    ## Fix missing AIRMASS
    #bdam, = np.where(str.airmass lt 0.9,nbdam)
    #for i=0,nbdam-1 do begin
    #  type = ['c4d','k4m','ksb']
    #  obs = ['ctio','kpno','kpno']
    #  MATCH,str[bdam[i]].instrument,type,ind1,ind2,/sort
    #  obsname = obs[ind2]
    #  OBSERVATORY,obsname,obstr
    #  lat = obstr.latitude
    #  lon = obstr.longitude
    #  jd = date2jd(str[bdam[i]].dateobs)
    #  ra = str[bdam[i]].ra
    #  dec = str[bdam[i]].dec
    #  str[bdam[i]].airmass = AIRMASS(jd,ra,dec,lat,lon)
    #endfor
    # THIS IS STILL RETURNING -1, IS ONE OF THE VALUES WRONG??


    # APPLY RELEASE-DATE CUTS
    list1 = fits.getdata(basedir+'lists/decam_instcal_list.fits',1)
    list2 = fits.getdata(basedir+'lists/mosaic3_instcal_list.fits',1)
    list3 = fits.getdata(basedir+'lists/bok90prime_instcal_list.fits',1)
    elist = np.hstack((list1,list2,list3))
    fluxfile = strmid(elist['fluxfile'],10)
    ind1,ind2 = utils.match(fluxfile,fil
    MATCH,fluxfile,file,ind1,ind2,/sort,count=nmatch
    # some don't match because they were from a previous version
    #  of the input list
    release_date = strarr(len(calstr))+'2020-01-01 00:00:00'
    release_date[ind2] = list[ind1].release_date
    release_year = int(strmid(release_date,0,4))
    release_month = int(strmid(release_date,5,2))
    release_day = int(strmid(release_date,8,2))
    release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
    #release_cutoff = [2017,4,24]  # v1 - April 24, 2017
    release_cutoff = [2017,10,11]  # v2 - Oct 11, 2017
    release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
    gdrelease, = np.where(release_mjd le release_cutoff_mjd,ngdrelease,comp=bdrelease,ncomp=nbdrelease)
    print(str(ngdrelease)+' exposures are PUBLIC')
    calstr = calstr[gdrelease]  # impose the public data cut

    # Zero-point structure
    zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},10)
    zpstr[0:6].instrument = 'c4d'
    zpstr[0:6].filter = ['u','g','r','i','z','Y','VR']
    zpstr[0].amcoef = [-1.60273, -0.375253]   # c4d-u
    zpstr[1].amcoef = [0.277124, -0.198037]   # c4d-g
    zpstr[2].amcoef = [0.516382, -0.115443]   # c4d-r
    zpstr[3].amcoef = [0.380338, -0.067439]   # c4d-i
    zpstr[4].amcoef = [0.123924, -0.096877]   # c4d-z
    zpstr[5].amcoef = [-1.06529, -0.051967]   # c4d-Y
    zpstr[6].amcoef = [1.004357, -0.081105]   # c4d-VR
    # Mosiac3 z-band
    zpstr[7].instrument = 'k4m'
    zpstr[7].filter = 'z'
    zpstr[7].amcoef = [-2.687201, -0.73573]   # k4m-z
    # Bok 90Prime, g and r
    zpstr[8].instrument = 'ksb'
    zpstr[8].filter = 'g'
    zpstr[8].amcoef = [-2.859646, -1.40837]   # ksb-g
    zpstr[9].instrument = 'ksb'
    zpstr[9].filter = 'r'
    zpstr[9].amcoef = [-4.008771, -0.25718]   # ksb-r
    nzpstr = len(zpstr)

    #STOP,'DOUBLE-CHECK THESE ZERO-POINTS!!!'

    # APPLY QA CUTS IN ZEROPOINT AND SEEING
    if ~nocuts:
        print('APPLYING QA CUTS')
        #fwhmthresh = 3.0  # arcsec, v1
        fwhmthresh = 2.0  # arcsec, v2
        #filters = ['u','g','r','i','z','Y','VR']
        #nfilters = len(filters)
        #zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
        #zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        badzpmask = bytarr(len(calstr)) + 1

        for i in range(nzpstr):
            ind, = np.where(calstr.instrument eq zpstr[i].instrument and calstr.filter eq zpstr[i].filter and calstr.success eq 1,nind)
            print(zpstr['instrument'][i]+'-'+zpstr['filter'][i]+' '+str(nind)+' exposures')
            if nind>0:
                calstr1 = calstr[ind]
                zpterm = calstr1.zpterm
                bdzp, = np.where(finite(zpterm) eq 0,nbdzp)  # fix Infinity/NAN
                if nbdzp gt 0 then zpterm[bdzp] = 999999.9
                am = calstr1.airmass
                mjd = calstr1.mjd
                bdam, = np.where(am lt 0.9,nbdam)
                if nbdam gt 0 then am[bdam] = median(am)
                glactc,calstr1.ra,calstr1.dec,2000.0,glon,glat,1,/deg

                # Measure airmass dependence
                gg0, = np.where(abs(zpterm) lt 50 and am lt 2.0,ngg0)
                coef0 = robust_poly_fitq(am[gg0],zpterm[gg0],1)
                zpf = poly(am,coef0)
                sig0 = mad(zpterm[gg0]-zpf[gg0])
                gg, = np.where(abs(zpterm-zpf) lt (3.5*sig0 > 0.2),ngg)
                coef = robust_poly_fitq(am[gg],zpterm[gg],1)
                print(zpstr['instrument'][i]+'-'+zpstr['filter'][i]+' '+str(coef))
                # Trim out bad exposures to determine the correlations and make figures
                gg, = np.where(np.abs(zpterm-zpf) lt (3.5*sig0 > 0.2) and calstr1.airmass lt 2.0 and calstr1.fwhm lt 2.0 and calstr1.rarms lt 0.15 and $
                           calstr1.decrms lt 0.15 and calstr1.success eq 1 and calstr1.wcscal eq 'Successful' and calstr1.zptermerr lt 0.05 and $
                           calstr1.zptermsig lt 0.08 and (calstr1.ngoodchipwcs eq calstr1.nchips) and $
                           (calstr1.instrument ne 'c4d' or calstr1.zpspatialvar_nccd le 5 or (calstr1.instrument eq 'c4d' and calstr1.zpspatialvar_nccd gt 5 and calstr1.zpspatialvar_rms lt 0.1)) and $
                           np.abs(glat) gt 10 and calstr1.nrefmatch gt 100 and calstr1.exptime ge 30,ngg)

                # Zpterm with airmass dependence removed
                relzpterm = zpterm + 25   # 25 to get "absolute" zpterm
                relzpterm -= zpstr[i].amcoef[1]*(am-1)

                # CURRENTLY K4M/KSB HAVE EXPTIME-DEPENDENCE IN THE ZEROPOINTS!!
                if zpstr[i].instrument eq 'k4m' or zpstr[i].instrument eq 'ksb':
                    print('REMOVING EXPTIME-DEPENDENCE IN K4M/KSB ZEROPOINTS!!!')
                    relzpterm += 2.5*alog10(calstr1.exptime)

                # Fit temporal variation in zpterm
                mjd0 = 56200L
                xx = calstr1[gg].mjd-mjd0
                yy = relzpterm[gg]
                invvar = 1.0/calstr1[gg].zptermerr^2
                nord = 3
                bkspace = 200 #20
                sset1 = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=yfit1)
                sig1 = mad(yy-yfit1)
                gd, = np.where(yy-yfit1 gt -3*sig1,ngd)      
                # refit
                sset = bspline_iterfit(xx[gd],yy[gd],invvar=invvar[gd],nord=nord,bkspace=bkspace)
                yfit = bspline_valu(xx,sset)
                allzpfit = bspline_valu(calstr1.mjd-mjd0,sset)

                # Make some figures
                # ZPterm vs. airmass
                file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_airmass'
                ps_open,file,/color,thick=4,/encap
                hess,am[gg],relzpterm[gg],dx=0.01,dy=0.02,xr=[0.9,2.5],yr=[-0.5,0.5]+median(relzpterm[gg]),xtit='Airmass',ytit='Zero-point',$
                     tit=zpstr[i].instrument+'-'+zpstr[i].filter
                x = scale_vector(findgen(100),0.5,2.0)
                oplot,x,poly(x,coef),co=250
                ps_close
                ps2png,file+'.eps',/eps
                # ZPterm vs. time (density)
                file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_time_density'
                ps_open,file,/color,thick=4,/encap
                hess,calstr1[gg].mjd-mjd0,relzpterm[gg],dx=2,dy=0.02,yr=[-0.5,0.5]+median(relzpterm[gg]),xtit='Time (days)',ytit='Zero-point',$
                     tit=zpstr[i].instrument+'-'+zpstr[i].filter
                oplot,calstr1[gg].mjd-mjd0,allzpfit[gg],ps=1,sym=0.3,co=250
                xyouts,50,-0.45+median(relzpterm[gg]),'MJD!d0!n = '+strtrim(mjd0,2),align=0,charsize=1.2
                ps_close
                ps2png,file+'.eps',/eps
                # ZPterm vs. time (points)
                file = plotsdir+zpstr[i].instrument+'-'+zpstr[i].filter+'_zpterm_time'
                ps_open,file,/color,thick=4,/encap
                plot,calstr1[gg].mjd-mjd0,relzpterm[gg],ps=1,sym=0.5,yr=[-0.5,0.5]+median(relzpterm[gg]),xs=1,ys=1,xtit='Time (days)',ytit='Zero-point',$
                     tit=zpstr[i].instrument+'-'+zpstr[i].filter,thick=1
                oplot,calstr1[gg].mjd-mjd0,allzpfit[gg],ps=1,sym=0.3,co=250
                xyouts,50,-0.45+median(relzpterm[gg]),'MJD!d0!n = '+strtrim(mjd0,2),align=0,charsize=1.2
                ps_close
                ps2png,file+'.eps',/eps

                # Remove temporal variations to get residual values
                relzpterm -= allzpfit

                # Find the GOOD exposures
                #------------------------
                # We are using ADDITIVE zpterm 
                #  calmag = instmag + zpterm
                # if there are clouds then instmag is larger/fainter
                #  and zpterm is smaller (more negative)
                #bdind, = np.where(calstr[ind].zpterm-medzp lt -zpthresh[i],nbdind)
                gdind, = np.where(relzpterm ge -zpstr[i].thresh and relzpterm le zpstr[i].thresh,ngdind,comp=bdind,ncomp=nbdind)
                print('  '+str(nbdind)+' exposures with ZPTERM below the threshold')
                if ngdind>0: badzpmask[ind[gdind]] = 0

        # Get bad DECaLS and SMASH exposures
        badexp = bytarr(len(calstr))
        READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/smash_badexposures.txt',smashexpnum,format='A',comment='#',/silent
        MATCH,int(calstr.expnum),int(smashexpnum),ind1,ind2,/sort,count=nmatch
        if nmatch>0:
            badexp[ind1] = 1
            badexp[ind1] = badexp[ind1] AND (calstr[ind1].instrument eq 'c4d')   # make sure they are DECam exposures
        READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/decals_bad_expid.txt',decalsexpnum,format='A',comment='#',/silent
        MATCH,int(calstr.expnum),int(decalsexpnum),ind1,ind2,/sort,count=nmatch
        if nmatch>0:
            badexp[ind1] = 1
            badexp[ind1] = badexp[ind1] AND (calstr[ind1].instrument eq 'c4d')   # make sure they are DECam exposures
        READCOL,'/home/dnidever/projects/noaosourcecatalog/obslog/mzls_bad_expid.txt',mzlsexpnum,format='A',comment='#',/silent
        MATCH,int(calstr.expnum),int(mzlsexpnum),ind1,ind2,/sort,count=nmatch
        if nmatch>0:
            badexp[ind1] = 1
            badexp[ind1] = badexp[ind1] AND (calstr[ind1].instrument eq 'k4m')   # make sure they are Mosaic3 exposures

        # Final QA cuts
        #  Many of the short u-band exposures have weird ZPTERMs, not sure why
        #  There are a few exposures with BAD WCS, RA>360!
        bdexp, = np.where(calstr.success eq 0 or $                              # SE failure
                      calstr.wcscal ne 'Successful' or $                    # CP WCS failure
                      calstr.fwhm gt fwhmthresh or $                        # bad seeing
                      calstr.ra gt 360 or $                                 # bad WCS/coords
                      calstr.rarms gt 0.15 or calstr.decrms gt 0.15 or $       # bad WCS
                      badzpmask eq 1 or $                                # bad ZPTERM
                      calstr.zptermerr gt 0.05 or $                         # bad ZPTERMERR
                      calstr.nrefmatch lt 5 or $                            # few phot ref match
                      badexp eq 1 or $                                   # bad SMASH/LS exposure
                      #calstr.ngoodchipwcs lt calstr.nchips or $                # not all chips astrom calibrated
                      (calstr.instrument eq 'c4d' and calstr.zpspatialvar_nccd gt 5 and calstr.zpspatialvar_rms gt 0.1),nbdexp)  # bad spatial zpterm
        # rarms/decrms, nrefmatch
        print('QA cuts remove '+str(nbdexp)+' exposures')
        # Remove
        torem = bytarr(nchstr)
        for i=0,nbdexp-1 do torem[calstr[bdexp[i]].chipindx:calstr[bdexp[i]].chipindx+calstr[bdexp[i]].nchips-1]=1
        bdchstr, = np.where(torem eq 1,nbdchstr)
        REMOVE,bdchstr,chstr
        REMOVE,bdexp,calstr
        # Get new CHIPINDEX values
        #   make two arrays of old and new indices to transfer 
        #   the new index values into an array with the size of
        #   the old CHSTR
        trimoldindex = lindgen(nchstr)                    # index into original array, but "bad" ones removed/trimed
        remove,bdchstr,trimoldindex
        trimnewindex = lindgen(len(trimoldindex))  # new index of trimmed array
        newindex = lonarr(nchstr)-1
        newindex[trimoldindex] = trimnewindex             # new index in original array
        newchipindex = newindex[calstr.chipindx]
        str.chipindx = newchipindex
        ncalstr = len(calstr)

        # SHOULD INCLUDE CUTS ON ZTERMERR OR NPHOTMATCH
        #STOP,'SHOULD INCLUDE CUTS ON ZTERMERR OR NPHOTMATCH'


    #STARTRUNNING:

    # CREATE LIST OF HEALPIX AND OVERLAPPING EXPOSURES
    # Which healpix pixels have data
    listfile = basedir+'lists/nsc_healpix_list.fits'
    if keyword_set(makelist) or file_test(listfile) eq 0:
        print('Finding the Healpix pixels with data')
        radius = 1.1
        healstr = replicate({file:'',base:'',pix:0L},1e5)
        nhealstr = len(healstr)
        cnt = 0LL
        for i in range(ncalstr):
            if i mod 1e3 eq 0 then print(str(i))
            theta = (90-calstr[i].dec)/radeg
            phi = calstr[i].ra/radeg
            ANG2VEC,theta,phi,vec
            QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

            # Use the chip corners to figure out which ones actually overlap
            chstr1 = chstr[calstr[i].chipindx:calstr[i].chipindx+calstr[i].nchips-1]
            #  rotate to tangent plane so it can handle RA=0/360 and poles properly
            ROTSPHCEN,chstr1.vra,chstr1.vdec,calstr[i].ra,calstr[i].dec,vlon,vlat,/gnomic
            #  loop over healpix
            overlap = bytarr(nlistpix)
            for j in range(nlistpix):
                PIX2VEC_RING,nside,listpix[j],vec,vertex
                vertex = transpose(reform(vertex))  # [1,3,4] -> [4,3]
                VEC2ANG,vertex,hdec,hra,/astro
                ROTSPHCEN,hra,hdec,calstr[i].ra,calstr[i].dec,hlon,hlat,/gnomic
                #  loop over chips
                for k=0,calstr[i].nchips-1 do overlap[j] >= DOPOLYGONSOVERLAP(hlon,hlat,vlon[*,k],vlat[*,k])
            # Only keep the healpix with real overlaps
            gdlistpix, = np.where(overlap eq 1,ngdlistpix)
            if ngdlistpix>0:
                listpix = listpix[gdlistpix]
                nlistpix = ngdlistpix
            else:
                del(listpix)
                nlistpix = 0

            if nlistpix eq 0 then stop,'No healpix for this exposure.  Something is wrong!'

            # Add new elements to array
            if (cnt+nlistpix)>nhealstr:
                old = healstr
                healstr = replicate({file:'',base:'',pix:0L},nhealstr+1e4)
                healstr[0:nhealstr-1] = old
                nhealstr += 1e4
                undefine,old

            # Add to the structure
            healstr[cnt:cnt+nlistpix-1].file = calstr[i].expdir+'/'+calstr[i].base+'_cat.fits'
            healstr[cnt:cnt+nlistpix-1].base = calstr[i].base
            healstr[cnt:cnt+nlistpix-1].pix = listpix
            cnt += nlistpix

        # Trim extra elements
        healstr = healstr[0:cnt-1]
        nhealstr = len(healstr)

        # Get uniq pixels
        ui = uniq(healstr.pix,sort(healstr.pix))
        upix = healstr[ui].pix
        nupix = len(upix)
        print(calstr(nupix)+' Healpix pixels have overlapping data')

        # Get start/stop indices for each pixel
        idx = sort(healstr.pix)
        healstr = healstr[idx]
        q = healstr.pix
        lo, = np.where(q ne shift(q,1),nlo)
        #hi, = np.where(q ne shift(q,-1))
        hi = [lo[1:nlo-1]-1,nhealstr-1]
        nexp = hi-lo+1
        index = replicate({pix:0L,lo:0L,hi:0L,nexp:0L},nupix)
        index.pix = upix
        index.lo = lo
        index.hi = hi
        index.nexp = nexp
        npix = len(index)

        # Replace /net/dl1/ with /dl1/ so it will work on all machines
        healstr.file = repstr(healstr.file,'/net/dl1/','/dl1/')

        # Write the full list plus an index
        print('Writing list to '+listfile)
        MWRFITS,healstr,listfile,/create
        MWRFITS,index,listfile,/silent
        # Copy to local directory for faster reading speed
        file_copy,listfile,localdir+'dnidever/nsc/instcal/'+version+'/',/over
        # PUT NSIDE IN HEADER!!

    # Using existing list
    else:
        print('Reading list from '+listfile)
        healstr = fits.getdata(listfile,1)
        index = fits.getdata(listfile,2)
        upix = index['pix']
        npix = len(index)
        # Copy to local directory for faster reading speed
        file_copy,listfile,localdir+'dnidever/nsc/instcal/'+version+'/',/over


    # Load the list of healpix pixels for this server to be run LOCALLY
    pixfile = basedir+'lists/combine_pix_'+host+'.txt'
    READLINE,pixfile,pixlist,count=npixlist
    rnd = sort(randomu(1,npixlist))   # RANDOMIZE!!
    pixlist = int(pixlist[rnd])
    print('Running '+str(npixlist)+' jobs on '+host+' with nmult='+str(nmulti))
    cmd = "nsc_instcal_combine,"+strtrim(pixlist,2)+",nside="+strtrim(nside,2)+",version='"+version+"',/local,/filesexist"
    if keyword_set(redo) then cmd+=',/redo'
    cmddir = strarr(npixlist)+localdir+'dnidever/nsc/instcal/'+version+'/tmp/'

    # Now run the combination program on each healpix pixel
    a = '' & read,a,prompt='Press RETURN to start'
    PBS_DAEMON,cmd,cmddir,jobs=jobs,/hyperthread,/idle,prefix='nsccmb',nmulti=nmulti,wait=1


    ## Make the commands
    #cmd = "nsc_instcal_combine,"+strtrim(index.pix,2)+",nside="+strtrim(nside,2)+",version='"+version+"'"
    #if keyword_set(redo) then cmd+=',/redo'
    #cmddir = strarr(npix)+localdir+'dnidever/nsc/instcal/'+version+'/tmp/'

    ## Check if the output file exists
    #if not keyword_set(redo) then begin
    #  outfiles = dir+'combine/'+strtrim(upix/1000,2)+'/'+strtrim(upix,2)+'.fits.gz'
    #  test = file_test(outfiles)
    #  gd, = np.where(test eq 0,ngd,comp=bd,ncomp=nbd)
    #  if nbd gt 0 then begin
    #    print,strtrim(nbd,2),' files already exist and /redo not set.'
    #  endif 
    #  if ngd eq 0 then begin
    #    print,'No files to process'
    #    return
    #  endif
    #  print,strtrim(ngd,2),' files left to process'
    #  cmd = cmd[gd]
    #  cmddir = cmddir[gd]
    #endif

    ## Prioritize longest-running jobs FIRST
    ## Use prediction program
    #PIX2ANG_RING,nside,index.pix,theta,phi
    #ra = phi*radeg
    #dec = 90-theta*radeg
    #glactc,ra,dec,2000.0,glon,glat,1,/deg
    #dt = predictcombtime(glon,glat,index.nexp)
    ## Do the sorting
    #hsi = reverse(sort(dt))
    #cmd = cmd[hsi]
    #cmddir = cmddir[hsi]
    #dt = dt[hsi]
    #index = index[hsi]

    # Divide into three using total times
    #tot = total(dt>10)
    #totcum = total(dt>10,/cum)
    #print,min(where(totcum ge tot/3))
    #print,min(where(totcum ge 2*tot/3))

    #ncmd = len(cmd)
    #nhalf = ncmd/2

    ## Randomize 1st half for hulk/thing/gp09
    #cmd1 = cmd[0:(nhalf-1)]
    #cmdadir1 = cmddir[0:(nhalf-1)]
    #pix1 = index[0:(nhalf-1)].pix
    #index1 = index[0:(nhalf-1)]
    ## now randomize
    #rnd = sort(randomu(1,len(cmd1)))
    #cmd1 = cmd1[rnd]
    #cmddir1 = cmddir1[rnd]
    #pix1 = pix1[rnd]
    #index1 = index1[rnd]

    # Slice it up
    ## hulk, 1st
    ##cmd = cmd[0:(nhalf-1):3]
    ##cmddir = cmddir[0:(nhalf-1):3]   
    ##pix = index[0:(nhalf-1):3].pix
    #cmd = cmd1[0:(nhalf/3)-1]
    #cmddir = cmddir1[0:(nhalf/3)-1]
    #pix = pix1[0:(nhalf/3)-1]

    # thing, 2nd
    ##cmd = cmd[1:(nhalf-1):3]
    ##cmddir = cmddir[1:(nhalf-1):3]
    ##pix = index[1:(nhalf-1):3].pix
    #cmd = cmd1[(nhalf/3):(2*nhalf/3)-1]
    #cmddir = cmddir1[(nhalf/3):(2*nhalf/3)-1]
    #pix = pix1[(nhalf/3):(2*nhalf/3)-1]

    # gp09, 3rd
    ##cmd = cmd[2:(nhalf-1):3]
    ##cmddir = cmddir[2:(nhalf-1):3]
    ##pix = index[2:(nhalf-1):3].pix
    #cmd = cmd1[(2*nhalf/3):*]
    #cmddir = cmddir1[(2*nhalf/3):*]
    #pix = pix1[(2*nhalf/3):*]

    # gp05
    #cmd = cmd[nhalf:*:4]
    #cmddir = cmddir[nhalf:*:4]
    #pix = index[nhalf:*:4].pix

    # gp06
    #cmd = cmd[nhalf+1:*:4]
    #cmddir = cmddir[nhalf+1:*:4]
    #pix = index[nhalf+1:*:4].pix

    # gp07
    #cmd = cmd[nhalf+2:*:4]
    #cmddir = cmddir[nhalf+2:*:4]
    #pix = index[nhalf+2:*:4].pix

    # gp08
    #cmd = cmd[nhalf+3:*:4]
    #cmddir = cmddir[nhalf+3:*:4]
    #pix = index[nhalf+3:*:4].pix

    ## Prioritize longest-running jobs FIRST
    ## Load the DECam run times
    #sum1 = mrdfits(dir+'nsccmb_summary_hulk.fits',1)
    #sum2 = mrdfits(dir+'nsccmb_summary_thing.fits',1)
    #sum3 = mrdfits(dir+'nsccmb_summary_gp09.fits',1)
    #sum = [sum1,sum2,sum3]
    #si = sort(sum.mtime)
    #sum = sum[si]
    ## only keep fairly recent ones
    #gd, = np.where(sum.mtime gt 1.4897704e+09,ngd)
    #sum = sum[gd]
    ## Deal with duplicates
    #dbl = doubles(sum.pix,count=ndbl)
    #alldbl = doubles(sum.pix,/all,count=nalldbl)
    #torem = bytarr(nalldbl)
    #for i=0,ndbl-1 do begin
    #  MATCH,sum[alldbl].pix,sum[dbl[i]].pix,ind1,ind2,/sort,count=nmatch
    #  torem[ind1[0:nmatch-2]] = 1
    #endfor
    #bd=where(torem eq 1,nbd)
    #remove,alldbl[bd],sum
    #dt = lonarr(len(index))-1
    #MATCH,index.pix,sum.pix,ind1,ind2,/sort,count=nmatch
    #dt[ind1] = sum[ind2].dt
    ## Do the sorting
    #hsi = reverse(sort(dt))
    #cmd = cmd[hsi]
    #cmddir = cmddir[hsi]
    #dt = dt[hsi]
    #
    ## Divide into three using total times
    #tot = total(dt>10)
    #totcum = total(dt>10,/cum)
    #print,min(where(totcum ge tot/3))
    #print,min(where(totcum ge 2*tot/3))

    ## Start with healpix with low NEXP and far from MW midplane, LMC/SMC
    #pix2ang_ring,nside,index.pix,theta,phi
    #pixra = phi*radeg
    #pixdec = 90-theta*radeg
    #glactc,pixra,pixdec,2000.0,pixgl,pixgb,1,/deg
    #cel2lmc,pixra,pixdec,palmc,radlmc
    #cel2smc,pixra,pixdec,rasmc,radsmc
    #gdpix, = np.where(index.nexp lt 50 and np.abs(pixgb) gt 10 and radlmc gt 5 and radsmc gt 5,ngdpix)
    #
    #outfile = dldir+'users/dnidever/nsc/instcal/combine/'+strtrim(index.pix,2)+'.fits'

    # Now run the combination program on each healpix pixel
    PBS_DAEMON,cmd,cmddir,jobs=jobs,/hyperthread,/idle,prefix='nsccmb',nmulti=nmulti,wait=1

    # RUN NSC_COMBINE_SUMMARY WHEN IT'S DONE!!!

    ## Load all the summary/metadata files
    #print,'Creating Healpix summary file'
    #sumstr = replicate({pix:0L,nexposures:0L,nobjects:0L,success:0},nupix)
    #sumstr.pix = upix
    #for i=0,nupix-1 do begin
    #  if (i+1) mod 5000 eq 0 then print,i+1
    #  file = dir+'combine/'+strtrim(upix[i],2)+'.fits'
    #  if file_test(file) eq 1 then begin
    #    meta = MRDFITS(file,1,/silent)
    #    sumstr[i].nexposures = len(meta)
    #    hd = headfits(file,exten=2)
    #    sumstr[i].nobjects = sxpar(hd,'naxis2')
    #    sumstr[i].success = 1
    #  endif else begin
    #    sumstr[i].success = 0
    #  endelse
    #endfor
    #gd, = np.where(sumstr.success eq 1,ngd)
    #print,strtrim(ngd,2),' Healpix successfully processed'
    #print,'Writing summary file to ',dir+'combine/nsc_instcal_combine.fits'
    #MWRFITS,sumstr,dir+'combine/nsc_instcal_combine.fits',/create

    # End logfile
    #------------
    #JOURNAL
