import os
import numpy as np
from glob import glob
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from dlnpyutils import utils as dln,coords
import healpy as hp
import time
import traceback
import subprocess
from . import utils

def measure(version='v4',nosources=False,quick=False):
    """ Make the nsc_measure_summary.fits summary file """

    dldir,mssdir,localdir = utils.rootdirs()
    basedir = os.path.join(dldir,'instcal/',version)

    t0 = time.time()

    # Find all of the directories
    print('Getting the exposure directories')
    c4d_expdirs = dln.readlines(basedir+'/c4d/allmeas_expdir.txt')
    k4m_expdirs = []
    ksb_expdirs = []
    #c4d_expdirs = glob(os.path.join(basedir,'c4d/20??/20??????/*'))
    #c4d_expdirs = [d for d in c4d_expdirs if os.path.isdir(d)]
    #k4m_expdirs = glob(os.path.join(basedir,'k4m/20??/20??????/*'))
    #k4m_expdirs = [d for d in k4m_expdirs if os.path.isdir(d)]
    #ksb_expdirs = glob(os.path.join(basedir,'ksb/20??/20??????/*'))
    #ksb_expdirs = [d for d in ksb_expdirs if os.path.isdir(d)]
    expdirs = c4d_expdirs+k4m_expdirs+ksb_expdirs

    nexpdirs = len(expdirs)
    print(nexpdirs,'exposure directories')

    # Create table
    dtyp = [('dir',str,200),('instrument',str,3),('base',str,50),('measfile',str,200),('measexists',bool),
            ('nchips',int),('nmeas',int),('logexists',bool),('success',bool),
            ('runtime',float),('logdate',int)]
    exptab = np.zeros(nexpdirs,dtype=np.dtype(dtyp))
    exptab['dir'] = expdirs

    instrument = []
    if len(c4d_expdirs)>0:
        instrument += len(c4d_expdirs)*['c4d']
    if len(k4m_expdirs)>0:
        instrument += len(c4d_expdirs)*['k4m']
    if len(ksb_expdirs)>0:
        instrument += len(c4d_expdirs)*['ksb']
    exptab['instrument'] = instrument
    exptab['base'] = [os.path.basename(d) for d in expdirs]

    # Loop through the exposure directories
    for i in range(len(expdirs)):
        if i % 5000 == 0:
            print(i)
        dir1 = exptab['dir'][i]
        base1 = exptab['base'][i]
        measfile = os.path.join(dir1,base1+'_meas.fits')
        exptab['measfile'][i] = measfile
        exptab['measexists'][i] = os.path.exists(measfile)

        ## Get chip files
        #if quick:
        #    chipfiles1 = glob(dir1+'/'+base1+'_[1-9].fits')
        #    chipfiles2 = glob(dir1+'/'+base1+'_[1-9][0-9].fits')
        #    chipfiles = chipfiles1+chipfiles2
        #    nchipfiles = len(chipfiles)
        #    exptab['nchips'][i] = nchipfiles
        #else:
        #    nchipfiles = 0
        # First chip date
        #if nchipfiles > 0 and quick==False:
        #    exptab['chip1date'][i] = os.path.getmtime(chipfiles[0])
            
        # It succeeded if the final log file exists
        logfile = os.path.join(dir1,base1+'.log')
        exptab['logexists'][i] = os.path.exists(logfile)
        if exptab['logexists'][i] and quick==False:
            exptab['logdate'][i] = os.path.getmtime(logfile)
        if exptab['measexists'][i] and exptab['logexists'][i]:
            exptab['success'][i] = True
        if exptab['logexists'][i]:
            loglines = dln.readlines(logfile)
            try:
                date1 = Time(loglines[0][:19])
                date2 = Time(loglines[-1][:19])
                dtime = (date2-date1).sec
                exptab['runtime'][i] = dtime
            except:
                pass
        ## Success, have logfile and chip files
        #if quick==False:
        #    if nchipfiles > 0 and exptab['logfile_success'][i]:
        #        exptab['success'][i] = True
        #else:
        #    # just check if logfile and first fits catalog exist
        #    if os.path.exists(logfile) and os.path.exists(dir1+'/'+base1+'_1.fits'):
        #        exptab['success'][i] = True
        # dt, need chip files
        #if nchipfiles > 0 and quick==False:
        #    mtime1 = os.path.getmtime(chipfiles[0])
        #    mtime2 = os.path.getmtime(chipfiles[-1])
        #    dt = mtime2-mtime1
        #    # this is the time between the end of first chip and last chip
        #    #  correct for that
        #    exptab['dt'][i] = dt * nchipfiles / (nchipfiles-1.0)
        if exptab['success'][i] and nosources==False and quick==False:
            hdu = fits.open(measfile)
            nobj = np.zeros(len(hdu),int)
            for j in range(len(hdu)-1):
                nobj[j] = hdu[j+1].header['naxis2']
            hdu.close()
            exptab['nchips'][i] = len(hdu)-1
            nsrc = np.sum(nobj)
            exptab['nmeas'][i] = nsrc
            #lines = dln.readlines(logfile)
            #glines = dln.grep(lines,'sextracted',index=True)
            #for j in range(len(glines)):
            #    line1 = lines[glines[j]]
            #    arr = line1.split('\\')
            #    g2 = dln.grep(arr,'Objects: detected',index=True)
            #    if len(g2)==0:
            #        continue
            #    line2 = arr[g2[0]]
            #    pos = line2.find('sextracted')
            #    if pos == -1:
            #        continue
            #    nsrc = int(line2[pos+10:])
            #    exptab['nsources'][i] += nsrc

        #import pdb; pdb.set_trace()

    # Write out the file
    outfile = os.path.join(basedir,'lists','nsc_measure_summary.fits')
    print('Writing summary file to ',outfile)
    exptab = Table(exptab)
    exptab.write(outfile,overwrite=True)

    print('dt = {:.1f} sec.'.format(time.time()-t0))

def calibrate(version='v4'):
    """ Make the nsc_calibrate_summary.fits summary file """

    dldir,mssdir,localdir = utils.rootdirs()
    basedir = os.path.join(dldir,'instcal/',version)

    t0 = time.time()

    # Find all of the directories
    print('Getting the exposure directories')    

    c4d_expdirs = dln.readlines(basedir+'/c4d/allmeas_expdir.txt')
    k4m_expdirs = []
    ksb_expdirs = []
    #c4d_expdirs = glob(os.path.join(basedir,'c4d/20??/20??????/*'))
    #c4d_expdirs = [d for d in c4d_expdirs if os.path.isdir(d)]
    #k4m_expdirs = glob(os.path.join(basedir,'k4m/20??/20??????/*'))
    #k4m_expdirs = [d for d in k4m_expdirs if os.path.isdir(d)]
    #ksb_expdirs = glob(os.path.join(basedir,'ksb/20??/20??????/*'))
    #ksb_expdirs = [d for d in ksb_expdirs if os.path.isdir(d)]
    expdirs = c4d_expdirs+k4m_expdirs+ksb_expdirs

    nexpdirs = len(expdirs)
    print(nexpdirs,'exposure directories')

    ## Load the list
    #listtab = Table.read(basedir+'/lists/nsc_calibrate_healpix_list.fits')
    #nlist = len(listtab)
    #print(nlist,'exposures to check')
    #listtab['expdir'] = [str(e).strip() for e in listtab['expdir']]
    #listtab['instrument'] = [str(e).strip() for e in listtab['instrument']]
    #listtab['filter'] = [str(f).strip() for f in listtab['filter']]

    # Chip offsets
    #chipoff = Table.read(dldir+'dnidever/nsc/instcal/decam_chip_xyoff.fits')

    #meassumfile = ''
    #msumtab = Table.read(meassumfile)

    #dtyp = [(),(),()]
    #exptab = np.zeros(len(files,dtype=np.dtype(dtyp)))

    outfile = os.path.join(basedir,'lists','nsc_calibrate_summary.fits')

    # Load the exposure and chip meta files
    calbasedir = '/home1/09970/dnidever/scratch1/nsc/instcal/v4'
    expdata = []
    chipdata = []
    for i in range(len(expdirs)):
        dir1 = expdirs[i]
        arr = dir1.split('/')
        base = arr[-1]
        night = arr[-2]
        year = arr[-3]
        instrument = arr[-4]

        measfile = os.path.join(calbasedir,instrument,year,night,base,base+'_meas.fits')
        metafile = os.path.join(calbasedir,instrument,year,night,base,base+'_meta.fits')
        if os.path.exists(metafile):
            hdu = fits.open(metafile)
            exptab1 = Table(hdu[1].data)
            # fix names
            for c in ['file','maskfile','wtfile']:
                ffile = exptab1[c][0]
                if ffile.find('INFO')>-1:
                    exptab1[c] = ffile[ffile.find(']')+1:].strip()
            exptab1['measfile'] = measfile
            chiptab1 = []
            for j in range(len(hdu)-2):
                chiptab1.append(Table(hdu[j+2].data))
            chiptab1 = vstack(chiptab1)
            hdu.close()
            chiptab1['file'] = exptab1['file'][0]
            chiptab1['wtfile'] = exptab1['wtfile'][0]
            chiptab1['maskfile'] = exptab1['maskfile'][0]
            chiptab1['measfile'] = measfile
            chiptab1['base'] = base
            expdata.append(exptab1)
            chipdata.append(chiptab1)
            print(i+1,base,len(chiptab1),exptab1['nsources'][0])
        else:
            print(i+1,base)

        # Save
        if (i % 500 == 0 and i > 0) or (i==len(expdirs)-1):
            print('Writing to',outfile)
            exptab = vstack(expdata)
            for c in ['file','wtfile','maskfile','base']:
                exptab[c] = np.array([f.strip() for f in exptab[c]])
            chiptab = vstack(chipdata)
            for c in ['expdir','filename','measfile']:
                chiptab[c] = np.array([f.strip() for f in chiptab[c]])
            ohdu = fits.HDUList()
            ohdu.append(fits.table_to_hdu(exptab))
            ohdu.append(fits.table_to_hdu(chiptab))
            ohdu.writeto(outfile,overwrite=True)
            ohdu.close()

    # Write out the file
    #print('Writing summary file to ',outfile)
    #exptab = Table(exptab)
    #exptab.write(outfile,overwrite=True)

    print('dt = {:.1f} sec.'.format(time.time()-t0))

    import pdb; pdb.set_trace()

def calibratedirs(expdirs,outfile,version='v4'):
    """ Make a calibrate summary file for a subset of exposures """

    dldir,mssdir,localdir = utils.rootdirs()
    basedir = os.path.join(dldir,'instcal/',version)

    t0 = time.time()

    print('Gathering summary data for',len(expdirs),'exposures')

    # Load the exposure and chip meta files
    calbasedir = '/home1/09970/dnidever/scratch1/nsc/instcal/v4'
    expdata = []
    chipdata = []
    for i in range(len(expdirs)):
        dir1 = expdirs[i]
        arr = dir1.split('/')
        base = arr[-1]
        night = arr[-2]
        year = arr[-3]
        instrument = arr[-4]

        measfile = os.path.join(calbasedir,instrument,year,night,base,base+'_meas.fits')
        metafile = os.path.join(calbasedir,instrument,year,night,base,base+'_meta.fits')
        if os.path.exists(metafile)==False:
            print(i+1,base,'no meta file')
            continue
        
        try:
            hdu = fits.open(metafile)
            exptab1 = Table(hdu[1].data)
            # fix names
            for c in ['file','maskfile','wtfile']:
                ffile = exptab1[c][0]
                if ffile.find('INFO')>-1:
                    exptab1[c] = ffile[ffile.find(']')+1:].strip()
            exptab1['measfile'] = measfile
            chiptab1 = []
            for j in range(len(hdu)-2):
                chiptab1.append(Table(hdu[j+2].data))
            chiptab1 = vstack(chiptab1)
            hdu.close()
            chiptab1['file'] = exptab1['file'][0]
            chiptab1['wtfile'] = exptab1['wtfile'][0]
            chiptab1['maskfile'] = exptab1['maskfile'][0]
            chiptab1['measfile'] = measfile
            chiptab1['base'] = base
            expdata.append(exptab1)
            chipdata.append(chiptab1)
            print(i+1,base,len(chiptab1),exptab1['nsources'][0])
        except:
            traceback.print_exc()

    # Save
    print('Writing to',outfile)
    exptab = vstack(expdata)
    for c in ['file','wtfile','maskfile','base']:
        exptab[c] = np.array([f.strip() for f in exptab[c]])
    chiptab = vstack(chipdata)
    for c in ['expdir','filename','measfile']:
        chiptab[c] = np.array([f.strip() for f in chiptab[c]])
    ohdu = fits.HDUList()
    ohdu.append(fits.table_to_hdu(exptab))
    ohdu.append(fits.table_to_hdu(chiptab))
    ohdu.writeto(outfile,overwrite=True)
    ohdu.close()

def calibratechunkscombine():
    """ Combine the chunks of calibration summary information """
    basedir = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/summary/calibrate'
    files = glob(basedir+'/calibrate_summary*.fits')
    files.sort()
    num = [int(os.path.basename(f)[17:-5]) for f in files]
    si = np.argsort(num)
    files = np.array(files)[si]
    print('found',len(files),'calibration chunk summary files')
    expdata = []
    chipdata = []
    for i in range(len(files)):
        print(i+1,files[i])
        hdu = fits.open(files[i])
        # exposure table
        exptab = Table(hdu[1].data)
        # Trim down the sizes
        for c in ['file','wtfile','maskfile','dateobs','instrument','base','filter','wcscal','measfile']:
            exptab[c] = np.array([f for f in exptab[c]])
        for c in ['exptime','airmass','fwhm','rarms','decrms','ebv','zpterm','zptermerr',
                  'zptermsig','zpspatialvar_rms','zpspatialvar_range','depth95','depth10sig']:
            exptab[c] = exptab[c].astype(np.float32)
        for c in ['zptype','zpspatialvar_nccd']:
            exptab[c] = exptab[c].astype(np.int16)
        expdata.append(exptab)
        # chip table
        chtab = Table(hdu[2].data)
        chtab['instrument'] = np.array([f.astype(str) for f in chtab['instrument']])
        del chtab[['expdir','file','wtfile','maskfile','filename','nmeas']]
        for c in ['rarms','rastderr','decrms','decstderr','zpterm','zptermerr','depth95','depth10sig']:
            chtab[c] = chtab[c].astype(np.float32)
        for c in ['ccdnum','zptype']:
            chtab[c] = chtab[c].astype(np.int16)
        #chtab.write(files[i].replace('.fits','_chip.csv'),format='csv')
        chipdata.append(chtab)
        hdu.close()

    # Write to final output file
    exptab = vstack(expdata)
    del expdata
    chiptab = vstack(chipdata)
    del chipdata

    # Add healpix to each table
    nside = 128
    exptab['pix128'] = hp.ang2pix(nside,exptab['ra'],exptab['dec'],lonlat=True)
    chiptab['pix128'] = hp.ang2pix(nside,chiptab['cenra'],chiptab['cendec'],lonlat=True)

    #np.save('/home1/09970/dnidever/scratch1/nsc/instcal/v4/lists/nsc_calibrate_summary_chips.fits',chiptab)
    chiptab.write('/home1/09970/dnidever/scratch1/nsc/instcal/v4/lists/nsc_calibrate_summary_chips.fits')
    exptab.write('/home1/09970/dnidever/scratch1/nsc/instcal/v4/lists/nsc_calibrate_summary_exp.fits')

    #outfile = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/lists/nsc_calibrate_summary.fits'
    #print('Writing summary results to',outfile)
    #ohdu = fits.HDUList()
    #ohdu.append(fits.table_to_hdu(exptab))
    #ohdu.append(fits.table_to_hdu(chiptab))
    #ohdu.writeto(outfile,overwrite=True)
    #ohdu.close()

    import pdb; pdb.set_trace()

def combine():
    """ Make the nsc_combine_summary.fits summary file """
    pass


def combinehealpix():
    """ Apply QA cuts and make healpix file """
    # from nsc_instcal_combine_qacuts.pro


    # Combine all of the data
    dldir,mssdir,localdir = utils.rootdirs()
    basedir = os.path.join(dldir,'instcal/',version)
    #host = first_el(strsplit(longhost,'.',/extract))
    #basedir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
    if os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/')==False:
        os.makedirs(localdir+'dnidever/nsc/instcal/'+version+'/')
    plotsdir = basedir+'plots/'
    if os.path.exists(plotsdir)==False:
        os.makedirs(plotsdir)
    nside = 128
    time = time.time()

    # Restore the calibration summary file                                                                                       
    temp = Table.read(basedir+'lists/nsc_calibrate_summary.fits.gz',1)
    #schema = temp[0]
    #struct_assign,{dum:''},schema
    #schema = create_struct(schema,'chipindx',-1,'NGOODCHIPWCS',0,'wcscal','')
    #str = replicate(schema,n_elements(temp))
    #struct_assign,temp,tab,/nozero
    tab = temp.copy()
    tab['chipindex'] = -1
    tab['ngoodchipwcs'] = 0
    tab['wcscal'] = 50*' '
    tab['expdir'] = np.array([str(e) for e in tab['expdir']])
    tab['instrument'] = np.array([str(e) for e in tab['instrument']])
    tab['metafile'] = np.array([str(e) for e in tab['metafile']])
    tab['file'] = np.array([str(f) for f in tab['file']])
    tab['base'] = np.array([str(f) for f in tab['base']])
    tab['filter'] = np.array([str(f) for f in tab['filter']])
    # Add WCSCAL and TELSTAT information
    coords = Table.read(basedir+'lists/allcoords.fits.gz',1)
    coords['file'] = np.array([str(f) for f in coords['file']])
    coords['wcscal'] = np.array([str(f) for f in coords['wcscal']])
    coords['telstat'] = np.array([str(f) for f in coords['telstat']])
    fluxfile = tab['file']
    g, = np.where(fluxfile[:4] == '/net')
    if ng > 0:
        fluxfile[g] = fluxfile[g][4:]
    _,ind1,ind2 = np.intersect1d(fluxfile,coords['file'],return_indices=True)
    ## v3, 490617 out of 490623 matches, only 6 did not match                                                                    
    tab['wcscal'][ind1] = coords['wcscal'][ind2]    # Failed (37712), Poor (0), Successful (452905)                                    
    # Only want exposures with successful SE processing                                                                          
    gd, = np.where(tab['success']==True)
    print(len(gd),' successful exposures')
    tab = tab[gd]
    si = np.argsort(tab['expdir'])
    tab = tab[si]
    chtab = Table.read(basedir+'lists/nsc_calibrate_summary.fits.gz',2)
    chtab['expdir'] = strtrim(chtab['expdir'],2)
    chtab['instrument'] = strtrim(chtab['instrument'],2)
    nchtab = len(chtab)
    # Get indices for CHSTR
    chindex = dln.create_index(chtab['expdir'])
    #siexp = sort(chtab[]expdir)
    #chstr = chstr[siexp]
    #expdir = chtab[]expdir
    #brklo = where(expdir ne shift(expdir,1),nbrk)
    #brkhi = [brklo[1:nbrk-1]-1,n_elements(expdir)-1]
    #nchexp = brkhi-brklo+1
    #if nstr ne n_elements(brklo) then stop,'number of exposures in STR and CHSTR do not match'
    tab['chipindx'] = brklo
    tab['nchips'] = nchexp
    # Getting number of good chip WCS for each exposures
    for i in range(len(tab)):
        tab['ngoodchipwcs'][i] = np.sum(chtab['ngaiamatch'][brklo[i]:brkhi[i]] > 0)
    # Fixing absolute paths of flux filename
    filename = tab['file']
    filename = [f.replace('/net/mss1/','/') for f in filename]
    filename = [f.replace('/mss1/','/') for f in filename]
    #g1, = np.where(stregex(filename,'/net/mss1/',/boolean) == True)
    #if len(g1) > 0:
    #    filename[g1] = strmid(filename[g1],10)
    #g2, = np.where(stregex(filename,'/mss1/',/boolean) == True)
    #if len(g2) > 0:
    #    filename[g2] = strmid(filename[g2],6)
    # Fixing very negative RAs
    print('FIXING NEGATIVE RAs in TAB and CHTAB')
    bdra, = np.where(chtab['cenra'] < 0)
    _,uibd = np.unique(chtab['expdir'][bdra],return_index=True)
    #MATCH,tab['expdir'],chstr[bdra[uibd]].expdir,ind1,ind2,/sort,count=nmatch
    _,ind1,ind2 = np.intersect1d(tab['expdir'],chtab['expdir'][bdra[uibd]],return_indices=True)
    nmatch = len(ind1)
    for i in range(nmatch):
        _,ind3,ind4 = np.intersect1d(chtab['expdir'][bdra],tab['expdir'][ind1[i]],return_indices=True)
        #MATCH,chstr[bdra].expdir,str[ind1[i]].expdir,ind3,ind4,/sort
        # Fix TAB RA
        chra = chtab['cenra'][bdra[ind3]]
        bd1, = np.where(chra < -180)
        if len(bd1) > 0:
            chra[bd1] += 360
        cenra = np.mean([np.min(chra),np.max(chra)])
        if cenra < 0:
            cenra += 360
        tab['ra'][ind1[i]] = cenra
        # Fix CHSTR CENRA
        bd2, = np.where(chra < 0)
        if len(bd2) > 0:
            chra[bd2] += 360
        chtab['CENRA'][bdra[ind3]] = chra
        # Fix CHSTR VRA
        vra = chtab['vra'][bdra[ind3]]
        bd3, = np.where(vra < 0)
        if len(bd3) > 0:
            vra[bd3] += 360
        chstr[bdra[ind3]].vra = vra

    # Zero-point structure                                                                                                       
    dtyp = [('instrument',str,3),('filter',str,2),('amcoef',float,2),('thresh',float)]
    zptab = np.zeros(10,dtype=np.dtype(dtyp))
    zptab['thresh'] = 0.5
    zptab['instrument'][:6] = 'c4d'
    zptab['filter'][:6] = ['u','g','r','i','z','Y','VR']
    zptab['amcoef'][0] = [-1.60273, -0.375253]   # c4d-u                                                                            
    zptab['amcoef'][1] = [0.277124, -0.198037]   # c4d-g                                                                            
    zptab['amcoef'][2] = [0.516382, -0.115443]   # c4d-r  changed a bit, fine                                                       
    zptab['amcoef'][3] = [0.380338, -0.067439]   # c4d-i                                                                            
    zptab['amcoef'][4] = [0.074517, -0.067031]   # c4d-z                                                                            
    zptab['amcoef'][5] = [-1.07800, -0.060014]   # c4d-Y                                                                            
    zptab['amcoef'][6] = [1.111859, -0.083630]   # c4d-VR                                                                           
    # Mosiac3 z-band                                                                                                             
    zptab['instrument'][7] = 'k4m'
    zptab['filter'][7] = 'z'
    zptab['amcoef'][7] = [2.232800, -0.73573]   # k4m-z                                                                             
    # Bok 90Prime, g and r                                                                                                       
    zptab['instrument'][8] = 'ksb'
    zptab['filter'][8] = 'g'
    zptab['amcoef'][8] = [1.055275, -0.30629]   # ksb-g                                                                             
    zptab['instrument'][9] = 'ksb'
    zptab['filter'][9] = 'r'
    zptab['amcoef'][9] = [0.836968, -0.19646]   # ksb-r                                                                             
    nzptab = len(zptab)

    # APPLY QA CUTS IN ZEROPOINT AND SEEING
    if nocuts==False:
        print('APPLYING QA CUTS')
        fwhmthresh = 2.0  # arcsec, v2
        #filters = ['u','g','r','i','z','Y','VR']
        #nfilters = n_elements(filters)
        #zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
        #zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        badzpmask = bytarr(n_elements(str)) + 1

        for i in range(nzptab):
            ind, = np.where((tab['instrument'] == zptab['instrument'][i]) & (tab['filter'] == zptab['filter'][i]) & (tab['success']==True))
            print(zptab['instrument'][i],'-',zptab['filter'][i],' ',len(ind),' exposures')
            if len(ind) > 0:
                tab1 = tab[ind]
                ## Fix Infinity/NAN values
                zpterm = tab1['zpterm']
                bdzp, = np.where(np.isfinite(zpterm) == 0)  # fix Infinity/NAN
                if len(bdzp)>0:
                    zpterm[bdzp] = 999999.9
                ## Correct "DES" zeropoints,  DES exposures are in electrons and
                ## CP are in ADU, so there's an offset of 2.5*log(gain)=2.5*log(4.41)=1.611
                gdes, = np.where(tab1['plver'][:3]=='DES')
                if ngdes > 0:
                    print('Offsetting ',len(gdes),' DES exposure zero-points')
                    zpterm[gdes] -= 1.611

                ## CORRECT K4M/KSB for exptime-dependence in the zero-points
                ##   this is because the image units are counts/sec.
                if zptab['instrument'][i] == 'k4m' or zptab['instrument'][i] == 'ksb':
                    print('REMOVING EXPTIME-DEPENDENCE IN K4M/KSB ZEROPOINTS!!!')
                    zpterm += 2.5*alog10(tab1['exptime'])
                am = tab1['airmass']
                mjd = tab1['mjd']
                bdam, = np.where(am < 0.9)
                if len(bdam) > 0:
                    am[bdam] = np.median(am)
                coo = SkyCoord(tab1['ra'],tab1['dec'],unit='degree',frame='icrs')
                glon = coo.galactic.l.degree
                glat = coo.galactic.b.egree

                # Measure airmass dependence
                gg0, = np.where((n.abs(zpterm) < 50) & (am < 2.0))
                coef0 = robust_poly_fitq(am[gg0],zpterm[gg0],1)
                zpf = np.polyval(coef0,am)
                sig0 = dln.mad(zpterm[gg0]-zpf[gg0])
                gg, = np.where(np.abs(zpterm-zpf) < np.maximum(3.5*sig0,0.2))
                coef = robust_poly_fitq(am[gg],zpterm[gg],1)
                print(zptab['instrument'][i]+'-'+zptab['filter'][i],' ',coef)
                # Trim out bad exposures to determine the correlations and make figures
                gg, = np.where((np.abs(zpterm-zpf) < np.maximum(3.5*sig0,0.2)) & (tab1['airmass'] < 2.0) & (tab1['fwhm'] < 2.0) & (tab1['rarms'] < 0.15) &
                               (tab1['decrms'] < 0.15) & (tab1['success']==True) & (tab1['wcscal']=='Successful') & (tab1['zptermerr'] < 0.05) &
                               (tab1['zptermsig'] < 0.08) &
                               ((tab1['instrument'] != 'c4d') | (tab1['zpspatialvar_nccd']<=5) | ((tab1['instrument']=='c4d') & (tab1['zpspatialvar_nccd']>5) & (tab1['zpspatialvar_rms']<0.1))) &
                               (np.abs(glat) > 10) & (tab1['nrefmatch'] > 100) & (tab1['exptime'] >= 30))
                ## I removed WCSCAL check because there are ~38k exposures with
                ## WCSCAL=Failed but my DECRMS and RARMS is small.
                ## and tab1.wcscal eq 'Successful'
                print(ngg)

                # Zpterm with airmass dependence removed
                relzpterm = zpterm + 25   # 25 to get "absolute" zpterm
                relzpterm -= zptab['zmcoef'][i][1]*(am-1)

                # Fit temporal variation in zpterm
                mjd0 = 56200
                xx = tab1['mjd'][gg]-mjd0
                yy = relzpterm[gg]
                invvar = 1.0/tab1['zptermerr'][gg]**2
                nord = 3
                bkspace = 200 #20
                sset1 = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=yfit1)
                sig1 = mad(yy-yfit1)
                gd = where(yy-yfit1 > -3*sig1,ngd)
                # refit
                sset = bspline_iterfit(xx[gd],yy[gd],invvar=invvar[gd],nord=nord,bkspace=bkspace)
                yfit = bspline_valu(xx,sset)
                allzpfit = bspline_valu(tab1['mjd']-mjd0,sset)
            
                # Remove temporal variations to get residual values
                relzpterm -= allzpfit


                # Find the GOOD exposures
                #------------------------
                # We are using ADDITIVE zpterm
                #  calmag = instmag + zpterm
                # if there are clouds then instmag is larger/fainter
                #  and zpterm is smaller (more negative)
                #bdind = where(str[ind].zpterm-medzp lt -zpthresh[i],nbdind)
                goodmask = ((relzpterm >= -zptab['thresh'][i]) & (relzpterm <= zptab['thresh'][i]))
                gdind, = np.where(goodmask)
                bdind, = np.where(~goodmask)
                print('  ',len(bdind),'exposures with ZPTERM below the threshold')
                if len(gdind) > 0:
                    badzpmask[ind[gdind]] = 0

        # Get bad DECaLS and SMASH exposures
        badexp = np.zeros(len(tab),bool)
        smashexpnum = dln.readlines('/home/dnidever/projects/noaosourcecatalog/obslog/'+version+'/smash_badexposures.txt')
        smashexpnum = [int(e) for e in smashexpnum]
        _,ind1,ind2 = np.intersect1d(tab['expnum'],smashexpnum,return_indices=True)
        if len(ind1) > 0:
            badexp[ind1] = True
            badexp[ind1] = (badexp[ind1] & (tab['instrument'][ind1] == 'c4d'))   # make sure they are DECam exposures
        decalsexpnum = dln.readlines('/home/dnidever/projects/noaosourcecatalog/obslog/'+version+'/decals_bad_expid.txt')
        decalsexpnum = [int(e) for e in decalsexpnum]
        _,ind1,ind2 = np.intersect1d(tab['expnum'],decalsexpnum,return_indices=True)
        if len(ind1) > 0:
            badexp[ind1] = True
            badexp[ind1] = (badexp[ind1] & (tab['instrument'][ind1].instrument == 'c4d'))   # make sure they are DECam exposures
        mzlsexpnum = dln.readlines('/home/dnidever/projects/noaosourcecatalog/obslog/'+version+'/mzls_bad_expid.txt')
        mzlsexpnum = [int(e) for e in mzlsexpnum]
        _,ind1,ind2 = np.intersect1d(tab['expnum'],mzlsexpnum,return_indices=True)
        if len(ind1) > 0:
            badexp[ind1] = True
            badexp[ind1] = (badexp[ind1] & (tab['instrument'][ind1] == 'k4m'))   # make sure they are Mosaic3 exposures

        ## Zero-point spatial variability threshold
        ##  varies with galactic latitude
        ##  |b|>10   0.15
        ##  |b|<=10  0.55
        coo = SkyCoord(tab['ra'],tab['dec'],unit='degree',frame='icrs')
        glon = coo.galactic.l.degree
        glat = coo.galactic.b.degree
        zpspvarthresh = (np.abs(glat) > 10)*0.15 + (np.abs(glat) <= 10)*0.55

        # Final QA cuts
        #  Many of the short u-band exposures have weird ZPTERMs, not sure why
        #  There are a few exposures with BAD WCS, RA>360!
        bdexp, = np.where((tab['success'] == False) |                          # SE failure
                          (tab['fwhm'] > fwhmthresh) |                         # bad seeing
                          (tab['ra'] > 360) |                                  # bad WCS/coords
                          (tab['rarms'] > 0.15) | (tab['decrms'] > 0.15) |     # bad WCS
                          (badzpmask == True) |                                # bad ZPTERM
                          (tab['zptermerr'] > 0.05) |                          # bad ZPTERMERR
                          (tab['nrefmatch'] < 5) |                             # few phot ref match
                          (badexp == 1) |                                      # bad SMASH/LS exposure
                          ((tab['instrument'] == 'c4d') & (tab['zpspatialvar_nccd'] > 5) & (tab['zpspatialvar_rms'] > zpspvarthresh)))  # bad spatial zpterm
                          ##tab[]wcscal ne 'Successful' or $                    ##CP WCS failure   TOO MANY FAILED
                          #tab[]ngoodchipwcs lt tab[]nchips or $                # not all chips astrom calibrated
        # rarms/decrms, nrefmatch
        print('QA cuts remove ',len(bdexp),' exposures')

        # Remove
        torem = np.zeros(len(chtab),bool)
        for i in range(len(bdexp)):
            torem[tab['chipindx'][bdexp[i]]:tab['chipindx'][bdexp[i]]+tab['nchips'][bdexp[i]]] = True
        bdchtab, = np.where(torem == True)
        chtab = np.delete(chtab,bdchtab)
        tab = np.delete(tab,bdexp)
        # Get new CHIPINDEX values
        #   make two arrays of old and new indices to transfer
        #   the new index values into an array with the size of
        #   the old CHSTR
        trimoldindex = lindgen(nchstr)                    # index into original array, but "bad" ones removed/trimed
        remove,bdchstr,trimoldindex
        trimnewindex = lindgen(n_elements(trimoldindex))  # new index of trimmed array
        newindex = lonarr(nchstr)-1
        newindex[trimoldindex] = trimnewindex             # new index in original array
        newchipindex = newindex[tab['chipindx']]
        tab['chipindx'] = newchipindex
        nstr = len(tab)
    else:
        print('SKIPPING QA CUTS')


    # CREATE LIST OF HEALPIX AND OVERLAPPING EXPOSURES
    # Which healpix pixels have data
    listfile = basedir+'lists/nsc_instcal_combine_healpix_list.fits'
    if os.path.exists(listfile)==False or redo:
        print('Finding the Healpix pixels with data')
        radius = 1.1
        dtyp = [('file',str,200),('base',str,50),('pix',int)]
        healtab = np.zeros(1000000,dtype=np.dtype(dtyp))
        nhealtab = len(healtab)
        cnt = 0
        for i in range(ntab):
            if i % 1e3 == 0:
                print(i)
            vec = hp.ang2vec(nside,tab['ra'][i],tab['dec'][i],lonlat=True)
            listpix = hp.query_disc(nside,vec,radius,inclusive=True,nest=False)
            nlistpix = len(listpix)
            #theta = (90-tab['dec'][i])/radeg
            #phi = tab['ra'][i]/radeg
            #ANG2VEC,theta,phi,vec
            #QUERY_DISC,nside,vec,radius,listpix,nlistpix,/deg,/inclusive

            # Use the chip corners to figure out which ones actually overlap
            chtab1 = chtab[tab['chipindx'][i]:tab['chipindx'][i]+tab['nchips'][i]]
            #  rotate to tangent plane so it can handle RA=0/360 and poles properly
            vlon,vlat = coords.rotsphcen(chtab1['vra'],chtab1['vdec'],tab['ra'][i],tab['dec'][i],gnomic=True)
            #  loop over healpix
            overlap = np.zeros(len(listpix),bool)
            for j in range(len(listpix)):
                vec,vertex = hp.pix2vec(nside,listpix[j],nest=False)
                hra,hdec = hp.vec2ang(vertex,lonlat=True)
                #PIX2VEC_RING,nside,listpix[j],vec,vertex
                #vertex = transpose(reform(vertex))  # [1,3,4] -> [4,3]                                                                 
                #VEC2ANG,vertex,hdec,hra,/astro
                hlon,hlat = coords.rotsphcen(hra,hdec,tab['ra'][i],tab['dec'][i],gnomic=True)
                #  loop over chips
                for k in range(tab['nchips'][i]):
                    overlap[j] >= coords.doPolygonsOverlap(hlon,hlat,vlon[:,k],vlat[:,k])
            # Only keep the healpix with real overlaps
            gdlistpix, = np.where(overlap==True)
            if len(gdlistpix) > 0:
                listpix = listpix[gdlistpix]
                nlistpix = len(gdlistpix)
            else:
                listpix = []
                nlistpix = 0

            # Add new elements to array
            if cnt+len(listpix) > nhealtab:
                old = healtab.copy()
                healtab = np.zeros(nhealtab+10000,dtype=np.dtype(dtyp))
                healtab[0:nhealtab] = old
                nhealtab += 1e4
                del old

            # Add to the structure
            healtab['file'][cnt:cnt+nlistpix] = tab['expdir'][i]+'/'+tab['base'][i]+'_cat.fits'
            healtab['base'][cnt:cnt+nlistpix] = tab['base'][i].base
            healtab['pix'][cnt:cnt+nlistpix] = listpix
            cnt += nlistpix

        # Trim extra elements
        healtab = healtab[:cnt]
        nhealtab = len(healtab)

        # Get uniq pixels
        _,ui = np.unique(healtab['pix'],return_index=True)
        upix = healtab['pix'][ui]
        nupix = n_elements(upix)
        print(nupix,'Healpix pixels have overlapping data')

        # Get start/stop indices for each pixel
        hindex = dln.create_index(healtab['pix'])
        #idx = sort(healtab.pix)
        #healtab = healtab[idx]
        #q = healtab.pix
        #lo = where(q ne shift(q,1),nlo)
        ##hi = where(q ne shift(q,-1))
        #hi = [lo[1:nlo-1]-1,nhealtab-1]
        #nexp = hi-lo+1
        #index = replicate({pix:0L,lo:0L,hi:0L,nexp:0L},nupix)
        #index.pix = upix
        #index.lo = lo
        #index.hi = hi
        #index.nexp = nexp
        npix = len(index['value'])
        
        # Replace /net/dl1/ with /dl1/ so it will work on all machines
        healtab['file'] = [f.replace('/net/dl1/','/dl1/') for f in healtab['file']]
        
        # Replace /net/dl1/ with /dl1/ so it will work on all machines
        healtab['file'] = [f.replace('/net/dl1/','/dl1/') for f in healtab['file']]
        
        # Write the full list plus an index
        print('Writing list to ',listfile)
        hdu = fits.HDUList()
        hdu.append(fits.table_to_hdu(healtab))
        hdu[1].header['nside'] = nside
        hdu.append(fits.table_to_hdu(index))
        hdu.writeto(listfile,overwrite=True)
        hdu.close()
        if os.path.exists(listfile+'.gz'):
            os.remove(listfile+'.gz')
        out = subprocess.call(['gzip',listfile],noshell=True)
        #MWRFITS,healtab,listfile,/create
        ## Add NSIDE to header
        #hd0 = headfits(listfile,exten=0)
        #sxaddpar,hd0,'nside',nside
        #modfits,listfile,0,hd0,exten_no=0
        #MWRFITS,index,listfile,/silent
        #if file_test(listfile+'.gz') eq 1 then file_delete,listfile+'.gz',/allow
        #spawn,['gzip',listfile],/noshell
            
    else:
        print(listfile,' EXISTS and redo NOT set')


    print('dt = {:.1f} sec.'.format(time.time()-t0))


