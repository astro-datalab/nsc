import os
import numpy as np
from glob import glob
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from dlnpyutils import utils as dln,coords
from scipy.stats import binned_statistic
import healpy as hp
import time
import traceback
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm as LN
from . import utils,db

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

def combine_qacuts(version='v4'):
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
    t0 = time.time()

    #basedir = '/Users/nidever/datalab/nsc/v4/'
    
    # Load the full exposure list
    decamlist1 = Table.read(os.path.join(basedir,'lists','decam_instcal_list.fits.gz'))
    for c in decamlist1.colnames: decamlist1[c].name=c.lower()
    decamlist1['base'] = [os.path.basename(f).replace('.fits.fz','').strip() for f in decamlist1['fluxfile']]
    decamlist1['plver'] = decamlist1['plver'].astype(str)
    decamlist2 = Table.read(os.path.join(basedir,'lists','decam_instcal_list_v3.fits.gz'))
    for c in decamlist2.colnames: decamlist2[c].name=c.lower()
    decamlist2['base'] = [os.path.basename(f).replace('.fits.fz','').strip() for f in decamlist2['fluxfile']]
    decamlist2['plver'] = decamlist2['plver'].astype(str)
    decamlist3 = Table.read(os.path.join(basedir,'lists','decam_instcal_list_20240727.fits.gz'))
    decamlist3['base'] = np.array([str(b).strip() for b in decamlist3['base']])
    decamlist3['plver'] = decamlist3['plver'].astype(str)


    # Restore the calibration summary file
    exptab = Table.read(os.path.join(basedir,'lists','nsc_calibrate_summary_exp.fits'))
    #schema = temp[0]
    #struct_assign,{dum:''},schema
    #schema = create_struct(schema,'chipindx',-1,'NGOODCHIPWCS',0,'wcscal','')
    #str = replicate(schema,n_elements(temp))
    #struct_assign,temp,tab,/nozero
    #tab = temp.copy()
    exptab['chipindex'] = -1
    exptab['ngoodchipwcs'] = 0
    #exptab['wcscal'] = 50*' '
    #exptab['expdir'] = np.array([str(e) for e in exptab['expdir']])
    exptab['instrument'] = np.array([str(e) for e in exptab['instrument']])
    #exptab['metafile'] = np.array([str(e) for e in exptab['metafile']])
    exptab['file'] = np.array([str(f) for f in exptab['file']])
    exptab['base'] = np.array([str(f) for f in exptab['base']])
    exptab['filter'] = np.array([str(f) for f in exptab['filter']])
    exptab['plver'] = 50*' '

    # Getting PLVER from the decam_instcal_list catalogs
    _,ind1,ind2 = np.intersect1d(exptab['base'],decamlist1['base'],return_indices=True)
    exptab['plver'][ind1] = decamlist1['plver'][ind2]
    _,ind1,ind2 = np.intersect1d(exptab['base'],decamlist2['base'],return_indices=True)
    exptab['plver'][ind1] = decamlist2['plver'][ind2]
    _,ind1,ind2 = np.intersect1d(exptab['base'],decamlist3['base'],return_indices=True)
    exptab['plver'][ind1] = decamlist3['plver'][ind2]
    # about 7000 rows are still missing plver values
    # tu
    base2 = np.array([b[:2] for b in exptab['base']])
    bb, = np.where((base2=='tu') & (exptab['plver']==50*' '))
    # 1157
    exptab['plver'][bb] = 'V3.1.0'
    bb, = np.where(exptab['plver']==50*' ')
    # v1
    basever = np.array([b.split('_')[-1] for b in exptab['base']])
    bb, = np.where((basever=='v1') & (exptab['plver']==50*' '))
    # V4.1 to V5.5
    # MT1, V4.8
    bb, = np.where((basever=='MT1') & (exptab['plver']==50*' '))    
    exptab['plver'][bb] = 'V4.8'

    # 1891

    # Deal with duplicate versions for exposures
    _,ui = np.unique(exptab['expnum'],return_index=True)
    expindex = dln.create_index(exptab['expnum'])

    #(Pdb) np.sum(expindex['num']>1)
    #46841
    dupind, = np.where(expindex['num']>1)
    keepind = []
    missingplver = 0
    for i in range(len(expindex['value'])):
        ind = expindex['index'][expindex['lo'][i]:expindex['hi'][i]+1]
        nind = len(ind)
        if nind==1:
            keepind.append(ind[0])
        else:
            plver = exptab['plver'][ind]
            problem = np.sum(plver==50*' ')
            if problem > 0:
                # v# and v#
                basever = np.array([b.split('_')[-1] for b in exptab['base'][ind]])
                basever1 = np.array([b[:1] for b in basever])
                if np.sum(basever1=='v')==len(ind):
                    si = np.argsort(basever)[::-1]
                    keepind.append(ind[si[0]])
                    continue
                # if one is lsXX, use it
                basever2 = np.array([b[:2] for b in basever])
                if np.sum(basever2=='ls')>0:
                    bestind, = np.where(basever2=='ls')
                    if len(bestind)>1:
                        print('multiple ls versions')
                        import pdb; pdb.set_trace()
                    keepind.append(ind[bestind[0]])
                    continue
                # Pick the version with the lowest RARMS
                si = np.argsort(exptab['rarms'][ind])
                keepind.append(ind[si[0]])
                continue

                #print('missing plver')
                #print(missingplver+1,exptab['base'][ind])
                #missingplver += 1
                #import pdb; pdb.set_trace()
            else:
                si = np.argsort(plver)[::-1]
            keepind.append(ind[si[0]])


    # 572204 unique exposures
    exptab = exptab[keepind]

    # Add galactic coordinates
    coo = SkyCoord(exptab['ra'],exptab['dec'],unit='degree',frame='icrs')
    exptab['glon'] = coo.galactic.l.degree
    exptab['glat'] = coo.galactic.b.degree
    
    # Add WCSCAL and TELSTAT information
    #coords = Table.read(basedir+'lists/allcoords.fits.gz',1)
    #coords['file'] = np.array([str(f) for f in coords['file']])
    #coords['wcscal'] = np.array([str(f) for f in coords['wcscal']])
    #coords['telstat'] = np.array([str(f) for f in coords['telstat']])
    #fluxfile = tab['file']
    #g, = np.where(fluxfile[:4] == '/net')
    #if ng > 0:
    #    fluxfile[g] = fluxfile[g][4:]
    #_,ind1,ind2 = np.intersect1d(fluxfile,coords['file'],return_indices=True)
    ## v3, 490617 out of 490623 matches, only 6 did not match
    #tab['wcscal'][ind1] = coords['wcscal'][ind2]    # Failed (37712), Poor (0), Successful (452905)

    ## Only want exposures with successful SE processing
    #gd, = np.where(exptab['success']==True)
    #print(len(gd),' successful exposures')
    #exptab = exptab[gd]
    #si = np.argsort(exptab['expdir'])
    #exptab = exptab[si]
    #chtab = Table.read(basedir+'lists/nsc_calibrate_summary.fits.gz',2)
    #chtab['expdir'] = strtrim(chtab['expdir'],2)
    #chtab['instrument'] = strtrim(chtab['instrument'],2)
    #nchtab = len(chtab)
    ## Get indices for CHSTR
    #chindex = dln.create_index(chtab['expdir'])
    ##siexp = sort(chtab[]expdir)
    ##chstr = chstr[siexp]
    ##expdir = chtab[]expdir
    ##brklo = where(expdir ne shift(expdir,1),nbrk)
    ##brkhi = [brklo[1:nbrk-1]-1,n_elements(expdir)-1]
    ##nchexp = brkhi-brklo+1
    ##if nstr ne n_elements(brklo) then stop,'number of exposures in STR and CHSTR do not match'
    #exptab['chipindx'] = brklo
    #exptab['nchips'] = nchexp
    ## Getting number of good chip WCS for each exposures
    #for i in range(len(exptab)):
    #    exptab['ngoodchipwcs'][i] = np.sum(chtab['ngaiamatch'][brklo[i]:brkhi[i]] > 0)
    ## Fixing absolute paths of flux filename
    #filename = exptab['file']
    #filename = [f.replace('/net/mss1/','/') for f in filename]
    #filename = [f.replace('/mss1/','/') for f in filename]
    ##g1, = np.where(stregex(filename,'/net/mss1/',/boolean) == True)
    ##if len(g1) > 0:
    ##    filename[g1] = strmid(filename[g1],10)
    ##g2, = np.where(stregex(filename,'/mss1/',/boolean) == True)
    ##if len(g2) > 0:
    ##    filename[g2] = strmid(filename[g2],6)
    ## Fixing very negative RAs
    #print('FIXING NEGATIVE RAs in EXPTAB and CHTAB')
    #bdra, = np.where(chtab['cenra'] < 0)
    #_,uibd = np.unique(chtab['expdir'][bdra],return_index=True)
    ##MATCH,exptab['expdir'],chstr[bdra[uibd]].expdir,ind1,ind2,/sort,count=nmatch
    #_,ind1,ind2 = np.intersect1d(exptab['expdir'],chtab['expdir'][bdra[uibd]],return_indices=True)
    #nmatch = len(ind1)
    #for i in range(nmatch):
    #    _,ind3,ind4 = np.intersect1d(chtab['expdir'][bdra],exptab['expdir'][ind1[i]],return_indices=True)
    #    #MATCH,chstr[bdra].expdir,str[ind1[i]].expdir,ind3,ind4,/sort
    #    # Fix TAB RA
    #    chra = chtab['cenra'][bdra[ind3]]
    #    bd1, = np.where(chra < -180)
    #    if len(bd1) > 0:
    #        chra[bd1] += 360
    #    cenra = np.mean([np.min(chra),np.max(chra)])
    #    if cenra < 0:
    #        cenra += 360
    #    exptab['ra'][ind1[i]] = cenra
    #    # Fix CHSTR CENRA
    #    bd2, = np.where(chra < 0)
    #    if len(bd2) > 0:
    #        chra[bd2] += 360
    #    chtab['CENRA'][bdra[ind3]] = chra
    #    # Fix CHSTR VRA
    #    vra = chtab['vra'][bdra[ind3]]
    #    bd3, = np.where(vra < 0)
    #    if len(bd3) > 0:
    #        vra[bd3] += 360
    #    chstr[bdra[ind3]].vra = vra

    # Zero-point structure
    dtyp = [('instrument',str,3),('filter',str,2),('amcoef',float,2),('thresh',float)]
    zptab = np.zeros(10,dtype=np.dtype(dtyp))
    zptab['thresh'] = 0.5
    zptab['instrument'][:7] = 'c4d'
    zptab['filter'][:7] = ['u','g','r','i','z','Y','VR']
    #zptab['amcoef'][0] = [-1.60273, -0.375253]   # c4d-u
    #zptab['amcoef'][1] = [0.277124, -0.198037]   # c4d-g
    #zptab['amcoef'][2] = [0.516382, -0.115443]   # c4d-r  changed a bit, fine
    #zptab['amcoef'][3] = [0.380338, -0.067439]   # c4d-i
    #zptab['amcoef'][4] = [0.074517, -0.067031]   # c4d-z
    #zptab['amcoef'][5] = [-1.07800, -0.060014]   # c4d-Y
    #zptab['amcoef'][6] = [1.111859, -0.083630]   # c4d-VR
    zptab['amcoef'][0] = [-1.34558, -0.449928]   # c4d-u
    zptab['amcoef'][1] = [0.324919, -0.192905]   # c4d-g
    zptab['amcoef'][2] = [0.514610, -0.106299]   # c4d-r  changed a bit, fine
    zptab['amcoef'][3] = [0.401896, -0.066105]   # c4d-i
    zptab['amcoef'][4] = [0.120886, -0.061114]   # c4d-z
    zptab['amcoef'][5] = [-1.03403, -0.043717]   # c4d-Y
    zptab['amcoef'][6] = [1.179283, -0.113482]   # c4d-VR
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
    print('APPLYING QA CUTS')
    fwhmthresh = 2.0  # arcsec, v2
    #filters = ['u','g','r','i','z','Y','VR']
    #nfilters = n_elements(filters)
    #zpthresh = [2.0,2.0,2.0,2.0,2.0,2.0,2.0]
    #zpthresh = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    badzpmask = np.ones(len(exptab),bool)
    exptab['bad'] = False
    exptab['badzpterm'] = False
    exptab['zpterm_corr'] = np.nan
    for i in range(nzptab):
        ind, = np.where((exptab['instrument'] == zptab['instrument'][i]) &
                        (exptab['filter'] == zptab['filter'][i]))
        print(' ')
        print('{:s}  {:d} exposures'.format(zptab['instrument'][i]+'-'+zptab['filter'][i],len(ind)))

        if len(ind)==0:
            print('No exposures for',zptab['instrument'][i]+'-'+zptab['filter'][i])
            continue

        exptab1 = exptab[ind]
        ## Fix Infinity/NAN values
        zpterm = exptab1['zpterm']
        bdzp, = np.where(np.isfinite(zpterm) == 0)  # fix Infinity/NAN
        if len(bdzp)>0:
            zpterm[bdzp] = 999999.9

        ## Correct "DES" zeropoints,  DES exposures are in electrons and
        ## CP are in ADU, so there's an offset of 2.5*log(gain)=2.5*log(4.41)=1.611
        plver3 = np.array([p[:3] for p in exptab1['plver']])
        basever = np.array([b.split('_')[-1] for b in exptab1['base']])
        gdes, = np.where((basever=='d1') | (basever=='d2') | (basever=='d3'))
        #gdes, = np.where(plver3=='DES')
        if len(gdes) > 0:
            print('  Offsetting',len(gdes),'DES exposure zero-points')
            zpterm[gdes] -= 1.611

        ## CORRECT K4M/KSB for exptime-dependence in the zero-points
        ##   this is because the image units are counts/sec.
        if zptab['instrument'][i] == 'k4m' or zptab['instrument'][i] == 'ksb':
            print('REMOVING EXPTIME-DEPENDENCE IN K4M/KSB ZEROPOINTS!!!')
            zpterm += 2.5*alog10(exptab1['exptime'])
        am = exptab1['airmass']
        mjd = exptab1['mjd']
        bdam, = np.where(am < 0.9)
        if len(bdam) > 0:
            am[bdam] = np.median(am)
        coo = SkyCoord(exptab1['ra'],exptab1['dec'],unit='degree',frame='icrs')
        glon = coo.galactic.l.degree
        glat = coo.galactic.b.degree

        # Measure airmass dependence
        gg0, = np.where((np.abs(zpterm) < 50) & (am < 2.0))
        ambins = np.arange(1.0,2.5,0.2)
        res1,bin_edges,binnumber = binned_statistic(am[gg0],zpterm[gg0],bins=ambins,statistic=np.nanmedian)
        xbins = bin_edges+0.5*(bin_edges[1]-bin_edges[0])
        gdbins, = np.where(np.isfinite(res1))
        coef0 = np.polyfit(xbins[gdbins],res1[gdbins],1)
        zpf = np.polyval(coef0,am)
        sig0 = dln.mad(zpterm[gg0]-zpf[gg0])
        # outlier rejection
        gg, = np.where(np.abs(zpterm-zpf) < np.maximum(3.5*sig0,0.2))
        res2,bin_edges,binnumber = binned_statistic(am[gg],zpterm[gg],bins=ambins,statistic=np.nanmedian)
        xbins = bin_edges+0.5*(bin_edges[1]-bin_edges[0])
        gdbins, = np.where(np.isfinite(res2))
        coef = np.polyfit(xbins[gdbins],res2[gdbins],1)
        zpf = np.polyval(coef,am)
        delta_zpterm = zpterm-zpf
        #coef0 = np.polyfit(am[gg0],zpterm[gg0],1)
        ##coef0 = robust_poly_fitq(am[gg0],zpterm[gg0],1)
        #zpf = np.polyval(coef0,am)
        #sig0 = dln.mad(zpterm[gg0]-zpf[gg0])
        #gg, = np.where(np.abs(zpterm-zpf) < np.maximum(3.5*sig0,0.2))
        #coef = np.polyfit(am[gg],zpterm[gg],1)
        ##coef = robust_poly_fitq(am[gg],zpterm[gg],1)
        print(' ',zptab['instrument'][i]+'-'+zptab['filter'][i],'airmass term:',coef)

        # Save figure of airmass dependence
        backend = matplotlib.rcParams['backend']
        matplotlib.use('Agg')
        plt.figure(1,figsize=(10,8))
        plt.clf()
        plt.hist2d(am[gg],zpterm[gg],norm=LN(),bins=50,cmap='Greys')
        plt.scatter(xbins[gdbins],res2[gdbins],s=50,marker='+',c='r')
        #plt.scatter(am[gg],zpterm[gg],s=20)
        plt.plot([0.9,np.max(am[gg])],np.polyval(coef,[0.9,np.max(am[gg])]),c='r')
        plt.xlabel('Airmass',fontsize=16)
        plt.ylabel('Relative '+zptab['filter'][i]+' zeropoint (mag)',fontsize=16)
        plt.title(zptab['filter'][i]+' zeropoint airmass dependence',fontsize=18)
        txt = '{:s} zpterm = {:.3f} * AM + {:.3f}'.format(zptab['filter'][i],coef[1],coef[0])
        plt.annotate(txt,xy=[0.9,0.9],xycoords='axes fraction',ha='right',c='r',fontsize=15)
        plt.savefig(plotsdir+'/'+zptab['filter'][i]+'_zpterm_airmass.png',bbox_inches='tight')
        plt.close()
        print('  Saving to',zptab['filter'][i]+'_zpterm_airmass.png')
        matplotlib.use(backend)
                
                
        # Trim out bad exposures to determine the correlations and make figures
        gg, = np.where((np.abs(delta_zpterm) < np.maximum(3.5*sig0,0.2)) &
                       (exptab1['airmass'] < 2.0) & (exptab1['fwhm'] < 2.0) & (exptab1['rarms'] < 0.15) &
                       (exptab1['decrms'] < 0.15) &
                       (exptab1['wcscal']=='Successful') & (exptab1['zptermerr'] < 0.05) &
                       (exptab1['zptermsig'] < 0.08) &
                       ((exptab1['instrument'] != 'c4d') | (exptab1['zpspatialvar_nccd']<=5) |
                        ((exptab1['instrument']=='c4d') & (exptab1['zpspatialvar_nccd']>5) & (exptab1['zpspatialvar_rms']<0.1))) &
                       (np.abs(glat) > 10) & (exptab1['nrefmatch'] > 100) & (exptab1['exptime'] >= 30))
        ## I removed WCSCAL check because there are ~38k exposures with
        ## WCSCAL=Failed but my DECRMS and RARMS is small.
        ## and exptab1.wcscal eq 'Successful'
        print('  {:d} exposures used for temporal analysis'.format(len(gg)))

        # Zpterm with airmass dependence removed
        relzpterm = zpterm + 25   # 25 to get "absolute" zpterm
        relzpterm -= zptab['amcoef'][i][1]*(am-1)

        # Fit temporal variation in zpterm
        mjd0 = 56200
        xx = exptab1['mjd'][gg]-mjd0
        yy = relzpterm[gg]
        invvar = 1.0/exptab1['zptermerr'][gg]**2
        nord = 3
        bkspace = 200 #20
        knots = np.arange(np.min(xx)+0.5*bkspace,np.max(xx),bkspace)
        wt = 1.0/exptab1['zptermerr'][gg]
        wt /= np.nansum(wt)
        spl = dln.bspline(xx,yy,w=wt,knots=knots,nord=nord)
        yfit1 = spl(xx)
        sig1 = dln.mad(yy-yfit1)
        gd, = np.where((yy-yfit1) > -3*sig1)
        #sset1 = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=yfit1)
        #sig1 = mad(yy-yfit1)
        #gd = where(yy-yfit1 > -3*sig1,ngd)
        # refit
        #sset = bspline_iterfit(xx[gd],yy[gd],invvar=invvar[gd],nord=nord,bkspace=bkspace)
        #yfit = bspline_valu(xx,sset)
        #allzpfit = bspline_valu(exptab1['mjd']-mjd0,sset)
        spl2 = dln.bspline(xx[gd],yy[gd],w=wt[gd],knots=knots,nord=nord)
        allzpfit = spl2(exptab1['mjd']-mjd0)
        lowoutofrange, = np.where(exptab1['mjd']-mjd0 < np.min(xx[gd]))
        if len(lowoutofrange)>0:
            allzpfit[lowoutofrange] = spl2(np.min(xx[gd]))
        hioutofrange, = np.where(exptab1['mjd']-mjd0 > np.max(xx[gd]))
        if len(hioutofrange)>0:
            allzpfit[hioutofrange] = spl2(np.max(xx[gd]))
                
        # Save figure of temporal dependence
        backend = matplotlib.rcParams['backend']
        matplotlib.use('Agg')
        plt.figure(1,figsize=(10,8))
        plt.hist2d(xx,yy,bins=100,norm=LN(),cmap='Greys')
        #plt.scatter(xx,yy,s=20,c='k')
        plt.scatter(xx[gd],yy[gd],s=2,c='blue')
        plt.scatter(xx,spl2(xx),s=2,c='r')
        allmjd = np.arange(0,np.max(exptab1['mjd']-mjd0))
        allfit = spl2(allmjd)
        plt.plot(allmjd,allfit,c='r')
        plt.xlabel('MJD-MJD0',fontsize=16)
        plt.ylabel(zptab['filter'][i]+' zero-point magnitude',fontsize=16)
        plt.title(zptab['filter'][i]+' temporal dependence',fontsize=18)
        #plt.annotate(txt,xy=[0.9,0.9],xycoords='axes fraction',ha='right',c='r',fontsize=15)
        plt.savefig(plotsdir+'/'+zptab['filter'][i]+'_zpterm_mjd.png',bbox_inches='tight')
        print('  Saving to',zptab['filter'][i]+'_zpterm_mjd.png')
        matplotlib.use(backend)

                
        # Remove temporal variations to get residual values
        relzpterm -= allzpfit
        exptab1['zpterm_corr'] = relzpterm
        exptab['zpterm_corr'][ind] = relzpterm

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
        print('  {:d} exposures with ZPTERM below the threshold ({:.1f}%)'.format(len(bdind),len(bdind)/len(ind)*100))
        if len(gdind) > 0:
            badzpmask[ind[gdind]] = 0
            exptab['bad'][ind[bdind]] = True
                    
    # Get bad DECaLS and SMASH exposures
    datadir = utils.datadir()
    obslogdir = os.path.abspath(datadir+'../../../obslog/v3/')
    badexp = np.zeros(len(exptab),bool)
    smashlines = dln.readlines(obslogdir+'/smash_badexposures.txt',comment='#',noblank=True)
    smashexpnum = [int(l.split()[0]) for l in smashlines]
    _,ind1,ind2 = np.intersect1d(exptab['expnum'],smashexpnum,return_indices=True)
    if len(ind1) > 0:
        badexp[ind1] = True
        badexp[ind1] = (badexp[ind1] & (exptab['instrument'][ind1] == 'c4d'))   # make sure they are DECam exposures
    decalslines = dln.readlines(obslogdir+'/decals_bad_expid.txt',comment='#',noblank=True)
    decalsexpnum = [int(l.split()[0]) for l in decalslines]
    _,ind1,ind2 = np.intersect1d(exptab['expnum'],decalsexpnum,return_indices=True)
    if len(ind1) > 0:
        badexp[ind1] = True
        badexp[ind1] = (badexp[ind1] & (exptab['instrument'][ind1] == 'c4d'))   # make sure they are DECam exposures
    mzlslines = dln.readlines(obslogdir+'/mzls_bad_expid.txt',comment='#',noblank=True)
    mzlsexpnum = [int(l.split()[0]) for l in mzlslines]
    _,ind1,ind2 = np.intersect1d(exptab['expnum'],mzlsexpnum,return_indices=True)
    if len(ind1) > 0:
        badexp[ind1] = True
        badexp[ind1] = (badexp[ind1] & (exptab['instrument'][ind1] == 'k4m'))   # make sure they are Mosaic3 exposures

    ## Zero-point spatial variability threshold
    ##  varies with galactic latitude
    ##  |b|>10   0.15
    ##  |b|<=10  0.55
    zpspvarthresh = (np.abs(exptab['glat']) > 10)*0.15 + (np.abs(exptab['glat']) <= 10)*0.55


    # Final QA cuts
    #  Many of the short u-band exposures have weird ZPTERMs, not sure why
    #  There are a few exposures with BAD WCS, RA>360!
    bdexp, = np.where((exptab['fwhm'] > fwhmthresh) |                         # bad seeing
                      (exptab['ra'] > 360) |                                  # bad WCS/coords
                      (exptab['rarms'] > 0.15) | (exptab['decrms'] > 0.15) |     # bad WCS
                      (badzpmask == True) |                                # bad ZPTERM
                      (exptab['zptermerr'] > 0.05) |                          # bad ZPTERMERR
                      (exptab['nrefmatch'] < 5) |                             # few phot ref match
                      (badexp == 1) |                                      # bad SMASH/LS exposure
                      ((exptab['instrument'] == 'c4d') & (exptab['zpspatialvar_nccd'] > 5) &
                       (exptab['zpspatialvar_rms'] > zpspvarthresh)))  # bad spatial zpterm
                      ##exptab[]wcscal ne 'Successful' or $                    ##CP WCS failure   TOO MANY FAILED
                          #exptab[]ngoodchipwcs lt exptab[]nchips or $                # not all chips astrom calibrated
        # rarms/decrms, nrefmatch
    print('QA cuts remove ',len(bdexp),' exposures')    
    exptab = np.delete(exptab,bdexp)
    exptab = Table(exptab)
    
    # Save the final exposures list
    outfile = os.path.join(basedir,'lists','nsc_instcal_combine_exposures.fits')
    print('Writing final list of exposures to',outfile)
    exptab.write(outfile,overwrite=True)
    
    import pdb; pdb.set_trace()
    

def combinehealpix(version='v4',nocuts=False):
    """ Make the healpix exposure and chip lists for the final list of exposures """

    dldir,mssdir,localdir = utils.rootdirs()
    basedir = os.path.join(dldir,'instcal/',version)
    if os.path.exists(localdir+'dnidever/nsc/instcal/'+version+'/')==False:
        os.makedirs(localdir+'dnidever/nsc/instcal/'+version+'/')
    plotsdir = basedir+'plots/'
    if os.path.exists(plotsdir)==False:
        os.makedirs(plotsdir)
    nside = 128
    t0 = time.time()

    # Load the final list of exposures
    expfile = os.path.join(basedir,'lists','nsc_instcal_combine_exposures.fits')
    if os.path.exists(expfile)==False:
        print(expfile,'not found.  Make sure to run summary.combine_qacuts() first')
        return
    exptab = Table.read(expfile)

    # Go through the calibration summary files
    localdir = '/home1/09970/dnidever/scratch1/nsc/instcal/'
    calibfiles = glob(os.path.join(localdir,version,'summary','calibrate','calibrate_summary*.fits'))
    calibfiles.sort()
    print(len(calibfiles),' calibration files')

    dbfile = os.path.join(basedir,'lists','nsc_instcal_combine_healpix_list.db')
    if os.path.exists(dbfile):
        os.remove(dbfile)

    for i in range(len(calibfiles)):
        ctab = Table.read(calibfiles[i],2)
        index = dln.create_index(ctab['base'])
        _,ind1,ind2 = np.intersect1d(index['value'],exptab['base'],return_indices=True)
        print('{:d} {:<25s} {:>8d} {:>8d}'.format(i+1,os.path.basename(calibfiles[i]),
                                                  len(ind1),len(index['value'])))
        # build index array of chips to keep
        keepind = []
        for j in range(len(ind1)):
            ind = index['index'][index['lo'][ind1[j]]:index['hi'][ind1[j]]+1]
            keepind.append(ind)
        keepind = np.hstack(keepind)
        chipinfo = ctab[keepind]
        chipinfo = chipinfo[['instrument','measfile','base','ccdnum','nsources',
                             'cenra','cendec','zpterm','depth95']]
        chipinfo['pix'] = hp.ang2pix(nside,chipinfo['cenra'],chipinfo['cendec'],lonlat=True)
        # Now write this to the database
        db.writecat(chipinfo,dbfile,'hlist')
        
    # Indexing the pix column
    print('Creating index on pix column')
    db.createindex(dbfile,col='pix',table='hlist',unique=False,verbose=True)

    import pdb; pdb.set_trace()


    #dbfile='/net/dl2/dnidever/nsc/instcal/v3/lists/nsc_instcal_combine_healpix_list.db'
    #dbutils.writecat(healstr,dbfile,'hlist')
    #dbutils.createindex(dbfile,col='pix',table='hlist',unique=False,verbose=True)
    #out=dbutils.query(dbfile,'hlist',where='PIX=65368')

    print('dt = {:.1f} sec.'.format(time.time()-t0))



def combine():
    """ Make the nsc_combine_summary.fits summary file """
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
    t0 = time.time()

    import pdb; pdb.set_trace()

