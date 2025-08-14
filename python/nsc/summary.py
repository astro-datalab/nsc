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
        chtab = Table(hdu[2].data)
        chtab.write(files[i].replace('.fits','_chip.csv'),format='csv')
        expdata.append(Table(hdu[1].data))
        #chipdata.append(Table(hdu[2].data))
        hdu.close()

    # Write to final output file
    exptab = vstack(expdata)
    #chiptab = vstack(chipdata)

    # Add healpix to each table
    nside = 128
    exptab['pix128'] = hp.ang2pix(nside,exptab['ra'],exptab['dec'],lonlat=True)
    #chiptab['pix128'] = hp.ang2pix(nside,chiptab['ra'],chiptab['dec'],lonlat=True)

    outfile = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/lists/nsc_calibrate_summary.fits'
    print('Writing summary results to',outfile)
    ohdu = fits.HDUList()
    ohdu.append(fits.table_to_hdu(exptab))
    #ohdu.append(fits.table_to_hdu(chiptab))
    ohdu.writeto(outfile,overwrite=True)
    ohdu.close()

def combine():
    """ Make the nsc_combine_summary.fits summary file """
    pass
