import os
import numpy as np
from glob import glob
from datetime import datetime
from astropy.io import fits,ascii
from astropy.table import Table,vstack

def readlines(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    return lines

def create_index(arr):
    narr = len(arr)
    si = np.argsort(arr)
    sarr = np.array(arr)[si]
    brklo, = np.where(sarr != np.roll(sarr,1))
    nbrk = len(brklo)
    brkhi = np.hstack((brklo[1:nbrk]-1,narr-1))
    num = brkhi-brklo+1
    index = {'index':np.atleast_1d(si),'value':np.atleast_1d(sarr[brklo]),
             'num':np.atleast_1d(num),'lo':np.atleast_1d(brklo),'hi':np.atleast_1d(brkhi)}
    return index

def checklogs():
    """ Check the logs for exceptions."""

    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_20241221_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_20250101_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250108_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250117_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250125_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250125b_left.fits.gz')
    #tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250130_left.fits.gz')
    tab = Table.read('../../lists/decam_instcal_list_exptime10sec_corrupted_fwhm2_20250211_left.fits.gz')
    tab['base'] = tab['base'].astype(str)

    # Remove exposures that have meas files
    #measlines = readlines('../allmeas_122924.txt')
    #measlines = readlines('../allmeas_010225.txt')
    #measlines = readlines('../allmeas_010425.txt')
    #measlines = readlines('../allmeas_010725.txt')
    measlines = readlines('../allmeas_021425.txt')
    measbase = [os.path.basename(l)[:-10] for l in measlines]
    _,ind1,ind2 = np.intersect1d(measbase,tab['base'],return_indices=True)
    #del tab[ind2]
    #base = np.delete(base,ind2)
    tab['done'] = False
    tab['done'][ind2] = True
    bd, = np.where(tab['done']==False)
    print(len(bd),' exposures do not have meas.fits files')

    # Find the latest log file for the ones that are left
    #loglines = readlines('../alllogs_122924.txt')
    #loglines = readlines('../alllogs_010225.txt')
    #loglines = readlines('../alllogs_010425.txt')
    #loglines = readlines('../alllogs_010725.txt')
    loglines = readlines('../alllogs_021425.txt')
    logbase = [os.path.basename(l).split('.')[0] for l in loglines]
    ulogbase = np.unique(logbase)
    index = create_index(logbase)
    #_,ind1,ind2 = np.intersect1d(logbase,base,return_indices=True)
    #_,ind1,ind2 = np.intersect1d(ulogbase,base,return_indices=True)
    _,ind1,ind2 = np.intersect1d(index['value'],tab['base'],return_indices=True)
    # 62002 of 63595 have log files
    print(len(ind1),' exposures have log files')

    #import pdb; pdb.set_trace()

    loglines = np.array(loglines)

    #outfile = 'logtab_122924.fits'
    #outlst = 'logtab_122924.lst'
    outfile = 'logtab_021425.fits'
    outlst = 'logtab_021425.lst'
    fout = open(outlst,'w')
    #fout = open(outlst,'a')


    #dt = [('base',str,50),('nlogfiles',int),('logfile',str,100),('timestamp',float),
    #      ('timespan',float),('fwhm',float),('corrupted',bool),('exception',bool),
    #      ('exceptline',str,100),('nospace',bool),('psferror',bool),('npsfstars',int),
    #      ('lastchip',int)]
    #logtab = np.zeros(len(ind1),dtype=np.dtype(dt))
    #logtab['fwhm'] = np.nan
    #logtab['npsfstars'] = -1
    #logtab['lastchip'] = -1
    #logtab['timespan'] = -1.0
    tab['nlogfiles'] = -1
    tab['logfile'] = np.zeros(len(tab),(str,100))
    tab['timestamp'] = np.nan
    tab['timespan'] = -1.0
    tab['fwhm'] = np.nan
    tab['corrupted'] = False
    tab['exception'] = False
    tab['exceptline'] = np.zeros(len(tab),(str,100))
    tab['nospace'] = False
    tab['nsources'] = -1
    tab['psferror'] = False
    tab['npsfstars'] = -1
    tab['lastchip'] = -1
    start = 0

    # Start from where we left off
    #logtab = Table.read(outfile)
    #start = 50813

    for i in range(start,len(ind1)):
        tindx = ind2[i]
        ind = index['index'][index['lo'][ind1[i]]:index['hi'][ind1[i]]+1]
        nind = len(ind)
        #logtab['base'][i] = base[ind2[i]]
        tab['nlogfiles'][tindx] = nind
        loglines1 = loglines[ind]
        loglines1 = np.array(['.'+l for l in loglines1])
        timestamp = np.zeros(nind,float)
        for j in range(nind):
            timestamp[j] = os.path.getmtime(loglines1[j])
        latest = np.argmax(timestamp)
        tab['timestamp'][tindx] = timestamp[latest]
        tab['logfile'][tindx] = loglines1[latest]
        if os.path.exists(loglines1[latest]):
            llines = readlines(loglines1[latest])
        else:
            llines = []
        if os.path.exists(loglines1[latest].replace('.log','.err')):
            elines = readlines(loglines1[latest].replace('.log','.err'))
        else:
            elines = []
        # get processing time
        #if os.path.exists(logtab['logfile'][i]):
        #    ctime = os.path.getctime(logtab['logfile'][i])
        #    mtime = os.path.getmtime(logtab['logfile'][i])
        #    logdt = mtime-ctime
        #    logtab['timespan'][i] = logdt
        if len(llines)>0:
            # at end of .log it often has "Total time = "
            timelines = [l for l in llines if l.find('Total time')>-1]
            if len(timelines)>0:
                dum = timelines[-1].split()
                dt = float(dum[3])
                tab['timespan'][tindx] = dt
        if len(elines)>0 and tab['timespan'][tindx]<=0:
            # otherwise .err has timestamp on each line
            timelines = [e for e in elines if (e.find('2')>-1 and e.find('[')>-1 and e.find(']')>-1)]
            if len(timelines)>0:
                try:
                    tl1 = timelines[0]
                    t1 = datetime.fromisoformat(tl1[:19])
                    tl2 = timelines[-1]
                    t2 = datetime.fromisoformat(tl2[:19])
                    dt = t2-t1
                    dt = dt.seconds
                    tab['timespan'][tindx] = dt
                except:
                    pass
        lexception = [e for e in llines if e.find('Exception:')>-1]
        lexception += [e for e in llines if e.find('Error: ')>-1]
        if len(lexception)>0:
            tab['exception'][tindx] = True
            tab['exceptline'][tindx] = ','.join(lexception)
        eexception = [e for e in elines if e.find('Exception:')>-1]
        eexception += [e for e in elines if e.find('Error: ')>-1]
        if len(eexception)>0:
            tab['exception'][tindx] = True
            tab['exceptline'][tindx] = ','.join(eexception)
        fwlines = [e for e in elines if e.find('FWHM = ')>-1]
        if len(fwlines)>0:
            fwhmarr = len(fwlines)*[None]
            for f in range(len(fwlines)):
                dum = fwlines[f].split('FWHM = ')
                dum2 = dum[-1].split('arcsec')
                fwhmarr[f] = float(dum2[0])
            fwhm = np.mean(np.array(fwhmarr))
            tab['fwhm'][tindx] = fwhm
        lbad = [e for e in llines if e.find('Some files have problems')>-1]
        if len(lbad)>0:
            tab['corrupted'][tindx] = True
        enospace = [e for e in elines if e.find('No space left')>-1]
        if len(enospace)>0:
            tab['nospace'][tindx] = True
        # Nsource detected
        ensources = [e for e in elines if e.find('sources detected')>-1]
        ensources = [e for e in ensources if e.find('new sources detected')==-1]
        if len(ensources)>0:
            dum = ensources[-1].split('sources detected')
            nsources = int(dum[0].split()[-1])
            tab['nsources'][tindx] = nsources
        # Check psf problems
        epsfexception = [e for e in elines if e.find('no psf success')>-1]
        if len(epsfexception)>0:
            tab['psferror'][tindx] = True
        # Get number of psf stars
        candlines = [e for e in elines if e.find('suitable candidates were found')>-1]
        if len(candlines)>0:
            dum = candlines[-1].split('suitable')
            npsfstars = int(dum[0].split()[-1])
            tab['npsfstars'][tindx] = npsfstars
        if len(candlines)==0:
            candlines = [e for e in elines if e.find('PSF stars')>-1 and e.find('INFO')>-1]
            candlines = [e for e in candlines if e.find('PSF stars')==len(e)-9]
            candlines = [e for e in candlines if e.find('sexpickpsf')==-1]
            candlines = [e for e in candlines if e.find('Too few PSF stars')==-1]
            if len(candlines)>0:
                dum = candlines[-1].split('PSF stars')
                try:
                    npsfstars = int(dum[0].split()[-1])
                except:
                    print('npsfstars problem')
                    import pdb; pdb.set_trace()
                tab['npsfstars'][tindx] = npsfstars
        # Check what chip it to go before it crashed
        chlines = [e for e in elines if e.find('Processing subimage')>-1]
        if len(chlines)>0:
            dum = chlines[-1].split('subimage')
            lastchip = int(dum[1].split()[0])
            tab['lastchip'][tindx] = lastchip
        cmt = '{:d} {:} {:d} {:} {:} {:}'.format(i+1,tab['base'][tindx],nind,tab['exception'][tindx],tab['corrupted'][tindx],tab['nospace'][tindx])
        if tab['npsfstars'][tindx]>=0:
            cmt += '  Npsfstars='+str(tab['npsfstars'][tindx])
        if np.isfinite(tab['fwhm'][tindx]):
            cmt += '  FWHM='+str(tab['fwhm'][tindx])+' arcsec'
        if tab['lastchip'][tindx] > 0:
            cmt += '  lastchip='+str(tab['lastchip'][tindx])
        if tab['timespan'][tindx] > -1:
            cmt += '  dt={:.1f}'.format(tab['timespan'][tindx])
        print(cmt)
        fout.write(cmt+'\n')
        fout.flush()

        if i % 1000 == 0 and i>0:
            Table(tab).write(outfile,overwrite=True)

        #import pdb; pdb.set_trace()

    print(np.sum(tab['corrupted']),'corrupted')

    Table(tab).write(outfile,overwrite=True)

if __name__ == '__main__':
    checklogs()
