#!/usr/bin/env python
import os
import sys
import numpy as np
from glob import glob
import subprocess
import shutil
from astropy.io import fits
from astropy.table import Table,vstack
import traceback

def taccifyname(filename):
    """ Modify mss filenames for tacc"""
    imagedir = '/scratch1/09970/dnidever/nsc/instcal/v4/images'
    if filename[:9]=='/net/mss1':
        base = os.path.basename(filename)
        if base[:3]=='c4d':
            instrument = base[:3]
            night = '20'+base[4:10]
        else:
            instrument = 'c4d'   # assume all decam for now 
            head = fits.getheader(filename,0)
            dateobs = head['DATE-OBS']
            night = dateobs[:4]+dateobs[5:7]+dateobs[8:10]
        newfilename = os.path.join(imagedir,instrument,night[:4],night,base)
    else:
        newfilename = filename
    return newfilename

def readlines(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    return lines

def checkfile(filename):
    """ run nsccheck and parse the results."""
    # Check that the InstCal files are not corrupted
    cmd = ['python3','/home1/09970/dnidever/projects/nsc/bin/nsccheck',filename]
    res = subprocess.run(cmd,shell=False,capture_output=True)
    txt = res.stdout
    if type(txt)==bytes:
        txt = txt.decode()
    txt = txt.split('\n')
    if 'Checking' in txt[0]:
        txt.pop(0)
    if txt[-1]=='':
        txt.pop(len(txt)-1)
    if len(txt)==0:
        okay = False
    else:
        okay = txt[0][:2]=='OK'
    return okay

def checkcorrupted(number,nprocs):
    """
    Check if nsc c4d files are corrupted
    """
    basedir = '/scratch1/09970/dnidever/nsc/instcal/v4/'

    tab = Table.read('/scratch1/09970/dnidever/nsc/instcal/v4/lists/decam_instcal_list_exptime10sec_20241221_left.fits.gz')
    nperproc = len(tab)//int(nprocs) + 1
    lo = (int(number)-1)*nperproc
    hi = lo + nperproc
    print('Number: {:d}'.format(int(number)))
    print('Nprocs: {:d}'.format(int(nprocs)))
    print('Nperproc: {:d}'.format(nperproc))
    print('lo: {:d}'.format(lo))
    print('hi: {:d}'.format(hi))
    tab = tab[lo:hi+1]
    ntab = len(tab)
    tab['okay'] = False
    tab['corrupted'] = np.zeros((ntab,3),bool)

    # Make sure the filename columns are long enough
    tab['fluxfile'] = tab['fluxfile'].astype((str,150))
    tab['wtfile'] = tab['wtfile'].astype((str,150))
    tab['maskfile'] = tab['maskfile'].astype((str,150))

    imlistfile = basedir+'images/allimages.lst'
    imlines = readlines(basedir+'images/allimages.lst')
    imlines = [basedir+'images/'+f[1:] for f in imlines]  # make absolute
    imlines = np.array(imlines)
    imbase = np.array([os.path.basename(f) for f in imlines])
    fluxfiles = tab['fluxfile']
    fluxbase = np.array([os.path.basename(f) for f in fluxfiles])
    _,ind1,ind2 = np.intersect1d(fluxbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' fluxfile names')
        tab['fluxfile'][ind1] = imlines[ind2]
    wtfiles = tab['wtfile']
    wtbase = np.array([os.path.basename(w) for w in wtfiles])
    _,ind1,ind2 = np.intersect1d(wtbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' wtfile names')
        tab['wtfile'][ind1] = imlines[ind2]
    maskfiles = tab['maskfile']
    maskbase = np.array([os.path.basename(m) for m in maskfiles])
    _,ind1,ind2 = np.intersect1d(maskbase,imbase,return_indices=True)
    if len(ind1)>0:
        print(len(ind1),' maskfile names')
        tab['maskfile'][ind1] = imlines[ind2]

    outfile = '/scratch1/09970/dnidever/nsc/instcal/v4/images/c4d/checkcorrupted_'+str(number)+'.fits'
    outlog = '/scratch1/09970/dnidever/nsc/instcal/v4/images/c4d/checkcorrupted_'+str(number)+'.lst'

    # If the output file exists already, then load it and start where we left off
    if os.path.exists(outfile):
        print('loading previously saved table')
        tab = Table.read(outfile)
        start = np.max(np.where(tab['okay']==True)[0])+1
        print('starting with index',start)
        fout = open(outlog,'a')  # append
        fout.write('\n\n')
        fout.write('starting again with index '+str(start)+'\n')
        fout.flush()
    else:
        start = 0
        fout = open(outlog,'w')

    for i in range(start,ntab):
        #fluxfile = taccifyname(tab['fluxfile'][i])
        #wtfile = taccifyname(tab['wtfile'][i])
        #maskfile = taccifyname(tab['maskfile'][i])
        fluxfile = tab['fluxfile'][i]
        wtfile = tab['wtfile'][i]
        maskfile = tab['maskfile'][i]

        # Check that the InstCal files are not corrupted
        fchk = checkfile(fluxfile)
        wchk = checkfile(wtfile)
        mchk = checkfile(maskfile)

        okay = [fchk,wchk,mchk]
        corrupt = [not o for o in okay]
        tab['corrupted'][i] = corrupt
        tab['okay'][i] = np.sum(okay)==3

        cmt = '{:5d} {:30s}  {:}  {:}'.format(i+1,tab['base'][i],tab['okay'][i],str(tab['corrupted'][i]))
        print(cmt)
        fout.write(cmt+'\n')
        fout.flush()

        if i % 100 == 0 and i>0:
            print('saving to',outfile)
            tab.write(outfile,overwrite=True)

    tab.write(outfile,overwrite=True)
    fout.close()


if __name__ == '__main__':
    if len(sys.argv)==3:
        checkcorrupted(sys.argv[1],sys.argv[2])
    else:
        print('syntax: checkcorrupted.py number nprocs')
