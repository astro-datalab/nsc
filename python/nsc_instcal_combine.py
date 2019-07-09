#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
#from astropy.wcs import WCS
from astropy.table import Table
#import time
#import shutil
#import re
#import subprocess
#import glob
#import logging
#import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve

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

    # on thing/hulk use
    if (host == "thing") | (host == "hulk"):
        dir = "/dl1/users/dnidever/nsc/instcal/"+verdir
        mssdir = "/mss1/"
        tmproot = "/d0/dnidever/nsc/instcal/"+verdir+"tmp/"
    # on gp09 use
    if (host == "gp09") | (host == "gp08") | (host == "gp07") | (host == "gp06") | (host == "gp05"):
        dir = "/net/dl1/users/dnidever/nsc/instcal/"+verdir
        mssdir = "/net/mss1/"
        tmproot = "/data0/dnidever/nsc/instcal/"+verdir+"tmp/"

    t0 = time.time()


    # Not enough inputs
    n = len(sys.argv)
    if n < 1:
        print "Syntax - nsc_instcal_combine pix --version v# --nside ### --redo YES/NO --outdir OUTDIR"
        sys.exit()

    # Inputs
    pix = args.pix
    version = args.version
    nside = args.nside
    redo = args.redo
    if (redo=='True' | redo=='TRUE' | redo=='Yes' | redo='YES' | redo=='Y' | redo='y'):
        redo = True
    else:
        redo = False
    outdir = args.outdir
    if outdir == '': outdir=dir+'combine/'

    # Check that the directory exist
    if os.path.exists(expdir) == False:
        print(expdir+" NOT FOUND")
        sys.exit()

    # Check if output file already exists
    subdir = str(int(pix)/1000)    # use the thousands to create subdirectory grouping
    outfile = outdir+'/'+subdir+'/'+str(pix)+'.fits'
    if (os.path.exists(outfile) | os.path.exists(outfile+'.gz')) & ~redo:
        print(outfile+' EXISTS already and REDO not set')
        sys.exit()

    print("Combining InstCal SExtractor catalogs for Healpix pixel = "+str(pix,2))

    # Load the list
    listfile = localdir+'dnidever/nsc/instcal/'+version+'/nsc_healpix_list.fits'
    if os.path.exists(list) is False:
        print(listfile+" NOT FOUND")
        sys.exist()
    healstr = fits.getdata(listfile,1)
    index = fits.getdata(listfile,2)
    # Find our pixel
    ind, = np.where(index['pix'] == pix)
    nind = len(ind)
    if nind == 0:
        print("No entries for Healpix pixel '"+str(pix)+"' in the list")
        sys.exit()
    ind = ind[0]
    list = healstr[index[ind].lo:index[ind].hi]
    nlist = n_elements(list)
    # GET EXPOSURES FOR NEIGHBORING PIXELS AS WELL
    #  so we can deal with the edge cases
    NEIGHBOURS_RING,nside,pix,neipix,nneipix
    for i=0,nneipix-1 do begin
        ind = where(index.pix eq neipix[i],nind)
        if nind gt 0 then begin
            ind = ind[0]
            list1 = healstr[index[ind].lo:index[ind].hi]
            push,list,list1
    # Use entire exposure files
    # Get unique values
    ui = np.uniq(list.file,sort(list.file))
    list = list[ui]
    nlist = len(list)
    print(str(nlist)+' exposures that overlap this pixel and neighbors')


