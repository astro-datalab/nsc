#!/usr/bin/env python

# download files to outdir

# Imports
from astropy.table import Table,vstack
import numpy as np
import os
from phot_tempest import *
import sys

# Main Code
if __name__=="__main__":

    # Arguments:
    # 1. file name
    # 2. file type (f, w, m, r) for (flux, weight, mask, or raw)
    # 3. input list (c4d, ksb, k4m)
    # 4. outdir for where files will be downloaded

    if len(sys.argv)<5:
        input_name = sys.argv[1]
        input_type = sys.argv[2]
        outdir = sys.argv[3]
        if input_type=="r": filecol = "RAWNAME"
        elif input_type=="f": filecol = "FLUXFILE"
        elif input_type=="w": filecol = "WTFILE"
        elif input_type=="m": filecol = "MASKFILE"

        listfile1 = "decam_instcal_list.fits.gz"
        listfile2 = "bok90prime_instcal_list.fits.gz"
        listfile3 = "mosaic3_instcal_list.fits.gz"

        tab1 = Table.read("/home/x25h971/nsc/instcal/v4/lists/"+listfile1)
        tab2 = Table.read("/home/x25h971/nsc/instcal/v4/lists/"+listfile2)
        tab3 = Table.read("/home/x25h971/nsc/instcal/v4/lists/"+listfile3)
        tab=vstack([tab1,tab2,tab3])
        filenames = np.array([tab[filecol][i].strip().split("/")[-1] for i in range(0,len(tab))])
        expo = tab[filenames==input_name]

        rawname = expo['RAWNAME'][0].strip()
        fluxfile = expo['FLUXFILE'][0].strip().split("/")[-1]
        wtfile = expo['WTFILE'][0].strip().split("/")[-1]
        maskfile = expo['MASKFILE'][0].strip().split("/")[-1]


    else:
        rawname = sys.argv[1]
        fluxfile = sys.argv[2]
        wtfile = sys.argv[3]
        maskfile = sys.argv[4]
        outdir = sys.argv[5]


    # check if files already exist!
    if os.path.exists(outdir+fluxfile) and os.path.exists(outdir+wtfile) and os.path.exists(outdir+maskfile):
        print("all files exist!")
    else:
        print("file names = ",rawname,fluxfile,wtfile,maskfile,outdir)
        getdata(rawname,fluxfile,wtfile,maskfile,outdir)
