#!/usr/bin/env python

# exposure_copy.py copies exposures from the NSC DR3 decam_instcal_list
# to a common repository in your current working directory (wherever this
# script is), to allow for easy but manual Globus transfers of exposures to
# the CCA.

# Arguments:
#    explist - fits file with table of exposures to copy/transfer, with
#              fluxfile, wtfile, and maskfile columns
#    transfer_lim - maximum amount of data [MB] per transfer (sets the number
#                   of exposures per sub-repository)
#    nside - HEALPix nside to sort healpix into; each has 

# The script will create a directory called "exp_copies/YYYYMMDDHHMMSS/" in
# your current working directory, and in within each healpix subdir there
# will be X number of "transfer repositories" with Y exposures each.  Each "transfer repository"
# contains all the files for one manual Globus transfer of <transfer_lim> MB

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from argparse import ArgumentParser
from astropy.time import Time,TimezoneInfo
from astropy.table import *
import astropy.units as u
from datetime import datetime,timezone
import healpy as hp
import numpy as np
import os
import socket
import shutil
import sys
import time

#-----------------------------------------------------------------------------
# Funcs
#-----------------------------------------------------------------------------

def makedir(dir):
    '''makes a directory with name "dir" if it does not exist
    Arguments:
    ----------
    dir (str)
            directory name
    Returns:
    --------
        None; directory "dir" is created if not already there
    '''
    if not os.path.exists(dir):
        os.mkdir(dir)


#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # -----------
    # -- Setup --
    # -----------
    
    # Initiate input arguments
    parser = ArgumentParser(description='Copy exposure filess from "explist" to several common repositories (of a set size "transfer_lim") to transfer with Globus from DataLab to CCA servers')
    parser.add_argument('--explist',type=str, nargs=1, default="decam_exposures_tomeas.fits.gz",help='Filename with list of exposures with rawname,flux/wt/maskfile columns (must be in same directory as script)')
    parser.add_argument('--transfer_lim',type=str,nargs=1,default="12000",help="max amount of data per transfer, in MBytes (default 4000, determines number of exposures (3 files each) per repo)") # each file is ~250-400 MB, so 1 exposure <1200MB of data
    parser.add_argument('--nside',type=str,nargs=1,default="32",help="HEALPix nside to split exposures into") 
    args = parser.parse_args()

    # Get input arguments
    explistfile = args.explist[0]               # name of file with exposure list
    print("explistfile = ",explistfile)
    transfer_limit = int(args.transfer_lim)  # max amount of data per transfer [MBytes]
    print("transfer limit = ",transfer_limit)
    ns = args.nside

    # Start time, get host & local directory
    t00 = time.time()
    utc_bzn = TimezoneInfo(utc_offset=-7*u.hour)
    ut = Time(datetime.now(),scale='utc')
    t0 = ut.to_datetime(timezone=utc_bzn)
    timestamp = str(t0.year)+str(t0.month).zfill(2)+str(t0.day).zfill(2)+str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2)
    localdir = os.getcwd()                      # where this script is
    mssdir = "/net/mss1/"                       # where exposures are stored on dl machines
    ####mssdir = "/home/group/davidnidever/nsc/instcal/v4/exposures/"
    copydir = localdir+"/exposure_copies/"+timestamp+"/"          # where exposures will be copied to, divided
    makedir(copydir)                                            # into one repo per planned transfer

    # ----------------------------------------
    # -- Get list of exposures to transfer --
    # ----------------------------------------
    
    # Read in list, check for missing columns
    exposures = Table.read(explistfile)         # list of exposures
    expcols = exposures.colnames                # columns
    for col in ["RAWNAME","FLUXFILE","WTFILE","MASKFILE","RA","DEC"]:
        if col not in expcols:
            print("missing colum ",col," from exposure input list ",explistfile,"!!!")
            sys.exit()
    exposures['ring'+str(ns)] = hp.ang2pix(nside=ns,exposures['RA'],exposures['DEC'],lonlat=True)
    nhpix = hp.nside2npix(ns)


    # Group exposures into subdirs (1 per hpix) and within that, into further subdirs (1 per transfer batch)
    nexposures = len(exposures)                 # Number of exposures
    nfiles = nexposures*3                       # Number of exposure files
    filesize = 400                              # max exposure file size [MB]
    exposuresize = filesize*3
    nexposures_per_repo = transfer_limit//exposuresize  # Nmber of files per repository
    ####nrepos = (nexposures//nexposures_per_repo)+1        # Number of repositories/transfers
    ####subdir_list = np.array([copydir+"expgroup_"+str(i).zfill(len(str(nrepos))) for i in range(nrepos)]) # list of subdirectories in copydir
    subdir_list = np.array([copydir+"ring32_"+str(i).zfill(len(str(nhpix))) for i in range(nhpix)]) # list of subdirectories in copydir
    print("splitting ",nexposures," exposures into ",np.unique(subdir_list)," HPix repositories")
    #exposures['copydir'] = [subdir_list[(i-nrepos*(i//nrepos))] for i in range(nexposures)]
    #print(exposures['copydir'])


    # Copy exposure files, grouped in directories by transfer batch
    for sbdr in subdir_list:
        makedir(sbdr)
        #transnum = exposures['copydir'] == sbdr
        hpexposures = exposures['ring'+str(ns)] == sbdr
        nrepos = (len(hpexposures)//nexposures_per_repo)+1
        rcount = 0
        #for exp in exposures[transnum]:
        for exp in exposures[hpexposures]:
            # Change the root directory name
            #  /net/mss1/blah/blah/
            fluxfile = exp['FLUXFILE'].strip()
            #fluxfile = mssdir+fluxfile.split("/")[-1]
            lo = fluxfile.find('/mss1/')
            fluxfile = mssdir+fluxfile[(lo+6):]
            wtfile = exp['WTFILE'].strip()
            #wtfile = mssdir+wtfile.split("/")[-1]
            lo = wtfile.find('/mss1/')
            wtfile = mssdir+wtfile[(lo+6):]
            maskfile = exp['MASKFILE'].strip()
            #maskfile = mssdir+maskfile.split("/")[-1]
            lo = maskfile.find('/mss1/')
            maskfile = mssdir+maskfile[(lo+6):]
            # Copy flux,wt,maskfiles from /net/mss1/blah to new common repo
            for fl in [fluxfile,wtfile,maskfile]:
                src = fl
                dst = sbdr+"/transfer_"+str(rcount)+"/"+fl.split("/")[-1]
                if os.path.exists(fl):
                    shutil.copy(src,dst)
                else: print(fl," missing!")
            rcount+=1
            rcount = rcount-(rcount//nrepos)*nrepos

    print(time.time()-t00," seconds")
    efile = copydir+explistfile
    print("saving exposure file to ",efile)
    exposures.write(efile)


#    if explistfile!="decam_exposures_tomeas.fits.gz": 
        
            
            
# Do the transfers
#    for sbdr in subdir_list:
#        hpdir = "hgroup_"+str(sbdr)+"/"
#        user = "x25h971" # username for file source machine (dl) ***************************************************
#        host = source+"datalab.noirlab.edu" #"hyalite.rci.montana.edu" # hostname for file source machine (dl)*************************************************************************
#        cmd = "scp -l "+str(bw_limit)+" -r "+user+"@"+host+":"+hpdir+" "+dest_user+"@"+destination+":"+sbdr
#        os.system(cmd)    
