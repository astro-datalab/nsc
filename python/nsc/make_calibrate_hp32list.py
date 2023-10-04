#!/usr/bin/env python 


#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from argparse import ArgumentParser
from astropy.table import Table,Column
from astropy.io import fits
from astropy.wcs import WCS
from dlnpyutils import utils as dln, coords
import healpy as hp
import logging
import numpy as np
import os
import socket
import subprocess
import sys
import time


#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Make list of healpix ')
    parser.add_argument('--version', type=str, nargs=1, help='Version number')
    parser.add_argument('--nside', type=str, nargs=1, help='HEALPix NSIDE')
    parser.add_argument('--partition',type=str,nargs=1,help='Delimited list of partitions to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposures that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=1, help='Max number of exposures to process at any given time')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list healpix to calibrate')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)     # NSC version
    nside = dln.first_el(args.nside)         # HEALPix NSIDE
    redo = args.redo                         # if called, redo = True
    partitions=args.partition[0].split(',')  # the slurm partitions to submit jobs to
    npar=len(partitions)                     # number of slurm partitions
    maxjobs = int(args.maxjobs[0])           # maximum number of jobs to maintain running at any time
    nchan = maxjobs//npar                    # number of job channels per partition cpar -> nchan
    sleep_time=10                            # seconds to sleep between checking job status
    inputlist = args.list                    # list of exposures to analyze
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories 
    # -tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"   # location of operations
    expdir = basedir+"exposures/"                        # where exposures will be
    outfiledir = basedir+"outfiles/"                     # a place for the job files
    makedir(expdir)
    makedir(outfiledir)

    # Load the list of exposures
    list1 = fits.open(os.path.join(basedir, 'lists/decam_instcal_list.fits.gz'))[1].data
    list2 = fits.open(os.path.join(basedir, 'lists/mosaic3_instcal_list.fits.gz'))[1].data
    list3 = fits.open(os.path.join(basedir, 'lists/bok90prime_instcal_list.fits.gz'))[1].data
    str_data = [list1, list2, list3]
    del list1, list2, list3
    # Setup exposure data structure
    nstr = len(str_data)
    str_data['filter'] = str_data['filter'].strip()
    str_data['fluxfile'] = str_data['fluxfile'].strip()
    str_data['maskfile'] = str_data['maskfile'].strip()
    str_data['wtfile'] = str_data['wtfile'].strip()
    print(f'{nstr} InstCal images')


    # Get good RA/DEC
    bd = [i for i, (ra, dec) in enumerate(zip(str_data['ra'], str_data['dec'])) if not (isfinite(ra) and isfinite(dec))]
    print(f'Getting RA/DEC for {len(bd)} exposures')
    # loop through 'bad' exposures
    for i in range(len(bd)):
        if i % 50 == 0: print(i)
        file = str_data['fluxfile'][bd[i]].strip()
        if file.startswith('/net'):file = file[4:]
        # get RA,DEC
        if os.path.exists(file):
            try:
                with fits.open(file) as hdul:
                    wcs = WCS(hdul[1].header)
                    ra, dec = wcs.all_pix2world(1, 1, 1)
                    str_data['ra'][bd[i]] = ra
                    str_data['dec'][bd[i]] = dec
            except Exception as e:
                pass

    # Start output file
    list_data = np.empty(nstr, dtype=[('expdir', 'U100'), ('pix', np.int64), ('instrument', 'U10'), ('filter', 'U10')])
    list_data['instrument'] = np.char.strip(str_data['instrument'])
    list_data['filter'] = np.char.strip(np.char.upper(np.char.strip(str_data['filter'], 2)))
    # rename filter "bok-r"" to "r"
    g = np.where(np.char.upper(np.char.strip(str_data['filter'], 4)) == 'BOKR')[0]
    if len(g) > 0:
        list_data['filter'][g] = 'R'
    # rename filter "k4m-zd" to "z"
    g = np.where((str_data['instrument'] == 'k4m') & (np.char.upper(np.char.strip(str_data['filter'], 2)) == 'ZD'))[0]
    if len(g) > 0:
        list_data['filter'][g] = 'Z'

    # Calculate the output directory
    for i in range(nstr):
        if i % 1000 == 0: print(i)
        dateobs = str_data['date_obs'][i].strip()
        night = dateobs[:4] + dateobs[5:7] + dateobs[8:10]
        baseroot = os.path.splitext(os.path.basename(str_data['fluxfile'][i].strip()))[0]
        expdir = os.path.join(dir, str_data['instrument'][i].strip(), night, baseroot, '')
        list_data['expdir'][i] = expdir

    # Calculate the healpix coordinates
    theta = (90 - str_data['dec']) / radeg
    phi = str_data['ra'] / radeg
    ipring = hp.ang2pix(nside, theta, phi, lonlat=True, nest=False)
    list_data['pix'] = ipring

    # Save the list
    output_fits_file = os.path.join(dir, 'lists/nsc_calibrate_healpix_list.fits')
    print(f'Writing list to {output_fits_file}')
    fits.writeto(output_fits_file, list_data, overwrite=True)