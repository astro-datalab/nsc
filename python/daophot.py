#!/usr/bin/env python

import os
import sys
import numpy as np
import scipy
from astropy.io import fits
import photutils
from skimage import measure, morphology
from scipy.cluster import vq
#import gaps
import matplotlib.pyplot as plt
import pylab
from scipy.signal import argrelmin
import scipy.ndimage.filters as filters
import time
from dlnpyutils import utils as dln

def rdpsf(psffile):
    """ Load a DAOPHOT .psf file"""
    # Check if the file exists
    if os.path.exists(psffile) is False:
        raise ValueError(psffile+" NOT FOUND")

    # Check RDPSF in mathsubs.f

    #      INTEGER FUNCTION  RDPSF  (PSFFIL, IPSTYP, PAR, MAXPAR, NPAR,
    #     .     PSF, MAXPSF, MAXEXP, NPSF, NEXP, NFRAC, 
    #     .     PSFMAG, BRIGHT, XPSF, YPSF)
    #
    # Read in the point-spread function
    #
    #      IMPLICIT NONE
    #      INTEGER MAXPSF, MAXTYP, MAXPAR, MAXEXP
    #      PARAMETER (MAXTYP=7)
    #
    #      REAL PAR(MAXPAR), PSF(MAXPSF,MAXPSF,MAXEXP)
    #
    #      CHARACTER*30 PSFFIL
    #      CHARACTER*8 LABEL, CHECK
    #      REAL PSFMAG, BRIGHT, XPSF, YPSF
    #      INTEGER I, J, K, IPSTYP, NPSF, NPAR, NEXP, NFRAC, ISTAT
    #      INTEGER NTERM, NPARAM
    #
    #      CALL INFILE (3, PSFFIL, ISTAT)
    #      IF (ISTAT .NE. 0) THEN
    #         RDPSF = -1
    #         RETURN
    #      END IF

    # Read header line
    #      READ (3,302,IOSTAT=ISTAT) LABEL, NPSF, NPAR, NEXP, NFRAC, PSFMAG, 
    #     .     BRIGHT, XPSF, YPSF
    #  302 FORMAT (1X, A8, 4I5, F9.3, F15.3, 2F9.1)    
    lines = dln.readlines(psffile)
    nlines = len(lines)
    line1 = lines[0]
    fmt = '(1X, A8, 4I5, F9.3, F15.3, 2F9.1)'
    label,npsf,npar,nexp,nfrac,psfmag,bright,xpsf,ypsf = dln.fread(line1,fmt)
    label = label.strip()
    header = {'label':label, 'npsf':npsf, 'npar':npar, 'nexp':nexp, 'nfrac':nfrac,
              'psfmag':psfmag, 'bright':bright, 'xpsf':xpsf, 'ypsf':ypsf}
    
    # Checking something here, maybe that the parameters are okay
    #      DO IPSTYP=1,MAXTYP
    #         I = NPARAM(IPSTYP, 1., CHECK, PAR, MAXPAR)
    #         IF ((LABEL .EQ. CHECK) .AND. (I .EQ. NPAR)) GO TO 1100
    #      END DO
    #      CALL STUPID ('Inappropriate PSF: '//LABEL)


    # Read in the parameters
    # 1100 READ (3,301,IOSTAT=ISTAT) (PAR(I), I=1,NPAR)
    #  301 FORMAT (1X, 6E13.6)
    line1 = lines[1]
    par = np.zeros(npar,float)
    for i in range(npar):
        par[i] = float(line1[1+i*13:1+(i+1)*13])
    
    # Read in the data
    #      NTERM = NEXP+NFRAC
    #      IF (NTERM .GE. 1) THEN
    #         DO K=1,NTERM
    #            READ (3,311,IOSTAT=ISTAT) ((PSF(I,J,K), I=1,NPSF), J=1,NPSF)
    #  311       FORMAT (1X, 6E13.6)
    #         END DO
    #      END IF
    nterm = nexp + nfrac
    header['nterm'] = nterm
    if nterm>=1:
        # Put all of the data into one long array
        bigline = ''
        for i in range(nlines-2):
            bigline += lines[i+2][1:]
        nnum = int(len(bigline)/13)
        numbers = np.zeros(nnum,float)
        for i in range(nnum):
            numbers[i] = float(bigline[i*13:(i+1)*13])
        # Now put in 3D array
        psf = np.zeros((npsf,npsf,nterm),float)
        n2 = npsf*npsf
        for k in range(nterm):
            numbers1 = numbers[k*n2:(k+1)*n2]
            psf[:,:,k] = numbers1.reshape(npsf,npsf)

        # plt.imshow(psf[:,:,0],origin='lower',aspect='auto')   

    # header, par, psf
    return (header, par, psf)
        
    #return (ipstyp, par, maxpar, npar, psf, maxpsf, maxexp, npsf, nexp, nfrac, psfmag, bright, xpsf, ypsf)
    
    
class PSF:
    """ DAOPHOT PSF class."""
    
    def __init__(self,header,par,psf):
        # Initalize the psf object
        self.header = header
        self.par = par
        self.psf = psf

    def call(self,x,y,mag,full=False,deriv=False,origin=0):
        """ Create a PSF image."""

        # CHECK ADDSTAR.F

        # have option to return the derivative
        
    def __str__(self):
        pass

    @classmethod
    def read(self,filename):
        """ Read in a PSF from a .psf file."""
        # Check if the file exists
        if os.path.exists(filename) is False:
            raise ValueError(filename+" NOT FOUND")
        # Load the file
        header, par, psf = rdpsf(filename)
        # Initalize the psf object
        return PSF(header, par, psf)
        
    def write(self,outfile):
        """ Write the PSF to a file."""
        pass


