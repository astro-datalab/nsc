#!/usr/bin/env python
#
# UTILS.PY - Utility functions.
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180823'  # yyyymmdd

import re
import logging
#import os
#import sys
#import numpy as np
#import warnings
#from astropy.io import fits
#from astropy.utils.exceptions import AstropyWarning
#from astropy.table import Table, Column
#import time
#import shutil
#import glob
#import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve
#import astropy.stats
#import struct
#import tempfile

# Ignore these warnings, it's a bug
#warnings.filterwarnings("ignore", message="numpy.dtype size changed")
#warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


# Standard grep function that works on string list
def grep(lines,expr,index=False):
    '''
    Similar to the standard unit "grep" but run on a list of strings.
    Returns a list of the matching lines unless index=True is set,
    then it returns the indices.

    Parameters
    ----------
    lines : list
          The list of string lines to check.
    expr  : str
          Scalar string expression to search for.
    index : bool, optional
          If this is ``True`` then the indices of matching lines will be
          returned instead of the actual lines.  index is ``False`` by default.

    Returns
    -------
    out : list
        The list of matching lines or indices.

    Example
    -------

    Search for a string and return the matching lines:

    .. code-block:: python

        mlines = grep(lines,"hello")

    Search for a string and return the indices of the matching lines:

    .. code-block:: python

        index = grep(lines,"hello",index=True)

    '''
    out = []
    cnt = 0L
    for l in lines:
        m = re.search(expr,l)
        if m != None:
            if index is False:
                out.append(l)
            else:
                out.append(cnt)
        cnt = cnt+1
    return out


# Read in all lines of files
def readlines(fil=None):
    '''
    Read in all lines of a file.
    
    Parameters
    ----------
    file : str
         The name of the file to load.
   
    Returns
    -------
    lines : list
          The list of lines from the file

    Example
    -------

    .. code-block:: python

       lines = readlines("file.txt")

    '''
    if fil is None:
        print("File not input")
        return
    f = open(fil,'r')
    lines = f.readlines()
    f.close()
    return lines


# Write all lines to file
def writelines(lines=None,fil=None):
    '''
    Write a list of lines to a file.
    
    Parameters
    ----------
    lines : list
         The list of lines to write to a file.
    fil : str
        The filename to write the lines to.
   
    Returns
    -------
    Nothing is returned.  The lines are written to `fil`.

    Example
    -------

    .. code-block:: python

       writelines(lines,"file.txt")

    '''
    if lines is None:
        print("No lines input")
        return
    if fil is None:
        print("No file name input")
        return
    f = open(fil,'w')
    f.writelines(lines)
    f.close()


# Remove indices from a list
def remove_indices(lst,index):
    '''
    This will remove elements from a list given their indices.

    Parameters
    ----------
    lst : list
          The list from which to remove elements.
    index : list or array
          The list or array of indices to remove.

    Returns
    -------
    newlst : list
           The new list with indices removed.

    Example
    -------

    Remove indices 1 and 5 from array `arr`.

    .. code-block:: python

        index = [1,5]
        arr  = range(10)
        arr2 = remove_indices(arr,index)
        print(arr)
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    '''
    newlst = []
    for i in range(len(lst)):
       if i not in index: newlst.append(lst[i])
    return newlst


# Little function used by numlines
def blocks(files, size=65536):
    '''
    This is a small utility function used by numlines()
    '''
    while True:
        b = files.read(size)
        if not b: break
        yield b


# Read number of lines in a file
def numlines(fil):
    '''
    This function quickly counts the number of lines in a file.

    Parameters
    ----------
    fil : str
          The filename to check the number of lines.

    Returns
    -------
    nlines : int
           The number of lines in `fil`.

    Example
    -------

    .. code-block:: python

        n = numlines("file.txt")

    '''
    with open(fil, "r") as f:
        return (sum(bl.count("\n") for bl in blocks(f)))

    # Could also use this
    #count=0
    #for line in open(fil): count += 1


# Set up basic logging to screen
def basiclogger():
    '''
    This sets up a basic logger that write just to the screen.
    '''
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)-2s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    return logger






