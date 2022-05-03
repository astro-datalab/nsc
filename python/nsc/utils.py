#!/usr/bin/env python

import os
import time
import socket
import numpy as np
from astropy.io import fits
from astropy.table import Table
import subprocess
from dlnpyutils import utils as dln

def rootdirs():
    # Return the NSC root directories for various machines

    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    if host.find('thing') > -1 or host.find('hulk') > -1:
        #dldir = '/dl1/users/'
        dldir = '/net/dl2/'
        mssdir = '/mss1/'
        localdir = '/d0/'
    elif host.find('gp09') > -1 or host.find('gp08') > -1 or host.find('gp07') > -1 or \
         host.find('gp06') > -1 or host.find('gp05') > -1:
        #dldir = '/net/dl1/users/'
        dldir = '/net/dl2/'
        mssdir = '/net/mss1/'
        localdir = '/data0/'
    else:
        raise ValueError(host+' UNKNOWN')

    return dldir,mssdir,localdir

def func_poly2d_wrap(x,*args):
    """ thin wrapper for curve_fit"""
    xx = x[0]
    yy = x[1]
    return func_poly2d(xx,yy,*args)

def func_poly2d(x,y,*args):
    """ 2D polynomial surface"""

    p = args
    np = len(p)
    if np==0:
        a = p[0]
    elif np==3:
        a = p[0] + p[1]*x + p[2]*y
    elif np==4:
        a = p[0] + p[1]*x + p[2]*x*y + p[3]*y
    elif np==6:
        a = p[0] + p[1]*x + p[2]*x**2 + p[3]*x*y + p[4]*y + p[5]*y**2
    elif np==8:
        a = p[0] + p[1]*x + p[2]*x**2 + p[3]*x*y + p[4]*(x**2)*y + p[5]*x*y**2 + p[6]*y + p[7]*y**2
    elif np==11:
        a = p[0] + p[1]*x + p[2]*x**2.0 + p[3]*x**3.0 + p[4]*x*y + p[5]*(x**2.0)*y + \
            p[6]*x*y**2.0 + p[7]*(x**2.0)*(y**2.0) + p[8]*y + p[9]*y**2.0 + p[10]*y**3.0
    elif np==15:
        a = p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4 + p[5]*y + p[6]*x*y + \
            p[7]*(x**2)*y + p[8]*(x**3)*y + p[9]*y**2 + p[10]*x*y**2 + p[11]*(x**2)*y**2 + \
            p[12]*y**3 + p[13]*x*y**3 + p[14]*y**4
    elif np==21:
        a = p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4 + p[5]*x**5 + p[6]*y + p[7]*x*y + \
            p[8]*(x**2)*y + p[9]*(x**3)*y + p[10]*(x**4)*y + p[11]*y**2 + p[12]*x*y**2 + \
            p[13]*(x**2)*y**2 + p[14]*(x**3)*y**2 + p[15]*y**3 + p[16]*x*y**3 + p[17]*(x**2)*y**3 + \
            p[18]*y**4 + p[19]*x*y**4 + p[20]*y**5
    else:
        raise Exception('Only 3, 4, 6, 8, 11 amd 15 parameters supported')

    return a
