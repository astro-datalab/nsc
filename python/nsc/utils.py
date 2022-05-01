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
