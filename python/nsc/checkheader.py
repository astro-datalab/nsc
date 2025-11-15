import os
import numpy as np
from glob import glob
from astropy.table import Table,vstack
from . import utils

def checkheader(expdir):

    dldir,mssdir,localdir = utils.rootdirs()
    instrument = 'c4d'

    base = os.path.basename(expdir)
    night = expdir.split('/')[-2]

    # v4+ use separate header file
    version = 'v4'
    # _header.fits
    headfile = os.path.join(expdir,base+'_header.fits')
    if os.path.exists(headfile):
        return True,headfile,1
    # .hdr
    headfile = os.path.join(expdir,base+'.hdr')
    if os.path.exists(headfile):
        return True,headfile,2
    # header directory .hdr
    headfile = os.path.join(dldir,'instcal',version,
                            'header',instrument,night,base+'.hdr')
    if os.path.exists(headfile):
        return True,headfile,3
    # Instcal files on tacc
    #/home1/09970/dnidever/scratch1/nsc/instcal/v4/images/c4d/2020/20200130
    headfile = '/home1/09970/dnidever/scratch1/nsc/instcal/v4/images/'
    headfile += '/'.join(expdir.split('/')[-4:])+'.fits.fz'
    if os.path.exists(headfile):
        return True,headfile,4
    # different version .hdr file
    # sometimes there's a different version, i.e. _d2 instead of _ls11 
    base2 = '_'.join(base.split('_')[:-1])
    headfile = glob(os.path.join(dldir,'instcal',version,
                                 'header',instrument,night,base2+'_*.hdr'))
    if len(headfile)>0:
        headfile = headfile[0]
        return True,headfile,5
    else:
        headfile = ''
    if os.path.exists(headfile)==False:
        print(headfile+' not found')

    return False,'',-1

def checkall(filename):
    tab = Table.read(filename)
    tab['headexists'] = False
    tab['headfile'] = np.zeros(len(tab),dtype=(str,200))
    tab['headfiletype'] = -1
    for i in range(len(tab)):
        out = checkheader(tab['expdir'][i])
        print(i+1,out[0],out[2])
        tab['headexists'][i] = out[0]
        if out[0]:
            tab['headfile'][i] = out[1]
            tab['headfiletype'][i] = out[2]
    return tab
