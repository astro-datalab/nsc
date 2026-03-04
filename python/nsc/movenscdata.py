import os
import numpy as np
from glob import glob
import shutil
from astropy.table import Table
from dlnpyutils import utils as dln


def movenscdata():
    """ Move NSC measurement data that combine has already used to a special tomove/ directory """

    rootdir = '/home/group/davidnidever/nsc/instcal/v4/'
    combdir = rootdir+'combine/'
    
    #print('Loading large chips/healpix table')
    #tab = Table.read(rootdir+'lists/nsc_calibrate_summary_chips_final.fits')
    #
    ## unique exposures-healpix combinations
    #exphpix = np.char.array(tab['base'].astype(str))+'-'+np.char.array(tab['pix128'].astype(str))
    #_,ui = np.unique(exphpix,return_index=True)
    ## len(ui)
    ## 9477840
    #tab2 = tab[ui]
    #tab2.write(rootdir+'lists/nsc_calibrate_summary_unique_exposure_healpix.fits')

    print('Loading large exposure/healpix table')
    tab = Table.read(rootdir+'lists/nsc_calibrate_summary_unique_exposure_healpix.fits')

    # Unique healpix
    upix = np.unique(tab['pix128'])
    print(len(upix),'unique healpix')
    hpixindex = dln.create_index(tab['pix128'])  # healpix index
    
    # Unique exposures
    ubase = np.unique(tab['base'])
    print(len(ubase),'unique exposures')

    # Create exposure index and exposure table
    expindex = dln.create_index(tab['base'])
    nexp = len(expindex['value'])
    eind = np.zeros(nexp,int)
    for i in range(nexp):
        ind = expindex['index'][expindex['lo'][i]:expindex['hi'][i]+1]
        nind = len(ind)
        eind[i] = ind[0]
    exptab = tab[eind]

    
    # Figure out which healpix are done
    print('Figuring out which healpix are done')
    exists = np.zeros(len(upix),bool)
    for i in range(len(upix)):
        objfile = os.path.join(combdir,str(upix[i]//1000),str(upix[i])+'.fits.gz')
        exists[i] = os.path.exists(objfile)
    done, = np.where(exists)
    upixdone = upix[done]
    print(len(upixdone),'healpix are done')

    # Healpix left to do
    notdone, = np.where(exists==False)
    upixleft = upix[notdone]
    print(len(upixleft),'healpix left')

    # Get all of the exposures for healpix that are DONE
    #  then make sure that they are not needed for any other healpix
    _,_,ind2 = np.intersect1d(upixdone,hpixindex['value'],return_indices=True)
    ind = []
    for i in range(len(ind2)):
        j = ind2[i]
        ind1 = hpixindex['index'][hpixindex['lo'][j]:hpixindex['hi'][j]+1]
        ind.append(ind1)
    ind = np.concatenate(ind)
    hdonebase = np.unique(tab['base'][ind])
    print(len(hdonebase),'exposures for healpix that are DONE')

    # Get all of the exposures for healpix that are LEFT
    _,_,ind2 = np.intersect1d(upixleft,hpixindex['value'],return_indices=True)
    leftind = []
    for i in range(len(ind2)):
        j = ind2[i]
        ind1 = hpixindex['index'][hpixindex['lo'][j]:hpixindex['hi'][j]+1]
        leftind.append(ind1)
    leftind = np.concatenate(leftind)
    hleftbase = np.unique(tab['base'][leftind])
    print(len(hleftbase),'exposures for healpix that are LEFT')
    
    # The ones that are NOT needed can be moved
    # Exposures that are needed for the healpix that are left
    _,ind1,ind2 = np.intersect1d(hleftbase,exptab['base'],return_indices=True)
    basetomove = exptab['base'].copy()
    basetomove = np.delete(basetomove,ind2)
    print(len(basetomove),'exposures are not needed anymore')
    
    _,ind1,ind2 = np.intersect1d(basetomove,exptab['base'],return_indices=True)
    tomoveexptab = exptab[ind2]
    dirstomove = ['/home/group/davidnidever'+os.path.dirname(f[f.find('/nsc'):]) for f in tomoveexptab['measfile']]
    
    # NOTE, some directories were ALREADY moved

    # Move the directories to "tomove/"
    moved = []
    for i in range(len(dirstomove)):
        indir = dirstomove[i]
        if os.path.exists(indir)==False: continue
        outdir = indir.replace('/c4d/','/c4d/tomove/')
        basedir = os.path.dirname(outdir)
        os.makedirs(basedir,exist_ok=True)
        shutil.move(indir,outdir)
        moved.append(indir)

    print(len(moved),'exposure directories moved')
        
    import pdb; pdb.set_trace()
