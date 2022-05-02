#!/usr/bin/env python

import os
import re
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from dlnpyutils import utils as dln

def modelmag(tab,instfilt,dec,eqnfile):
    """
    This calculates the model magnitudes for the NSC catalog
    given a catalog with the appropriate information

    Parameters
    ----------
    tab : table
       Catalog of sources with appropriate magnitude
         columns.
    instfilt : str
       Short instrument and filter name, e.g. 'c4d-g'.
    dec : float
       The declination of the exposure.
    eqnfile : str
       File with the model magnitude equations.

    Returns
    -------
    model_mag : numpy array
       An [Nsource,3] array with model magnitudes, errors and color.

    Example
    -------

    model_mag = modelmag(cat,'c4d-g',-50.0,'modelmag_equations.txt')

    By D. Nidever  Feb 2019
    Translated to Python by D. Nidever, April 2022
    """

    # This calculates the model magnitude for stars given the
    # the magnitudes in reference catalogs
    # NUV - Galex NUV magnitude
    # GMAG - Gaia G magnitude
    # JMAG - 2MASS J magnitude
    # KMAG - 2MASS Ks magnitude
    # APASS_GMAG - APASS g magnitue
    # APASS_RMAG - APASS r magnitude
    # EBV  - E(B-V) reddening

    ntab = len(tab)
    for n in tab.colnames:
        tab[n].name = n.upper()
    tabcols = np.char.array(tab.colnames)

    ## Load the model magnitude equation information
    ## band, dec range, color equation, color min/max range, quality cuts, model mag equation
    if os.path.exists(eqnfile)==False:
        raise ValueError(eqnfile+' NOT FOUND')
    eqnstr = Table.read(eqnfile,format='ascii')
    for c in eqnstr.colnames:
        eqnstr[c].name = c.lower()
    neqn = len(eqnstr)
    ## Get COLOR and DEC ranges
    eqnstr['colorlim'] = np.zeros((len(eqnstr),2),float)
    eqnstr['declim'] = np.zeros((len(eqnstr),2),float)
    for i in range(len(eqnstr)):
        # Color range
        cr = np.char.array(eqnstr['colorange'][i].split(','))
        cr = cr.replace('[','').replace(']','')
        eqnstr['colorlim'][i] = np.array(cr).astype(float)
        # DEC range
        dr = np.char.array(eqnstr['decrange'][i].split(','))
        dr = dr.replace('[','').replace(']','')
        eqnstr['declim'][i] = np.array(dr).astype(float)

    ## Get the band for this INSTRUMENT-FILTER and DEC.
    gd, = np.where((np.char.array(eqnstr['instrument'])+'-'+np.char.array(eqnstr['band']) == instfilt) & 
                   (dec >= eqnstr['declim'][:,0]) & (dec <= eqnstr['declim'][:,1]))
    if len(gd)==0:
        raise ValueError('No model magnitude equation for INSTRUMENT-FILTER='+instfilt+' and DEC=%.2f' % dec)
    if len(gd) > 1:
        print('Found multiple magnitude equation for INSTRUMENT-FILTER='+instfilt+' and DEC=%.2f. Using the first one' % dec)
        gd = gd[0]
    eqnstr1 = eqnstr[gd]
    eqnstr1 = dict(zip(eqnstr1.colnames, eqnstr1[0]))  # convert to dictionary

    ## No parentheses allowed
    if eqnstr1['coloreqn'].find('(') != -1 or eqnstr1['coloreqn'].find(')') != -1 or \
       eqnstr1['modelmageqn'].find('(') != -1 or eqnstr1['modelmageqn'].find(')') != -1:
        raise ValueError('No parentheses allowed in the model magnitude equations')
    
    ## Are we using color?
    if eqnstr1['colorlim'][0] < -10 and eqnstr1['colorlim'][1] > 10 and eqnstr1['modelmageqn'].upper().find('COLOR') == -1:
        usecolor = False
    else:
        usecolor = True

    ## Get all columns that we need
    coloreqn = eqnstr1['coloreqn']
    qualitycuts = eqnstr1['qualitycuts']
    modelmageqn = eqnstr1['modelmageqn']
    coloreqn_cols = re.split('[-+*^]',coloreqn)
    modelmageqn_cols = re.split('[-+*^]',modelmageqn)

    if usecolor:
        cols = np.char.array(coloreqn_cols+modelmageqn_cols).upper()
    else:
        cols = np.char.array(modelmageqn_cols).upper()

    ## Remove numbers and "COLOR"
    bd, = np.where(cols.isnumeric() | (cols.upper() == 'COLOR') | (cols == '??'))
    if len(bd)>0:
        if len(bd) < len(cols):
            cols = np.delete(cols,bd)
        else:
            cols = None
    ncols = len(cols)
    ## No columns left
    if ncols==0:
        raise ValueError('No columns to use.')
    ## Only unique columns
    cols = np.unique(cols)
    ncols = len(cols)
    ind1,ind2 = dln.match(tabcols,cols)
    ntagmatch = len(ind1)
    if ntagmatch < ncols:
        leftind = np.arange(ncols)
        if ntagmatch>0:
            leftind = np.delete(leftind,ind2)
        print('Needed columns missing. '+' '.join(cols[leftind]))

    ## Make the color
    ##  replace the columns by TAB[GD].COLUMN
    if usecolor:
        coloreqn_cols = np.char.array(re.split('[-+*^]',coloreqn)).upper()
        coloreqn_cols = np.unique(coloreqn_cols)  # unique ones
        bd, = np.where(coloreqn_cols.isnumeric())  ## Remove numbers
        if len(bd)>0:
            coloreqn_cols = np.delete(coloreqn_cols,bd)
        colcmd = coloreqn.upper()
        for i in range(len(coloreqn_cols)):
            colcmd = colcmd.replace(coloreqn_cols[i],"tab['"+coloreqn_cols[i]+"']")
        color = eval(colcmd)
        color = np.array(color)
    else:
        color = np.zeros(ntab,float)

    ## Make quality cuts
    magcolsind, = np.where((dln.find(cols,'MAG$') > -1) & (cols.find('^E_') == -1) | (cols=='NUV'))
    ## make sure all magnitudes are good (<50) and finite
    goodmask = np.ones(ntab,bool)
    for i in range(len(magcolsind)):
        magind, = np.where(tabcols.upper() == cols[magcolsind[i]].upper())
        if len(magind) == 0:
            print(cols[magcolsind[i]].upper()+' column NOT found')
            return -999999.
        goodmask &= ((tab[tabcols[magind[0]]] < 50) & (tab[tabcols[magind[0]]] > 0) & np.isfinite(tab[tabcols[magind[0]]]))

    ## input quality cuts
    ##  replace <=, <, >, >=, =, &, |
    qualitycuts = qualitycuts.replace('<=',' <= ')
    qualitycuts = qualitycuts.replace('>=',' >= ')
    qualitycuts = qualitycuts.replace('>',' > ')
    qualitycuts = qualitycuts.replace('<',' < ')
    qualitycuts = qualitycuts.replace('=',' == ')
    qualitycuts = qualitycuts.replace('&',' & ')
    qualitycuts = qualitycuts.replace('|',' | ')
    ## fix column names
    qualitycuts_cols = qualitycuts.split()
    for i in range(len(qualitycuts_cols)):
        col = qualitycuts_cols[i]
        colind, = np.where(tabcols==col.upper())
        if len(colind)>0:
            qualitycuts = qualitycuts.replace(tabcols[colind[0]],"tab['"+tabcols[colind[0]]+"']")
    goodmask &= eval(qualitycuts)

    ## Apply the color range
    if usecolor:
        goodmask &= ((color >= eqnstr1['colorlim'][0]) & (color <= eqnstr1['colorlim'][1]))
    ## Get the sources that pass all cuts
    gd, = np.where(goodmask==True)
    if len(gd)==0:
        print('No good sources left')
        return -999999.

    # Make the model magnitude
    ##  replace the columns by TAB[GD].COLUMN
    modelmageqn_cols = np.char.array(re.split('[-+*^]',modelmageqn)).upper()
    bd, = np.where(modelmageqn_cols.isnumeric() | (modelmageqn_cols.upper() == 'COLOR'))  ## Remove numbers and "COLOR"
    if len(bd)>0:
        modelmageqn_cols = np.delete(modelmageqn_cols,bd)
    modelmageqn_cols = np.unique(modelmageqn_cols)  # unique ones
    magcmd = modelmageqn.upper()
    for i in range(len(modelmageqn_cols)):
        magcmd = magcmd.replace(modelmageqn_cols[i],"tab['"+modelmageqn_cols[i]+"'][gd]")
    magcmd = magcmd.replace('COLOR','COLOR[gd]')
    modelmag_gd = eval(magcmd)
    modelmag = np.zeros(ntab,float)+99.99
    modelmag[gd] = modelmag_gd

    ## Make the error structure
    ##  Each magnitude has an E_MAG error except for PS and Gaia GMAG
    ## If we are using PS or GMAG then add the errors for the
    adderrcols = []
    psmagind, = np.where((dln.find(cols,'^PS_') > -1) & (dln.find(cols,'MAG$') > -1))
    if len(psmagind) > 0:
        adderrcols += ['E_'+cols[psmagind]]
    nadderrcols = len(adderrcols)
    ## Making error structure
    errcolind = np.where(dln.find(tabcols,'^E_') > -1)
    errcols = tabcols[errcolind]
    errdt = []
    for i in range(len(errcols)):
        errdt += [(errcols[i],float)]
    if nadderrcols > 0:
        for i in range(nadderrtags):
            errdt += [(adderrcols[i],float)]
    err = np.zeros(ntab,dtype=np.dtype(errdt))
    err = Table(err)
    for c in errcols:
        err[c] = 0.001
    for n in err.colnames:
        err[n] = tab[n]
    ## leave the PS errors at 0.001
    ## convert NAN or 99.99 to 9.99 to be consistent
    for c in err.colnames:
        bd = np.where((err[c] > 10.0) | (np.isfinite(err[c]) == False))
        if len(bd) > 0:
            err[c][bd] = 9.99

    ## Calculate the color errors
    ## get the columns
    if usecolor:
        colorerr_cols = np.char.array(re.split('[-+*^]',coloreqn)).upper()
        colorerr_cols = np.unique(colorerr_cols)  # unique ones
        bd, = np.where(colorerr_cols.isnumeric() | (colorerr_cols.upper() == 'EBV'))  ## Remove numbers and "EBV"
        if len(bd)>0:
            colorerr_cols = np.delete(colorerr_cols,bd)
        ## use - and + signs to break apart the components that need to be squared
        coloreqn_terms = np.char.array(re.split('[-+]',coloreqn)).upper()
        ## remove any terms that don't have a COLORERR_COLS in them
        okay = np.zeros(len(coloreqn_terms),bool)
        for i in range(len(coloreqn_terms)):
            for j in range(len(colorerr_cols)):
                okay[i] |= (coloreqn_terms[i].find(colorerr_cols[j]) > -1)
        bd, = np.where(okay == False)
        if len(bd) > 0:
            coloreqn_terms = np.delete(coloreqn_terms,bd)
        ## Now create the equation, add in quadrature
        colorerrcmd = 'np.sqrt( '+ '+'.join('('+coloreqn_terms.upper()+')**2') +' )'
        #colorerrcmd = colorerrcmd.upper()
        for i in range(len(colorerr_cols)):
            colorerrcmd = colorerrcmd.replace(colorerr_cols[i],"err['E_"+colorerr_cols[i]+"'][gd]")
        colorerr_gd = eval(colorerrcmd)
        colorerr = np.zeros(ntab,float)+9.99
        colorerr[gd] = colorerr_gd
    else:
        colorerr = np.zeros(ntab,float)

    ## The modelmag errors
    ## get the columns
    modelmagerr_cols = np.char.array(re.split('[-+*^]',modelmageqn)).upper()
    modelmagerr_cols = np.unique(modelmagerr_cols)  # unique ones
    bd, = np.where(modelmagerr_cols.isnumeric() | (modelmagerr_cols.upper() == 'EBV'))  ## Remove numbers and "EBV"
    if len(bd)>0:
        modelmagerr_cols = np.delete(modelmagerr_cols,bd)
    ##   use - and + signs to break apart the components that need to be  squared
    modelmageqn_terms = np.char.array(re.split('[-+]',modelmageqn)).upper()
    ## remove any terms that don't have a COLORERR_COLS in them
    okay = np.zeros(len(modelmageqn_terms),bool)
    for i in range(len(modelmageqn_terms)):
        for j in range(len(modelmagerr_cols)):
            okay[i] |= (modelmageqn_terms[i].find(modelmagerr_cols[j]) > -1)
    bd, = np.where(okay == False)
    if len(bd) > 0:
        modelmageqn_terms = np.delete(modelmageqn_terms,bd)
    ## Now create the equation, add in quadrature 
    modelmagerrcmd = 'np.sqrt( '+ '+'.join('('+modelmageqn_terms+')**2') +' )'
    for i in range(len(modelmageqn_cols)):
        modelmagerrcmd = modelmagerrcmd.replace(modelmageqn_cols[i],"err['E_"+modelmageqn_cols[i].upper()+"'][gd]")
    modelmagerrcmd = modelmagerrcmd.replace('COLOR','COLORERR[gd]')
    modelmagerr_gd = eval(modelmagerrcmd)
    modelmagerr = np.zeros(ntab,float)+99.90
    modelmagerr[gd] = modelmagerr_gd

    ## combine modelmag and modelmagerr
    ## Combine mags and errors
    mags = np.zeros((len(tab),3),float)
    mags[:,0] = modelmag
    mags[:,1] = modelmagerr
    mags[:,2] = color

    return mags
