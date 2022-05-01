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
       Structure of source with appropriate magnitude
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

    cols = np.char.array(tab.colnames)

    ## Load the model magnitude equation information
    ## band, dec range, color equation, color min/max range, quality cuts, model mag equation
    if os.path.exists(eqnfile)==False:
        raise ValueError(eqnfile+' NOT FOUND')
    eqnstr = Table.read(eqnfile,format='ascii')
    #eqnstr = IMPORTASCII(eqnfile,/header,/silent)
    neqn = len(eqnstr)
    ## Get COLOR range
    dum = strsplitter(eqnstr['colorange'],',',/extract)
    dum = dum.replace('[','').replace(']','')
    add_tag,eqnstr,'colorlim',[0.0,0.0],eqnstr
    eqnstr['colorlim'][0] = float(reform(dum[0,:]))
    eqnstr['colorlim'][1] = float(reform(dum[1,:]))
    ## Get DEC range
    dum = strsplitter(eqnstr['decrange'],',',/extract)
    dum = dum.replace('[','').replace(']','')
    add_tag,eqnstr,'declim',[0.0,0.0],eqnstr
    eqnstr['declim'][0] = float(reform(dum[0,:]))
    eqnstr['declim'][1] = float(reform(dum[1,:]))


    ## Get the band for this INSTRUMENT-FILTER and DEC.
    gd, = np.where((str(eqnstr['instrument'])+'-'+str(eqnstr['band']) == instfilt) & (dec >= eqnstr['declim'][0]) & (dec <= eqnstr['declim'][1]))
    if len(gd)==0:
        raise ValueError('No model magnitude equation for INSTRUMENT-FILTER='+instfilt+' and DEC='+stringize(dec,ndec=2))
    if len(gd) > 1:
        print('Found multiple magnitude equation for INSTRUMENT-FILTER='+instfilt+' and DEC=%.2f. Using the first one' % dec)
        gd = gd[0]
    eqnstr1 = eqnstr[gd]

    ## No parentheses allowed
    if eqnstr1['coloreqn'].find('(') != -1 or eqnstr1['coloreqn'].find(')') != -1 or
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
    #coloreqn_cols = strsplit(coloreqn,'-+*^',/extract)
    #modelmageqn_cols = strsplit(modelmageqn,'-+*^',/extract)    
    if usecolor:
        cols = np.append(coloreqn_cols, modelmageqn_cols).upper()
    else:
        cols = modelmageqn_cols.upper()
    ## Remove numbers and "COLOR"
    bd, = np.where((valid_num(cols) == 1) | (cols.upper() == 'COLOR') | (cols == '??'))
    if len(bd)>0:
        if len(bd) < len(cols):
            cols = np.delete(cols,bd)
            #REMOVE,bd,cols
        else:
            cols = None
    ncols = len(cols)
    ## No columns left
    if ncols==0:
        print('No columns to use.')
        return -999999.
    ## Only unique columns
    ucols,ui = np.unique(cols,return_index=True)
    cols = cols[ui]
    ncols = len(cols)
    ind1,ind2 = dln.match(cols,cols)
    ntagmatch = len(ind1)
    #MATCH,tags,cols,ind1,ind2,/sort,count=ntagmatch  
    if ntagmatch < ncols:
        leftind = np.arange(ncols)
        if ntagmatch>0:
            leftind = np.delete(leftind,ind2)
            #remove,ind2,leftind
        print('Needed columns missing. '+' '.join(cols[leftind]))

    ## Make the color
    ##  replace the columns by TAB[GD].COLUMN
    if usecolor:
        coloreqn_cols = re.split('[-+*^]',coloreqn).upper()
        #coloreqn_cols = strsplit(coloreqn,'-+*^',/extract).upper()        
        coloreqn_cols = coloreqn_cols[uniq(coloreqn_cols,sort(coloreqn_cols))]  # unique ones
        coloreqn_cols = np.char.array(coloreqn_cols)
        bd, = np.where(valid_num(coloreqn_cols) == 1)  ## Remove numbers
        if len(bd)>0:
            coloreqn_cols = np.delete(coloreqn_cols,bd)
            #REMOVE,bd,coloreqn_cols
        colcmd = [c.upper() for c in coloreqn]
        for i in range(len(coloreqn_cols)):
            colcmd = colcmd.replace(coloreqn_cols[i],"tab['"+coloreqn_cols[i]+"']")
        dum = exec('color='+colcmd)
    else:
        color = np.zeros(ntab,float)

    ## Make quality cuts
    magcolsind, = np.where((cols.find('mag$') > -1) & (cols.find('^e_') == -1) | (cols=='NUV'))
    ## make sure all magnitudes are good (<50) and finite
    goodmask = np.ones(ntab,bool)
    for i in range(len(magcolsind)):
        magind = np.where(cols.upper() == cols[magcolsind[i]].upper())
        if len(magind) == 0:
            print(cols[magcolsind[i]].upper()+' column NOT found')
            return -999999.
        goodmask &= ((tab[cols[magind[0]]] < 50) & (tab[cols[magind[0]]] > 0) & np.isfinite(tab[cols[magind[0]]]))

    ## input quality cuts
    ##  replace <=, <, >, >=, =, &, |
    qualitycuts = qualitycuts.replace('<=',' le ')
    qualitycuts = qualitycuts.replace('>=',' ge ')
    qualitycuts = qualitycuts.replace('>',' gt ')
    qualitycuts = qualitycuts.replace('<',' lt ')
    qualitycuts = qualitycuts.replace('=',' eq ')
    qualitycuts = qualitycuts.replace(,'&',' and ')
    qualitycuts = qualitycuts.replace(,'|',' or ')
    ## fix columns names
    qualitycuts_cols = qualitycuts.split(' ')
    for i in range(len(qualitycuts_cols)):
        col = qualitycuts_cols[i]
        colind = np.where(cols == col.upper())
        if len(colind)>0:
            qualitycuts = qualitycuts.replace(col,'tab.'+col)
    dum = exec('goodmask &= '+qualitycuts)
    ## Apply the color range
    if usecolor:
        goodmask &= (color >= eqnstr1['colorlim'][0] and color <= eqnstr1['colorlim'][1])
    ## Get the sources that pass all cuts
    gd, = np.where(goodmask==1)
    if len(gd)==0:
        print('No good sources left')
        return -999999.

    # Make the model magnitude
    ##  replace the columns by TAB[GD].COLUMN
    modelmageqn_cols = re.split('[-+*^]',modelmageqn).upper()
    #modelmageqn_cols = strsplit(modelmageqn,'-+*^',/extract).upper()    
    bd, = np.where((valid_num(modelmageqn_cols) == 1) | (modelmageqn_cols.upper() == 'COLOR'))  ## Remove numbers and "COLOR"
    if len(bd)>0:
        modelmageqn_cols = np.delete(modelmageqn_cols,bd)
        #REMOVE,bd,modelmageqn_cols
    modelmageqn_cols = modelmageqn_cols[uniq(modelmageqn_cols,sort(modelmageqn_cols))]  # unique ones
    magcmd = strupcase(modelmageqn)
    for i in range(len(modelmageqn_cols)):
        magcmd = magcmd.replace(modelmageqn_cols[i],"tab['"+modelmageqn_cols[i]+"'][gd]")
    magcmd = magcmd.replace('COLOR','COLOR[gd]')
    dum = exec('modelmag_gd='+magcmd)
    modelmag = np.zeros(ntab,float)+99.99
    modelmag[gd] = modelmag_gd

    ## Make the error structure
    ##  Each magnitude has an E_MAG error except for PS and Gaia GMAG
    ## If we are using PS or GMAG then add the errors for the
    adderrtags = []
    #if (np.where(cols == 'GMAG'))[0] ne -1 then push,adderrtags,'E_GMAG'
    psmagind, = np.where((cols.find('^PS_') > -1) & (cols.find('MAG$') > -1)
    if len(psmagind) > 0:
        adderrtags += ['E_'+cols[psmagind]]
    nadderrtags = len(adderrtags)
    ## Making error structure
    errtagind = np.where(cols.find('^E_') > -1)
    errtags = cols[errtagind]
    errschema = create_struct(errtags[0],0.001)
    if nerrtags > 1:
        for i in np.arange(1,nerrtags):
            errschema = create_struct(errschema,errtags[i],0.001)
    if nadderrtags > 0:
        for i in range(nadderrtags):
            errschema = create_struct(errschema,adderrtags[i],0.001)
    err = np.zeros(ntab,dtype=np.dtype(edt))
    #struct_assign,tab,err,/nozero
    for n in tab.colnames:
        err[n] = tab[n]
    #if (np.where(cols == 'GMAG'))[0] ne -1 then err.e_gmag = 2.5*alog10(1.0+tab.e_fg/tab.fg)
    ## leave the PS errors at 0.001
    ## convert NAN or 99.99 to 9.99 to be consistent
    for i=0,n_tags(err)-1:
        bd = np.where((err.(i) > 10.0) | (np.isfinite(err.(i)) == False))
        if len(bd) > 0:
            err[bd].(i) = 9.99

    ## Calculate the color errors
    ## get the columns
    if usecolor:
        colorerr_cols = re.split('[-+*^]',coloreqn).upper()
        #colorerr_cols = strsplit(coloreqn,'-+*^',/extract).upper()                         
        colorerr_cols = colorerr_cols[uniq(colorerr_cols,sort(colorerr_cols))]  # unique ones
        bd, = np.where((valid_num(colorerr_cols) == 1) | (colorerr_cols.upper() == 'EBV'))  ## Remove numbers and "EBV"
        if len(bd)>0:
            colorerr_cols = np.delete(colorerr_cols,bd)
            #REMOVE,bd,colorerr_cols
        ## use - and + signs to break apart the components that need to be  squared
        coloreqn_terms = re.split('[-+]',coloreqn).upper()
        #coloreqn_terms = strupcase(strsplit(coloreqn,'-+',/extract))                         
        ## remove any terms that don't have a COLORERR_COLS in them
        okay = np.zeros(len(coloreqn_terms),bool)
        for i in range(len(coloreqn_terms)):
            for j in range(len(colorerr_cols)):
                okay[i] |= (coloreqn_terms[i].find(colorerr_cols[j]) < -1)
        bd, = np.where(okay == False)
        if len(bd) > 0:
            coloreqn_terms = np.delete(coloreqn_terms,bd)
            #REMOVE,bd,coloreqn_terms
        ## Now create the equation, add in quadrature
        colorerrcmd = 'np.sqrt( '+ '+'.join('('+coloreqn_terms+')**2') +' )'
        colorerrcmd = strupcase(colorerrcmd)
        for i in range(len(colorerr_cols)):
            colorerrcmd = colorerrcmd.replace(colorerr_cols[i],"err['e_"+colorerr_cols[i]+"'][gd]")
        dum = exec('colorerr_gd='+colorerrcmd)
        colorerr = np.zeros(ntab,float)+9.99
        colorerr[gd] = colorerr_gd
    else:
        colorerr = np.zeros(ntab,float)

    ## The modelmag errors
    ## get the columns
    modelmagerr_cols = re.split('[-+*^]',modelmageqn).upper()
    #modelmagerr_cols = strupcase(strsplit(modelmageqn,'-+*^',/extract))                         
    modelmagerr_cols = modelmagerr_cols[uniq(modelmagerr_cols,sort(modelmagerr_cols))]  # unique ones
    bd, = np.where((valid_num(modelmagerr_cols) == 1) | (modelmagerr_cols.upper() == 'EBV'))  ## Remove numbers and "EBV"
    if len(bd)>0:
        modelmagerr_cols = np.delete(modelmagerr_cols,bd)
        #REMOVE,bd,modelmagerr_cols
    ##   use - and + signs to break apart the components that need to be  squared
    modelmageqn_terms = re.split('[-+]',modelmageqn).upper()
    #modelmageqn_terms = strupcase(strsplit(modelmageqn,'-+',/extract))                         
    ## remove any terms that don't have a COLORERR_COLS in them
    okay = np.zeros(len(modelmageqn_terms),bool)
    for i in range(len(modelmageqn_terms)):
        for j in range(len(modelmagerr_cols)):
            okay[i] |= (modelmageqn_terms[i].find(modelmagerr_cols[j]) > -1)
    bd, = np.where(okay == False)
    if len(bd) > 0:
        modelmageqn_terms = np.delete(modelmageqn_terms,bd)
        #REMOVE,bd,modelmageqn_terms
    ## Now create the equation, add in quadrature 
    modelmagerrcmd = 'np.sqrt( '+ '+'.join('('+modelmageqn_terms+')**2') +' )'
    modelmagerrcmd = modelmagerrcmd.upper()
    for i in range(len(modelmageqn_cols)):
        modelmagerrcmd = modelmagerrcmd.replace(modelmageqn_cols[i],"err['e_"+modelmageqn_cols[i]+"'][gd]")
    modelmagerrcmd = modelmagerrcmd.replace('COLOR','COLORERR[gd]')
    dum = exec('modelmagerr_gd='+modelmagerrcmd)
    modelmagerr = np.zeros(ntab,float)+99.90
    modelmagerr[gd] = modelmagerr_gd

    ## combine modelmag and modelmagerr
    ## Combine mags and errors
    mags = np.zeros((len(tab),3),float)
    mags[:,0] = modelmag
    mags[:,1] = modelmagerr
    mags[:,2] = color

    return mags
