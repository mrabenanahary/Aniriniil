#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 14:09:36 2020

@author: mialy94
=================================================================
= Auteur : Mialy RABENANAHARY
= 
= globals.py
= 
= Module contenant toutes les variables globales utilisées dans les scripts
= qui font appel à ce module
= 
=
=================================================================
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from scipy.optimize import curve_fit
import csv
import os
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
import matplotlib.text as mtext
#from durhampy.data.globals import *

def tangent_slope(y1,y2,x1,x2):
    return (y2-y1)/(x2-x1)

def list_tangent_slope(xarray,yarray):
    for chaine in [xarray,yarray]:
        test = True
        try:
            test = (len(np.shape(chaine))==1)
            assert test
        except:
            raise AssertionError('Autre chose qu une liste 1D a ete inseree')
        test = True
        try : 
                if(np.shape(chaine)[0]<1):
                    test = False
                assert test
        except:
                raise AssertionError('Cette liste est vide')        
                
    try:
        assert len(xarray)==len(yarray)
    except:
        raise AssertionError('xarray et yarray n\'ont pas la même taille')
    
    outi = []
    for i,y_i in enumerate(yarray):
        outj = []
        for j,y_j in enumerate(yarray):
            if(i!=j) : outj.append((y_i-y_j)/(xarray[i]-xarray[j]))
        outi.append(outj)
    return outi
                