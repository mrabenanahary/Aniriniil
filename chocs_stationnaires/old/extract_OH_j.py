#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 16:22:37 2019

@author: rrabena
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from scipy.optimize import curve_fit
import csv
import os
from astropy.io import ascii
plt.style.use('classic')

#import constants

"""Pour le fit linéaire """
def func(x, a, b):
    return a * x + b


"""Calcul d'une fonction de partition"""
"""Attention: Eu ici est pris en unités K de température"""
def partition(Tex,gu,Eu):
    a = np.multiply(gu,np.exp(np.divide(-Eu,Tex))) #plus de KB car Eu déja en K
    return a,np.sum(a)

GrapheOut = True

mpl.rcdefaults()

Latex=True
font = {'family' : 'serif',
        'size'   : 17,
        'weight' : 'extra bold'}
plt.rc('text', usetex=Latex)
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage[squaren,Gray]{SIunits}',
    r'\usepackage{amsmath}',
    r'\usepackage[version=4]{mhchem}']
plt.rcParams['axes.linewidth']=2


plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2

plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 2

global KB, clight

#in CS units
KB = 1.38064852e-16 #erg/K
clight = 2.99792458e10 #cm/s
hPlanck = 6.62606885e-27 #erg.s


prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
labelTextColor= colors

"""========================================================================"""
"""Modeles de chocs C/J: input"""
originalDir=os.getcwd()
dirFiles = ['j-n4-b0.1','j-n5-b0.1','j-n6-b0.1']
labeldirFiles = {}
vitDirFiles = [[5,10,15,19,20,25,30],[5,10,15,19,20,25,30],[5,10,15,19,20,25,30]]
ageDirFiles = [[1460,390,670,1000,745,780,1050],[440,60,160,200,120,170,160],[40,10,15,15,10,30,30]]

axTableXScale='linear'
axTableYScale='log'
ActiverYang=False
ActiverMyObs=True

moleculeName = ['OI(63microns)','OI(145microns)']
imoleculeName = [18,19]
iTimeN = [1,1]
OutputName = ['OH_int_int.out','OH_int_int_145.out']
molLambdaMicrons = [63.2,145.5]
test = ascii.read('intensity.out',header_start=0,data_start=0)

for l,molecule in enumerate(moleculeName):
    lambdaMicrons=molLambdaMicrons[l]
    originalDir=os.getcwd()
    for i,directory in enumerate(dirFiles):
        os.chdir(directory)
        currentDirLvl1=os.getcwd()
        for j,vitesse in enumerate(vitDirFiles[i]):        
            #for k,tempsJ in enumerate((ageDirFiles[i])[j]):
                vitDirName="v{:1}".format(str(vitesse))
                os.chdir(vitDirName)
                Molfound=False
                iMol=0
                TableMol=ascii.read('intensity.out')
                while(Molfound==False):
                    if (iMol==len(TableMol)):
                        iMol=iMol-1
                        """Conversion de la brillance de surface à une intensité intégrée en K.km/s
                            TdV=InuDnu*1.e-12*lambdaMicrons**3/(2.*KB)/1.e5"""
                        TdV=(TableMol[iMol])[imoleculeName[l]]*1.e-12*lambdaMicrons**3/(2.*KB)/1.e5
                        with open(OutputName[l],'w') as OutputFile:
                            OutputFile.write('name lambda(microns) InuDnu(erg/s/cm2/sr) TdV(K.km.s-1)') #Header
                            OutputFile.write('\n')
                            OutputFile.write(moleculeName[l]+' '+str(molLambdaMicrons[l])+' '+str((TableMol[iMol])[imoleculeName[l]])+' '+str(TdV)+'\n') 
                        Molfound=True                        
                    else:
                        if ((((TableMol[iMol])[iTimeN[l]])>=(ageDirFiles[i])[j])):
                            """Conversion de la brillance de surface à une intensité intégrée en K.km/s
                                TdV=InuDnu*1.e-12*lambdaMicrons**3/(2.*KB)/1.e5"""
                            TdV=(TableMol[iMol])[imoleculeName[l]]*1.e-12*lambdaMicrons**3/(2.*KB)/1.e5
                            with open(OutputName[l],'w') as OutputFile:
                                OutputFile.write('name lambda(microns) InuDnu(erg/s/cm2/sr) TdV(K.km.s-1)') #Header
                                OutputFile.write('\n')
                                OutputFile.write(moleculeName[l]+' '+str(molLambdaMicrons[l])+' '+str((TableMol[iMol])[imoleculeName[l]])+' '+str(TdV)+'\n') 
                            Molfound=True                        
                        else:
                            iMol=iMol+1
                            Molfound=False
                os.chdir(currentDirLvl1)
        os.chdir(originalDir)
