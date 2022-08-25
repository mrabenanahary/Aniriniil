# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
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
        'size'   : 34,
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
dirFiles = ['c-n4-b1','c-n5-b1','c-n6-b1']
labeldirFiles = {'c-n4-b1':r'C $n_H=10^4 cm^{-3}$ $b=1$'
                ,'c-n5-b1':r'C $n_H=10^5 cm^{-3}$ $b=1$'
                ,'c-n6-b1':r'C $n_H=10^6 cm^{-3}$ $b=1$'
                }
vitDirFiles = [[5,10,15,20,25,30],[5,10,15,20,25,30],[5,10,15,20,25,30]]
axTableXScale='log'
axTableYScale='log'



originalDir=os.getcwd()
for i,directory in enumerate(dirFiles):
    os.chdir(directory)
    currentDirLvl1=os.getcwd()
    for vitesse in vitDirFiles[i]:
        figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(14,14))
        figRatio.subplots_adjust(wspace=0) 
        ax.grid(True)
        ax.set_xlim(10,1E6)
        ax.set_xscale(axTableXScale)
        ax.set_yscale(axTableYScale)
        locmin = mtick.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=6)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(mtick.NullFormatter())
        ax.set_ylabel(r'$T_{n}$ (K) ')
        ax.set_xlabel(r'$t_i$ (annees)')
        
        #ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
        X,Y=([],[])
        vitDirName="v{}".format(str(vitesse))
        if os.path.exists(vitDirName):
            os.chdir(vitDirName)
            currentDirLvl2=os.getcwd()
            """Numerateur"""
            asciiTemp = ascii.read('mhd_phys.out', header_start=0,data_start=1)
            for lignes in asciiTemp:
                X.append(lignes[2])
                Y.append(lignes[14])
            ax.plot(X,Y,linewidth=5,label=r'$T_n=f(t_i)$'+' {}'.format(directory+'-'+vitDirName))
            lgd = plt.legend(fontsize=18,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.3))
            """Enregistrement des figures"""
            fileTitle="Tn=f(timeI)_{}".format(directory+'-'+vitDirName)
            fileTitle=fileTitle.replace('\ce{','')
            fileTitle=fileTitle.replace('}','')
            fileTitle=fileTitle.replace('{','')
            fileTitle=fileTitle.replace('$','')
            fileTitle=fileTitle.replace(' ','_')
            plt.savefig(fileTitle+".pdf", bbox_inches='tight')    
            #plt.savefig(fileTitle+".png", bbox_inches='tight')    
            #plt.show()
        plt.close()  
        if os.path.exists(vitDirName):
            with open('comparaison_CO.dat','w') as outputFile:
                outputFile.write('timeI TdV') #Header
                outputFile.write('\n')
                Molfound=False
                iMol=0
                TableMol=ascii.read('CO_sum_em.out',header_start=1,data_start=2)
                Header=ascii.read('CO_sum_em.out',header_start=1,data_start=1)
                Header=Header[0]
                while(Molfound==False):
                    if (((TableMol[iMol])[2])<1E5):
                        iMol=iMol+1
                        Molfound=False
                    else:
                        TdV= (TableMol[iMol])[39]
                        outputFile.write(str(((TableMol[iMol])[2]))+' '+str(TdV)) 
                        outputFile.write('\n')
                        Molfound=True
                TdV= (TableMol[-1])[39]
                outputFile.write(str(((TableMol[-1])[2]))+' '+str(TdV)) 
                outputFile.write('\n')
        os.chdir(currentDirLvl1)
    os.chdir(originalDir)
   
