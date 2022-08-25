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
from labellines import labelLine, labelLines
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
        'size'   : 22,
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



    
figRatio, [ax,ax2,ax3,ax4,ax5,ax6] = plt.subplots(1,6, sharey=True,figsize=(35,7))
figRatio.subplots_adjust(wspace=0) 
ax.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax5.grid(True)
ax6.grid(True)
ax.set_xscale('linear')
ax2.set_xscale('linear')
ax3.set_xscale('linear')
ax4.set_xscale('linear')
ax5.set_xscale('linear')
ax6.set_xscale('linear')
ax.set_yscale('linear')
ax.set_ylabel(r'$log(\frac{N_u}{g_u})$')
ax.set_xlabel(r'$E_{up}$ (K)')
ax2.set_xlabel(r'$E_{up}$ (K)')
ax3.set_xlabel(r'$E_{up}$ (K)')
ax4.set_xlabel(r'$E_{up}$ (K)')
ax5.set_xlabel(r'$E_{up}$ (K)')
ax6.set_xlabel(r'$E_{up}$ (K)')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
ax2.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
ax3.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
ax4.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
ax5.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
ax6.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))


densiteVitesseAge=['cj-n5-b1/v19-t50','cj-n5-b1/v19-t100','cj-n5-b1/v19-t150','cj-n5-b1/v19-t200','cj-n5-b1/v19-t250','cj-n5-b1/v19-t500']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax.plot(X, Y,'-D', markersize=5, alpha=0.6,markeredgecolor='black',label=element)
    

densiteVitesseAge=['cj-n5-b2/v19-t50','cj-n5-b2/v19-t100','cj-n5-b2/v19-t150','cj-n5-b2/v19-t200','cj-n5-b2/v19-t250','cj-n5-b2/v19-t500']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax2.plot(X, Y,'--o', markersize=5, alpha=0.3,markeredgecolor='black',label=element)
    

densiteVitesseAge=['cj-n5-b2/v25-t50','cj-n5-b2/v25-t100','cj-n5-b2/v25-t150','cj-n5-b2/v25-t200','cj-n5-b2/v25-t250','cj-n5-b2/v25-t500']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax3.plot(X, Y,'-D', markersize=5, alpha=0.3,markeredgecolor='black',label=element)
    


densiteVitesseAge=['c-n5-b1/v5','c-n5-b1/v10','c-n5-b1/v15','c-n5-b1/v19']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax4.plot(X, Y,'--D', markersize=5, alpha=0.6,markeredgecolor='black',label=element)


densiteVitesseAge=['j-n5-b0.1/v5','j-n5-b0.1/v10','j-n5-b0.1/v15','j-n5-b0.1/v19','j-n5-b0.1/v20','j-n5-b0.1/v25','j-n5-b0.1/v30']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax5.plot(X, Y,'-o', markersize=5, alpha=0.6,markeredgecolor='black',label=element)

densiteVitesseAge=['j-n5-b1/v5','j-n5-b1/v10','j-n5-b1/v15','j-n5-b1/v20','j-n5-b1/v30']
"""dirFiles = ['cj-n4-b1','cj-n5-b1']
vitDirFiles = [[5,10],[25]]
ageDirFiles = [[[100,250,500],[100,250,500]],[[100,150,250,500]]]    
densiteVitesseAge=[]
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):        
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            densiteVitesseAge.append(directory+'/v{:1}-t{:2}'.format(vitesse,tempsJ))
"""            
figObs=ascii.read('H2_int_int_avec_error.out')
figCJ=[]
for i,element in enumerate(densiteVitesseAge):
    X=[]
    Y=[]
    asciiTemp=ascii.read(element+'/excit.out')
    for l in range(2,6):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    for l in range(7,9):
        X.append(asciiTemp[l][2])
        Y.append(asciiTemp[l][3])
    X.append(asciiTemp[13][2])
    Y.append(asciiTemp[13][3])
    ax6.plot(X, Y,'-o', markersize=5, alpha=0.6,markeredgecolor='black',label=element)
    #ax5.plot(X, Y,'-o', markersize=5, alpha=0.6,markeredgecolor='black',label=element)
    #ax6.plot(X, Y,'-o', markersize=5, alpha=0.6,markeredgecolor='black',label=element)


ax.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')
ax2.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')
ax3.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')
ax4.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')
ax5.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')
ax6.scatter([row[0] for row in figObs], [row[-2] for row in figObs], zorder = 3, s=50, marker='D',facecolor='yellow',edgecolor='black',label=r'Observations')

ax.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')
ax2.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')
ax3.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')
ax4.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')
ax5.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')
ax6.errorbar([row[0] for row in figObs], [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')

lgd = ax.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.62))
lgd = ax2.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.42))
lgd = ax3.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.42))
lgd = ax4.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.42))
lgd = ax5.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.42))
lgd = ax6.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.42))

plt.savefig('H2_DEXT_CJ_C_J.pdf', bbox_inches='tight')

plt.show()
plt.close()  
