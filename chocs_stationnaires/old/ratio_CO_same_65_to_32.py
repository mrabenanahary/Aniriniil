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

fig4data = ascii.read("ratio_CO_same.dat")
Header = (  ascii.read("ratio_CO_same.dat",header_start=0,data_start=0))
Header = Header[0]
print(fig4data)

fig, ax = plt.subplots(1,1, sharex=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
ax.grid(False)

xmin=0
xmax=35
ymin=1e-1
ymax=1e3
ax.set_xlim(0,xmax)
#ax.set_ylim(ymin,ymax)

#ax.set_title(r'\textbf{Observed $CO$ $\int T\mathrm{d}V$ against shock waves models}',fontsize=22)
ax.set_yscale('linear')
ax.set_xscale('linear')
#ax.xaxis.set_major_locator(mtick.AutoMinorLocator(5))
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5)) # <<<<<<<<<<<<<< tres important a adapter!! 
#ax.tick_params(bottom="on", top="on",left="on",right="on")

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='both', direction='in')
#ax.tick_params(axis='both', direction='in',labelbottom="on", labeltop="on", labelleft="on", labelright="on")

#ax.set_xlim(1e9,1e20)
ax.set_xlabel(r'$v_{s}$ (K.km/s)',fontsize=23,fontweight='extra bold')
ax.set_ylabel(r'rapport de $\sum T\mathrm{d}v$ CO(6-5)/CO(3-2)',fontsize=23,fontweight='extra bold')

labelTextColor=['red','blue','green','brown','black']
Size=np.size(fig4data)


#J_up=np.array([0,1,2,3,5,6,7])+2
ratio1X = np.arange(xmin,xmax,0.1)
ratio2X = np.arange(xmin,xmax,0.1)
ratio1Y = np.full(ratio1X.shape,(fig4data[0])[1])
ratio2Y = np.full(ratio1X.shape,(fig4data[1])[1])


#liste de tous les labels utilisés
#liste de tous les labels utilisés
labels=[r'$n_H=10^4 cm^{-3}$ $b=1$',
        r'$n_H=10^5 cm^{-3}$ $b=1$',
        r'$n_H=10^6 cm^{-3}$ $b=1$']

        
#dirlist = os.popen('ls -1d c-*b1').read().split()
dirlist = ['c-n4-b1']
vlist = [5, 10, 15, 20, 25, 30]



#ulist = os.popen('ls -1d '+cheminDensiteChamp+"/"+"v*").read().split()
ulist = ['c-n4-b1/v5', 'c-n4-b1/v10', 'c-n4-b1/v15', 'c-n4-b1/v20', 'c-n4-b1/v25', 'c-n4-b1/v30']
J_up = []
Y=[]
for kk,vitesse in enumerate(ulist):
    """CO(6-5)/CO(3-2)"""
    chemin= vitesse + "/" + "CO_int_int.out"
    fig4data = ascii.read(chemin,header_start=None,data_start=2)
    #Header = (ascii.read(chemin,header_start=None,data_start=2))
    print(fig4data)        
    J_up.append(vlist[kk])
    #print(Y)
    Y.append((fig4data[5])[1]/(fig4data[2])[1])
ax.plot(J_up,Y,'-o',markersize=10,label=labels[0])  

ulist = ['c-n5-b1/v5', 'c-n5-b1/v10', 'c-n5-b1/v15', 'c-n5-b1/v20', 'c-n5-b1/v25', 'c-n5-b1/v30']
J_up = []
Y=[]
for kk,vitesse in enumerate(ulist):
    """CO(6-5)/CO(3-2)"""
    chemin= vitesse + "/" + "CO_int_int.out"
    fig4data = ascii.read(chemin,header_start=None,data_start=2)
    #Header = (ascii.read(chemin,header_start=None,data_start=2))
    print(fig4data)        
    J_up.append(vlist[kk])
    #print(Y)
    Y.append((fig4data[5])[1]/(fig4data[2])[1])
ax.plot(J_up,Y,'-o',markersize=10,label=labels[1])  

ulist = ['c-n6-b1/v5', 'c-n6-b1/v10', 'c-n6-b1/v15', 'c-n6-b1/v20', 'c-n6-b1/v25', 'c-n6-b1/v30']
J_up = []
Y=[]
for kk,vitesse in enumerate(ulist):
    """CO(6-5)/CO(3-2)"""
    chemin= vitesse + "/" + "CO_int_int.out"
    fig4data = ascii.read(chemin,header_start=None,data_start=2)
    #Header = (ascii.read(chemin,header_start=None,data_start=2))
    print(fig4data)        
    J_up.append(vlist[kk])
    #print(Y)
    Y.append((fig4data[5])[1]/(fig4data[2])[1])
ax.plot(J_up,Y,'-o',markersize=10,label=labels[2])  



#ax.plot(ratio1X,ratio1Y,label=r'CO(6-5)/CO(3-2) Observ.')        
ax.plot(ratio2X,ratio2Y,label=r'CO(6-5)/CO(3-2) Observ.')        



#plt.text(15, 3e2,r'\textbf{$\ce{CO}$}',fontsize=36,weight='bold')        

lgd = plt.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.3))
#labelTextColor.append('black')
#for ii,text in enumerate(lgd.get_texts()):
#    text.set_color(labelTextColor[ii])

plt.savefig("ratio_CO_same_65_to_32.png", bbox_inches='tight')
plt.savefig("ratio_CO_same_65_to_32.pdf", bbox_inches='tight')
plt.show()    
plt.close()
