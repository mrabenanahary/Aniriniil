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

fig4data = ascii.read("K_CO_TdV.txt")
Header = (  ascii.read("K_CO_TdV.txt",header_start=0,data_start=0))
Header = Header[0]
print(fig4data)

fig, ax = plt.subplots(1,1, sharex=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
ax.grid(False)

xmax=23
ymin=1
ymax=1e3
ax.set_xlim(0,xmax)
ax.set_ylim(ymin,ymax)

#ax.set_title(r'\textbf{Observed $CO$ $\int T\mathrm{d}V$ against shock waves models}',fontsize=22)
ax.set_yscale('log')
ax.set_xscale('linear')
#ax.xaxis.set_major_locator(mtick.AutoMinorLocator(5))
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5)) # <<<<<<<<<<<<<< tres important a adapter!! 
#ax.tick_params(bottom="on", top="on",left="on",right="on")

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='both', direction='in')
#ax.tick_params(axis='both', direction='in',labelbottom="on", labeltop="on", labelleft="on", labelright="on")

#ax.set_xlim(1e9,1e20)
ax.set_xlabel(r'$J_{up}$',fontsize=23,fontweight='extra bold')
ax.set_ylabel(r'\textbf{$\int$ TdV} ($K.km.s^{-1}$)',fontsize=23,fontweight='extra bold')

labelTextColor=['red','blue','green','brown','black']
Size=np.size(fig4data)

""" Traitement de TOUS les modeles de chocs presents dans le dossier en cours"""
Ytot=[]
#J_up=np.array([0,1,2,3,5,6,7])+2
E_up=np.array([509.850, 
               1015.153,
               1681.678,
               2503.870,
               4586.377,
               5829.758,
               7196.995])

#liste de tous les labels utilisés
labels=[r'$n_H=10^6 cm^{-3}$ $b=1$ $V_s=5 km.s^{-1}$ (choc C)',
           r'$n_H=10^6 cm^{-3}$ $b=1$ $V_s=10 km.s^{-1}$ (choc C)',
           r'$n_H=10^6 cm^{-3}$ $b=1$ $V_s=20 km.s^{-1}$ (choc J)',
           r'$n_H=10^6 cm^{-3}$ $b=1$ $V_s=25 km.s^{-1}$ (choc J)',
           r'$n_H=10^6 cm^{-3}$ $b=1$ $V_s=30 km.s^{-1}$ (choc J)'
           ]

        
#dirlist = os.popen('ls -1d c-*b1').read().split()
dirlist = ['c-n6-b1']

for ii,densiteChamp in enumerate(dirlist):
    cheminDensiteChamp=densiteChamp
    #ulist = os.popen('ls -1d '+cheminDensiteChamp+"/"+"v*").read().split()
    ulist = ['c-n6-b1/v5', 'c-n6-b1/v10', 'c-n6-b1/v20', 'c-n6-b1/v25', 'c-n6-b1/v30']
    for kk,vitesse in enumerate(ulist):
        J_up = []
        Y=[]
        Ytot.append([])
        chemin= vitesse + "/" + "CO_int_int.out"
        fig4data = ascii.read(chemin,header_start=None,data_start=2)
        #Header = (ascii.read(chemin,header_start=None,data_start=2))
        print(fig4data)
        
        fig4index = range(0,np.size(fig4data))


        
        for levels in fig4index:
            J_up.append(fig4data[levels][0])
            Y.append(fig4data[levels][1])
        Ytot.append(0)
        Ytot[kk]=Y
        print("Y=",Y)
        #print(Y)
        ax.plot(J_up,Y,'-o',markersize=10,label=labels[kk])  
    

plt.text(15, 2e-3,r'\textbf{$\ce{CO}$}',fontsize=36,weight='bold')        
plt.text(1, 1e-2,r'Tous les modeles: $v_s=20$ $km.s^{-1}$,',fontsize=18,weight='bold')    
plt.text(1, 5e-3,r'$n_H=10^4$ $cm^{-3}$, $B=100$ $\mu G$',fontsize=18,weight='bold')    
"""    

#steady state
Y=[]
Ytot.append([])
labels.append(Header[7])
labels[ii+1]=labels[ii+1].replace('_',' ')
for levels in figIndex:
    Y.append(fig4data[levels][7])
Ytot[ii+1]=Y
#print(Y)
#ax.scatter(np.array(J_up),np.array(Y),'o',edgecolor='black',facecolor=None,label=labels[ii])
#ax.scatter(np.array(J_up),np.array(Y),s=100,marker='o',edgecolors='black',facecolors='none',label='eee')
ax.plot(np.array(J_up),np.array(Y),'--o',markersize=10,color=labelTextColor[ii+1],markeredgecolor=labelTextColor[ii+1],markerfacecolor='none',label=labels[ii+1])
"""

#observationnal data:
Xobs=np.array([3,
               4,
               6,
               7]
               )
Yobs=np.array([195.3,
               130.79,
               146.42,
               77.32])
ax.scatter(Xobs,Yobs,s=100,zorder=6,marker='D',edgecolors='black',facecolors='yellow',label=r'Observations')

#legend_properties = {'weight':'bold'}

#plt.plot([1,2,3], [4,5,6], label='Test')
#plt.legend(prop=legend_properties)

lgd = plt.legend(fontsize=12,loc='best')
labelTextColor.append('black')
for ii,text in enumerate(lgd.get_texts()):
    text.set_color(labelTextColor[ii])

plt.savefig("CO_TdV_c-n6.png", bbox_inches='tight')
plt.savefig("CO_TdV_c-n6.pdf", bbox_inches='tight')
plt.show()    
plt.close()
