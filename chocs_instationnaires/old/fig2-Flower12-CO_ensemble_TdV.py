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










prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
labelTextColor= colors

originalDir=os.getcwd()
dirFiles = ['./chocs_cj_18_06_19']
nameFiles = ['CO_int_int.out']
nameFilesXHeader = [0]
nameFilesYHeader = [[1]]
figData = {}
dataStart = [2]
HeaderStart = [None] 
dirlist = []
vitesseTime = []
linewdthsize=[4,3,2,1]

#figExcit, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
#ax.grid(False)
#figExcit.subplots_adjust(wspace=0)



for i,directory in enumerate(dirFiles):       
    figData[directory]={}
    os.chdir(directory)
    dirlist.append(os.popen('ls -1d cj-*b1').read().split())
    currentDirLvL1=os.getcwd()
    for k,dirdir in enumerate(dirlist[i]):
        (figData[directory])[dirdir]={}
        os.chdir(dirdir)
        #print(os.getcwd())
        vitesseTime.append(os.popen('ls -1d v*-t*').read().split())
        currentDirLvL2=os.getcwd()
        for l,vt in enumerate(vitesseTime[k]):
            ((figData[directory])[dirdir])[vt]={}
            os.chdir(vt)
            #print(vt)
            if os.path.exists(nameFiles[i]): 
                asciiTemp=ascii.read(nameFiles[i],header_start=HeaderStart[i],data_start=dataStart[i])
                (((figData[directory])[dirdir])[vt])['X']=[row[nameFilesXHeader[i]] for row in asciiTemp]
                (((figData[directory])[dirdir])[vt])['Y']=[]
                for m,Yheader in enumerate(nameFilesYHeader[i]):
                    #print(Yheader)
                    #if(iAx2%16==0):
                    #    iAx=0
                    #    figExcit, [axTable[0],axTable[1],axTable[2],axTable[3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                    #    figExcit.subplots_adjust(wspace=0) 
                    (((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                    #X=(((figData[directory])[dirdir])[vt])['X']
                    #Y=(((figData[directory])[dirdir])[vt])['Y']
                    #axTable[iAx].plot(X,Y,label=dirdir+vt,alpha=0.65,linewidth=linewdthsize[iAx2%4])
                    #iAx2=iAx2+1
                    #if(iAx2%4==0):
                    #    lgd = axTable[iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
                    #    iAx=iAx+1
                    #if(iAx2%16==0):
                    #    plt.show()
                    #    plt.close()                        
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)

iAx = 0
iAx2 = 0
axTable = [[],[],[],[]]

for i,directory in enumerate(dirFiles):       
    #figData[directory]={}
    os.chdir(directory)
    dirlist.append(os.popen('ls -1d cj-*b1').read().split())
    currentDirLvL1=os.getcwd()
    for k,dirdir in enumerate(dirlist[i]):
        #(figData[directory])[dirdir]={}
        os.chdir(dirdir)
        #print(os.getcwd())
        vitesseTime.append(os.popen('ls -1d v*-t*').read().split())
        currentDirLvL2=os.getcwd()
        for l,vt in enumerate(vitesseTime[k]):
            #((figData[directory])[dirdir])[vt]={}
            os.chdir(vt)
            currentDirLvL3=os.getcwd()
            #print(vt)
            if os.path.exists(nameFiles[i]): 
                #asciiTemp=ascii.read(nameFiles[i],header_start=HeaderStart[i],data_start=dataStart[i])
                #(((figData[directory])[dirdir])[vt])['X']=[row[nameFilesXHeader[i]] for row in asciiTemp]
                #(((figData[directory])[dirdir])[vt])['Y']=[]
                for m,Yheader in enumerate(nameFilesYHeader[i]):
                    #print(Yheader)
                    if(iAx2%16==0):
                        iAx=0
                        figExcit, [axTable[0],axTable[1],axTable[2],axTable[3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                        figExcit.subplots_adjust(wspace=0) 
                    #(((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                    X=(((figData[directory])[dirdir])[vt])['X']
                    Y=(((figData[directory])[dirdir])[vt])['Y']
                    axTable[iAx].plot(X,Y,label=dirdir+vt,alpha=0.65,linewidth=linewdthsize[iAx2%4])
                    iAx2=iAx2+1
                    if(iAx2%4==0):
                        lgd = axTable[iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
                        iAx=iAx+1
                    if(iAx2%16==0):
                        os.chdir(originalDir)
                        figExcit.savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
                        plt.show()
                        plt.close()                        
                        os.chdir(currentDirLvL3)
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)

if(iAx2%16!=0):
    lgd = axTable[iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
    figExcit.savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
    plt.show()
    plt.close()                        