p#!/usr/bin/env python2
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
from matplotlib.ticker import AutoMinorLocator
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
dirFiles = ['cj-n5-b2']
labeldirFiles = {}
vitDirFiles = [[30]]
ageDirFiles = [[[50,100,150,200,250,500]]]

axTableXScale='log'
axTableYScale='log'
ActiverYang=True
ActiverMyObs=False

"""OI(63$\mu m$)"""
Eup1=227.71
raieNomNum = 'OI(63$\mu m$) ($E_{up}=$'+str(Eup1)+' K)'
nameFilesNum = 'OH_int_int.out'
nameFilesXHeaderNum = 0 #inutile
iRaieNum = 0
nameFilesYHeaderNum = 2
figDataNum = {}
dataStartNum = 1
HeaderStartNum = 0
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=63.18 #a remplir si les données Yang 2017 sont activées



"""OI(145$\mu m$)"""
Eup2=326.58
raieNomDen = 'OI(145$\mu m$) ('+str(Eup2)+' K)'
nameFilesDen = 'OH_int_int_145.out'
nameFilesXHeaderDen = 0 #inutile
iRaieDen = 0
nameFilesYHeaderDen = 3
figDataDen = {}
dataStartDen = 1
HeaderStartDen = 0
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=145.54 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}

figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))
figRatio.subplots_adjust(wspace=0) 
ax.grid(True,axis='both',alpha=0.5,zorder=0)
ax.set_xlim(7,1000)
ax.set_xscale(axTableXScale)
ax.set_yscale(axTableYScale)
ax.set_ylabel(r'$T_e$ (K) ')
ax.set_xlabel(r'$t_i$ (ans)')
#ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
labelList = {}









"""Sauvegarde des TdV"""

originalDir=os.getcwd()
for i,directory in enumerate(dirFiles):
    os.chdir(directory)
    currentDirLvl1=os.getcwd()
    figDataNum[directory]={}
    figDataDen[directory]={}
    labeldirFiles[directory]={}
    Yratio[directory]={}
    for j,vitesse in enumerate(vitDirFiles[i]):        
        (figDataNum[directory])[vitesse]={}
        (figDataDen[directory])[vitesse]={}
        (Yratio[directory])[vitesse]={}
        (labeldirFiles[directory])[vitesse]={}
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            """Label"""
            ((labeldirFiles[directory])[vitesse])[tempsJ]=directory+"-v"+str(vitesse)+"-t"+str(tempsJ)
            vitDirName="v{:1}-t{:2}".format(str(vitesse),str(tempsJ))
            os.chdir(vitDirName)
            currentDirLvl2=os.getcwd()
            """Numerateur"""
            asciiTemp = ascii.read('mhd_phys.out',header_start=HeaderStartNum,data_start=dataStartNum)
            ((figDataNum[directory])[vitesse])[tempsJ]={'Tn':[],'Ti':[]}
            for line in asciiTemp:
                ((figDataNum[directory])[vitesse])[tempsJ]['Tn'].append(line[16])
                ((figDataNum[directory])[vitesse])[tempsJ]['Ti'].append(line[15])
            #temp = (figDataNum[directory])[vitesse][tempsJ]/(1.e-12*lambdaNumMicrons**3/(2.*KB)/1.e5)
            #((figDataNum[directory])[vitesse])[tempsJ]= temp
            #print(directory,vitDirName,(figDataNum[directory])[vitesse][tempsJ])
            """Dénominateur"""
            asciiTemp2 = ascii.read('mhd_phys.out',header_start=HeaderStartDen,data_start=dataStartDen)
            ((figDataDen[directory])[vitesse])[tempsJ]=[]
            for line in asciiTemp2:
                ((figDataDen[directory])[vitesse])[tempsJ].append(line[2])
            #print(directory,vitDirName,(figDataDen[directory])[vitesse][tempsJ]/(1.e-12*lambdaDenMicrons**3/(2.*KB)/1.e5))
            """Ratio"""
            #((Yratio[directory])[vitesse])[tempsJ]=((figDataNum[directory])[vitesse])[tempsJ]/((figDataDen[directory])[vitesse])[tempsJ]
            os.chdir(currentDirLvl1)
    os.chdir(originalDir)

"""Plot des modeles"""
X=[]
Y=[]
dirdir=[]
l=0
for i,directory in enumerate(dirFiles):
    dirdir.append(directory)
    for j,vitesse in enumerate(vitDirFiles[i]):
        X=[]
        Y=[]
        Y2=[]
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            X=((figDataDen[directory])[vitesse])[tempsJ]
            Y=((figDataNum[directory])[vitesse])[tempsJ]['Tn']
            #Y2=((figDataNum[directory])[vitesse])[tempsJ]['Ti']
            #Y.append(((Yratio[directory])[vitesse])[tempsJ])
            print(tempsJ)
            ax.plot(X,Y,'-',label=((labeldirFiles[directory])[vitesse])[tempsJ],zorder=len((ageDirFiles[i])[j])-k,color=labelTextColor[k])
            #ax.plot(X,Y2,'--',label=((labeldirFiles[directory])[vitesse])[tempsJ],zorder=len((ageDirFiles[i])[j])-k,color=labelTextColor[k])
            print(max(Y))
            Xline=np.full((100,),tempsJ)
            itempsJ=np.where(np.array(X)<tempsJ)[0][-1]
            Yline=np.linspace(0,Y[itempsJ],num=100)
            ax.plot(Xline,Yline,'--',color='black',linewidth=1.0,zorder=1,alpha=0.5)
            if(tempsJ<10000.): 
                plt.text(1.05*tempsJ,3.5*Y[itempsJ],r'\textbf{%s ans}' % (str(tempsJ)),color=labelTextColor[k],fontsize=12)
                #plt.text(tempsJ,1e-1,r'\textbf{%s}' % (str(tempsJ)),color='black',fontsize=12,rotation=90)
            #ax.axvline(x=tempsJ,ymin=0,ymax=1,ls='--',color='grey',linewidth=0.5)
        l=l+1
    

"""========================================================================"""        
xmin=50
xmax=500
locmin = mtick.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(mtick.NullFormatter())

lgd = plt.legend(fontsize=14,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.45))

"""Enregistrement des figures"""
fileTitle="Te_"+"{}".format(raieNomNum)
fileTitle=fileTitle.replace('\ce{','')
fileTitle=fileTitle.replace('}','')
fileTitle=fileTitle.replace('{','')
fileTitle=fileTitle.replace('$','')
fileTitle=fileTitle.replace(' ','_')
fileTitle=fileTitle.replace('.','v')
fileTitle=fileTitle.replace('$\mu m$','m')
l=0
for i,directory in enumerate(dirFiles):
    for j,vitesse in enumerate(vitDirFiles[i]):
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            if(tempsJ<10000.): 
                plt.text(tempsJ,1.5*ax.get_ylim()[0],r'%s' % (str(tempsJ)),color='black',fontsize=12,rotation=90)
            #ax.axvline(x=tempsJ,ymin=0,ymax=1,ls='--',color='grey',linewidth=0.5)
        l=l+1
plt.savefig(fileTitle+".pdf", bbox_inches='tight')    
plt.savefig(fileTitle+".png", bbox_inches='tight')    
plt.show()
plt.close()      
