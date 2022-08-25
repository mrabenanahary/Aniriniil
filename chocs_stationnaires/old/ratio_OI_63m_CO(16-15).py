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
dirFiles = ['c-n4-b1','c-n5-b1','c-n6-b1']
labeldirFiles = {'c-n4-b1':r'C $n_H=10^4 cm^{-3}$ $b=1$'
                ,'c-n5-b1':r'C $n_H=10^5 cm^{-3}$ $b=1$'
                ,'c-n6-b1':r'C $n_H=10^6 cm^{-3}$ $b=1$'
                }
vitDirFiles = [[5,10,15,19],[5,10,15,19],[5,10,15,19]]

axTableXScale='linear'
axTableYScale='log'
ActiverYang=True
ActiverMyObs=False

"""OI(63$\mu m$)"""
Eup1=227.71
raieNomNum = 'OI(63$\mu m$) ('+str(Eup1)+' K)'
nameFilesNum = 'OH_int_int.out'
nameFilesXHeaderNum = 0 #inutile
iRaieNum = 0
nameFilesYHeaderNum = 3
figDataNum = {}
dataStartNum = 1
HeaderStartNum = 0
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=63.18 #a remplir si les données Yang 2017 sont activées



"""CO(16-15)"""
Eup2=751.72
raieNomDen = 'CO(16-15) ('+str(Eup2)+' K)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0 #inutile
iRaieDen = 15
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=162.82 #a remplir si les données Yang 2017 sont activées


"""Ratio"""
Yratio = {}

figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))
figRatio.subplots_adjust(wspace=0) 
ax.grid(True)
ax.set_xscale(axTableXScale)
ax.set_yscale(axTableYScale)
ax.set_ylabel(r'$\sum T_{mb}\Delta v$ '+r'{:1}/{:2}'.format(raieNomNum,raieNomDen))
ax.set_xlabel(r'$v_s$ (km/s)')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))

originalDir=os.getcwd()
for i,directory in enumerate(dirFiles):
    os.chdir(directory)
    currentDirLvl1=os.getcwd()
    figDataNum[directory]={}
    figDataDen[directory]={}
    Yratio[directory]={}
    for vitesse in vitDirFiles[i]:
        vitDirName="v{}".format(str(vitesse))
        os.chdir(vitDirName)
        currentDirLvl2=os.getcwd()
        """Numerateur"""
        asciiTemp = ascii.read(nameFilesNum,header_start=HeaderStartNum,data_start=dataStartNum)
        
        (figDataNum[directory])[vitesse]=(asciiTemp[iRaieNum])[nameFilesYHeaderNum]
        print(directory,vitDirName,(figDataNum[directory])[vitesse]/(1.e-12*lambdaNumMicrons**3/(2.*KB)/1.e5))
        """Dénominateur"""
        asciiTemp2 = ascii.read(nameFilesDen,header_start=HeaderStartDen,data_start=dataStartDen)
        
        (figDataDen[directory])[vitesse]=(asciiTemp2[iRaieDen])[nameFilesYHeaderDen]
        print(directory,vitDirName,(figDataDen[directory])[vitesse]/(1.e-12*lambdaDenMicrons**3/(2.*KB)/1.e5))
        """Ratio"""
        (Yratio[directory])[vitesse]=(figDataNum[directory])[vitesse]/(figDataDen[directory])[vitesse]
        os.chdir(currentDirLvl1)
    os.chdir(originalDir)

"""Plot des modeles"""
X=[]
Y=[]
dirdir=[]
for i,directory in enumerate(Yratio.keys()):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in sorted((Yratio[directory]).keys()):
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'-o',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])
    
"""Numérateur"""
figDataNum = {}

"""Dénominateur"""
figDataDen = {}
"""Ratio"""
Yratio = {}

originalDir=os.getcwd()
dirFiles = ['c-n4-b1','c-n5-b1','c-n6-b1']
labeldirFiles = {'c-n4-b1':r'J $n_H=10^4 cm^{-3}$ $b=1$'
                ,'c-n5-b1':r'J $n_H=10^5 cm^{-3}$ $b=1$'
                ,'c-n6-b1':r'J $n_H=10^6 cm^{-3}$ $b=1$'}
vitDirFiles = [[20,25,30],[20,25,30],[20,25,30]]

for i,directory in enumerate(dirFiles):
    os.chdir(directory)
    currentDirLvl1=os.getcwd()
    figDataNum[directory]={}
    figDataDen[directory]={}
    Yratio[directory]={}
    for vitesse in vitDirFiles[i]:
        vitDirName="v{}".format(str(vitesse))
        os.chdir(vitDirName)
        currentDirLvl2=os.getcwd()
        """Numerateur"""
        asciiTemp = ascii.read(nameFilesNum,header_start=HeaderStartNum,data_start=dataStartNum)
        (figDataNum[directory])[vitesse]=(asciiTemp[iRaieNum])[nameFilesYHeaderNum]
        """Dénominateur"""
        asciiTemp2 = ascii.read(nameFilesDen,header_start=HeaderStartDen,data_start=dataStartDen)
        (figDataDen[directory])[vitesse]=(asciiTemp2[iRaieDen])[nameFilesYHeaderDen]
        """Ratio"""
        (Yratio[directory])[vitesse]=(figDataNum[directory])[vitesse]/(figDataDen[directory])[vitesse]
        os.chdir(currentDirLvl1)
    os.chdir(originalDir)

"""Plot des modeles"""
X=[]
Y=[]
dirdir=[]
for i,directory in enumerate(Yratio.keys()):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in sorted((Yratio[directory]).keys()):
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'-.x',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])

"""========================================================================"""
"""Modeles de chocs J: input"""
originalDir=os.getcwd()
dirFiles = ['j-n4-b0.1','j-n5-b0.1','j-n6-b0.1']
labeldirFiles = {'j-n4-b0.1':r'J $n_H=10^4 cm^{-3}$ $b=0.1$'
                ,'j-n5-b0.1':r'J $n_H=10^5 cm^{-3}$ $b=0.1$'
                ,'j-n6-b0.1':r'J $n_H=10^6 cm^{-3}$ $b=0.1$'
                }
vitDirFiles = [[5,10,15,20,25,30],[5,10,15,20,25,30],[5,10,15,20,25,30]]
axTableXScale='linear'
axTableYScale='log'
ActiverYang=True
ActiverMyObs=False

"""OI(63$\mu m$)"""
Eup1=227.71
raieNomNum = 'OI(63$\mu m$) ('+str(Eup1)+' K)'
nameFilesNum = 'OH_int_int.out'
nameFilesXHeaderNum = 0 #inutile
iRaieNum = 0
nameFilesYHeaderNum = 3
figDataNum = {}
dataStartNum = 1
HeaderStartNum = 0
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=63.18 #a remplir si les données Yang 2017 sont activées


"""CO(16-15)"""
Eup2=751.72
raieNomDen = 'CO(16-15) ('+str(Eup2)+' K)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0 #inutile
iRaieDen = 15
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=162.82 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}
"""
figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))
figRatio.subplots_adjust(wspace=0) 
ax.grid(True)
ax.set_xscale(axTableXScale)
ax.set_yscale(axTableYScale)
ax.set_ylabel(r'$\sum T_{mb}\Delta v$ '+r'{:1}/{:2}'.format(raieNomNum,raieNomDen))
ax.set_xlabel(r'$v_s$ (km/s)')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
"""
originalDir=os.getcwd()
for i,directory in enumerate(dirFiles):
    os.chdir(directory)
    currentDirLvl1=os.getcwd()
    figDataNum[directory]={}
    figDataDen[directory]={}
    Yratio[directory]={}
    for vitesse in vitDirFiles[i]:
        vitDirName="v{}".format(str(vitesse))
        os.chdir(vitDirName)
        currentDirLvl2=os.getcwd()
        """Numerateur"""
        asciiTemp = ascii.read(nameFilesNum,header_start=HeaderStartNum,data_start=dataStartNum)
        (figDataNum[directory])[vitesse]=(asciiTemp[iRaieNum])[nameFilesYHeaderNum]
        """Dénominateur"""
        asciiTemp = ascii.read(nameFilesDen,header_start=HeaderStartDen,data_start=dataStartDen)
        (figDataDen[directory])[vitesse]=(asciiTemp[iRaieDen])[nameFilesYHeaderDen]
        """Ratio"""
        (Yratio[directory])[vitesse]=(figDataNum[directory])[vitesse]/(figDataDen[directory])[vitesse]
        os.chdir(currentDirLvl1)
    os.chdir(originalDir)

"""Plot des modeles"""
X=[]
Y=[]
dirdir=[]
for i,directory in enumerate(Yratio.keys()):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in sorted((Yratio[directory]).keys()):
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'--d',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])



"""========================================================================"""        
"""Observations hors Yang"""
TelescopeNameNum='HIFI'
TelescopeNameDen='HIFI'
raieNomRatio = TelescopeNameNum+' '+raieNomNum+'/'+TelescopeNameDen+' '+raieNomDen
nameFilesRatio = "ratio_H2O_same.dat"
iRatio = 0 #numero de ligne Ascii
jRatio = 1 #numero de colonne Ascii
HeaderStartRatio=0
DataStartRatio=0
xmin=np.amin(X)
xmax=np.amax(X)

if (ActiverMyObs==True):
    fig4data = ascii.read(nameFilesRatio)
    Header = (ascii.read(nameFilesRatio,header_start=HeaderStartRatio,data_start=DataStartRatio))
    Header = Header[0]
    
    ratioX = np.arange(xmin,xmax,0.1)
    ratioY = np.full(ratioX.shape,(fig4data[iRatio])[jRatio])
    #print(fig4data)
    
    ax.plot(ratioX,ratioY,'-',label=raieNomRatio,color='red')
    


"""========================================================================"""        
"""Observations  Yang (si activées)"""
ratioXYang = np.arange(xmin,xmax,0.1)
FluxNumYang=516. #1e-18 W/m2
FluxDenYang=171. #1e-18 W/m2
errorNumYang = 14.6 #1e-18 W/m2
errorDenYang = 19.0 #1e-18 W/m2
lambdaNumerateur=lambdaNumMicrons #microns
lambdaDenominateur=lambdaDenMicrons #microns

if(ActiverYang==True):
    labelRatioYang='PACS'+' '+raieNomNum+'/'+raieNomDen
    Numerateur=FluxNumYang#CO1615Yang
    Denominateur=FluxDenYang#OI63Yang
    ratioYY = Numerateur/Denominateur
    correction = ((lambdaNumerateur/lambdaDenominateur)**3) #rapport des TdV
    errNum=errorNumYang
    errDen=errorDenYang
    ratioYYang = np.array(ratioYY)*correction # pour convertir le rapport de flux en rapport de TdV
    errorRatio = ratioYYang*((errNum/Numerateur)+(errDen/Denominateur)) 
    
    ax.fill_between(ratioXYang,ratioYYang-errorRatio,ratioYYang+errorRatio,label=labelRatioYang,color='pink',alpha=0.7)


lgd = plt.legend(fontsize=18,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.5))

"""Enregistrement des figures"""
fileTitle="ratio_"+"{}".format(raieNomNum)+"_"+"{}".format(raieNomDen)
fileTitle=fileTitle.replace('\ce{','')
fileTitle=fileTitle.replace('}','')
fileTitle=fileTitle.replace('{','')
fileTitle=fileTitle.replace('$','')
fileTitle=fileTitle.replace(' ','_')
plt.savefig(fileTitle+".pdf", bbox_inches='tight')    
plt.savefig(fileTitle+".png", bbox_inches='tight')    
plt.show()
plt.close()  
