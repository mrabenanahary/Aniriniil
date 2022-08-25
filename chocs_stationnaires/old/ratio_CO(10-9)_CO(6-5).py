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
        'size'   : 14,
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
ActiverMyObs=True

"""Numérateur"""
raieNomNum = 'CO(10-9)'
nameFilesNum = 'CO_int_int.out'
nameFilesXHeaderNum = 0
iRaieNum = 9
nameFilesYHeaderNum = 1
figDataNum = {}
dataStartNum = 2
HeaderStartNum = None
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=260.25 #a remplir si les données Yang 2017 sont activées



"""Dénominateur"""
raieNomDen = 'CO(6-5)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0
iRaieDen = 5
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=433.54 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}

figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(5,5))
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
for i,directory in enumerate(dirFiles):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in vitDirFiles[i]:
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'-o',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])
"""========================================================================"""
"""Modeles de chocs C/J: input"""
originalDir=os.getcwd()
dirFiles = ['c-n4-b1','c-n5-b1','c-n6-b1']
labeldirFiles = {'c-n4-b1':r'J $n_H=10^4 cm^{-3}$ $b=1$'
                ,'c-n5-b1':r'J $n_H=10^5 cm^{-3}$ $b=1$'
                ,'c-n6-b1':r'J $n_H=10^6 cm^{-3}$ $b=1$'
                }
vitDirFiles = [[20,25,30],[20,25,30],[20,25,30]]
axTableXScale='linear'
axTableYScale='log'
ActiverYang=True
ActiverMyObs=True

"""Numérateur"""
raieNomNum = 'CO(10-9)'
nameFilesNum = 'CO_int_int.out'
nameFilesXHeaderNum = 0
iRaieNum = 9
nameFilesYHeaderNum = 1
figDataNum = {}
dataStartNum = 2
HeaderStartNum = None
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=260.25 #a remplir si les données Yang 2017 sont activées



"""Dénominateur"""
raieNomDen = 'CO(6-5)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0
iRaieDen = 5
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=433.54 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}

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
for i,directory in enumerate(dirFiles):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in vitDirFiles[i]:
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'--x',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])

"""========================================================================"""
"""Modeles de chocs J: input"""
originalDir=os.getcwd()
dirFiles = ['j-n4-b0.1','j-n5-b0.1','j-n6-b0.1']
labeldirFiles = {'j-n4-b0.1':r'J $n_H=10^4 cm^{-3}$ $b=0.1$'
                ,'j-n5-b0.1':r'J $n_H=10^5 cm^{-3}$ $b=0.1$'
                ,'j-n6-b0.1':r'J $n_H=10^6 cm^{-3}$ $b=0.1$'
                }
vitDirFiles = [[5,10,15,20,25,30],[5,10,15,20,25,30],[5,10,15,20,25,30]]
ActiverYang=True
ActiverMyObs=True

"""Numérateur"""
raieNomNum = 'CO(10-9)'
nameFilesNum = 'CO_int_int.out'
nameFilesXHeaderNum = 0
iRaieNum = 9
nameFilesYHeaderNum = 1
figDataNum = {}
dataStartNum = 2
HeaderStartNum = None
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=260.25 #a remplir si les données Yang 2017 sont activées



"""Dénominateur"""
raieNomDen = 'CO(6-5)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0
iRaieDen = 5
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=433.54 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}
"""
figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(14,14))
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
for i,directory in enumerate(dirFiles):
    X.append([])
    Y.append([])
    dirdir.append(directory)
    for vitesse in vitDirFiles[i]:
        X[i].append(vitesse)
        Y[i].append((Yratio[directory])[vitesse])
        
for n,graphe in enumerate(dirdir):
    if(labeldirFiles!=None): ax.plot(X[n],Y[n],'--d',label=labeldirFiles[graphe],markersize=7,color=labelTextColor[n])


"""========================================================================"""        
"""Observations hors Yang"""
TelescopeNameNum='HIFI'
TelescopeNameDen='APEX'
raieNomRatio = TelescopeNameNum+' '+raieNomNum+'/'+TelescopeNameDen+' '+raieNomDen
nameFilesRatio = "ratio_CO_same.dat"
iRatio = 5 #numero de ligne Ascii
jRatio = 1 #numero de colonne Ascii
HeaderStartRatio=0
DataStartRatio=0
xmin=np.amin(X)
xmax=np.amax(X)

if(ActiverMyObs==True):
    fig4data = ascii.read(nameFilesRatio)
    Header = (ascii.read(nameFilesRatio,header_start=HeaderStartRatio,data_start=DataStartRatio))
    Header = Header[0]


    
    ratioX = np.arange(xmin,xmax,0.1)
    ratioY = (fig4data[iRatio])[jRatio]#np.full(ratioX.shape,(fig4data[iRatio])[jRatio])
    #print(fig4data)
    Numerateur=(fig4data[iRatio])[-4]
    Denominateur=(fig4data[iRatio])[-3]    
    #ax.plot(ratioX,ratioY,'-',label=raieNomRatio,color='red')
    errNum=(fig4data[iRatio])[-2]
    errDen=(fig4data[iRatio])[-1]
    ratioYYang = ratioY  # pour convertir le rapport de flux en rapport de TdV
    errorRatio = ratioYYang*((errNum/Numerateur)+(errDen/Denominateur)) 
    
    ax.fill_between(ratioX,ratioYYang-errorRatio,ratioYYang+errorRatio,label=raieNomRatio,color='red',alpha=0.7)


"""========================================================================"""        
"""Observations  Yang (si activées)"""
ratioXYang = np.arange(xmin,xmax,0.1)
FluxNumYang=467. #1e-18 W/m2
FluxDenYang=279. #1e-18 W/m2
errorNumYang = 17.2 #1e-18 W/m2
errorDenYang = 4.07 #1e-18 W/m2
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

lgd = plt.legend(fontsize=8,loc='best',ncol=1)


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
