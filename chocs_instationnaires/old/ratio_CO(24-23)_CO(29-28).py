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
dirFiles = ['cj-n4-b1','cj-n5-b1','cj-n6-b1']
labeldirFiles = {}
vitDirFiles = [[5,10,19],[5,10,19],[5,10,19]]
ageDirFiles = [[[50,100,200,500],
                [50,100,200,500],
                [50,100,200,500]],
                [[50,100,150,200,250,500],
                 [50,100,150,200,250,500],
                 [50,100,150,200,250,500]],
                 [[50,100,150,200,250,500],
                 [50,100,150,200,250,500],
                 [50,100,150,200,250,500]]]

axTableXScale='linear'
axTableYScale='log'
ActiverYang=True
ActiverMyObs=False
"""Numérateur"""
Eup1=1656.47
raieNomNum = 'CO(24-23) ('+str(Eup1)+' K)'
nameFilesNum = 'CO_int_int.out'
nameFilesXHeaderNum = 0
iRaieNum = 23
nameFilesYHeaderNum = 1
figDataNum = {}
dataStartNum = 2
HeaderStartNum = None
dirlistNum = []
vitesseTimeNum = []
lambdaNumMicrons=108.75 #a remplir si les données Yang 2017 sont activées



"""Dénominateur"""
Eup2=2399.82
raieNomDen = 'CO(29-28) ('+str(Eup2)+' K)'
nameFilesDen = 'CO_int_int.out'
nameFilesXHeaderDen = 0
iRaieDen = 28
nameFilesYHeaderDen = 1
figDataDen = {}
dataStartDen = 2
HeaderStartDen = None
dirlistDen = []
vitesseTimeDen = []
lambdaDenMicrons=90.16 #a remplir si les données Yang 2017 sont activées

"""Ratio"""
Yratio = {}

figRatio, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))
figRatio.subplots_adjust(wspace=0) 
ax.grid(True)
#ax.set_ylim(1,1e3)
ax.set_xscale(axTableXScale)
ax.set_yscale(axTableYScale)
ax.set_ylabel(r'$\sum T_{mb}\Delta v$ '+r'{:1}/{:2}'.format(raieNomNum,raieNomDen))
ax.set_xlabel(r'Age (ans)')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
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
        #(labeldirFiles[directory])[vitesse]={}
        """Label"""
        ((labeldirFiles[directory])[vitesse])=directory+"-v"+str(vitesse)
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            vitDirName="v{:1}-t{:2}".format(str(vitesse),str(tempsJ))
            os.chdir(vitDirName)
            currentDirLvl2=os.getcwd()
            """Numerateur"""
            asciiTemp = ascii.read(nameFilesNum,header_start=HeaderStartNum,data_start=dataStartNum)
            ((figDataNum[directory])[vitesse])[tempsJ]=(asciiTemp[iRaieNum])[nameFilesYHeaderNum]
            """Dénominateur"""
            asciiTemp2 = ascii.read(nameFilesDen,header_start=HeaderStartDen,data_start=dataStartDen)
            ((figDataDen[directory])[vitesse])[tempsJ]=(asciiTemp2[iRaieDen])[nameFilesYHeaderDen]
            """Ratio"""
            ((Yratio[directory])[vitesse])[tempsJ]=((figDataNum[directory])[vitesse])[tempsJ]/((figDataDen[directory])[vitesse])[tempsJ]
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
        for k,tempsJ in enumerate((ageDirFiles[i])[j]):
            X.append(tempsJ)
            Y.append(((Yratio[directory])[vitesse])[tempsJ])
        ax.plot(X,Y,'-o',label=((labeldirFiles[directory])[vitesse]),markersize=10,color=labelTextColor[l])
        l=l+1

"""========================================================================"""        
"""Modeles PdF IRAS4B"""   
Age=np.array([[250,500],[100,125,150,200,250],[50,100,200]])
    
"""Numérateur"""
nameFilesNumPdF = 'CO_IRSA4B_CJ_chocs.dat'
nameFilesXHeaderNumPdF = 0
iRaieNumPdF = 23
nameFilesYHeaderNumPdF = [range(5,7),range(7,12),range(12,15)]
figDataNumPdF = {}
dataStartNumPdF = 1
HeaderStartNumPdF = 0
#dirlistNumPdF = []
#vitesseTimeNumPdF = []


"""Dénominateur"""
nameFilesDenPdF = 'CO_IRSA4B_CJ_chocs.dat'
nameFilesXHeaderDenPdF = 0
iRaieDenPdF = 28
nameFilesYHeaderDenPdF = [range(5,7),range(7,12),range(12,15)]
figDataDenPdF = {}
dataStartDenPdF = 1
HeaderStartDenPdF = 0
#dirlistDenPdF = []
#vitesseTimeDenPdF = []

"""Ratio"""
YratioPdF = {}

"""Sauvegarde des TdV"""
figDataNumPdF = ascii.read(nameFilesNumPdF)
HeaderNumPdF = (  ascii.read(nameFilesNumPdF,header_start=HeaderStartNumPdF,data_start=HeaderStartNumPdF))
HeaderNumPdF = HeaderNumPdF[0]
#print(figDataNumPdF)

figDataDenPdF = ascii.read(nameFilesDenPdF)
HeaderDenPdF = (  ascii.read(nameFilesDenPdF,header_start=HeaderStartDenPdF,data_start=HeaderStartNumPdF))
HeaderDenPdF = HeaderDenPdF[0]
#print(figDataDenPdF)


YNum=[]
for i,vitesseRang in enumerate(nameFilesYHeaderNumPdF):
    YNum.append([])
    for j,numeroOfModelsNum in enumerate(vitesseRang):
    #    X=[]
        
        #Ytot.append([])
        #labels.append(HeaderNumPdF[numeroOfModels])
        #labels[ii]=labels[ii].replace('_',' ')
        #labels[ii]=labels[ii].replace('yr','ans')
        YNum[i].append(figDataNumPdF[iRaieNumPdF][numeroOfModelsNum])
        #print(Y)
        #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)
labels=[]
YDen=[]
for i,vitesseRang2 in enumerate(nameFilesYHeaderDenPdF):
    YDen.append([])
    labels.append(HeaderDenPdF[vitesseRang2[0]])
    labels[i]=labels[i].replace('_',' ')
    labels[i]=labels[i].replace('yr','ans')
    for k,numeroOfModelsDen in enumerate(vitesseRang2):
        #Ytot.append([])

        #print(fig4dataCO[iraieCO][numeroOfModels])
        YDen[i].append(figDataDenPdF[iRaieDenPdF][numeroOfModelsDen])
        #print(Y)
        #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)    

ModeleRatio = []#np.array(YNum)/np.array(YDen)
for i,vitesseRang in enumerate(nameFilesYHeaderNumPdF):
    ModeleRatio.append([])
    #for j,numeroOfModelsNum in enumerate(vitesseRang):
    #    X=[]        
    #Ytot.append([])
    #labels.append(HeaderNumPdF[numeroOfModels])
    #labels[ii]=labels[ii].replace('_',' ')
    #labels[ii]=labels[ii].replace('yr','ans')
    ModeleRatio[i].append(np.array(YNum[i])/np.array(YDen[i]))
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)
    ax.plot(Age[i],ModeleRatio[i][0],'--o',label=labels[i][:11]+" (PdF13)",markersize=5)    

   

"""========================================================================"""        
xmin=50
xmax=500
"""Observations  Yang (si activées)"""
ratioXYang = np.arange(xmin,xmax,0.1)
FluxDenYang=165. #1e-18 W/m2
FluxNumYang=508. #1e-18 W/m2
errorDenYang = 23.6 #1e-18 W/m2
errorNumYang = 43.9 #1e-18 W/m2
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
    
lgd = plt.legend(fontsize=14,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.45))

"""Enregistrement des figures"""
fileTitle="ratio_"+"{}".format(raieNomNum)+"_"+"{}".format(raieNomDen)
fileTitle=fileTitle.replace('\ce{','')
fileTitle=fileTitle.replace('}','')
fileTitle=fileTitle.replace('{','')
fileTitle=fileTitle.replace('$','')
fileTitle=fileTitle.replace(' ','_')
fileTitle=fileTitle.replace('.','v')
plt.savefig(fileTitle+".pdf", bbox_inches='tight')    
plt.savefig(fileTitle+".png", bbox_inches='tight')    
plt.show()
plt.close()      
