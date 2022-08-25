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

"""Modeles de chocs CJ: input"""
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

figExcit = []

#figExcit, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
#ax.grid(False)
#figExcit.subplots_adjust(wspace=0)

""" Sauvegarde des TdV """

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
                    (((figData[directory])[dirdir])[vt])['Y'].append([row[Yheader] for row in asciiTemp])
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
axTable = []

""" Plot des TdV """

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
                        figExcit.append([])
                        axTable.append([[],[],[],[]])
                        figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                        figExcit[-1].subplots_adjust(wspace=0) 
                    #(((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                    X=(((figData[directory])[dirdir])[vt])['X']
                    Y=((((figData[directory])[dirdir])[vt])['Y'])[m]
                    axTable[-1][iAx].plot(X,Y,'-o',markersize=3,label=dirdir+vt,alpha=0.65,linewidth=2)
                    iAx2=iAx2+1
                    if(iAx2%4==0):
                        lgd = axTable[-1][iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
                        iAx=iAx+1
                    if(iAx2%16==0):
                        os.chdir(originalDir)
                        #figExcit[-1].savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
                        #plt.show()
                        #plt.close()                        
                        os.chdir(currentDirLvL3)
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)


"""Observables: input"""
originalDir=os.getcwd()
dirFilesObs = ['./']
nameFilesObs = ['Table2_yang17_CO.txt']
LabelObs = ['PACS+SPIRE']
nameFilesXHeaderObs = [0]
nameFilesYHeaderObs = [[4]]
figObs = {}
dataStartObs = [1]
HeaderStartObs = [0] 
"""Corrections"""
FluxToTdV=[True]
iLambdaFluxToTdV=[3]
theta_SFluxToTdV=[10]

dirlist = []
vitesseTime = []

"""Sauvegarde des Observables"""

for i,directory in enumerate(dirFilesObs):       
    figObs[directory]={}
    os.chdir(directory)
    for h,nf in enumerate(nameFilesObs):            
        (figObs[directory])[nf]={}
        if os.path.exists(nameFilesObs[h]): 
            asciiTemp=ascii.read(nameFilesObs[h],header_start=HeaderStartObs[h],data_start=dataStartObs[h])
            ((figObs[directory])[nf])['X']=[row[nameFilesXHeaderObs[h]] for row in asciiTemp]
            ((figObs[directory])[nf])['Y']=[]
            if(FluxToTdV[i]==True) : 
                ((figObs[directory])[nf])['lambda_microns']=[]
            for m,Yheader in enumerate(nameFilesYHeaderObs[h]):
                #print(Yheader)
                #if(iAx2%16==0):
                #    iAx=0
                #    figExcit, [axTable[0],axTable[1],axTable[2],axTable[3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                #    figExcit.subplots_adjust(wspace=0) 
                ((figObs[directory])[nf])['Y'].append([row[Yheader] for row in asciiTemp])
                if(FluxToTdV[h]==True) : 
                    ((figObs[directory])[nf])['lambda_microns'].append([row[iLambdaFluxToTdV[i]] for row in asciiTemp])
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
    os.chdir(originalDir)

"""Plot des Observables"""
iAx = 0
iAx2 = 0
iAx3 = 0
iFig=0

for n,directoryObs in enumerate(nameFilesObs):       
    iAx = 0
    iAx2 = 0
    iFig=0
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
                if os.path.exists(nameFiles[n]): 
                    #asciiTemp=ascii.read(nameFiles[i],header_start=HeaderStart[i],data_start=dataStart[i])
                    #(((figData[directory])[dirdir])[vt])['X']=[row[nameFilesXHeader[i]] for row in asciiTemp]
                    #(((figData[directory])[dirdir])[vt])['Y']=[]
                    for m,Yheader in enumerate(nameFilesYHeader[n]):
                        #print(Yheader)
                        #if(iAx2%16==0):
                        #    iAx=0
                        #    figExcit.append([])
                        #    axTable.append([[],[],[],[]])
                        #    figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                        #    figExcit[-1].subplots_adjust(wspace=0) 
                        #(((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                        X=(figObs[directoryObs])['X']
                        Y=((figObs[directoryObs])['Y'])[m]
                        if(FluxToTdV[n]==True) : 
                            (figObs[directory])['lambda_microns'].append([row[iLambdaFluxToTdV[i]] for row in asciiTemp])
                        #axTable[iFig][iAx].scatter(np.array(X),np.array(Y),marker='D')#,label=LabelObs[i])
                        iAx2=iAx2+1
                        if(iAx2%4==0):
                            axTable[iFig][iAx3].scatter(np.array(X),np.array(Y),marker='D',label=LabelObs[n])#,label=LabelObs[i])
                            lgd = axTable[iFig][iAx3].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
                            iAx=iAx+1
                            iAx3=iAx3+1
                        if(iAx2%16==0):
                            os.chdir(originalDir)
                            #figExcit[iFig].savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
                            #plt.show()
                            #plt.close()                        
                            os.chdir(currentDirLvL3)
                            iFig=iFig+1
                            iAx3=0
                os.chdir(currentDirLvL2)
            os.chdir(currentDirLvL1)
        os.chdir(originalDir)
if(iAx2%16!=0):
    axTable[iFig][iAx3].scatter(np.array(X),np.array(Y),marker='D',label=LabelObs[n])#,label=LabelObs[i])
    #lgd = axTable[-1][iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
    #figExcit[-1].savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
    #plt.show()
    #plt.close()
    
""" Enregistrement des figures """
iAx = 0
iAx2 = 0
iFig=0
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
                    #    figExcit.append([])
                    #    axTable.append([[],[],[],[]])
                    #    figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                    #    figExcit[-1].subplots_adjust(wspace=0) 
                    #(((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                    #X=(((figData[directory])[dirdir])[vt])['X']
                    #Y=(((figData[directory])[dirdir])[vt])['Y']
                    #axTable[-1][iAx].plot(X,Y,label=dirdir+vt,alpha=0.65,linewidth=linewdthsize[iAx2%4])
                    iAx2=iAx2+1
                    if(iAx2%4==0):
                        #lgd = axTable[-1][iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
                        iAx=iAx+1
                    if(iAx2%16==0):
                        os.chdir(originalDir)
                        figExcit[iFig].savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
                        plt.show()
                        plt.close()                        
                        os.chdir(currentDirLvL3)
                        iFig=iFig+1
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)

if(iAx2%16!=0):
    lgd = axTable[-1][iAx].legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
    figExcit[-1].savefig("CO_TdV_{}.pdf".format(iAx2), bbox_inches='tight')
    plt.show()
    plt.close()                        