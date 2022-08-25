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







mydir='cj-n5-b1'
myvittmps='v*-t*'


prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
labelTextColor= colors

"""Modeles de chocs CJ: input"""
originalDir=os.getcwd()
dirFiles = ['./']
nameFiles = ['CO_int_int.out']
nameFilesXHeader = [0]
nameFilesYHeader = [[1]]
figData = {}
dataStart = [2]
HeaderStart = [None] 
dirlist = []
vitesseTime = []

ShareY=False
axTableXScale='linear'
axTableYScale='log'
linewdthsize=[4,3,2,1]
axTableXLabel=r'$J_{up}$'
axTableYLabel=r'$\sum T_{mb}\Delta v$ CO (K.km/s)'
figExcit = []

#figExcit, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
#ax.grid(False)
#figExcit.subplots_adjust(wspace=0)

""" Sauvegarde des TdV """

for i,directory in enumerate(dirFiles):       
    figData[directory]={}
    os.chdir(directory)
    dirlist.append(os.popen('ls -1d '+mydir).read().split())
    currentDirLvL1=os.getcwd()
    for k,dirdir in enumerate(dirlist[i]):
        (figData[directory])[dirdir]={}
        os.chdir(dirdir)
        vitesseTime.append(os.popen('ls -1d '+myvittmps).read().split())
        currentDirLvL2=os.getcwd()
        for l,vt in enumerate(vitesseTime[k]):
            ((figData[directory])[dirdir])[vt]={}
            os.chdir(vt)
            if os.path.exists(nameFiles[i]): 
                asciiTemp=ascii.read(nameFiles[i],header_start=HeaderStart[i],data_start=dataStart[i])
                (((figData[directory])[dirdir])[vt])['X']=[row[nameFilesXHeader[i]] for row in asciiTemp]
                (((figData[directory])[dirdir])[vt])['Y']=[]
                for m,Yheader in enumerate(nameFilesYHeader[i]):
                    (((figData[directory])[dirdir])[vt])['Y'].append([row[Yheader] for row in asciiTemp])                      
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)



iAx = 0
iAx2 = 0
axTable = []


""" Plot des TdV """

for i,directory in enumerate(dirFiles):       
    os.chdir(directory)
    dirlist.append(os.popen('ls -1d cj-*b1').read().split())
    currentDirLvL1=os.getcwd()
    for k,dirdir in enumerate(dirlist[i]):
        os.chdir(dirdir)
        vitesseTime.append(os.popen('ls -1d '+myvittmps).read().split())
        currentDirLvL2=os.getcwd()
        for l,vt in enumerate(vitesseTime[k]):
            os.chdir(vt)
            currentDirLvL3=os.getcwd()
            if os.path.exists(nameFiles[i]): 
                for m,Yheader in enumerate(nameFilesYHeader[i]):
                    if(iAx2%20==0):
                        iAx=0
                        figExcit.append([])
                        axTable.append([[],[],[],[]])
                        figExcit[-1], [axTable[-1][0],axTable[-1][1]] = plt.subplots(1,2, sharey=ShareY,figsize=(10,4))
                        axTable[-1][0].set_ylabel(axTableYLabel)
                        axTable[-1][0].set_xlabel(axTableXLabel)
                        axTable[-1][1].set_xlabel(axTableXLabel)
                        #axTable[-1][2].set_xlabel(axTableXLabel)
                        #axTable[-1][3].set_xlabel(axTableXLabel)
                        if(ShareY): figExcit[-1].subplots_adjust(wspace=0) 
                        for element in [axTable[-1][0],axTable[-1][1]]:
                            element.grid(True)
                            element.set_yscale(axTableYScale)
                            element.set_xscale(axTableXScale)
                            #element.set_ylim(1e-4,1e3)
                        axTable[-1][0].set_xscale(axTableXScale)
                        axTable[-1][0].set_yscale(axTableYScale)
                    X=(((figData[directory])[dirdir])[vt])['X']
                    Y=((((figData[directory])[dirdir])[vt])['Y'])[m]
                    axTable[-1][iAx].plot(X,Y,'-o',markersize=2,label=dirdir+'-'+vt,alpha=0.65,linewidth=1.0)
                    iAx2=iAx2+1
                    if(iAx2%10==0):
                        lgd = axTable[-1][iAx].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
                        iAx=iAx+1
                    if(iAx2%20==0):
                        os.chdir(originalDir)          
                        os.chdir(currentDirLvL3)
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)








"""Modeles de chocs CJ: PdF"""

""" Sauvegarde des TdV """

File='CO_IRSA4B_CJ_chocs.dat'
DirName='CO_IRSA4B_CJ_chocs.dat'
figData[DirName]={}
asciiTemp=ascii.read(File,header_start=0,data_start=1)
HeaderasciiTemp=ascii.read(File,header_start=0,data_start=0)
HeaderasciiTemp=HeaderasciiTemp[0]
(figData[DirName])['X']=[row[0] for row in asciiTemp]
(figData[DirName])['Y']=[]
(figData[DirName])['labels']=[]
for m,Yheader in enumerate(range(5,14)):
    (figData[DirName])['Y'].append([row[Yheader] for row in asciiTemp])
    (figData[DirName])['labels'].append((HeaderasciiTemp[Yheader]).replace('_',' '))
    



iAx = 0
iAx2 = 0
iAx3 = 0
iFig=0
#axTable = []

""" Plot des TdV """

for i,directory in enumerate(dirFiles):       
    os.chdir(directory)
    dirlist.append(os.popen('ls -1d '+mydir).read().split())
    currentDirLvL1=os.getcwd()
    for k,dirdir in enumerate(dirlist[i]):
        os.chdir(dirdir)
        vitesseTime.append(os.popen('ls -1d '+myvittmps).read().split())
        currentDirLvL2=os.getcwd()
        for l,vt in enumerate(vitesseTime[k]):
            os.chdir(vt)
            currentDirLvL3=os.getcwd()
            if os.path.exists(nameFiles[i]): 
                for m,Yheader in enumerate(nameFilesYHeader[i]):
                    if(iAx2%20==0):
                        iAx=0
                    #    figExcit.append([])
                    #    axTable.append([[],[],[],[]])
                    #    figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                    #    figExcit[-1].subplots_adjust(wspace=0) 
                    #X=((figObs[directoryobs])[nf])['X']
                    #Y=(((figObs[directoryobs])[nf])['Y'])[m]
                    iAx2=iAx2+1
                    if(iAx2%10==0): 
                        for n,graphe in enumerate((figData[DirName])['Y']):
                            X=(figData[DirName])['X']
                            Y=((figData[DirName])['Y'])[n]
                            #axTable[iFig][iAx3].plot(X,Y,'--D',label=r'(PdF) {}'.format(((figData[DirName])['labels'])[n]),linewidth=1,markersize=0.5,color=labelTextColor[n])
                        #lgd = axTable[iFig][iAx3].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
                        iAx=iAx+1
                        iAx3=iAx3+1
                    if(iAx2%20==0):
                        os.chdir(originalDir)          
                        os.chdir(currentDirLvL3)
                        iFig=iFig+1
                        iAx3=0
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)

if(iAx2%20!=0):
    for n,graphe in enumerate((figData[DirName])['Y']):
        X=(figData[DirName])['X']
        Y=((figData[DirName])['Y'])[n]
        #axTable[iFig][iAx3].plot(X,Y,'--D',label=r'(PdF) {}'.format(((figData[DirName])['labels'])[n]),linewidth=1,markersize=0.5,color=labelTextColor[n])
    #lgd = axTable[iFig][iAx3].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
























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

Xobs=np.array([3,
               4,
               6,
               7]
               )
Yobs=np.array([195.3,
               130.79,
               146.42,
               77.32])
ErrorBar=np.multiply(Yobs,np.array([0.16726027621644057,0.16815469068687913,0.17167483799322486,0.27837070607375336]))
#ax.scatter(Xobs,Yobs,s=10,zorder=6,marker='D',edgecolors='black',facecolors='yellow',label=r'Obs. APEX')
#ax.errorbar(Xobs,Yobs, yerr=ErrorBar,fmt='none',capsize=4,capthick=2,elinewidth=2,color='black')

"""Sauvegarde des Observables"""

for d,directory in enumerate(dirFilesObs):
    os.chdir(directory)
    figObs[directory]={}
    for n,nf in enumerate(nameFilesObs):
        if os.path.exists(nf): 
            (figObs[directory])[nf]={}
            asciiTemp=ascii.read(nf,header_start=HeaderStart[n],data_start=dataStart[n])
            ((figObs[directory])[nf])['X']=[row[nameFilesXHeaderObs[n]] for row in asciiTemp]
            ((figObs[directory])[nf])['Y']=[]
            ((figObs[directory])[nf])['Yerr']=[]
            if(FluxToTdV[n]==True) : ((figObs[directory])[nf])['lambda_micron']=[row[iLambdaFluxToTdV[n]] for row in asciiTemp]
            for m,Yheader in enumerate(nameFilesYHeaderObs[n]):
                ((figObs[directory])[nf])['Y'].append([row[Yheader] for row in asciiTemp])
                ((figObs[directory])[nf])['Yerr'].append([np.float((row[-1].replace('[','')).replace(']','')) for row in asciiTemp])
                if(FluxToTdV[n]==True) :
                    #(((figObs[directory])[nf])['Y'])[0]=((1.5e-6)/(theta_SFluxToTdV[n]**2))*np.array((((figObs[directory])[nf])['Y'])[0])*(np.array(((figObs[directory])[nf])['lambda_micron'])**3)
                    #print((((figObs[directory])[nf])['Y'])[0])
                    temp=[]
                    for iii,Flux1em18w_per_mm2 in enumerate((((figObs[directory])[nf])['Y'])[0]):
                        temp.append((Flux1em18w_per_mm2*(1.5e-6)*((np.array(((figObs[directory])[nf])['lambda_micron'][iii]))**3))/((theta_SFluxToTdV[n])**2))
                    #print(np.array(temp)-(((figObs[directory])[nf])['Y'])[0])
                    (((figObs[directory])[nf])['Y'])[0]=temp
                    temp=[]
                    for iii,Flux1em18w_per_mm2 in enumerate((((figObs[directory])[nf])['Yerr'])[0]):
                        temp.append((Flux1em18w_per_mm2*(1.5e-6)*((np.array(((figObs[directory])[nf])['lambda_micron'][iii]))**3))/((theta_SFluxToTdV[n])**2))
                    ((figObs[directory])[nf])['Yerr'][0]=temp
    os.chdir(originalDir)
    
""" Plot des TdV des observables """
iAx = 0
iAx2 = 0
iAx3 = 0
iFig=0
for d,directoryobs in enumerate(dirFilesObs):
    os.chdir(directoryobs)
    for n,nf in enumerate(nameFilesObs):
        if os.path.exists(nf): 
            for i,directory in enumerate(dirFiles):       
                os.chdir(directory)
                dirlist.append(os.popen('ls -1d '+mydir).read().split())
                currentDirLvL1=os.getcwd()
                for k,dirdir in enumerate(dirlist[i]):
                    os.chdir(dirdir)
                    vitesseTime.append(os.popen('ls -1d '+myvittmps).read().split())
                    currentDirLvL2=os.getcwd()
                    for l,vt in enumerate(vitesseTime[k]):
                        os.chdir(vt)
                        currentDirLvL3=os.getcwd()
                        if os.path.exists(nameFiles[i]): 
                            for m,Yheader in enumerate(nameFilesYHeader[i]):
                                if(iAx2%20==0):
                                    iAx=0
                                #    figExcit.append([])
                                #    axTable.append([[],[],[],[]])
                                #    figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                                #    figExcit[-1].subplots_adjust(wspace=0) 
                                X=((figObs[directoryobs])[nf])['X']
                                Y=(((figObs[directoryobs])[nf])['Y'])[m]
                                Yerr=(((figObs[directoryobs])[nf])['Yerr'])[m]
                                iAx2=iAx2+1
                                if(iAx2%10==0): 
                                    axTable[iFig][iAx3].scatter(X,Y,s=15,zorder=3,marker='D',label=LabelObs[n],facecolors='white',edgecolors='black')
                                    axTable[iFig][iAx3].errorbar(X,Y, yerr=Yerr,fmt='none',capsize=4,capthick=2,elinewidth=1,color='black')
                                    axTable[iFig][iAx3].scatter(Xobs,Yobs,s=15,zorder=3,marker='D',label='Obs APEX.',facecolors='yellow',edgecolors='black')
                                    axTable[iFig][iAx3].errorbar(Xobs,Yobs, yerr=ErrorBar,fmt='none',capsize=4,capthick=2,elinewidth=1,color='blue')
                                    #axTable[iFig][iAx3].yaxis.set_minor_locator(AutoMinorLocator(10))
                                    axTable[iFig][iAx3].xaxis.set_minor_locator(AutoMinorLocator(5))
                                    lgd = axTable[iFig][iAx3].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
                                    iAx=iAx+1
                                    iAx3=iAx3+1
                                if(iAx2%20==0):
                                    os.chdir(originalDir)          
                                    os.chdir(currentDirLvL3)
                                    iFig=iFig+1
                                    iAx3=0
                        os.chdir(currentDirLvL2)
                    os.chdir(currentDirLvL1)
                os.chdir(originalDir)            
    os.chdir(originalDir)

if(iAx2%20!=0):
    axTable[iFig][iAx3].scatter(X,Y,s=15,zorder=3,marker='D',label=LabelObs[n],facecolors='white',edgecolors='black')
    axTable[iFig][iAx3].errorbar(X,Y, yerr=Yerr,fmt='none',capsize=4,capthick=2,elinewidth=1,color='black')
    axTable[iFig][iAx3].scatter(Xobs,Yobs,s=15,zorder=3,marker='s',label='Obs APEX.',facecolors='yellow',edgecolors='black')
    axTable[iFig][iAx3].errorbar(Xobs,Yobs, yerr=ErrorBar,fmt='none',capsize=4,capthick=2,elinewidth=1,color='blue')
    #axTable[iFig][iAx3+1].scatter(X,Y,s=15,zorder=3,marker='D',label=LabelObs[n],facecolors='white',edgecolors='black')
    #axTable[iFig][iAx3+1].errorbar(X,Y, yerr=Yerr,fmt='none',capsize=4,capthick=2,elinewidth=1,color='black')
    #axTable[iFig][iAx3+1].scatter(Xobs,Yobs,s=15,zorder=3,marker='s',label='Obs APEX.',facecolors='yellow',edgecolors='black')
    #axTable[iFig][iAx3+1].errorbar(Xobs,Yobs, yerr=ErrorBar,fmt='none',capsize=4,capthick=2,elinewidth=1,color='blue')
    #axTable[iFig][iAx3].yaxis.set_minor_locator(AutoMinorLocator(10))
    axTable[iFig][iAx3].xaxis.set_minor_locator(AutoMinorLocator(5))
    #axTable[iFig][iAx3+1].xaxis.set_minor_locator(AutoMinorLocator(5))
    #lgd = axTable[iFig][iAx3+1].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
    plt.show()
    plt.close()  


    
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
                    if(iAx2%20==0):
                        iAx=0
                    #    figExcit.append([])
                    #    axTable.append([[],[],[],[]])
                    #    figExcit[-1], [axTable[-1][0],axTable[-1][1],axTable[-1][2],axTable[-1][3]] = plt.subplots(1,4, sharey=True,figsize=(28,7))
                    #    figExcit[-1].subplots_adjust(wspace=0) 
                    #(((figData[directory])[dirdir])[vt])['Y']=[row[Yheader] for row in asciiTemp]
                    #X=(((figData[directory])[dirdir])[vt])['X']
                    #Y=(((figData[directory])[dirdir])[vt])['Y']
                    #axTable[-1][iAx].plot(X,Y,label=dirdir+vt,alpha=0.65,linewidth=linewdthsize[iAx2%10])
                    iAx2=iAx2+1
                    if(iAx2%10==0):
                        #lgd = axTable[-1][iAx].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
                        iAx=iAx+1
                    if(iAx2%20==0):
                        os.chdir(originalDir)
                        plt.close(figExcit[iFig-1])
                        figExcit[iFig].savefig("CO_TdV_{}_noPdF.pdf".format(iAx2), bbox_inches='tight')
                        plt.show()
                        plt.close()                        
                        os.chdir(currentDirLvL3)
                        iFig=iFig+1
            os.chdir(currentDirLvL2)
        os.chdir(currentDirLvL1)
    os.chdir(originalDir)
    
    

if(iAx2%20!=0):
    lgd = axTable[iFig][iAx].legend(fontsize=10,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.6))
    figExcit[iFig].savefig("CO_TdV_{}_noPdF.pdf".format(iAx2), bbox_inches='tight')
    plt.show()
    plt.close()                        