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
labelTextColor=colors

fig4dataoH2O = ascii.read("oH2O_flower12_CJ_chocs.dat")
HeaderoH2O = (  ascii.read("oH2O_flower12_CJ_chocs.dat",header_start=0,data_start=0))
HeaderoH2O = HeaderoH2O[0]
print(fig4dataoH2O)

fig4dataCO = ascii.read("oH2O_flower12_CJ_chocs.dat")
HeaderCO = (  ascii.read("oH2O_flower12_CJ_chocs.dat",header_start=0,data_start=0))
HeaderCO = HeaderCO[0]
print(fig4dataCO)


fig4dataH2OHIFI = ascii.read("Table2_Yang17_oH2O.txt")
HeaderH2OHIFI = (  ascii.read("Table2_Yang17_oH2O.txt",header_start=0,data_start=0))
HeaderH2OHIFI = HeaderH2OHIFI[0]
print(fig4dataH2OHIFI)

fig4dataCOHIFI = ascii.read("Table2_Yang17_oH2O.txt")
HeaderCOHIFI = (  ascii.read("Table2_Yang17_oH2O.txt",header_start=0,data_start=0))
HeaderCOHIFI = HeaderCOHIFI[0]
print(fig4dataCOHIFI)

fig, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
ax.grid(False)
fig.subplots_adjust(wspace=0)

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
labelRatio=r'PACS $\ce{H2O}$ $2_{12}-1_{01}$/$6_{16}-5_{05}$'
ax.set_xlabel(r'Age (ans)',fontsize=23,fontweight='extra bold')
ax.set_ylabel(labelRatio,fontsize=23,fontweight='extra bold')
#J_up=np.array([0,1,2,3,5,6,7])+2
xmin=50
xmax=500
"""v20"""
iraieHIFINumerateur=-6
iraieHIFINumerateur2=5
iraieHIFIDenominateur=2
iraieHIFIDenominateur2=5
ratioX = np.arange(xmin,xmax,0.1)
ilambdaNumerateur = 4
ilambdaDenominateur = 4
lambdaNumerateur = (fig4dataH2OHIFI[iraieHIFINumerateur])[ilambdaNumerateur] #microns
lambdaDenominateur = (fig4dataCOHIFI[iraieHIFIDenominateur])[ilambdaDenominateur] #microns

#ratioYY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])
ratioYY = (fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2]
correction = ((lambdaNumerateur/lambdaDenominateur)**3)
erroraa = np.float((((fig4dataH2OHIFI[iraieHIFINumerateur])[-1]).replace('[','')).replace(']',''))#*((lambdaNumerateur/lambdaDenominateur)**3)
errorbb = np.float((((fig4dataCOHIFI[iraieHIFIDenominateur])[-1]).replace('[','')).replace(']',''))#*((lambdaNumerateur/lambdaDenominateur)**3)
errorcc = (erroraa/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])*correction+ratioYY*(errorbb/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])*correction
ratioY = np.array(ratioYY)*correction # pour convertir le rapport de flux en rapport de TdV
#ratioY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])

ax.fill_between(ratioX,ratioY-errorcc,ratioY+errorcc,label=labelRatio,color=labelTextColor[0],alpha=0.7)
#ax.plot(ratioX,ratioY,'--',label=labelRatio,color=labelTextColor[0])

#['red','blue','green','brown','black']
Size=np.size(fig4dataoH2O)
J_up = []
figIndex=[250,500]

for levels in figIndex:
    J_up.append(levels)


#J_up.append(fig4data[11][1]) #pour S(7)
#labels=[]
#Ytot=[]
YH2O=[]
iraieH2O=-6
for ii,numeroOfModels in enumerate(range(7,9)):
#for ii,numeroOfModels in enumerate(range(5,7)):
#    X=[]
    
    #Ytot.append([])
    #labels.append(HeaderoH2O[numeroOfModels])
    #labels[ii]=labels[ii].replace('_',' ')
    #labels[ii]=labels[ii].replace('yr','ans')
    YH2O.append(fig4dataoH2O[iraieH2O][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)

labels=[]
#Ytot=[]
YCO=[]
iraieCO=16
#for ii,numeroOfModels in enumerate(range(5,7)):
for ii,numeroOfModels in enumerate(range(7,9)):    
    #Ytot.append([])
    labels.append(HeaderCO[numeroOfModels])
    labels[ii]=labels[ii].replace('_',' ')
    labels[ii]=labels[ii].replace('yr','ans')
    print(fig4dataCO[iraieCO][numeroOfModels])
    YCO.append(fig4dataCO[iraieCO][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)    

ModeleRatio = np.array(YH2O)/np.array(YCO)


ax.plot(J_up,ModeleRatio,'-o',label=labels[0][:11],color=labelTextColor[1],markersize=5)    




"""v25"""
#iraieHIFINumerateur=-1
#iraieHIFINumerateur2=3
#iraieHIFIDenominateur=-2
#iraieHIFIDenominateur2=5
#ratioX = np.arange(xmin,xmax,0.1)
#ratioY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])

#ax.plot(ratioX,ratioY,'--',label=labelRatio,color=labelTextColor[0])

#['red','blue','green','brown','black']
Size=np.size(fig4dataoH2O)
J_up = []
figIndex=[100,125,150,200,250]

for levels in figIndex:
    J_up.append(levels)


#J_up.append(fig4data[11][1]) #pour S(7)
#labels=[]
#Ytot=[]
YH2O=[]
#for ii,numeroOfModels in enumerate(range(7,12)):
for ii,numeroOfModels in enumerate(range(9,14)):
#    X=[]
    
    #Ytot.append([])
    #labels.append(HeaderoH2O[numeroOfModels])
    #labels[ii]=labels[ii].replace('_',' ')
    #labels[ii]=labels[ii].replace('yr','ans')
    YH2O.append(fig4dataoH2O[iraieH2O][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)

labels=[]
#Ytot=[]
YCO=[]
#for ii,numeroOfModels in enumerate(range(7,12)):
for ii,numeroOfModels in enumerate(range(9,14)):
    #Ytot.append([])
    labels.append(HeaderCO[numeroOfModels])
    labels[ii]=labels[ii].replace('_',' ')
    labels[ii]=labels[ii].replace('yr','ans')
    print(fig4dataCO[iraieCO][numeroOfModels])
    YCO.append(fig4dataCO[iraieCO][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)    

ModeleRatio = np.array(YH2O)/np.array(YCO)


ax.plot(J_up,ModeleRatio,'-o',label=labels[0][:11],color=labelTextColor[2],markersize=5)    


"""v30"""
#iraieHIFINumerateur=-1
#iraieHIFINumerateur2=3
#iraieHIFIDenominateur=-2
#iraieHIFIDenominateur2=5
#ratioX = np.arange(xmin,xmax,0.1)
#ratioY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])

#ax.plot(ratioX,ratioY,'--',label=labelRatio,color=labelTextColor[0])

#['red','blue','green','brown','black']
Size=np.size(fig4dataoH2O)
J_up = []
figIndex=[50,100,200]

for levels in figIndex:
    J_up.append(levels)


#J_up.append(fig4data[11][1]) #pour S(7)
#labels=[]
#Ytot=[]
YH2O=[]
for ii,numeroOfModels in enumerate(range(14,17)):
#for ii,numeroOfModels in enumerate(range(12,15)):
#    X=[]
    
    #Ytot.append([])
    #labels.append(HeaderoH2O[numeroOfModels])
    #labels[ii]=labels[ii].replace('_',' ')
    #labels[ii]=labels[ii].replace('yr','ans')
    YH2O.append(fig4dataoH2O[iraieH2O][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)

labels=[]
#Ytot=[]
YCO=[]
for ii,numeroOfModels in enumerate(range(14,17)):
#for ii,numeroOfModels in enumerate(range(12,15)):
    #Ytot.append([])
    labels.append(HeaderCO[numeroOfModels])
    labels[ii]=labels[ii].replace('_',' ')
    labels[ii]=labels[ii].replace('yr','ans')
    print(fig4dataCO[iraieCO][numeroOfModels])
    YCO.append(fig4dataCO[iraieCO][numeroOfModels])
    #print(Y)
    #ax.plot(J_up,Y,'-o',alpha=0.75,label=labels[ii],color=labelTextColor[ii],markersize=5)    

ModeleRatio = np.array(YH2O)/np.array(YCO)


ax.plot(J_up,ModeleRatio,'-o',label=labels[0][:11],color=labelTextColor[3],markersize=5)    

    
"""L1157"""
file1=ascii.read("K5_oH2O_TdV_HIFI.txt")
file2=ascii.read("K5_oH2O_TdV_HIFI.txt")
i1=12
i2=9
i1b=7
i2b=7
ax.plot(500,(file1[i1])[i1b]/(file2[i2])[i2b],'-o',label="C-J v20n4b1",color=labelTextColor[4],markersize=5)    

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







lgd = ax.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
#labelTextColor.append('black')
#for ii,text in enumerate(lgd.get_texts()):
#    text.set_color(labelTextColor[ii])




plt.savefig("H2O_to_H2O_2_12-1_01_to_6_16-5_05_CJ.png", bbox_inches='tight')
plt.savefig("H2O_to_H2O_2_12-1_01_to_6_16-5_05_CJ.pdf", bbox_inches='tight')
plt.show()    
plt.close()

