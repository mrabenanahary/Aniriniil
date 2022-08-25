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



os.chdir('./chocs_cj_18_06_19/cj-n5-b1')
modelesCJ=['v19-t50','v19-t100','v19-t250']
ligneCOaextraire=15
lambdaCOmicrons=162.82
ligneOIaextraire=[2117,2275]
lambdaOImicrons=63.18
labelRatio=r'OI 63$\mu m/CO(16-15)$'
age=np.array([50,100,250])
line1=[]
line2=[]


fig, ax = plt.subplots(1,1, sharey=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
ax.grid(False)
fig.subplots_adjust(wspace=0)
ax.set_yscale('log')
ax.set_xscale('linear')
ax.set_xlabel(r'Age (ans)',fontsize=23,fontweight='extra bold')
ax.set_ylabel(labelRatio,fontsize=23,fontweight='extra bold')

"""Yang17"""
#ratioYY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])
xmin=50
xmax=500
"""v20"""
ratioX = np.arange(xmin,xmax,0.1)
CO1615Yang=667. #1e-18 W/m2
OI63Yang=2009. #1e-18 W/m2
errorOI = 74.0 #1e-18 W/m2
errorCO = 23.9 #1e-18 W/m2
lambdaNumerateur=lambdaOImicrons #MICRONS
lambdaDenominateur=lambdaCOmicrons
Numerateur=OI63Yang#CO1615Yang
Denominateur=CO1615Yang#OI63Yang
ratioYY = Numerateur/Denominateur
correction = ((lambdaNumerateur/lambdaDenominateur)**3) #rapport des TdV
errNum=errorOI
errDen=errorCO
ratioY = np.array(ratioYY)*correction # pour convertir le rapport de flux en rapport de TdV
errorRatio = ratioY*((errNum/Numerateur)+(errDen/Denominateur)) 
#(errNum/Denominateur)*correction+ratioYY*(errorbb/Denominateur)*correction

#ratioY = np.full(ratioX.shape,(fig4dataH2OHIFI[iraieHIFINumerateur])[iraieHIFINumerateur2]/(fig4dataCOHIFI[iraieHIFIDenominateur])[iraieHIFIDenominateur2])

ax.fill_between(ratioX,ratioY-errorRatio,ratioY+errorRatio,label=labelRatio,color=labelTextColor[0],alpha=0.7)


for n,direct in enumerate(modelesCJ):
    os.chdir(direct)
    """CO"""
    fig4dataCO = ascii.read("CO_int_int.out",header_start=None,data_start=2)
    line1.append((fig4dataCO[ligneCOaextraire])[1])
    print(line1)
    """OI"""
    """
    fig4dataOI63 = ascii.read("intensity.out",header_start=0,data_start=1)
    Headerline1 = ascii.read("intensity.out",header_start=0,data_start=0)
    Headerline1=Headerline1[0]
    line2.append((fig4dataOI63[ligneOIaextraire[n]])[18])
    print(line2)
    """
    os.chdir('..')

os.chdir('..')
"""OI"""   
line2=[4.454106484E-006,4.135676163E-006,4.152904260E-006] #erg/s/cm2/sr from intensity.out


line2=np.array(line2)
line2=1.e-12*line2*lambdaOImicrons**3/(2.*KB)/1.e5 #K.km/s

ax.plot(age,line2/np.array(line1),'-o',label='cj-n5-b1'+(modelesCJ[0])[0:3])
lgd = ax.legend(fontsize=12,loc='lower center',ncol=2,bbox_to_anchor=(0.5, -0.32))
plt.savefig("CO_1615_to_OH_63microns.png", bbox_inches='tight')
plt.savefig("CO_1615_to_OH_63microns.pdf", bbox_inches='tight')
plt.show()
plt.close()




