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
    r'\usepackage[french]{babel}',
    r'\usepackage[utf8]{inputenc}',
    r'\usepackage{textcomp}',
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

if os.path.exists('./H2_results.txt') : os.remove('./H2_results.txt')
#filenames=["bhr71-h2-flux_ellipse_2v5x5arcsecond_1st_position_sud.dat",
#           "bhr71-h2-flux_circle_5arcsec_sud_2e_position.dat",
filenames=[           "bhr71-h2-flux.dat"]#,
#           "bhr71-h2-flux_circle-petit_lobe_5arcs.dat"]

#FileLabel=[r'$a\times b=2.5{}^\prime{}^\prime \times 5{}^\prime{}^\prime$ pos 1',
          # r'$r=2.5{}^\prime{}^\prime$ pos 2',
#FileLabel=[r'After correction by'+'\n'+r'the small-to-larger $d_{beam}=5\arcsecond\longrightarrow 20\arcsecond$'+'\n'+r'beam filling factor $\simeq 5.5\times 10^{-2} $']#,
FileLabel=[r'Observations']#,          
           #r'$r=5{}^\prime{}^\prime$ pos 1']
           #           r'$\begin{cases} d_{\text{small lobe}}=5{}^\prime{}^\prime\\ d_{\text{main lobe}}=20{}^\prime{}^\prime \\ \Longleftrightarrow$ Filling Factor = $5.4788\times 10^{-2} \end{cases}$' 
#FileLabel=[r'$d_{\text{small lobe}}=5{}^\prime{}^\prime$'+'\n'+r'$d_{\text{main lobe}}=20{}^\prime{}^\prime$ $\Longleftrightarrow$ Filling Factor = $5.4788\times 10^{-2}$']#,
#FileColors=["blue",
#           "orange",
#           "green",
FileColors=[           "red"]

#FileMarkers=["v",
#           "D",
#           "s",
FileMarkers=[           "X"]

#FileFF=[0.098195,
#        0.054788,
FileFF=[        0.054788]#,

fitResults=[]
#        0.25111]

"""====== AFFICHAGE DE nH et x(H)""" 
fig, ax = plt.subplots(1,1, sharex=True,figsize=(7,7))# = plt.figure(figsize=(9,6))
ax.grid(True)

for fileIndex,lineFiles in enumerate(filenames):
    IRS1InuDnu = []
    IRS2InuDnu = []
    with open(lineFiles,'r') as csvfile:
        lines = csv.reader(csvfile,delimiter=' ')
        for row in lines:
            print(row)
            IRS1InuDnu.append(np.double(row[0]))
            IRS2InuDnu.append(np.double(row[1]))
    
    IRS1InuDnu= np.array(IRS1InuDnu)
    IRS2InuDnu= np.array(IRS2InuDnu)
    
    
    
    with open('bhr71-H2-extinction.txt','r') as Myfile:
        MyfileLines = Myfile.readlines()
        Extinction = np.double(MyfileLines) #A_lambda/A_V
        
    nsize1=np.size(IRS1InuDnu)
    
    A_V = 2
    alpha_nu = np.exp(-(A_V/2.5)*Extinction)    
    
    # correction d'extinction
    IRS1InuDnu = IRS1InuDnu * (1./alpha_nu) 
    IRS2InuDnu = IRS2InuDnu * (1./alpha_nu)    
    
    print('IRS1 petit :',IRS1InuDnu)
    print('IRS2 petit :',IRS2InuDnu)
    
    #IRS1InuDnu=FileFF[fileIndex]*IRS1InuDnu
    #IRS2InuDnu=FileFF[fileIndex]*IRS1InuDnu
    
    #print('IRS1 grand :',IRS1InuDnu)
    #print('IRS2 grand :',IRS2InuDnu)
    
    """ ===== Calcul de l'ordonnée ln(N_u/g_u) du diagramme d'excitation ===== """
    #nu : S(0),S(1),S(2)SH,S(2)SL,S(3),S(4),S(5),S(6),S(7)
    lambda_ij = np.array([2.8219e-05,
                          1.7035e-05,
                          1.2278e-05,
                          1.2278e-05,
                          9.6645e-06,
                          8.0255e-06,
                          6.9089e-06,
                          6.1086e-06,
                          5.5112e-06])
    
    #nu=np.array([10623,17598,24415,24415,31018,37357,43385,49146,54006])*1e9 #Hz
    nu=(clight*1E-2)/lambda_ij
    
    #Integrated Intensity (cgs: K.cm/s)
    W = np.multiply(((clight**3)/(2*KB*nu**3)),IRS1InuDnu)
    
    #coefficient d'Einstein pour l'émission spontanée A_ul : S(0),S(1),S(2)SH,S(2)SL,S(3),S(4),S(5),S(6),S(7)
    A_upTolow=np.array([0.294929E-10,
                       0.477160E-09,
                       0.276066E-08,
                       0.276066E-08,
                       0.985741e-08,
                       0.264890e-07,
                       0.589233e-07,
                       0.114443e-06,
                       0.200584e-06])
    
    
    A1 = ((8*np.pi*KB*nu**2)/(hPlanck*(clight**3)))
    
    #coefficient gamma_up : S(0),S(1),S(2)SH,S(2)SL,S(3),S(4),S(5),S(6),S(7)
    gamma_up = np.divide(A1,A_upTolow)
    
    #N_up^Thin (cm^-2) : S(0),S(1),S(2)SH,S(2)SL,S(3),S(4),S(5),S(6),S(7)
    Ncol_upThin = np.multiply(gamma_up,W)     
        
    #statistical weight : S(0),S(1),S(2)SH,S(2)SL,S(3),S(4),S(5),S(6),S(7)
    g_up=np.array([5,21,9,9,33,13,45,17,57])
    
    #N_up^Thin (cm^-2)/g_u
    DiagramOrdonnee_exp = np.divide(Ncol_upThin,g_up)
    DiagramOrdonnee_log = np.log(DiagramOrdonnee_exp)
    
    
        
    """ ===== Calcul de l'abscisse E_u/k (en K) du diagramme d'excitation ===== """
    
    DiagramAbscisse_lin=np.array([509.850,
                        1015.153,
                        1681.678,
                        1681.678,
                        2503.870,
                        3474.434,
                        4586.377,
                        5829.758,
                        7196.995])
        
    DiagramAbscisse_lin_temp=[]
    DiagramOrdonnee_log_temp=[]
    for iii,tested in enumerate(DiagramOrdonnee_log):
        if (not np.isnan(tested)):
            #print(iii,tested)
            DiagramAbscisse_lin_temp.append(DiagramAbscisse_lin[iii])
            DiagramOrdonnee_log_temp.append(DiagramOrdonnee_log[iii])
        
    DiagramAbscisse_lin=np.array(DiagramAbscisse_lin_temp)
    #DiagramOrdonnee_log=np.array(DiagramOrdonnee_log_temp)-np.log(FileFF[fileIndex])
    DiagramOrdonnee_log=np.array(DiagramOrdonnee_log_temp)
        
    nSize=np.size(DiagramOrdonnee_log)#np.size(DiagramOrdonnee_log[0])    
    Slope = np.divide((DiagramOrdonnee_log[1:nSize]-DiagramOrdonnee_log[0:nSize-1]),(DiagramAbscisse_lin[1:nSize]-DiagramAbscisse_lin[0:nSize-1]))
    
    Tex = -1./Slope    
    
    MeanTex = np.mean(Tex,axis=0)
    
    Eu_Xlist={}    
        
    
        
    """v in [-4.5,60] km/s"""
    
    
    """Fit linéaire :"""
    Eu_Xlist["-4.5 to 60 km/s"]=[]
    popt, pcov, xdata, ydata = [],[],[],[]
    Eu_gammes=[r'510 K to 2503.87 K',r'3474 K to 7196 K']
    XToFit=[((DiagramAbscisse_lin)[1:5]),
           ((DiagramAbscisse_lin)[6:])]
    YToFit=[((DiagramOrdonnee_log)[1:5]),
           ((DiagramOrdonnee_log)[6:])]
    for ii,dataToFit in enumerate(XToFit):
        popt.append(0)
        pcov.append(0)
        xdata.append(0)
        ydata.append(0)
        xdata[ii] = XToFit[ii]
        ydata[ii] = YToFit[ii]
        Eu_Xlist["-4.5 to 60 km/s"].append(xdata[ii])
    
        #Fit for the parameters a, b of the function func:
        popt[ii], pcov[ii] = curve_fit(func, xdata[ii], ydata[ii])
        print('Fit (a,b): ', popt[ii])
    """
    for ii in range(0,nSize-1):
    #ii=0
        popt.append(0)
        pcov.append(0)
        xdata.append(0)
        ydata.append(0)
        
        xdata[ii] = ((DiagramAbscisse_lin)[ii:ii+2])
        ydata[ii] = ((DiagramOrdonnee_log)[ii:ii+2])
        Eu_Xlist["-4.5 to 60 km/s"].append(xdata[ii])
    
        #Fit for the parameters a, b of the function func:
        popt[ii], pcov[ii] = curve_fit(func, xdata[ii], ydata[ii])
        print('Fit (a,b): ', popt[ii])
    
    """
    
    AVec=[]
    BVec=[]
    XVec=[]
    
    a,b,x=[],[],[]
    Npts=10.
    alphamin=0.8
    alphamax=1.2
    for k in range(0,np.shape(popt)[0]):
        a.append(popt[k][0])
        b.append(popt[k][1])
        x.append(np.arange(alphamin*xdata[k][0],alphamax*xdata[k][-1],np.abs(alphamin*xdata[k][0]-alphamax*xdata[k][-1])/Npts))
        AVec.append(a[k])
        BVec.append(b[k])
        XVec.append(x[k])
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        #Extinction = np.double(MyfileLines) #A_lambda/A_V
    
        #for row in lines:
            #print(row)
            #IRS1InuDnu.append(np.double(row[0]))
            #IRS2InuDnu.append(np.double(row[1]))
    
    
    """ Etape 2: Traitement des diagrammes, calcul de N_i et N_tot"""
    
    weight_g_u=[]
    Energy_E_u=[]
    with open('H2_levels_Evj.txt','r') as Myfile:
        MyfileLines = Myfile.readlines()
        MyfileLines = MyfileLines[6:]
    
    for lines in MyfileLines:
        CLine=lines.split()
        if(np.double(CLine[1])==0.):
            weight_g_u.append(np.double(CLine[-2]))
            Energy_E_u.append(np.double(CLine[-1]))
    
    weight_g_u=np.array(weight_g_u)
    Energy_E_u=np.array(Energy_E_u)
    
    TexVec=[]
    logNi_over_Z_Tex={}
    logNi_over_Z_Tex["T_ex"]=[]
    logNi_over_Z_Tex["logNZ"]=[]
    logNi_over_Z_Tex["Z(T_ex)"]=[]
    logNi_over_Z_Tex["N_i"]=[]
    
    for k in range(0,np.shape(popt)[0]):
        TexVec.append(-1./np.array(AVec[k]))
        logNi_over_Z_Tex["T_ex"].append(TexVec[k])
        logNi_over_Z_Tex["logNZ"].append(BVec[k])
        CCC = partition(TexVec[k],weight_g_u,Energy_E_u)
        logNi_over_Z_Tex["Z(T_ex)"].append(CCC)
        ZTex = (logNi_over_Z_Tex["Z(T_ex)"])[k][1]
        Bi = (logNi_over_Z_Tex["logNZ"])[k]
        Ni = ZTex*np.exp(Bi)
        (logNi_over_Z_Tex["N_i"]).append(np.float(Ni))
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if GrapheOut:       
    
        #fig = plt.figure(figsize=(10,10))
        #ax3=fig.add_subplot(3, 1, 3)
        #ax2=fig.add_subplot(3, 1, 2)
        #ax1=fig.add_subplot(3, 1, 1)
        
        plt.subplots_adjust(wspace=0, hspace=0.1)
        #ax3 = plt.axes([.525, .18, .15, .5], facecolor='white')
        #n, bins, patches = plt.hist(s, 400, normed=1)
        #plt.title('Probability')
        #plt.xticks([])
        #plt.yticks([])
        #ax3.yaxis.tick_right()
        #ax3.yaxis.set_label_position("right")
        
        bshrink = 1.
        
        " === axe 1 "
        ax.tick_params(axis='x', direction='in')
        ax.tick_params(axis='y', direction='in')
        ax.tick_params(bottom="on", top="on",left="on",right="on")
        ax.set_yscale('linear')
        ax.set_xscale('linear')
        #ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='both',which='both', direction='in')
        ax.set_yscale('linear')
        ax.set_xscale('linear')
        #ax.set_xlim(1e9,1e20)
        ax.set_xlabel(r'$\frac{E_u}{k_B}$ (en K)',fontsize=30,fontweight='extra bold')
        ax.set_ylabel(r'$\log\left(\frac{N_{u}^{thin}}{g_u}\right)$',fontsize=30,fontweight='extra bold')
        
        
        ax.set_title(r'Diagramme d\textquotesingle excitation du $H_2$')
        #fig.suptitle(r'Vitesse du son $c$ et $\mathcal{M}(V_n,V_i)$')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * bshrink, box.height])
        
        for k in range(0,np.shape(popt)[0]):
            #if(k!=2 and k!=5 and k!=6): 
            if(k!=2):
                a1=AVec[k]
                N1=(logNi_over_Z_Tex["N_i"])[k]
                ZTex1=(logNi_over_Z_Tex["Z(T_ex)"])[k][1]
                ax.plot(XVec[k],func(XVec[k],AVec[k],BVec[k]),'--',color='black',linewidth=2)
                fitResults.append(r'$y=$-1/{0:.3f}$\times x + \log(${1:.2e}$/{2:.2e})$'.format(-1./a1,N1,ZTex1))
        #        print(k+1,' N_i(cm-2)= ', (logNi_over_Z_Tex["N_i"])[k])
        ax.scatter(np.append(DiagramAbscisse_lin[0:5],DiagramAbscisse_lin[6:]),np.append(DiagramOrdonnee_log[0:5],DiagramOrdonnee_log[6:]),s=150, facecolors='none', edgecolors=FileColors[fileIndex], label=FileLabel[fileIndex], marker=FileMarkers[fileIndex])
        
            #rect = [0.7,0.7,-0.3,0.2]
            # this is an inset axes over the main axes
            #ax = fig.add_subplot(111)
            #ax1 = fig.add_subplot(111)
        
            
        
        #ax1.grid(True)
        
        #ax2.yaxis.tick_right()
        #ax2.yaxis.set_label_position("right")
        #plt.tight_layout()
        lgd = plt.legend(fontsize=24,loc='best',framealpha=1)
        #lgd = ax.legend(fontsize=16,loc='best',ncol=1,framealpha=1)
        #ax1.legend(fontsize=11,loc='right',bbox_to_anchor=(1.3, 0.5),ncol=1)
        #ax2.legend(fontsize=11,loc='right',bbox_to_anchor=(1.3, 0.5),ncol=1)
        #ax3.legend(fontsize=11,loc='right',bbox_to_anchor=(1.3, 0.5),ncol=1)
        
        
            
    eee=[]
    for k in range(0,np.shape(popt)[0]):
            if(k!=2): 
                print(k+1,' N_i(cm-2)= ', (logNi_over_Z_Tex["N_i"])[k],
                      'T_ex(K)=',(logNi_over_Z_Tex["T_ex"])[k],
                      'log(N/Z)=',(logNi_over_Z_Tex["logNZ"])[k],
                      'Z(T_ex)=',(logNi_over_Z_Tex["Z(T_ex)"])[k][1])
                eee.append([k+1,'E_up1,E_up2=',DiagramAbscisse_lin[k], DiagramAbscisse_lin[k+1],
                            'N_i(cm-2)= ', (logNi_over_Z_Tex["N_i"])[k],
                      'T_ex(K)=',(logNi_over_Z_Tex["T_ex"])[k],
                      'log(N/Z)=',(logNi_over_Z_Tex["logNZ"])[k],
                      'Z(T_ex)=',(logNi_over_Z_Tex["Z(T_ex)"])[k][1]])
        
    with open('H2_results.txt','a') as output:
        output.write(lineFiles+'\n')
        for k in range(0,np.shape(eee)[0]):
            for mn,elements in enumerate(eee[k]):
                if(mn==np.size(eee[k])): output.write(np.str(elements))
                else: output.write(np.str(elements)+' ')
            output.write('\n')
        
plt.savefig("H2_Excit_Diag.png", bbox_inches='tight')
plt.savefig("H2_Excit_Diag.pdf", bbox_inches='tight')
plt.show()
plt.close()

fig, ax = plt.subplots(1,1, sharex=True,figsize=(10,7))# = plt.figure(figsize=(9,6))
ax.grid(False)

ax.tick_params(axis='x', direction='in')
ax.tick_params(axis='y', direction='in')
ax.tick_params(bottom="on", top="on",left="on",right="on")
#ax.set_title(r'\textbf{Observed $H_2$ diagram against shock waves models}',fontsize=18)
ax.set_yscale('linear')
ax.set_xscale('linear')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
#ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='both', direction='in')

#ax.set_xlim(1e9,1e20)
ax.set_xlabel(r'$\frac{E_{up}}{k_B} (K)$',fontsize=23,fontweight='extra bold')
ax.set_ylabel(r'\textbf{$\log\left(\frac{N_{up}}{g_u}\right)$}',fontsize=23,fontweight='extra bold')

labelTextColor=['red','blue','green','brown','black']

Ytot=[]

""" Traitement de TOUS les modeles de chocs presents dans le dossier en cours"""

J_up=np.array([0,1,2,3,5,6,7])+2
E_up=np.array([509.850, 
               1015.153,
               1681.678,
               2503.870,
               4586.377,
               5829.758,
               7196.995])

#liste de tous les labels utilisés
labels={'c-n6-b1/v5' : r'$V_s=5 km.s^{-1}$ (choc C)',
        'c-n6-b1/v10': r'$V_s=10 km.s^{-1}$ (choc C)',
        'c-n6-b1/v19': r'$V_s=19 km.s^{-1}$ (choc C)',
        'c-n6-b1/v20': r'$V_s=20 km.s^{-1}$ (choc J)',
        'c-n6-b1/v25': r'$V_s=25 km.s^{-1}$ (choc J)',
        'c-n6-b1/v30': r'$V_s=30 km.s^{-1}$ (choc J)'
           }

        
#dirlist = os.popen('ls -1d c-*b1').read().split()
dirlist = ['c-n6-b1']

for ii,densiteChamp in enumerate(dirlist):
    cheminDensiteChamp=densiteChamp
    #ulist = os.popen('ls -1d '+cheminDensiteChamp+"/"+"v*").read().split()
    ulist = ['c-n6-b1/v5', 'c-n6-b1/v10','c-n6-b1/v19', 'c-n6-b1/v20', 'c-n6-b1/v25', 'c-n6-b1/v30']
    for kk,vitesse in enumerate(ulist):
        Y=[]
        Ytot.append([])
        chemin= vitesse + "/" + "excit.out"
        fig4data = ascii.read(chemin)
        Header = (  ascii.read(chemin,header_start=0,data_start=0))
        Header = Header[0]
        print(fig4data)


        
        for levels in [2,3,4,5,7,8,9]:
            Y.append(fig4data[levels][3])
        Ytot.append(0)
        Ytot[kk]=Y
        print("Y=",Y)
        #print(Y)
        ax.plot(E_up,Y,'-o',markersize=12,label=labels[vitesse])               
"""

for ii,numeroOfModels in enumerate(range(3,7)):
#    X=[]
    Y=[]
    Ytot.append([])
    labels.append(Header[numeroOfModels])
    labels[ii]=labels[ii].replace('_',' ')
    labels[ii]=labels[ii].replace('yr','ans')
    for levels in figIndex:
        Y.append(fig4data[levels][numeroOfModels])
    Ytot[ii]=Y
    #print(Y)
    if ii==4 : labels[ii]=r'etat stationnaire'
    ax.plot(E_up,Y,'-o',label=labels[ii],color=labelTextColor[ii],markersize=12)
    
"""

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
ax.plot(np.array(E_up),np.array(Y),'--o',markersize=17,color=labelTextColor[ii+1],markeredgecolor=labelTextColor[ii+1],markerfacecolor='none',label=labels[ii+1])

plt.text(3200, 44,r'\textbf{$\ce{H2}$}',fontsize=36,weight='bold')    
plt.text(500, 34,r'Tous les modeles: $v_s=20$ $km.s^{-1}$,',fontsize=20,weight='bold')    
plt.text(500, 33,r'$n_H=10^6$ $cm^{-3}$, $B=100$ $\mu G$',fontsize=20,weight='bold')    
"""
#observationnal data:
#Xobs=np.array([2,3,4,4,5,7,8,9])
Xobs=np.array([E_up[0],
               E_up[1],
               E_up[2],
               E_up[2],
               E_up[3],
               E_up[4],
               E_up[5],
               E_up[6]])
figObs=ascii.read('H2_int_int_avec_error.out')
Yobs=np.array([row[-2] for row in figObs])
Yerr=np.array([row[-1] for row in figObs])
ax.scatter(Xobs,Yobs,zorder=6,s=50,marker='D',edgecolors='black',facecolors='yellow',label='Observations')
ax.errorbar(Xobs, [row[-2] for row in figObs], yerr=[row[-1] for row in figObs],zorder=5,fmt='none',capsize=7,capthick=2,elinewidth=2,color='black')

#legend_properties = {'weight':'bold'}

#plt.plot([1,2,3], [4,5,6], label='Test')
#plt.legend(prop=legend_properties)

lgd = plt.legend(fontsize=12,loc='best',ncol=2)
lgd.set_title(r'\underline{$n_H=10^6~cm^{-3}$ $b=1$ :}',prop={'size':12})
"""
labelTextColor.append('black')
for ii,text in enumerate(lgd.get_texts()):
    text.set_color(labelTextColor[ii])
"""

plt.savefig("H2_Excit_Diag_c-n6.png", bbox_inches='tight')
plt.savefig("H2_Excit_Diag_c-n6.pdf", bbox_inches='tight')
plt.show()    
plt.close()