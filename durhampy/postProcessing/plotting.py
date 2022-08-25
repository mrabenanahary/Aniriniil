#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 19:02:31 2020

@author: mialy94
=================================================================
= Auteur : Mialy RABENANAHARY
= 
= globals.py
= 
= Module contenant toutes les variables globales utilisées dans les scripts
= qui font appel à ce module
= 
=
=================================================================
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
import matplotlib.text as mtext
#from durhampy.data.globals import *

#variables globales
#global input_dir
def du_plot(
XTab,
YTab,
figSize=(5,5),
shareX=True,
subplotsArg=[1,1],
axGrid=True,
XScaleType = 'linear',
YScaleType = 'linear',
fontSize=8,
legendLocalisation='best',
Xlabel=r'',Ylabel=r'',
xMin=0,xMax=0,yMin=0,yMax=0,
typeOfPlot = 'scatter',
plotMarker = 'o',
plotLabel = r' out',
lineWidth = 1,
markerSize = 5,
moreCommands = """plt.tight_layout()""",
saveFigSwitch = False,
saveFigName="default_output.pdf",
saveFigBboxInches='tight',
supTitle=r'',
closePlot = True,
retakePlot = False,
inFig=['',''],
color='',
XFontsize='',
YFontsize=''
):
    if(retakePlot == False):    
        fig, ax = plt.subplots(subplotsArg[0],subplotsArg[1], sharex=shareX,figsize=figSize)# = plt.figure(figsize=(9,6))
    else:
        if(inFig==['','']):
            raise AssertionError('Il faut entrer une variable figure et un axe a traiter sous la forme inFig=[fig,ax]')
        else:
            fig,ax =inFig[0],inFig[1]           
    ax.grid(axGrid)
    plt.subplots_adjust(wspace=0, hspace=0.1)
    
    
    changeLimits=True
    changeLimits = changeLimits and (xMin==0)
    changeLimits = changeLimits and (xMax==0)
    if(not changeLimits):
        ax.set_xlim(xMin,xMax)
    changeLimits=True
    changeLimits = changeLimits and (yMin==0)
    changeLimits = changeLimits and (yMax==0)
    if(not changeLimits):
        ax.set_ylim(yMin,yMax)
    
    if(retakePlot == False):
        bshrink = 1.
        if(XFontsize==''): ax.set_xlabel(Xlabel)
        else: ax.set_xlabel(Xlabel,fontsize=XFontsize)
        if(YFontsize==''): ax.set_ylabel(Ylabel)
        else: ax.set_ylabel(Ylabel,fontsize=YFontsize)
        fig.suptitle(supTitle)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * bshrink, box.height])
        ax.set_yscale(YScaleType)
        ax.set_xscale(XScaleType)   
    
    
    
    if(typeOfPlot=='plot'): 
        if(color==''): ax.plot(XTab,YTab,plotMarker,label=plotLabel,linewidth=lineWidth,markersize=markerSize)
        else: ax.plot(XTab,YTab,plotMarker,label=plotLabel,linewidth=lineWidth,markersize=markerSize,color=color)
    elif(typeOfPlot=='scatter'):
        if(color==''): ax.scatter(XTab,YTab,marker=plotMarker,label=plotLabel,linewidth=lineWidth,s=markerSize)
        else: ax.scatter(XTab,YTab,marker=plotMarker,label=plotLabel,linewidth=lineWidth,s=markerSize,color=color)
    else:
        raise AssertionError('typeOfPlot=\"{:s}\" n est pas un type de plot reconnu par cette version du code'.format(typeOfPlot))
    
    lgd = ax.legend(fontsize=fontSize,loc=legendLocalisation)
    
    exec(moreCommands,locals(),globals())
    
    
    if(saveFigSwitch==True) : plt.savefig(saveFigName,bbox_inches=saveFigBboxInches )
    outfig,outax=fig,ax

    if(closePlot):
        plt.show()
        plt.close()    
    return outfig,outax
