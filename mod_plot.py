#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 23:24:13 2019

@author: mialy94
=================================================================
= Auteur : Mialy RABENANAHARY
= 
= globals.py
= 
= Module pour lire le fichier parametres entrer par l utilisateur
= 
=
=================================================================
"""

#global global_variables


from cycler import cycler
import numpy as np
import vtk as v
import csv as csv
import os
os.environ["PATH"] += os.pathsep + '/obs/mrabenanahary/usr/local/texlive/bin/x86_64-linux'
#print(os.getenv("PATH"))
from math import ceil
from astropy import units as u
import time
import struct
import fileinput
import yaml
from mod_exceptions import *
from mod_global_parameters import Variable,create_box_file_dict_name #,global_variables

from vtk.util.numpy_support import vtk_to_numpy  # thats what you need
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from astropy.io import ascii
import vtktonumpy as ah_vtk
import numpy_support as ah
import read
from matplotlib import cm
from pyevtk.hl import gridToVTK
import seaborn as sns
import pylatex




from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import combinations
from itertools import cycle
import string
from matplotlib.ticker import LogLocator
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
import matplotlib.patches as mpatches

import sys
# sys.path.append('./ParaView-5.8.0-RC1-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages')
#from paraview.simple import *
#from math import *
import matplotlib as mpl
import matplotlib.ticker as mtick
import vtk as v
#from vtk import *
from matplotlib.ticker import AutoMinorLocator
import matplotlib.text as mtext
from math import ceil

#import yt
# print(os.getcwd())

# importation des packages et modules durhampy personnalisées crées
import time
import struct
import vtktonumpy as ah_vtk
import numpy_support as ah
import read
#import xarray as xa
from matplotlib import cm
from pyevtk.hl import gridToVTK
import shutil

#importation des packages et modules durhampy personnalisées crées
#from src.postProcessing.plotting import *
#from src.postProcessing.calculus import *

#global global_variables

from mod_global_parameters import * #"""DONE!"""
from mod_read_params import * #"""DONE!"""
from mod_bazooka import *  #"""DONE!"""

def get_cycle(cmap, N=None, use_index="auto"):
    if isinstance(cmap, str):
        if use_index == "auto":
            if cmap in ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                        'Dark2', 'Set1', 'Set2', 'Set3',
                        'tab10', 'tab20', 'tab20b', 'tab20c']:
                use_index=True
            else:
                use_index=False
        cmap = mpl.cm.get_cmap(cmap)
    if not N:
        N = cmap.N
    if use_index=="auto":
        if cmap.N > 100:
            use_index=False
        elif isinstance(cmap, LinearSegmentedColormap):
            use_index=False
        elif isinstance(cmap, ListedColormap):
            use_index=True
    if use_index:
        ind = np.arange(int(N)) % cmap.N
        return cycler("color",cmap(ind))
    else:
        colors = cmap(np.linspace(0,1,N))
        return colors

my_cmap = np.array(
[
[5, 97, 255, 1],
[5, 99, 253, 1],
[5, 101, 252, 1],
[5, 103, 251, 1],
[5, 104, 250, 1],
[5, 106, 248, 1],
[5, 108, 247, 1],
[5, 110, 246, 1],
[5, 112, 245, 1],
[5, 113, 244, 1],
[5, 115, 242, 1],
[5, 117, 241, 1],
[5, 119, 240, 1],
[5, 121, 239, 1],
[5, 122, 237, 1],
[5, 124, 236, 1],
[5, 126, 235, 1],
[5, 128, 234, 1],
[5, 130, 233, 1],
[5, 131, 231, 1],
[5, 133, 229, 1],
[5, 134, 228, 1],
[5, 136, 226, 1],
[5, 137, 225, 1],
[5, 139, 223, 1],
[5, 140, 221, 1],
[5, 142, 220, 1],
[5, 143, 218, 1],
[5, 145, 217, 1],
[5, 146, 215, 1],
[5, 148, 213, 1],
[5, 150, 212, 1],
[5, 151, 210, 1],
[5, 153, 208, 1],
[5, 154, 207, 1],
[5, 156, 205, 1],
[5, 157, 204, 1],
[5, 159, 202, 1],
[5, 160, 200, 1],
[5, 161, 198, 1],
[5, 163, 196, 1],
[5, 164, 194, 1],
[5, 166, 192, 1],
[5, 167, 190, 1],
[5, 169, 188, 1],
[5, 170, 186, 1],
[5, 171, 184, 1],
[5, 173, 182, 1],
[5, 174, 181, 1],
[5, 176, 179, 1],
[5, 177, 177, 1],
[5, 178, 175, 1],
[5, 180, 173, 1],
[5, 181, 171, 1],
[5, 183, 169, 1],
[5, 184, 167, 1],
[5, 186, 165, 1],
[5, 187, 162, 1],
[5, 189, 160, 1],
[5, 190, 158, 1],
[5, 192, 156, 1],
[5, 193, 153, 1],
[5, 195, 151, 1],
[5, 197, 149, 1],
[5, 198, 147, 1],
[5, 200, 144, 1],
[5, 201, 142, 1],
[5, 203, 140, 1],
[5, 204, 138, 1],
[5, 206, 135, 1],
[5, 207, 133, 1],
[5, 209, 131, 1],
[5, 211, 129, 1],
[5, 212, 126, 1],
[5, 214, 123, 1],
[5, 215, 120, 1],
[5, 216, 118, 1],
[5, 218, 115, 1],
[5, 219, 112, 1],
[5, 221, 109, 1],
[6, 222, 106, 1],
[6, 223, 103, 1],
[6, 225, 101, 1],
[6, 226, 98, 1],
[6, 228, 95, 1],
[6, 229, 92, 1],
[6, 230, 89, 1],
[5, 232, 86, 1],
[5, 233, 83, 1],
[5, 235, 81, 1],
[4, 236, 78, 1],
[4, 237, 75, 1],
[14, 238, 69, 1],
[25, 239, 63, 1],
[36, 240, 58, 1],
[47, 241, 52, 1],
[58, 242, 46, 1],
[68, 243, 40, 1],
[78, 243, 38, 1],
[87, 244, 36, 1],
[96, 244, 34, 1],
[105, 244, 32, 1],
[115, 245, 31, 1],
[124, 245, 29, 1],
[131, 246, 26, 1],
[137, 246, 24, 1],
[143, 247, 21, 1],
[150, 248, 18, 1],
[156, 248, 15, 1],
[162, 249, 13, 1],
[168, 249, 12, 1],
[173, 250, 11, 1],
[178, 250, 11, 1],
[183, 250, 10, 1],
[188, 251, 10, 1],
[193, 251, 9, 1],
[198, 251, 8, 1],
[203, 252, 8, 1],
[208, 252, 7, 1],
[213, 252, 7, 1],
[218, 253, 6, 1],
[222, 253, 6, 1],
[227, 253, 5, 1],
[232, 253, 5, 1],
[237, 254, 5, 1],
[242, 254, 4, 1],
[247, 254, 4, 1],
[252, 255, 4, 1],
[255, 254, 5, 1],
[255, 252, 8, 1],
[255, 250, 10, 1],
[255, 248, 13, 1],
[255, 246, 16, 1],
[255, 245, 19, 1],
[255, 243, 21, 1],
[255, 241, 24, 1],
[255, 239, 27, 1],
[255, 237, 30, 1],
[255, 235, 33, 1],
[255, 233, 36, 1],
[255, 232, 39, 1],
[255, 230, 42, 1],
[255, 228, 45, 1],
[255, 226, 48, 1],
[255, 224, 50, 1],
[255, 222, 53, 1],
[255, 220, 55, 1],
[255, 218, 55, 1],
[255, 216, 55, 1],
[255, 214, 55, 1],
[255, 212, 55, 1],
[255, 210, 55, 1],
[255, 208, 55, 1],
[255, 206, 55, 1],
[255, 204, 55, 1],
[255, 202, 55, 1],
[255, 200, 55, 1],
[255, 198, 55, 1],
[255, 196, 55, 1],
[255, 194, 55, 1],
[255, 192, 55, 1],
[255, 190, 55, 1],
[255, 188, 55, 1],
[255, 186, 55, 1],
[255, 184, 55, 1],
[255, 182, 55, 1],
[255, 180, 55, 1],
[255, 178, 55, 1],
[255, 176, 55, 1],
[255, 174, 55, 1],
[255, 172, 55, 1],
[255, 170, 55, 1],
[255, 168, 55, 1],
[255, 166, 55, 1],
[255, 164, 55, 1],
[255, 162, 55, 1],
[255, 160, 55, 1],
[255, 158, 55, 1],
[255, 156, 55, 1],
[255, 154, 55, 1],
[255, 152, 55, 1],
[255, 150, 55, 1],
[255, 147, 55, 1],
[255, 145, 55, 1],
[255, 143, 55, 1],
[255, 140, 55, 1],
[255, 138, 55, 1],
[255, 136, 55, 1],
[255, 133, 55, 1],
[255, 131, 55, 1],
[255, 129, 55, 1],
[255, 126, 55, 1],
[255, 124, 55, 1],
[255, 122, 55, 1],
[255, 119, 55, 1],
[255, 117, 55, 1],
[255, 114, 55, 1],
[255, 112, 55, 1],
[255, 110, 55, 1],
[255, 107, 55, 1],
[255, 105, 55, 1],
[255, 102, 55, 1],
[255, 99, 55, 1],
[254, 96, 55, 1],
[254, 92, 54, 1],
[254, 89, 54, 1],
[254, 86, 54, 1],
[253, 83, 53, 1],
[253, 80, 52, 1],
[253, 77, 52, 1],
[252, 74, 51, 1],
[252, 71, 50, 1],
[252, 68, 49, 1],
[252, 64, 49, 1],
[252, 59, 50, 1],
[252, 54, 51, 1],
[252, 50, 52, 1],
[252, 45, 52, 1],
[253, 40, 53, 1],
[252, 37, 55, 1],
[250, 36, 56, 1],
[249, 35, 58, 1],
[247, 33, 60, 1],
[245, 32, 62, 1],
[243, 31, 63, 1],
[242, 29, 65, 1],
[240, 28, 67, 1],
[238, 26, 68, 1],
[236, 24, 70, 1],
[234, 23, 72, 1],
[232, 21, 73, 1],
[230, 20, 75, 1],
[228, 18, 77, 1],
[226, 16, 78, 1],
[224, 15, 80, 1],
[222, 13, 81, 1],
[220, 11, 83, 1],
[218, 10, 85, 1],
[215, 10, 86, 1],
[213, 11, 87, 1],
[211, 11, 88, 1],
[208, 11, 89, 1],
[206, 11, 90, 1],
[203, 11, 92, 1],
[201, 11, 93, 1],
[199, 11, 94, 1],
[196, 11, 95, 1],
[194, 12, 96, 1],
[191, 12, 97, 1],
[189, 12, 99, 1],
[186, 12, 100, 1],
[184, 12, 101, 1],
[182, 12, 102, 1],
[179, 12, 103, 1],
[177, 13, 104, 1],
[174, 13, 106, 1]
]
)

my_cmap=np.float64(my_cmap)
my_cmap[:,0:3]=my_cmap[:,0:3]/255

def set_shared_xlabel(a, xlabel, labelpad = 0.01,**kwargs):
    """Set a x label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    xlabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""
        
    varplus={}
    
    for k,v in kwargs.items():
        varplus[k]=v

    f = a[0].get_figure()
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    right = a[0].get_position().x1
    left = a[-1].get_position().x0

    # get the coordinates of the bottom side of the tick labels 
    y0 = 1
    for at in a:
        #at.set_xlabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.xaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.transformed(f.transFigure.inverted())#inverse_transformed(f.transFigure)
        yt = bboxes.y1
        if yt > y0:
            y0 = yt
    tick_label_bottom = y0

    # set position of label
    a[-1].set_xlabel(xlabel,**kwargs)
    a[-1].xaxis.set_label_coords((right + left)/2,tick_label_bottom - labelpad, transform=f.transFigure)


def set_shared_xlabel_2(a, xlabel, labelpad = 0.01,**kwargs):
    """Set a x label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    xlabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""
        
    varplus={}
    
    for k,v in kwargs.items():
        varplus[k]=v

    f = a[0].get_figure()
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    right = a[0].get_position().x1
    left = a[-1].get_position().x0

    # get the coordinates of the bottom side of the tick labels 
    y0 = 1
    for at in a:
        #at.set_xlabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.xaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.transformed(f.transFigure.inverted())#inverse_transformed(f.transFigure)
        yt = bboxes.y0
        if yt < y0:
            y0 = yt
    tick_label_bottom = y0

    # set position of label
    a[-1].set_xlabel(xlabel,**kwargs)
    a[-1].xaxis.set_label_coords((right + left)/2,tick_label_bottom - labelpad, transform=f.transFigure)

def set_shared_ylabel(a, ylabel, labelpad = 0.01,**kwargs):
    """Set a y label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    ylabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""

    f = a[0].get_figure()
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    top = a[0].get_position().y1
    bottom = a[-1].get_position().y0

    # get the coordinates of the left side of the tick labels 
    x0 = 1
    for at in a:
        at.set_ylabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.yaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.transformed(f.transFigure.inverted())
        xt = bboxes.x0
        if xt < x0:
            x0 = xt
    tick_label_left = x0

    # set position of label
    a[-1].set_ylabel(ylabel,**kwargs)
    a[-1].yaxis.set_label_coords(tick_label_left - labelpad,(bottom + top)/2, transform=f.transFigure)



def points2D_builder(CellsData,my_variables,my_bazooka):
    temp_CellsPoint = np.empty((my_bazooka.numberOfCells,my_bazooka.ndimProblem*2,))
    {CellsData.GetCellBounds(i,temp_CellsPoint[i]) for i in range(0,my_bazooka.numberOfCells)}
    for idim in range(my_bazooka.ndimProblem):
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,idim]=np.tile(temp_CellsPoint[:,2*idim:2*idim+2],my_bazooka.ndimProblem-1)
    r_offset = np.amin(np.tile(temp_CellsPoint[:,2*my_variables.r_:2*my_variables.r_+2],my_bazooka.ndimProblem-1))
    my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,my_variables.phi_]=my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:, :, my_variables.phi_]-0.5*my_bazooka.geo_dangle_phi
    if(r_offset<0): my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,my_variables.r_]=my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:, :, my_variables.r_]-r_offset
    
    r_z_phi_max = np.array([list(np.amax(
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, my_bazooka.numberOfCells)])#+r_offset
    r_z_phi_min = np.array([list(np.amin(
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, my_bazooka.numberOfCells)])#+r_offset


    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.r_] = np.abs(r_z_phi_max[:, my_variables.r_]*r_z_phi_max[:, my_variables.r_] -
                                                     r_z_phi_min[:, my_variables.r_]*r_z_phi_min[:, my_variables.r_])
    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.z_cyl_] = np.abs(r_z_phi_max[:, my_variables.z_cyl_] -
                                                         r_z_phi_min[:, my_variables.z_cyl_])
    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.phi_] = my_bazooka.geo_dangle_phi
    
    points2D = 0.5*(r_z_phi_max+r_z_phi_min)
    
    return points2D







"""============================================================================="""

def plot_maps(my_variables, CellsData, numberOfCells, ndimProblem, numberOfVertixes, numberOfVertixesUnderDimension):

    os.write(1,b"PLotting figures: Begins...\n")
    
    my_bazooka = Bazooka()
    my_bazooka.numberOfCells,\
    my_bazooka.ndimProblem,\
    my_bazooka.numbeOfVertixes,\
    my_bazooka.numberOfVertixesUnderDimension = \
        (numberOfCells,
         ndimProblem, 
         numberOfVertixes, 
         numberOfVertixesUnderDimension)
    #free some unnecessarily used memory
    my_bazooka.set_attributes(CellsData,my_variables)
    
    cell_positions=points2D_builder(CellsData,my_variables,my_bazooka) 
    
    for iplot,varnameplot in enumerate(my_variables.variable_to_plot):
        os.write(1,(">>> Fig. {0} of {1}: Plotting {2}...\n".format(str(iplot+1),str(np.size(my_variables.variable_to_plot)),my_variables.variable_to_plot[iplot])).encode('ascii'))
        
        

        LINE_STYLES = ['solid', 'dashed', 'dashdot', 'dotted']
        NUM_STYLES = len(LINE_STYLES)
        my_ls=LINE_STYLES[0]
        LINE_STYLES = cycle(LINE_STYLES)
        iLINE_STYLES = 0 
        
        if(my_variables.plot_use_sns):
            sns.reset_orig()  # get default matplotlib styles back
            clrs = sns.color_palette(my_variables.plot_sns_palette[iplot], n_colors=my_variables.plot_sns_palette_n[iplot])  # a list of RGB tuples
            NUM_COLORS = len(clrs)
            palette = cycle(clrs)

        else:
            #Default matplotlib palette
            #clrs=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
            clrs=get_cycle(my_variables.plot_mpl_palette[iplot], N=my_variables.plot_mpl_palette_n[iplot])
            NUM_COLORS = len(clrs)
            palette = cycle(clrs)
        
        data = read.extract(CellsData, my_variables.variable_to_plot[iplot], attribute_mode=my_variables.geo_2D_attribute_mode)
        if(my_variables.variable_to_plot[iplot] in ['v1','v2','v3']):
            data = data / 1e5
        
        nrow=1
        ncol=1
        
        h, w = my_variables.plot_height[iplot],my_variables.plot_width[iplot]
        xlabel = my_variables.variable_label[iplot] #r'$\mathrm{v_R}$ ($\mathrm{km~s^{-1}}$)'
        
        titles = my_variables.plot_title[iplot]
        title_pad = my_variables.plot_title_pad[iplot]
        xmin=my_variables.plot_xmin[iplot]
        xmax=my_variables.plot_xmax[iplot]
        dr = my_variables.plot_dr[iplot]
        ymin=my_variables.plot_ymin[iplot]
        ymax=my_variables.plot_ymax[iplot]
        vmin=my_variables.plot_vmin[iplot]
        vmax=my_variables.plot_vmax[iplot]
    
        
        graphName='{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/',''))+'-{0}'.format((my_variables.variable_to_plot[iplot]).replace(' ','_'))+'.pdf'
        cbar_x = 1.0
        cbar_y = 0.0
        cbar_width=0.02
        cbar_height=1.0
        majorLength=2*10
        majorWidth=2*0.8
        majorLengthAxes=15
        majorWidthAxes=1.6
        minorLength=majorLength*0.5
        minorWidth=majorWidth
        deltaX5=0.99
        deltaY5=0.9
        deltaY5b=0.975
        plot_fontsize = my_variables.plot_fontsize[iplot]
        plot_fontsize_cbr = my_variables.plot_fontsize_cbr[iplot]
        plot_fontsize_title = my_variables.plot_fontsize_title[iplot]
        plot_fontsize_global = my_variables.plot_fontsize_global[iplot]
        # Get the colormap colors, multiply them with the factor "a", and create new colormap
        
        
        darken_factor = my_variables.plot_darken_factor[iplot]
        vZcolormap = my_variables.plot_colormap[iplot]
        
        if(vZcolormap=='my_cmap'):
            listed_vZ_cmap = my_cmap
            listed_vZ_cmap[:,0:3] *= darken_factor
            listed_vZ_cmap = ListedColormap(listed_vZ_cmap)
        else:
            vZ_cmap = cm.get_cmap(vZcolormap, 2*256)
            newcolors = vZ_cmap(np.linspace(0,1,2*256))
            #white = np.array([255/255, 255/255, 255/255, 1])
            #gray = np.array([128/255, 128/255, 128/255, 1])
            #newcolors[0:2, :] = white
            #newcolors[254:256, :] = gray
            newcolors[0:5, :] = newcolors[0, :]
            newcolors[:,0:3] *= darken_factor
            listed_vZ_cmap = ListedColormap(newcolors)
            
            vZ_cmap2 = cm.get_cmap(vZcolormap, 2*256)
            newcolors = vZ_cmap(np.linspace(0,1,2*256))
            #white = np.array([255/255, 255/255, 255/255, 1])
            #gray = np.array([128/255, 128/255, 128/255, 1])
            #newcolors[0:2, :] = white
            #newcolors[254:256, :] = gray
            newcolors[0:7, :] = newcolors[0, :]
            newcolors[:,0:3] *= darken_factor
            listed_vZ_cmap2 = ListedColormap(newcolors)
        
        
        
        dark_turbo=listed_vZ_cmap
        
        """Parameters for matplotlib"""
        dpi=300
        mpl.rcParams['figure.dpi']= dpi
        mpl.rc("savefig", dpi=dpi)
        font = {'family' : 'serif',
                'size'   : plot_fontsize_global}
        mpl.rc('font', **font)
        
        #usetex = mpl.checkdep_usetex(True)
        #plt.rc('text', usetex=usetex)
        #params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
        #plt.rcParams.update(params)

        
        
        sum_max=max(ncol*abs(xmax-xmin),nrow*abs(ymax-ymin))
        sum_min=min(ncol*abs(xmax-xmin),nrow*abs(ymax-ymin))
        
        if(sum_max==ncol*abs(xmax-xmin)):
            h=h*(sum_max/sum_min)
        else:
            w=w*(sum_max/sum_min)
            
        title_pad = title_pad*(h/8.0)
        
        col=my_variables.variable_to_plot[iplot]
        
        fig, axs = plt.subplots(nrow, ncol,
                               subplot_kw={
                                   'xlim': [-xmax,xmax],
                                   'ylim': [ymin,ymax]
                                   },
                               figsize=(2*h,w))
        
        old_axs = axs
        axs = []
        axs.append([old_axs])
        del old_axs
        axs=np.array(axs)
        
        run_img = {}
        refcol={}
        
        i=0
        j=0
        xx, yy= np.meshgrid(np.unique(cell_positions[:,0]), np.unique(cell_positions[:,1]),indexing='xy')
        run_img[col]=[xx/1e17,yy/1e17]
        del xx,yy
        x_vals, x_idx = np.unique(cell_positions[:,0], return_inverse=True)
        y_vals, y_idx = np.unique(cell_positions[:,1], return_inverse=True)
        vals_array = np.zeros(x_vals.shape + y_vals.shape)
        vals_array[x_idx, y_idx] = data
        zz = vals_array.T
        del vals_array
        run_img[col].append(zz)
        refcol[col]=[]
        if((my_variables.variable_to_plot[iplot] in ['v1','v2','v3'])or(vmin<=0)):
            refcol[col].append(axs[i][j].pcolormesh((run_img[col])[0],(run_img[col])[1],(run_img[col])[2], shading='auto', cmap=dark_turbo, vmin=vmin,vmax=vmax,linewidth=0,rasterized=True))
            refcol[col].append(axs[i][j].pcolormesh(-(run_img[col])[0],(run_img[col])[1],(run_img[col])[2], shading='auto', cmap=dark_turbo, vmin=vmin,vmax=vmax,linewidth=0,rasterized=True))
        else:
            refcol[col].append(axs[i][j].pcolormesh((run_img[col])[0],(run_img[col])[1],(run_img[col])[2], shading='auto', cmap=dark_turbo, norm=mpl.colors.LogNorm(vmin=vmin,vmax=vmax),linewidth=0,rasterized=True))
            refcol[col].append(axs[i][j].pcolormesh(-(run_img[col])[0],(run_img[col])[1],(run_img[col])[2], shading='auto', cmap=dark_turbo, norm=mpl.colors.LogNorm(vmin=vmin,vmax=vmax),linewidth=0,rasterized=True))
        if(i<nrow-1): 
                axs[i][j].tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=True,      # ticks along the bottom edge are off
                top=True,         # ticks along the top edge are on
                labeltop=False,
                labelbottom=False,
                direction='inout',
                length=majorLengthAxes,
                width=majorWidthAxes
                )   
        else:
            axs[i][j].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=True,      # ticks along the bottom edge are off
            top=True,         # ticks along the top edge are on
            labeltop=False,
            labelbottom=True,
            direction='inout',
            length=majorLengthAxes,
            width=majorWidthAxes
            )
        if(j>0):
            axs[i][j].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left=True,      # ticks along the bottom edge are off
            right=True,         # ticks along the top edge are on
            labelleft=False,
            labelright=False,
            direction='inout',
            length=majorLengthAxes,
            width=majorWidthAxes
            )          
        else:
            axs[i][j].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left=True,      # ticks along the bottom edge are off
            right=True,         # ticks along the top edge are on
            labelleft=True,
            labelright=False,
            direction='inout',
            length=majorLengthAxes,
            width=majorWidthAxes
            )
        
        axs[i][j].set_xlim([-xmax,xmax])
        axs[i][j].set_xticks(np.arange(-xmax,xmax+dr,dr))
        axs[i][j].set_xticklabels(np.abs(np.round(np.arange(-xmax,xmax+dr,dr),2)))
        #axs[i][2*j+1].set_xlim([xmin,xmax])
        #axs[i][2*j+1].set_xlim([xmin,xmax])
        if(j!=0): axs[i][j].set_yticklabels([])
        if(j<ncol-1) : axs[i][j].set_xticks(axs[i][j].get_xticks()[:-1])
        
        axs[i][j].set_aspect("equal")
        
        axs[i][j].set_xlabel(r'$R$ ($\times 10^{17}$ cm)',fontsize=plot_fontsize)
        axs[i][j].set_ylabel(r'$z$ ($\times 10^{17}$ cm)',fontsize=plot_fontsize,labelpad=-0.0)
        axs[i][j].set_title(titles,fontsize=plot_fontsize_title,pad=title_pad)
        
        """ ====== Drawing multi-boxes when needed!! ======="""
        
        if((my_variables.frame_the_box)&(my_variables.plotting)&(my_variables.plotting_box)):
            for irctgel,rctgel in enumerate(my_variables.x_frame_min):
                box_lw = my_variables.plot_rect_lw[irctgel]
                
                left, bottom, width, height = (my_variables.x_frame_min[irctgel],
                                               my_variables.y_frame_min[irctgel],
                                               abs(my_variables.x_frame_max[irctgel]-my_variables.x_frame_min[irctgel]),
                                               abs(my_variables.y_frame_max[irctgel]-my_variables.y_frame_min[irctgel]))
                
                global_to_local_length = 1e17 #cm
                left, bottom, width, height = left/global_to_local_length,\
                                              bottom/global_to_local_length,\
                                              width/global_to_local_length,\
                                              height/global_to_local_length
                                              
                if((irctgel+1)%NUM_COLORS==0): my_ls=next(LINE_STYLES)
                rect=mpatches.Rectangle((left,bottom),width,height, 
                                fill=False,
                                linestyle=my_ls,
                                color=next(palette),
                                linewidth=box_lw)
                                #facecolor="red")
                axs[i][j].add_patch(rect)        
         
        box1 = axs[nrow-1][0].get_position()
        box2 = axs[0][ncol-1].get_position()
        #axi.set_position([box.x0*1.05, box.y0, box.width, box.height])
        
        minorticks = []
        majorticks = []
        if((my_variables.variable_to_plot[iplot] in ['v1','v2','v3'])or(vmin<=0)):
            logvmin=vmin
            logvmax=vmax
            #majorticks.append(vmin)
            for i in range(logvmin,logvmax):
                tt=i*10
                if(tt>=vmin and tt<=vmax) : majorticks.append(tt)
                for j in range(1,10):
                    t=5*j
                    if(t>=vmin and t<=vmax): minorticks.append(t)
            #majorticks.append(vmax)        
        else:
            logvmin=int(np.floor(np.log10(vmin)))
            logvmax=int(np.ceil(np.log10(vmax)))
            for i in range(logvmin,logvmax):
                tt=10**i
                if(tt>=vmin and tt<=vmax) : majorticks.append(tt)
                for j in range(1,10):
                    t=j*10**i
                    if(t>=vmin and t<=vmax): minorticks.append(t)
        
        
                
        # create color bar
        axColor = plt.axes([box1.x0+cbar_x*(box2.x1-box1.x0), box1.y0+cbar_y*(box2.y1-box1.y0), cbar_width*(box2.x1-box1.x0), cbar_height*(box2.y1-box1.y0)])        
        cbr = fig.colorbar(refcol[col][0], cax=axColor, orientation="vertical",ticks=majorticks)
        cbr.set_label(xlabel,fontsize=plot_fontsize_cbr,verticalalignment='baseline', rotation=270)
        #cbr.ax.set_title(r'$\mathrm{n_H}$ ($\mathrm{cm^{-3}}$)',fontsize=35,verticalalignment='baseline',loc='left',rotation=0)
        
        cbr.ax.yaxis.set_ticks(minorticks, minor=True)
        cbr.ax.tick_params(which='major', length=majorLength, width=majorWidth, direction='inout' ,top=True,bottom=False,labelbottom=False,labeltop=True,labelsize=plot_fontsize_cbr)
        cbr.ax.tick_params(which='minor', length=minorLength, width=minorWidth, direction='inout' ,top=True,bottom=False,labelbottom=False,labeltop=True)
        

        
        #plt.tight_layout()
        #plt.show()        
        fig.savefig(graphName)
        plt.close('all')
        
        if os.path.exists('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}'.format(graphName))):
            os.remove('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}'.format(graphName)))
        shutil.move('./{0}'.format(graphName),'{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/',''))+'/'+graphName)

    os.write(1,b"Plotting figures : Done...\n")
