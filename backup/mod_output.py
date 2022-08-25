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
import csv
import os
from math import ceil
from astropy import units as u
import time
import struct
import fileinput
import yaml
from mod_exceptions import *
from mod_global_parameters import Variable #,global_variables

from vtk.util.numpy_support import vtk_to_numpy  # thats what you need
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from astropy.io import ascii
import vtktonumpy as ah_vtk
import numpy_support as ah
import read
from matplotlib import cm
from pyevtk.hl import gridToVTK
from astropy.io import ascii

    
def code_output_to_vtk(my_variable,my_sky_grid):
    my_sky_grid.masse_arrays = my_sky_grid.masse_arrays.reshape((my_sky_grid.nx, my_sky_grid.ny, my_sky_grid.nz))
    my_sky_grid.masse_arrays_ism = my_sky_grid.masse_arrays_ism.reshape((my_sky_grid.nx, my_sky_grid.ny, my_sky_grid.nz))
    my_sky_grid.masse_arrays_jet = my_sky_grid.masse_arrays_jet.reshape((my_sky_grid.nx, my_sky_grid.ny, my_sky_grid.nz))
    if(my_variable.frame_the_box):
        my_sky_grid.masse_arrays_box_inc = my_sky_grid.masse_arrays_box_inc.reshape((my_sky_grid.nx, my_sky_grid.ny, my_sky_grid.nz))
        my_sky_grid.masse_arrays_box_exc = my_sky_grid.masse_arrays_box_exc.reshape((my_sky_grid.nx, my_sky_grid.ny, my_sky_grid.nz))
    
    to_be_output = {}
    need_to_replace=['.','-','e']
    replace_by=['v','m','d']
    
    to_be_output["masse_totale"] = my_sky_grid.masse_arrays
    del my_sky_grid.masse_arrays
    to_be_output["masse_ism"] = my_sky_grid.masse_arrays_ism
    del my_sky_grid.masse_arrays_ism
    to_be_output["masse_jet"] = my_sky_grid.masse_arrays_jet
    del my_sky_grid.masse_arrays_jet
    
    if(my_variable.frame_the_box):
        to_be_output["masse_box_(included)"] = my_sky_grid.masse_arrays_box_inc
        del my_sky_grid.masse_arrays_box_inc
        to_be_output["masse_box_(excluded)"] = my_sky_grid.masse_arrays_box_exc
        del my_sky_grid.masse_arrays_box_exc

        
    if(my_variable.floor_ism): 
        msg_floor_ism = str(my_variable.threshold_f_ism)
        for i in range(len(need_to_replace)): msg_floor_ism = msg_floor_ism.replace(need_to_replace[i],replace_by[i])
        to_be_output["masse_ism_(fism_supeq)"] = my_sky_grid.masse_arrays_ism_inc
        del my_sky_grid.masse_arrays_ism_inc
        to_be_output["masse_ism_(fism_inf)"] = my_sky_grid.masse_arrays_ism_exc
        del my_sky_grid.masse_arrays_ism_exc
        
        
    if(my_variable.floor_jet):
        msg_floor_jet = str(my_variable.threshold_f_jet)
        for i in range(len(need_to_replace)): msg_floor_jet = msg_floor_jet.replace(need_to_replace[i],replace_by[i])
        to_be_output["masse_jet_(fjet_supeq)"] = my_sky_grid.masse_arrays_jet_inc
        del my_sky_grid.masse_arrays_jet_inc
        to_be_output["masse_jet_(fjet_inf)"] = my_sky_grid.masse_arrays_jet_exc
        del my_sky_grid.masse_arrays_jet_exc

    gridToVTK(my_variable.output_file_vtk, my_sky_grid.x, my_sky_grid.y, (my_sky_grid.z *
                  u.cm/u.s).to_value(my_variable.unit_output_v), cellData=to_be_output)
        
    
    
    #gridToVTK(my_variable.output_file_vtk, my_sky_grid.x, my_sky_grid.y, (my_sky_grid.z *
    #              u.cm/u.s).to_value(my_variable.unit_output_v), cellData={"masse totale":my_sky_grid.masse_arrays, 
    #                                                                       "masse ism":my_sky_grid.masse_arrays_ism,
    
    del to_be_output
    
def code_output_impulsion(my_variable,my_transverse_layer): 
    if(my_variable.unit_momentum_flux_output!=u.g*u.cm/(u.s*u.s)): my_transverse_layer.fluxes[:,0] = list((my_transverse_layer.fluxes[:,0] * u.g*u.cm/(u.s*u.s)).to_value(my_variable.unit_momentum_flux_output))
    if(my_variable.unit_mass_flux_output!=u.g/(u.s)): my_transverse_layer.fluxes[:,1] = list((my_transverse_layer.fluxes[:,1] * u.g/u.s).to_value(my_variable.unit_mass_flux_output))
        
    
    output=""
    output+="center_r center_z center_phi radius inclination cells_r cells_z cells_phi cells_surface_radius cells_surface_phi flux_momentum flux_momentum_unit flux_mass flux_mass_unit\n"
    for iline,line in enumerate(my_transverse_layer.fluxes):
        output+=str(my_transverse_layer.fluxes_surface_center[iline][my_variable.r_])+" "
        output+=str(my_transverse_layer.fluxes_surface_center[iline][my_variable.z_cyl_])+" "
        output+=str(my_transverse_layer.fluxes_surface_center[iline][my_variable.phi_])+" "
        
        output+=str(my_transverse_layer.fluxes_surface_inclination[iline][my_variable.r_])+" "
        output+=str(my_transverse_layer.fluxes_surface_inclination[iline][-1])+" "
        
        output+=str(my_transverse_layer.fluxes_surface_r_z[iline][my_variable.r_])+" "
        output+=str(my_transverse_layer.fluxes_surface_r_z[iline][my_variable.z_cyl_])+" "
        output+=str(my_transverse_layer.fluxes_surface_r_z[iline][my_variable.phi_])+" "

        output+=str(my_transverse_layer.fluxes_surface_radius_z_inclination[iline][my_variable.r_])+" "
        output+=str(my_transverse_layer.fluxes_surface_radius_z_inclination[iline][-1])+" "
        
        output+=str(my_transverse_layer.fluxes[iline,0])+" "
        output+=str(my_variable.unit_momentum_flux_output).replace(' ','_')+" "
        output+=str(my_transverse_layer.fluxes[iline,1])+" "
        output+=str(my_variable.unit_mass_flux_output).replace(' ','_')+"\n"
    #print(output)
    dat_output = ascii.read(output,guess=False)
    ascii.write(dat_output,my_variable.output_impulsion_ascii,overwrite=True)
    