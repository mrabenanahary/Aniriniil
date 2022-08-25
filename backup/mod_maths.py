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
    
"""brute scheme methods"""
def layer_cells_surfaces(CellsData,my_variables,my_transverse_layer):
    temp_CellsPoint = np.empty((my_transverse_layer.numberOfCells,my_transverse_layer.ndimProblem*2,))
    {CellsData.GetCellBounds(i,temp_CellsPoint[i]) for i in range(0,my_transverse_layer.numberOfCells)}
    for idim in range(my_transverse_layer.ndimProblem):
        my_transverse_layer.cells_r_z_phi[:,:,idim]=np.tile(temp_CellsPoint[:,2*idim:2*idim+2],my_transverse_layer.ndimProblem-1)
    r_offset = np.amin(np.tile(temp_CellsPoint[:,2*my_variables.r_:2*my_variables.r_+2],my_transverse_layer.ndimProblem-1))
    my_transverse_layer.cells_r_z_phi[:,:,my_variables.phi_]=my_transverse_layer.cells_r_z_phi[:, :, my_variables.phi_]#-0.5*my_transverse_layer.geo_dangle_phi
    if(r_offset<0): my_transverse_layer.cells_r_z_phi[:,:,my_variables.r_]=my_transverse_layer.cells_r_z_phi[:, :, my_variables.r_]-r_offset
    
    
    r_z_phi_max = np.array([list(np.amax(
        my_transverse_layer.cells_r_z_phi[i], axis=0)) for i in range(0, my_transverse_layer.numberOfCells)])#+r_offset
    r_z_phi_min = np.array([list(np.amin(
        my_transverse_layer.cells_r_z_phi[i], axis=0)) for i in range(0, my_transverse_layer.numberOfCells)])#+r_offset
    
    #print(np.shape(r_z_phi_max),np.shape(r_z_phi_min))
    
    my_transverse_layer.cells_drr = np.abs(r_z_phi_max[:, my_variables.r_]*r_z_phi_max[:, my_variables.r_] -
                                                     r_z_phi_min[:, my_variables.r_]*r_z_phi_min[:, my_variables.r_]) # (r2^2-r1^2)

    points2D = 0.5*(r_z_phi_max+r_z_phi_min)
    
    #points2D[:,my_variables.phi_] = points2D[:,my_variables.phi_] + my_transverse_layer.geo_dangle_phi
    
    
    del temp_CellsPoint
    
    return my_transverse_layer,points2D,r_z_phi_max

def compute_momentum_flux_through_circle(my_variables,points2D,rmax,my_transverse_layer,center=[[0,5.0e16,0]],r_z_inclination=[[2.5e15,0.0]]):
    output = []
    output_r_z = [] 
    output_r_z_inclination = []
    #indexes = []
    for i,centre in enumerate(center):
        
        nearest_z = np.min(np.abs(points2D[:,my_variables.z_cyl_]-centre[my_variables.z_cyl_]))
        nearest_z = np.max(points2D[:,1][np.where(np.abs(points2D[:,1]-centre[my_variables.z_cyl_])==nearest_z)])
        nearest_r = np.min(np.abs(points2D[:,my_variables.r_]-centre[my_variables.r_]))
        nearest_r = np.max(points2D[:,my_variables.r_][np.where(np.abs(points2D[:,my_variables.r_]-centre[my_variables.r_])==nearest_r)])
        
        
        output_r_z.append([nearest_r,nearest_z,0.0])
        
        # int(rho*sign(v_z)*|v_z|*v.dS) = int(rho*sign(v_z)*|v_z|*v_z*dA) = int(rho*sign(v_z)*v_z^2*rdrdt)
        # = sum_over_each_cell(rho_i * v_{i,z}^2 * int(rdr,r=rmax_i,r=rmin_i)*int(dt,t=2pi,t=0))
        # = sum_over_each_cell(rho_i * v_{i,z}^2 * 0.5 * [rmax_i^2-rmin_i^2] * 2pi)
        # = pi * sum_over_each_cell(rho_i * v_{i,z}^2 * [rmax_i^2-rmin_i^2])
        
        concerned_indexes = np.where((points2D[:, my_variables.r_]<=r_z_inclination[i][my_variables.r_])&(np.abs(points2D[:,1]-nearest_z)==np.min(np.abs(points2D[:,1]-nearest_z))))[0]
        
        output_r_z_inclination.append([np.max((rmax[:,my_variables.r_])[concerned_indexes]),0.0])
        
        layer_drr = my_transverse_layer.cells_drr[concerned_indexes]
        layer_rho = my_variables.geo_rho_2D[concerned_indexes]
        layer_vz  = (my_variables.velocityVector['v2'])[concerned_indexes]
        
        
        #output.append(np.sum(layer_rho * (layer_vz**2) * np.pi * layer_drr))
        output.append([np.sum(layer_rho * np.sign(layer_vz) * (layer_vz**2) * np.pi * layer_drr),\
                       np.sum(layer_rho * layer_vz * np.pi * layer_drr)])
        
        
        """<<< Succesfully checked until this point!"""
        #print("M flow rate min",np.min(layer_rho),np.min(layer_vz),np.min(layer_drr))
        #print("M flow rate max",np.max(layer_rho),np.max(layer_vz),np.max(layer_drr))
        #print("M flow rate sum",np.sum(layer_rho),np.sum(layer_vz),np.sum(layer_drr))
        #indexes.append(concerned_indexes)
    return center,r_z_inclination,output_r_z,output_r_z_inclination,np.array(output)

        
    