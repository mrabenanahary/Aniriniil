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
def construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka,numberOfCells,ndimProblem):
    temp_CellsPoint = np.empty((numberOfCells,ndimProblem*2,))
    {CellsData.GetCellBounds(i,temp_CellsPoint[i]) for i in range(0,numberOfCells)}
    for idim in range(ndimProblem):
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,idim]=np.tile(temp_CellsPoint[:,2*idim:2*idim+2],ndimProblem-1)
    r_offset = np.amin(np.tile(temp_CellsPoint[:,2*my_variables.r_:2*my_variables.r_+2],ndimProblem-1))
    my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,my_variables.phi_]=my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:, :, my_variables.phi_]-0.5*my_bazooka.geo_dangle_phi
    if(r_offset<0): my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:,:,my_variables.r_]=my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[:, :, my_variables.r_]-r_offset
    
    r_z_phi_max = np.array([list(np.amax(
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, numberOfCells)])#+r_offset
    r_z_phi_min = np.array([list(np.amin(
        my_bazooka.cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, numberOfCells)])#+r_offset


    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.r_] = np.abs(r_z_phi_max[:, my_variables.r_]*r_z_phi_max[:, my_variables.r_] -
                                                     r_z_phi_min[:, my_variables.r_]*r_z_phi_min[:, my_variables.r_])
    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.z_cyl_] = np.abs(r_z_phi_max[:, my_variables.z_cyl_] -
                                                         r_z_phi_min[:, my_variables.z_cyl_])
    my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.phi_] = my_bazooka.geo_dangle_phi
    
    points2D = 0.5*(r_z_phi_max+r_z_phi_min)
    if(my_variables.frame_the_box):

        f_box_2D_inc = np.where(((points2D[:,my_variables.x_frame]>=my_variables.x_frame_min)&(points2D[:,my_variables.x_frame]<=my_variables.x_frame_max))&((points2D[:,my_variables.y_frame]>=my_variables.y_frame_min)&(points2D[:,my_variables.y_frame]<=my_variables.y_frame_max)),True,False)

        f_box_2D_exc = ~f_box_2D_inc



    
    
        f_box_2D_inc = np.where(f_box_2D_inc,1.0,0.0)
        f_box_2D_exc = np.where(f_box_2D_exc,1.0,0.0)
    
    points2D[:,my_variables.phi_] = points2D[:,my_variables.phi_] + 0.5*my_bazooka.geo_dangle_phi
    
    temp_dangle_t = my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.phi_]
    temp_drr = my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.r_]  # r2**2-r1**2
    temp_dz = my_bazooka.cellsPoints_drr_dz_dphi_2D_frame[:, my_variables.z_cyl_]
    my_bazooka.cellsPoints_volume_2D_frame = 0.5 * temp_drr * temp_dangle_t * \
        temp_dz  # 0.5 * (r2**2-r1**2) * (phi2-phi1) * (z2 - z1)
    del temp_dangle_t, temp_drr, temp_dz, temp_CellsPoint
    if(my_variables.frame_the_box):
        return my_bazooka,points2D,f_box_2D_inc,f_box_2D_exc
    else:
        return my_bazooka,points2D
    
def load_the_bazooka(my_bazooka,my_variables,points2D,numberOfCells,f_box_2D_inc=False,f_box_2D_exc=False):
    for n in range(0, Variable.nphi):  # phi \in [0,2\pi[
        # instance n
        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] = n*numberOfCells
        # instance phi:
        geo_angle_phi = n*my_bazooka.geo_dangle_phi
    
        geo_angle_phi_array = np.full((numberOfCells), geo_angle_phi)
        del geo_angle_phi
    
        # arrays components for which we assign the same values in this loop level
    
        my_bazooka.geo_position_3D_cell_index_in_2D_frame[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                               my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                               numberOfCells] = np.arange(0, numberOfCells)
    
        my_bazooka.geo_position_r_z_phi[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                             my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                             numberOfCells, my_variables.r_] = points2D[:, my_variables.r_]  # r in 3D
        my_bazooka.geo_position_r_z_phi[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                             my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                             numberOfCells, my_variables.z_cyl_] = points2D[:, my_variables.z_cyl_]  # z in 3D
        my_bazooka.geo_position_r_z_phi[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                             my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                             numberOfCells, my_variables.phi_] = geo_angle_phi_array  # phi in 3D
        my_bazooka.geo_position_cellsPoints_volume[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells] = my_bazooka.cellsPoints_volume_2D_frame 
        my_bazooka.geo_position_f_ism[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells] = my_variables.rhoTracers['ISM']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])
        my_bazooka.geo_position_f_jet[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells] = my_variables.rhoTracers['jet']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])
        if(my_variables.frame_the_box):
            my_bazooka.geo_position_f_box_inc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells] = f_box_2D_inc
            my_bazooka.geo_position_f_box_exc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells] = f_box_2D_exc
        # radial velocity
        my_bazooka.geo_cells_r_z_phi_v[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                            my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                            numberOfCells, my_variables.r_] = my_variables.velocityVector['v1']  # r in 3D
        my_bazooka.geo_cells_r_z_phi_v[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                            my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                            numberOfCells, my_variables.z_cyl_] = my_variables.velocityVector['v2']  # z in 3D
        my_bazooka.geo_cells_r_z_phi_v[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                            my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                            numberOfCells, my_variables.phi_] = my_variables.velocityVector['v3']  # z in 3D
        del geo_angle_phi_array
        # arrays components for which we use intermediate arrays
        geo_r_array = np.array(my_bazooka.geo_position_r_z_phi[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                    my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                    numberOfCells, my_variables.r_])
        geo_phi_array = np.array(my_bazooka.geo_position_r_z_phi[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                      my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                      numberOfCells, my_variables.phi_])
        # 3D cartesian coordinates
        my_bazooka.geo_position_x_y_z[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                           my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                           numberOfCells, my_variables.x_] = list(geo_r_array*np.cos(geo_phi_array))  # x in 3D
        my_bazooka.geo_position_x_y_z[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                           my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                           numberOfCells, my_variables.y_] = list(geo_r_array*np.sin(geo_phi_array))  # y in 3D
        my_bazooka.geo_position_x_y_z[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                           my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                           numberOfCells, my_variables.z_cart_] = points2D[:, my_variables.z_cyl_]  # z in 3D
        my_bazooka.geo_rho[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                numberOfCells] = my_variables.geo_rho_2D
        my_bazooka.geo_position_cellsPoints_masse[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                       my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                       numberOfCells] = my_bazooka.geo_rho[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                                 my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                                 numberOfCells] * \
            my_bazooka.geo_position_cellsPoints_volume[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                            my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                            numberOfCells] #raw mass
        my_bazooka.geo_position_cellsPoints_masse_ism[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                       my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                       numberOfCells] = my_bazooka.geo_position_cellsPoints_masse[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                                 my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                                 numberOfCells] *\
        my_bazooka.geo_position_f_ism[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells]
        my_bazooka.geo_position_cellsPoints_masse_jet[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                       my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                       numberOfCells] = my_bazooka.geo_position_cellsPoints_masse[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                                 my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                                 numberOfCells] *\
        my_bazooka.geo_position_f_jet[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells]
        if(my_variables.frame_the_box):
            my_bazooka.geo_position_cellsPoints_masse_box_inc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                       my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                       numberOfCells] = my_bazooka.geo_position_cellsPoints_masse[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                                 my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                                 numberOfCells] *\
        my_bazooka.geo_position_f_box_inc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells]
            my_bazooka.geo_position_cellsPoints_masse_box_exc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                       my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                       numberOfCells] = my_bazooka.geo_position_cellsPoints_masse[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                                                 my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                                                 numberOfCells] *\
        my_bazooka.geo_position_f_box_exc[my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n]:
                                        my_bazooka.geo_index_first_cell_in_2D_frame_iphi[n] +
                                        numberOfCells]
                        
        message = "Progression : {0:.2f} % --> point n° : {1} ".format(100*my_bazooka.geo_3D_index/my_bazooka.numberOf3Delem,my_bazooka.geo_3D_index)
        print(message)
        my_bazooka.geo_3D_index += numberOfCells
        
    print(''.ljust(len(message),' '))
    print("Progression : Finished.... 100 %")
    print("number of 3D points = {0}".format(len(my_bazooka.geo_rho)))
    
    
    # convert, if needed, the mass from gram to unit_mass :
    if(my_variables.unit_masse != u.g):
        my_bazooka.geo_position_cellsPoints_masse = (
            my_bazooka.geo_position_cellsPoints_masse*u.g).to_value(my_variables.unit_masse)
        my_bazooka.geo_position_cellsPoints_masse_ism = (
            my_bazooka.geo_position_cellsPoints_masse_ism*u.g).to_value(my_variables.unit_masse)
        my_bazooka.geo_position_cellsPoints_masse_jet = (
            my_bazooka.geo_position_cellsPoints_masse_jet*u.g).to_value(my_variables.unit_masse)
        if(my_variables.frame_the_box):
            my_bazooka.geo_position_cellsPoints_masse_box_inc = (
                my_bazooka.geo_position_cellsPoints_masse_box_inc*u.g).to_value(my_variables.unit_masse)
            my_bazooka.geo_position_cellsPoints_masse_box_exc = (
                my_bazooka.geo_position_cellsPoints_masse_box_exc*u.g).to_value(my_variables.unit_masse)
    
    
    
    
    print("\nMinimal mass found : {0} {1}".format(
        np.str(np.amin(my_bazooka.geo_position_cellsPoints_masse)), str(my_variables.unit_masse)))
    print("Maximal mass found : {0} {1}".format(
        np.str(np.amax(my_bazooka.geo_position_cellsPoints_masse)), str(my_variables.unit_masse)))
    
    return my_bazooka

def full_shoot_with_the_bazooka(my_bazooka,my_variables):
    # méthode 2 plus compliquée mais plus rapide
    my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.Gamma_] = my_bazooka.geo_position_r_z_phi[:, my_variables.r_]*np.sin(my_variables.geo_angle_i) *\
        np.cos(my_bazooka.geo_position_r_z_phi[:, my_variables.phi_]) +\
        my_bazooka.geo_position_r_z_phi[:, my_variables.z_cyl_]*np.cos(my_variables.geo_angle_i)  # r in 3D
    my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.Delta_] = -my_bazooka.geo_position_r_z_phi[:, my_variables.r_]*np.cos(my_variables.geo_angle_i) *\
        np.cos(my_bazooka.geo_position_r_z_phi[:, my_variables.phi_]) + \
        my_bazooka.geo_position_r_z_phi[:, my_variables.z_cyl_]*np.sin(my_variables.geo_angle_i)  # z in 3D
    my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.alpha_] = my_bazooka.geo_position_r_z_phi[:,
                                                                     my_variables.r_]*np.sin(my_bazooka.geo_position_r_z_phi[:, my_variables.phi_])  # phi in 3D
    if(my_variables.voxelling_mode == 'Voxel'):
        my_bazooka.geo_vobs = my_bazooka.geo_cells_r_z_phi_v[:, my_variables.phi_] *\
            np.sin(my_variables.geo_angle_i)*np.sin(my_bazooka.geo_position_r_z_phi[:, my_variables.phi_]) -\
            -my_bazooka.geo_cells_r_z_phi_v[:, my_variables.r_]*np.sin(my_variables.geo_angle_i)*np.cos(my_bazooka.geo_position_r_z_phi[:, my_variables.phi_]) -\
            my_bazooka.geo_cells_r_z_phi_v[:, my_variables.z_cyl_]*np.cos(my_variables.geo_angle_i)
    
    return my_bazooka