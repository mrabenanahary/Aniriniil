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
from mod_global_parameters import Variable,Sky_Grid,create_box_file_dict_name #,global_variables
    
"""railgun scheme methods"""

def compute_2D_grid_nodes_positions_volumes(CellsData,my_variables,my_railgun):
        CellsNodes = np.empty((my_railgun.numberOfCells,my_railgun.ndimProblem*2,))
        {CellsData.GetCellBounds(i,CellsNodes[i]) for i in range(0,my_railgun.numberOfCells)}
        for idim in range(my_railgun.ndimProblem):
            my_railgun.nodes_r_z_phi[:,:,idim]=np.tile(CellsNodes[:,2*idim:2*idim+2],my_railgun.ndimProblem-1)
            if(idim==my_variables.r_): r_offset = np.amin(my_railgun.nodes_r_z_phi[:,:,idim])
        #r_offset = np.amin(np.tile(CellsNodes[:,2*my_variables.r_:2*my_variables.r_+2],my_railgun.ndimProblem-1))
        del CellsNodes
        
        my_railgun.nodes_r_z_phi[:,:,my_variables.phi_]=my_railgun.nodes_r_z_phi[:, :, my_variables.phi_]-0.5*my_railgun.geo_dangle_phi
        
        # il arrive que pour des raisons de convéniences numériques, les noeuds de la grille ne soit pas délimité par min(r)=0 mais une petite valeur min(r)=dr<0:
        # ici, on translate les valeurs en r de toute la grille 2D de -dr pour que le domaine soit délimité par min(r)=0 
        if(r_offset<0): my_railgun.nodes_r_z_phi[:,:,my_variables.r_]=my_railgun.nodes_r_z_phi[:, :, my_variables.r_]-r_offset
        
        # construire les éléments positions en 3D de références, dans le repère 3D cylindrique, des noeuds de chaque cellule du plan 2D de la simu 2D d'AMRVAC
        # plongées dans une grille 3D d'épaisseur 1 celulle en phi
        r_z_phi_max = np.array([list(np.amax(
            my_railgun.nodes_r_z_phi[i], axis=0)) for i in range(0, my_railgun.numberOfCells)])#+r_offset
        r_z_phi_min = np.array([list(np.amin(
            my_railgun.nodes_r_z_phi[i], axis=0)) for i in range(0, my_railgun.numberOfCells)])#+r_offset
        del my_railgun.nodes_r_z_phi
        cells = 0.5*(r_z_phi_max+r_z_phi_min)
        my_railgun.bound_radius = np.max(r_z_phi_max[:,my_variables.r_])
        my_railgun.bound_z = np.abs(np.max(r_z_phi_max[:,my_variables.z_cyl_]))
        
        # construire de même les éléments dr^2, dz, dphi pour calculer la masse de chaque cellule de cette référence plongée en 3D
        #
        # - int(rdr,r_min,rmax) = 0.5 * [r_max^2-r_min^2] ==> 2*my_railgun.cells_drr_dz_dphi[:, my_variables.r_]
        # - int(dz,z_min,z_max) = [z_max-z_min] ==> my_railgun.cells_drr_dz_dphi[:, my_variables.z_cyl_]
        # - int(dphi,phi_min,phi_max) = [phi_max-phi_min] = dphi = 2*pi/nphi ==> my_railgun.cells_drr_dz_dphi[:, my_variables.phi_]
        #
        my_railgun.cells_drr_dz_dphi[:, my_variables.r_] = np.abs(r_z_phi_max[:, my_variables.r_]*r_z_phi_max[:, my_variables.r_] -
                                                         r_z_phi_min[:, my_variables.r_]*r_z_phi_min[:, my_variables.r_])
        my_railgun.cells_drr_dz_dphi[:, my_variables.z_cyl_] = np.abs(r_z_phi_max[:, my_variables.z_cyl_] -
                                                             r_z_phi_min[:, my_variables.z_cyl_])
        my_railgun.cells_drr_dz_dphi[:, my_variables.phi_] = my_railgun.geo_dangle_phi
        
        
        #si on cherche à sélectionner des zones à clip
        if(my_variables.frame_the_box):
            my_railgun.f_box_2D_inc={}
            my_railgun.f_box_2D_exc={}
            if(my_variables.frame_ism):
                my_railgun.f_box_2D_ism_inc={}
                my_railgun.f_box_2D_ism_exc={}
            if(my_variables.frame_jet):
                my_railgun.f_box_2D_jet_inc={}
                my_railgun.f_box_2D_jet_exc={}                
            for box_i,box_el in enumerate(my_variables.x_frame_min):
                my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                         my_variables.x_frame_max[box_i],
                                                         my_variables.y_frame_min[box_i],
                                                         my_variables.y_frame_max[box_i])
                my_railgun.f_box_2D_inc[my_dict_name] = np.where(((cells[:,my_variables.x_frame]>=my_variables.x_frame_min[box_i] )&\
                                                                  (cells[:,my_variables.x_frame]<=my_variables.x_frame_max[box_i]))&\
                                                                 ((cells[:,my_variables.y_frame]>=my_variables.y_frame_min[box_i] )&\
                                                                  (cells[:,my_variables.y_frame]<=my_variables.y_frame_max[box_i])),True,False)
                my_railgun.f_box_2D_exc[my_dict_name] = ~my_railgun.f_box_2D_inc[my_dict_name]
                my_railgun.f_box_2D_inc[my_dict_name] = np.where(my_railgun.f_box_2D_inc[my_dict_name],1.0,0.0)
                my_railgun.f_box_2D_exc[my_dict_name] = np.where(my_railgun.f_box_2D_exc[my_dict_name],1.0,0.0)
                if(my_variables.frame_ism):
                    my_railgun.f_box_2D_ism_inc[my_dict_name] = np.where(((cells[:,my_variables.x_frame]>=my_variables.x_frame_min[box_i] )&\
                                                                      (cells[:,my_variables.x_frame]<=my_variables.x_frame_max[box_i]))&\
                                                                     ((cells[:,my_variables.y_frame]>=my_variables.y_frame_min[box_i] )&\
                                                                      (cells[:,my_variables.y_frame]<=my_variables.y_frame_max[box_i])),True,False)
                    my_railgun.f_box_2D_ism_exc[my_dict_name] = ~my_railgun.f_box_2D_ism_inc[my_dict_name]
                    my_railgun.f_box_2D_ism_inc[my_dict_name] = np.where(my_railgun.f_box_2D_ism_inc[my_dict_name],1.0,0.0)
                    my_railgun.f_box_2D_ism_exc[my_dict_name] = np.where(my_railgun.f_box_2D_ism_exc[my_dict_name],1.0,0.0)
                if(my_variables.frame_jet):
                    my_railgun.f_box_2D_jet_inc[my_dict_name] = np.where(((cells[:,my_variables.x_frame]>=my_variables.x_frame_min[box_i] )&\
                                                                      (cells[:,my_variables.x_frame]<=my_variables.x_frame_max[box_i]))&\
                                                                     ((cells[:,my_variables.y_frame]>=my_variables.y_frame_min[box_i] )&\
                                                                      (cells[:,my_variables.y_frame]<=my_variables.y_frame_max[box_i])),True,False)
                    my_railgun.f_box_2D_jet_exc[my_dict_name] = ~my_railgun.f_box_2D_jet_inc[my_dict_name]
                    my_railgun.f_box_2D_jet_inc[my_dict_name] = np.where(my_railgun.f_box_2D_jet_inc[my_dict_name],1.0,0.0)
                    my_railgun.f_box_2D_jet_exc[my_dict_name] = np.where(my_railgun.f_box_2D_jet_exc[my_dict_name],1.0,0.0)                    
                
        # centrer les cellules (i.e. leur centre) de références en phi=0,
        # dans l'image de placer les noeuds des cellules en phi=+/-0.5*dphi=+/-0.5*2*pi/nphi
        cells[:,my_variables.phi_] = cells[:,my_variables.phi_] + my_railgun.geo_dangle_phi
        
        # calculer le volume desdites cellules :
        # int(r*dphi*dr*dz,{phi=-0.5*2*pi/nphi,phi=0.5*2*pi/nphi},{r=r_min,r=r_max},{z=z_min,z=z_max})
        # = 0.5 * (r_max^2-r_min^2) * (phi_max-phi_min) * (z_max - z_min)
        my_railgun.cells_Volume = 0.5 * \
            my_railgun.cells_drr_dz_dphi[:, my_variables.r_] * \
            my_railgun.cells_drr_dz_dphi[:, my_variables.phi_] * \
            my_railgun.cells_drr_dz_dphi[:, my_variables.z_cyl_]  
        del my_railgun.cells_drr_dz_dphi
        
        os.write(1,b"Done...\n")
        
        
        os.write(1,b"Computing fractions of ISM and jet and extracting density...\n")
        #fractions nécessaires pour le calcul des 
        my_railgun.cells_f_ISM = my_variables.rhoTracers['ISM']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])
        my_railgun.cells_f_jet = my_variables.rhoTracers['jet']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])



        del my_variables.rhoTracers['ISM'],my_variables.rhoTracers['jet']
        
        return my_variables,my_railgun,cells        


def building_empty_voxel_mesh(my_variables,my_railgun,my_sky_grid):
        v_x_max = my_variables.bound_offset_frac*np.max(my_railgun.bound_radius)
        v_x_min = -v_x_max
        
        
        v_y_max = my_variables.bound_offset_frac*np.sqrt(my_railgun.bound_radius**2+my_railgun.bound_z**2)
        v_y_min = -v_y_max
        
        del my_railgun.bound_radius,my_railgun.bound_z

        my_sky_grid.set_attributes(my_variables) #don't forget to update the SKy_Grid attributes according to the variables current settings
            

        #Construire la grille d'échantillonnage sur le plan du ciel
        
        # Dimensions
        #nx, ny, nz = my_variables.voxels_nGamma_cells, my_variables.voxels_nGamma_cells, 1
        if(v_x_max > v_y_max):
            nx = int(my_variables.voxels_nGamma_cells)
            if(v_x_min<0): X = np.linspace(0, v_x_max, num=nx, dtype='float64')
            else: X = np.linspace(v_x_min, v_x_max+v_x_min, num=nx, dtype='float64')
            nx = len(X)
            # Y = np.linspace(-v_y_max, v_y_max, , dtype='float64')
            dx = 0.5*(np.max(X[1:]-X[:-1])+np.min(X[1:]-X[:-1]))
            if(v_x_min<0):
                nxc = int(ceil(np.abs(v_x_min)/dx))
                v_x_min = -nxc*dx
                Xc = np.linspace(v_x_min, 0, num=nxc+1, dtype='float64')
                X = np.concatenate((Xc[:-1],X))
                nx = len(X)
                del Xc
            #X = np.array([-v_x_max+i*dx for i in range(0, nx+1)])
            #nx = len(X)
        
            ny = int(ceil(v_y_max/dx))
            v_y_max = ny*dx
            if(v_y_min<0): Y = np.linspace(0, v_y_max, num=ny+1, dtype='float64')
            else: Y = np.linspace(v_y_min, v_y_max+v_y_min, num=ny+1, dtype='float64')
            ny = len(Y)
            #ny = len(Y)
            dy = 0.5*(np.max(Y[1:]-Y[:-1])+np.min(Y[1:]-Y[:-1]))
            #Y = np.array([-v_y_max+i*dx for i in range(0, ny+1)])
            if(v_y_min<0):
                nyc = int(ceil(np.abs(v_y_min)/dy))
                v_y_min = -nyc*dy
                Yc = np.linspace(v_y_min, 0, num=nyc+1, dtype='float64')
                Y = np.concatenate((Yc[:-1],Y)) #to center in zero
                ny = len(Y)
                del Yc
            
            
        else:
            ny = int(my_variables.voxels_nGamma_cells)
            if(v_y_min<0) : Y = np.linspace(0, v_y_max, num=ny, dtype='float64')
            else: Y = np.linspace(v_y_min, v_y_max+v_y_min, num=ny, dtype='float64')
            ny = len(Y)
            # Y = np.linspace(-v_y_max, v_y_max, , dtype='float64')
            dy = 0.5*(np.max(Y[1:]-Y[:-1])+np.min(Y[1:]-Y[:-1]))
            #Y = np.array([-v_y_max+i*dy for i in range(0, ny+1)])
            #ny = len(Y)
            if(v_y_min<0):
                nyc = int(ceil(np.abs(v_y_min)/dy))
                v_y_min = -nyc*dy
                Yc = np.linspace(v_y_min, 0, num=nyc+1, dtype='float64')
                Y = np.concatenate((Yc[:-1],Y)) #to center in zero
                ny = len(Y)
                del Yc
        
            nx = int(ceil(v_x_max/dy))
            v_x_max = nx*dy
            if(v_x_min<0) : X = np.linspace(0, v_x_max, num=nx+1, dtype='float64')
            else : X = np.linspace(v_x_min, v_x_max+v_x_min, num=nx+1, dtype='float64')
            nx = len(X)
            #nx = len(X)
            dx = 0.5*(np.max(X[1:]-X[:-1])+np.min(X[1:]-X[:-1]))#dy
            #X = np.array([-v_x_max+i*dx for i in range(0, nx+1)])
            if(v_x_min<0):
                nxc = int(ceil(np.abs(v_x_min)/dx))
                v_x_min = -nxc*dx
                Xc = np.linspace(v_x_min, 0, num=nxc+1, dtype='float64')
                X = np.concatenate((Xc[:-1],X)) #to center in zero
                nx = len(X)
                del Xc
        if(my_variables.voxelling_mode == 'Masse'):
            nz = 1
            lz = 1
            dz = lz/nz
            Z = np.arange(0, lz+0.1*dz, dz, dtype='float64')
        elif(my_variables.voxelling_mode == 'Voxel'):
            v_max = my_variables.voxels_z_max #max(np.abs(np.amin(geo_vobs)), np.abs(np.amax(geo_vobs)))
            v_min = my_variables.voxels_z_min
            nz = int(ceil(v_max/my_variables.voxels_dz))
            #v_max = nz*my_variables.voxels_dz
            #dz = my_variables.voxels_dz
            #Z = np.arange(-lz, lz+0.1*dz, dz, dtype='float64')
            if(v_min<0): Z = np.linspace(0, v_max, num=nz+1, dtype='float64')
            else: Z = np.linspace(v_min, v_max+v_min, num=nz+1, dtype='float64')
            nz = len(Z)
            # Y = np.linspace(-v_y_max, v_y_max, , dtype='float64')
            dz = 0.5*(np.max(Z[1:]-Z[:-1])+np.min(Z[1:]-Z[:-1]))
            if(v_min<0):
                nzc = int(ceil(np.abs(v_min)/dz))
                v_min = -nzc*dz
                Zc = np.linspace(v_min, 0, num=nzc+1, dtype='float64')
                Z = np.concatenate((Zc[:-1],Z))
                nz = len(Z)
                del Zc
        
        
        print("v_X/Y/Z_min = ", v_x_min, v_y_min, v_min)
        print("X/Y/Z_min = ", np.min(X), np.min(Y), np.min(Z))
        print("v_X/Y/Z_max = ", v_x_max, v_y_max, v_max)
        print("X/Y/Z_max = ", np.max(X), np.max(Y), np.max(Z))
        print("dx,dy,dz = ",dx,dy,dz)
        
        #ncells = nx * ny * nz
        #npoints = (nx + 1) * (ny + 1) * (nz + 1)
        # Coordinates
        #Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')
        my_sky_grid.set_attributes(my_variables) 
        my_sky_grid.x = np.zeros((nx + 1, ny + 1, nz + 1))
        my_sky_grid.y = np.zeros((nx + 1, ny + 1, nz + 1))
        my_sky_grid.z = np.zeros((nx + 1, ny + 1, nz + 1))
        my_sky_grid.dx = dx
        my_sky_grid.dy = dy
        my_sky_grid.dz = dz
        my_sky_grid.nx = nx
        my_sky_grid.ny = ny
        my_sky_grid.nz = nz
        my_sky_grid.X = X
        my_sky_grid.Y = Y
        my_sky_grid.Z = Z
        del dx,dy,dz
        my_sky_grid.set_attributes(my_variables) 
        # We add some random fluctuation to make the grid more interesting
        #for k in range(nz+1):
        #    for j in range(ny+1):
        #        for i in range(nx+1):
        #            if(i==nx): my_sky_grid.x[i, j, k] = X[nx-1] + 0.5 * my_sky_grid.dx
        #            else: my_sky_grid.x[i, j, k] = X[i] - 0.5 * my_sky_grid.dx
        #            if(j==ny): my_sky_grid.y[i, j, k] = Y[ny-1] + 0.5 * my_sky_grid.dy
        #            else: my_sky_grid.y[i, j, k] = Y[j] - 0.5 * my_sky_grid.dy
        #            if(k==nz): my_sky_grid.z[i, j, k] = Z[nz-1] + 0.5 * my_sky_grid.dz
        #            else: my_sky_grid.z[i, j, k] = Z[k] - 0.5 * my_sky_grid.dz

        
        my_sky_grid.x,my_sky_grid.y,my_sky_grid.z=np.meshgrid(np.append(X - 0.5 * my_sky_grid.dx,X[nx-1] + 0.5 * my_sky_grid.dx),
                                                              np.append(Y - 0.5 * my_sky_grid.dy,Y[ny-1] + 0.5 * my_sky_grid.dy),
                                                              np.append(Z - 0.5 * my_sky_grid.dz,Z[nz-1] + 0.5 * my_sky_grid.dz),indexing='ij')
                    
        del X,Y,Z,nx,ny,nz
        
        os.write(1,b"Done...\n")
        
        os.write(1,b"Setting masses and converting their units if necessary...\n")
                    
            #2) projeter sur ces projections les masses correspondantes
        my_railgun.cells_Masse = my_variables.geo_rho_2D * my_railgun.cells_Volume #raw mass
        my_railgun.cells_Masse_ISM = my_railgun.cells_Masse * my_railgun.cells_f_ISM
        my_railgun.cells_Masse_jet = my_railgun.cells_Masse * my_railgun.cells_f_jet
        del my_variables.geo_rho_2D
        
        if(my_variables.frame_the_box):
            my_railgun.cells_Masse_inc={}
            my_railgun.cells_Masse_exc={}
            if(my_variables.frame_ism):
                  my_railgun.cells_Masse_box_ism_inc={}
                  my_railgun.cells_Masse_box_ism_exc={}
            if(my_variables.frame_jet):
                  my_railgun.cells_Masse_box_jet_inc={}
                  my_railgun.cells_Masse_box_jet_exc={}                  
            for box_i,box_el in enumerate(my_variables.x_frame_min):
                my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                         my_variables.x_frame_max[box_i],
                                                         my_variables.y_frame_min[box_i],
                                                         my_variables.y_frame_max[box_i])
                
                my_railgun.cells_Masse_inc[my_dict_name] = my_railgun.cells_Masse * my_railgun.f_box_2D_inc[my_dict_name]
                my_railgun.cells_Masse_exc[my_dict_name] = my_railgun.cells_Masse * my_railgun.f_box_2D_exc[my_dict_name]  
                if(my_variables.frame_ism):
                    my_railgun.cells_Masse_box_ism_inc[my_dict_name] = my_railgun.cells_Masse_ISM * my_railgun.f_box_2D_ism_inc[my_dict_name] 
                    my_railgun.cells_Masse_box_ism_exc[my_dict_name] = my_railgun.cells_Masse_ISM * my_railgun.f_box_2D_ism_exc[my_dict_name]                      
                if(my_variables.frame_jet):
                    my_railgun.cells_Masse_box_jet_inc[my_dict_name] = my_railgun.cells_Masse_jet * my_railgun.f_box_2D_jet_inc[my_dict_name]
                    my_railgun.cells_Masse_box_jet_exc[my_dict_name] = my_railgun.cells_Masse_jet * my_railgun.f_box_2D_jet_exc[my_dict_name]                                          
        if(my_variables.floor_ism):
            my_railgun.cells_Masse_ISM_inc = my_railgun.cells_Masse * np.where(my_railgun.cells_f_ISM<my_variables.threshold_f_ism,0,my_railgun.cells_f_ISM)
            my_railgun.cells_Masse_ISM_exc = my_railgun.cells_Masse * np.where(my_railgun.cells_f_ISM>=my_variables.threshold_f_ism,0,my_railgun.cells_f_ISM)
        if(my_variables.floor_jet):
            my_railgun.cells_Masse_jet_inc = my_railgun.cells_Masse * np.where(my_railgun.cells_f_jet<my_variables.threshold_f_jet,0,my_railgun.cells_f_jet)
            my_railgun.cells_Masse_jet_exc = my_railgun.cells_Masse * np.where(my_railgun.cells_f_jet>=my_variables.threshold_f_jet,0,my_railgun.cells_f_jet)
        
        del my_railgun.cells_f_ISM,my_railgun.cells_f_jet
        
        # convert, if needed, the mass from gram to unit_mass :
        if(my_variables.unit_masse != u.g):
            my_railgun.cells_Masse = (
                my_railgun.cells_Masse*u.g).to_value(my_variables.unit_masse)
            my_railgun.cells_Masse_ISM = (
                my_railgun.cells_Masse_ISM*u.g).to_value(my_variables.unit_masse)
            my_railgun.cells_Masse_jet = (
                my_railgun.cells_Masse_jet*u.g).to_value(my_variables.unit_masse)
            if(my_variables.frame_the_box):     
                for box_i,box_el in enumerate(my_variables.x_frame_min):
                    my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                             my_variables.x_frame_max[box_i],
                                                             my_variables.y_frame_min[box_i],
                                                             my_variables.y_frame_max[box_i])
    
                    my_railgun.cells_Masse_inc[my_dict_name] = (
                        my_railgun.cells_Masse_inc[my_dict_name]*u.g).to_value(my_variables.unit_masse)
                    my_railgun.cells_Masse_exc[my_dict_name] = (
                        my_railgun.cells_Masse_exc[my_dict_name]*u.g).to_value(my_variables.unit_masse)                  
                    if(my_variables.frame_ism):
                        my_railgun.cells_Masse_box_ism_inc[my_dict_name] = (
                            my_railgun.cells_Masse_box_ism_inc[my_dict_name]*u.g).to_value(my_variables.unit_masse)
                        my_railgun.cells_Masse_box_ism_exc[my_dict_name] = (
                            my_railgun.cells_Masse_box_ism_exc[my_dict_name]*u.g).to_value(my_variables.unit_masse)                                          
                    if(my_variables.frame_jet):
                        my_railgun.cells_Masse_box_jet_inc[my_dict_name] = (
                            my_railgun.cells_Masse_box_jet_inc[my_dict_name]*u.g).to_value(my_variables.unit_masse)
                        my_railgun.cells_Masse_box_jet_exc[my_dict_name] = (
                            my_railgun.cells_Masse_box_jet_exc[my_dict_name]*u.g).to_value(my_variables.unit_masse)                                                                  
            if(my_variables.floor_ism):
                my_railgun.cells_Masse_ISM_inc = (my_railgun.cells_Masse_ISM_inc * u.g).to_value(my_variables.unit_masse)
                my_railgun.cells_Masse_ISM_exc = (my_railgun.cells_Masse_ISM_exc * u.g).to_value(my_variables.unit_masse) 
            if(my_variables.floor_jet):
                my_railgun.cells_Masse_jet_inc = (my_railgun.cells_Masse_jet_inc * u.g).to_value(my_variables.unit_masse)
                my_railgun.cells_Masse_jet_exc = (my_railgun.cells_Masse_jet_exc * u.g).to_value(my_variables.unit_masse)
        
        # convert, if needed, the volume from gram to unit_volume :
        if(my_variables.unit_volume != u.cm**3):        
            my_railgun.cells_Volume = (my_railgun.cells_Volume*u.cm**3).to_value(my_variables.unit_volume)
        
        os.write(1,b"Done...\n")
        
        os.write(1,b"Building the 2D masses arrays...\n")
        are_mass_variables = ["masse arrays" in el.replace("_"," ") for el in Sky_Grid.list_of_variables]                
        print("Number of available masses PV output : {0}".format(str(np.sum(are_mass_variables))))
        listtotale = ""
        for iel,el in enumerate(np.array(Sky_Grid.list_of_variables)[are_mass_variables]):
            listtotale+="{}) {}\n".format(str(iel+1),el)
        print("List of available masses PV output:\n{0}\n".format(listtotale))
        del listtotale

        
        return my_variables,my_railgun,my_sky_grid