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
from mod_global_parameters import Variable,Sky_Grid,MV_spectrum,create_box_file_dict_name #,global_variables



""" Bazooka voxelling functions"""    
def voxelling_sky_plane_meshing_bazooka(my_sky_grid,my_bazooka,my_variables):
    # maximum bounds
    v_x_max = np.max(my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.alpha_])
    v_y_max = np.max(my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.Delta_])
    
    v_x_min = np.min(my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.alpha_])
    v_y_min = np.min(my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.Delta_])
    
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
    #my_sky_grid.set_attributes(my_variables) #don't forget to update the SKy_Grid attributes according to the variables current settings
    return my_sky_grid

def put_back(hist,masse,X,Y,h,nx,ny):
    for k in range(nx):
        for j in range(ny):
            pass

def sample_sky_plane_grid_bazooka(my_variables,my_bazooka,my_sky_grid):
    
    my_sky_grid.redefine_nx_ny_nz_arrays(my_variables,my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz)
    extenti = (np.min(my_sky_grid.x),
               np.max(my_sky_grid.x))
    extentj = (np.min(my_sky_grid.y),
               np.max(my_sky_grid.y))
    extentk = (np.min(my_sky_grid.z),
               np.max(my_sky_grid.z))
    
    if(my_variables.voxelling_scheme=='brute'):
            

        x_y_z_to_hist = np.stack((my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.alpha_],
                                  my_bazooka.geo_position_Gamma_Delta_alpha[:, my_variables.Delta_],
                                  my_bazooka.geo_vobs),axis=1)
        os.write(1,b"Computing volume PV...\n")
        hist,edges = np.histogramdd(x_y_z_to_hist,
                                      bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                      range=(extenti, extentj, extentk),
                                      weights=my_bazooka.geo_position_cellsPoints_volume)
        my_sky_grid.volume_arrays[:,:,:] += list(hist)            
        os.write(1,b"Computing total mass PV...\n")
        hist, edges = np.histogramdd(x_y_z_to_hist,
                                     bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                     range=(extenti, extentj, extentk),
                                     weights=my_bazooka.geo_position_cellsPoints_masse)
        #print(np.shape(xedges),np.shape(yedges))
        my_sky_grid.masse_arrays[:,:,:] = list(hist) 
        os.write(1,b"Computing ism mass PV...\n")
        hist, edges = np.histogramdd(x_y_z_to_hist,
                                     bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                     range=(extenti, extentj, extentk),
                                     weights=my_bazooka.geo_position_cellsPoints_masse_ism)
        my_sky_grid.masse_arrays_ism[:,:, :] = list(hist)
        os.write(1,b"Computing jet mass PV...\n")
        hist, edges = np.histogramdd(x_y_z_to_hist,
                                     bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                     range=(extenti, extentj, extentk),
                                     weights=my_bazooka.geo_position_cellsPoints_masse_jet)
        #print(np.shape(my_sky_grid.masse_arrays_jet[:, :, :]),np.shape(hist))
        my_sky_grid.masse_arrays_jet[:, :, :] = list(hist)
        #my_sky_grid.masse_arrays[k, j, :] = np.sum(masse_c)
        #my_sky_grid.masse_arrays_ism[k, j, :] = np.sum(masse_ism_c)
        
        #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)
        if(my_variables.frame_the_box):
            for box_i,box_el in enumerate(my_variables.x_frame_min):
                my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                         my_variables.x_frame_max[box_i],
                                                         my_variables.y_frame_min[box_i],
                                                         my_variables.y_frame_max[box_i])            
                os.write(1,("Computing box masses PV for "+my_dict_name+" ...\n").encode('ascii'))
                hist, edges = np.histogramdd(x_y_z_to_hist,
                                             bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                             range=(extenti, extentj, extentk),
                                             weights=my_bazooka.geo_position_cellsPoints_masse_box_inc[my_dict_name])
                (my_sky_grid.masse_arrays_box_inc[my_dict_name])[:,:, :] = list(hist)
                
                hist, edges = np.histogramdd(x_y_z_to_hist,
                                             bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                             range=(extenti, extentj, extentk),
                                             weights=my_bazooka.geo_position_cellsPoints_masse_box_exc[my_dict_name])
                (my_sky_grid.masse_arrays_box_exc[my_dict_name])[:, :, :] = list(hist)
                        #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)
                        
                if(my_variables.frame_ism):
                    my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                             my_variables.x_frame_max[box_i],
                                                             my_variables.y_frame_min[box_i],
                                                             my_variables.y_frame_max[box_i])            
                    os.write(1,("Computing box ISM masses PV for "+my_dict_name+" ...\n").encode('ascii'))
                    hist, edges = np.histogramdd(x_y_z_to_hist,
                                                 bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                 range=(extenti, extentj, extentk),
                                                 weights=my_bazooka.geo_position_cellsPoints_masse_box_ism_inc[my_dict_name])
                    (my_sky_grid.masse_arrays_box_ism_inc[my_dict_name])[:,:, :] = list(hist)
                    
                    hist, edges = np.histogramdd(x_y_z_to_hist,
                                                 bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                 range=(extenti, extentj, extentk),
                                                 weights=my_bazooka.geo_position_cellsPoints_masse_box_ism_exc[my_dict_name])
                    (my_sky_grid.masse_arrays_box_ism_exc[my_dict_name])[:, :, :] = list(hist)
                            #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)
                if(my_variables.frame_jet):
                    my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                             my_variables.x_frame_max[box_i],
                                                             my_variables.y_frame_min[box_i],
                                                             my_variables.y_frame_max[box_i])            
                    os.write(1,("Computing box JET masses PV for "+my_dict_name+" ...\n").encode('ascii'))
                    hist, edges = np.histogramdd(x_y_z_to_hist,
                                                 bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                 range=(extenti, extentj, extentk),
                                                 weights=my_bazooka.geo_position_cellsPoints_masse_box_jet_inc[my_dict_name])
                    (my_sky_grid.masse_arrays_box_jet_inc[my_dict_name])[:,:, :] = list(hist)
                    
                    hist, edges = np.histogramdd(x_y_z_to_hist,
                                                 bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                 range=(extenti, extentj, extentk),
                                                 weights=my_bazooka.geo_position_cellsPoints_masse_box_jet_exc[my_dict_name])
                    (my_sky_grid.masse_arrays_box_jet_exc[my_dict_name])[:, :, :] = list(hist)
                            #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)                            
                    
        if(my_variables.floor_ism): 
            os.write(1,b"Computing ism floored mass PV...\n")
            hist, edges = np.histogramdd(x_y_z_to_hist,
                                         bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                         range=(extenti, extentj, extentk), 
                                         weights=my_bazooka.geo_position_cellsPoints_masse_ism_inc)
            #print(">>>>>>>>>>",np.shape(my_sky_grid.masse_arrays_ism_inc[:, :, :]),np.shape(hist))
            my_sky_grid.masse_arrays_ism_inc[:,:, :] = list(hist)
            hist, edges = np.histogramdd(x_y_z_to_hist,
                                         bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                         range=(extenti, extentj, extentk),
                                         weights=my_bazooka.geo_position_cellsPoints_masse_ism_exc)
            my_sky_grid.masse_arrays_ism_exc[:, :, :] = list(hist)
        if(my_variables.floor_jet): 
            os.write(1,b"Computing jet floored mass PV...\n")
            hist, edges = np.histogramdd(x_y_z_to_hist,
                                         bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                         range=(extenti, extentj, extentk),
                                         weights=my_bazooka.geo_position_cellsPoints_masse_jet_inc)
            my_sky_grid.masse_arrays_jet_inc[:,:, :] = list(hist)
            hist, edges = np.histogramdd(x_y_z_to_hist,
                                         bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                         range=(extenti, extentj, extentk),
                                         weights=my_bazooka.geo_position_cellsPoints_masse_jet_exc)
            my_sky_grid.masse_arrays_jet_exc[:, :, :] = list(hist)


            
            
        print("Progression (in %): Finishing... 100 % ")
    else:
        raise NotImplementedError("No voxelling_scheme {0} implemented yet! Possibilities are : {1}".format(my_variables.voxelling_scheme,my_variables.voxelling_scheme_list[my_variables.voxelling_method]))
    
    return my_sky_grid


""" Railgun voxelling functions"""
  
def sample_sky_plane_grid_railgun(my_variables,my_railgun,my_sky_grid,my_mv_spectrum,cells):
    
        my_sky_grid.redefine_nx_ny_nz_arrays(my_variables,my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz)
        extenti = (np.min(my_sky_grid.x),
                   np.max(my_sky_grid.x))
        extentj = (np.min(my_sky_grid.y),
                   np.max(my_sky_grid.y))
        extentk = (np.min(my_sky_grid.z),
                   np.max(my_sky_grid.z))
        
        my_mv_spectrum.redefine_nz_arrays(my_variables,my_sky_grid.nz)
        
        #itérer l'opération sur l'ensemble des tranches de phi : d'où l'appelation railgun
        for n in range(0, my_variables.nphi):
            
            progression_bar="Progression (en %): {0:.2f}\n".format(100*n/my_variables.nphi)
            print(progression_bar)
            os.write(1,str.encode(progression_bar))            
            
            #1) projeter les positions des celulles d'une tranche de phi sur le plan du ciel
            geo_angle_phi = n*my_railgun.geo_dangle_phi
            my_railgun.cells_Delta_Alpha[:, my_variables.Delta_-1] = -cells[:, my_variables.r_]*np.cos(my_variables.geo_angle_i) *\
                np.cos(geo_angle_phi) + \
                cells[:, my_variables.z_cyl_]*np.sin(my_variables.geo_angle_i)  # z in 3D
            my_railgun.cells_Delta_Alpha[:, my_variables.alpha_-1] = cells[:,my_variables.r_]*np.sin(geo_angle_phi)  # phi in 3D
            
            
            if(my_variables.voxelling_mode == 'Voxel'):
                my_railgun.geo_vobs = my_variables.velocityVector['v3'] *\
                    np.sin(my_variables.geo_angle_i)*np.sin(geo_angle_phi) -\
                    -my_variables.velocityVector['v1']*np.sin(my_variables.geo_angle_i)*np.cos(geo_angle_phi) -\
                    my_variables.velocityVector['v2']*np.cos(my_variables.geo_angle_i)

            
            
            if(my_variables.voxelling_scheme=='progressive'):
                
                #for h in range(my_sky_grid.nz):
                #d = np.abs(my_railgun.geo_vobs-my_sky_grid.Z[:]) <= my_variables.fac_dz*my_sky_grid.dz
                #anti_d = ~d
                
                x_y_z_to_hist = np.stack((my_railgun.cells_Delta_Alpha[:, my_variables.alpha_-1],
                                          my_railgun.cells_Delta_Alpha[:, my_variables.Delta_-1],
                                          my_railgun.geo_vobs),axis=1)
                del my_railgun.geo_vobs
    
                hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk),
                                              weights=my_railgun.cells_Volume)
                my_sky_grid.volume_arrays[:,:,:] += list(hist) 
                hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk),
                                              weights=my_railgun.cells_Masse)
                #print(np.shape(xedges),np.shape(yedges))
                my_sky_grid.masse_arrays[:,:,:] += list(hist)
                my_mv_spectrum.masse_arrays += list(np.sum(np.sum(my_sky_grid.masse_arrays[:,:,:],axis=0),axis=0))
                hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk),
                                              weights=my_railgun.cells_Masse_ISM)
                my_sky_grid.masse_arrays_ism[:,:,:] += list(hist)
                my_mv_spectrum.masse_arrays_ism += list(np.sum(np.sum(my_sky_grid.masse_arrays_ism[:,:,:],axis=0),axis=0))
                hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk), 
                                              weights=my_railgun.cells_Masse_jet)
                #print(np.shape(my_sky_grid.masse_arrays_jet[:, :, :]),np.shape(hist))
                my_sky_grid.masse_arrays_jet[:, :,:] += list(hist)
                my_mv_spectrum.masse_arrays_jet += list(np.sum(np.sum(my_sky_grid.masse_arrays_jet[:,:,:],axis=0),axis=0))
                #my_sky_grid.masse_arrays[k, j, :] = np.sum(masse_c)
                #my_sky_grid.masse_arrays_ism[k, j, :] = np.sum(masse_ism_c)
                
                #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)
                if(my_variables.frame_the_box): 
                    for box_i,box_el in enumerate(my_variables.x_frame_min):
                        my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                                 my_variables.x_frame_max[box_i],
                                                                 my_variables.y_frame_min[box_i],
                                                                 my_variables.y_frame_max[box_i])                                
                        hist,edges = np.histogramdd(x_y_z_to_hist,
                                                  bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                  range=(extenti, extentj, extentk), 
                                                  weights=my_railgun.cells_Masse_inc[my_dict_name])
                        (my_sky_grid.masse_arrays_box_inc[my_dict_name])[:,:,:] += list(hist)
                        my_mv_spectrum.masse_arrays_box_inc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_inc[my_dict_name])[:,:,:],axis=0),axis=0))
                        
                        hist,edges = np.histogramdd(x_y_z_to_hist,
                                                  bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                  range=(extenti, extentj, extentk), 
                                                  weights=my_railgun.cells_Masse_exc[my_dict_name])
                        (my_sky_grid.masse_arrays_box_exc[my_dict_name])[:, :,:] += list(hist)
                        my_mv_spectrum.masse_arrays_box_exc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_exc[my_dict_name])[:,:,:],axis=0),axis=0))
                            #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)
                        if(my_variables.frame_ism):
                            my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                                     my_variables.x_frame_max[box_i],
                                                                     my_variables.y_frame_min[box_i],
                                                                     my_variables.y_frame_max[box_i])                              
                            hist,edges = np.histogramdd(x_y_z_to_hist,
                                                      bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                      range=(extenti, extentj, extentk), 
                                                      weights=my_railgun.cells_Masse_box_ism_inc[my_dict_name])
                            (my_sky_grid.masse_arrays_box_ism_inc[my_dict_name])[:,:,:] += list(hist)
                            my_mv_spectrum.masse_arrays_box_ism_inc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_ism_inc[my_dict_name])[:,:,:],axis=0),axis=0))
                            
                            hist,edges = np.histogramdd(x_y_z_to_hist,
                                                      bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                      range=(extenti, extentj, extentk), 
                                                      weights=my_railgun.cells_Masse_box_ism_exc[my_dict_name])
                            (my_sky_grid.masse_arrays_box_ism_exc[my_dict_name])[:, :,:] += list(hist)
                            my_mv_spectrum.masse_arrays_box_ism_exc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_ism_exc[my_dict_name])[:,:,:],axis=0),axis=0))
                                #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c)                            
                        if(my_variables.frame_jet):
                            my_dict_name = create_box_file_dict_name(my_variables.x_frame_min[box_i],
                                                                     my_variables.x_frame_max[box_i],
                                                                     my_variables.y_frame_min[box_i],
                                                                     my_variables.y_frame_max[box_i])                              
                            hist,edges = np.histogramdd(x_y_z_to_hist,
                                                      bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                      range=(extenti, extentj, extentk), 
                                                      weights=my_railgun.cells_Masse_box_jet_inc[my_dict_name])
                            (my_sky_grid.masse_arrays_box_jet_inc[my_dict_name])[:,:,:] += list(hist)
                            my_mv_spectrum.masse_arrays_box_jet_inc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_jet_inc[my_dict_name])[:,:,:],axis=0),axis=0))
                            
                            hist,edges = np.histogramdd(x_y_z_to_hist,
                                                      bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                                      range=(extenti, extentj, extentk), 
                                                      weights=my_railgun.cells_Masse_box_jet_exc[my_dict_name])
                            (my_sky_grid.masse_arrays_box_jet_exc[my_dict_name])[:, :,:] += list(hist)
                            my_mv_spectrum.masse_arrays_box_jet_exc[my_dict_name] += list(np.sum(np.sum((my_sky_grid.masse_arrays_box_jet_exc[my_dict_name])[:,:,:],axis=0),axis=0))
                                #my_sky_grid.masse_arrays_jet[k, j, :] = np.sum(masse_jet_c) 
        
                            
                if(my_variables.floor_ism): 
                    hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk), 
                                              weights=my_railgun.cells_Masse_ISM_inc)
                    #print(">>>>>>>>>>",np.shape(my_sky_grid.masse_arrays_ism_inc[:, :, :]),np.shape(hist))
                    my_sky_grid.masse_arrays_ism_inc[:,:,:] += list(hist)
                    my_mv_spectrum.masse_arrays_ism_inc += list(np.sum(np.sum(my_sky_grid.masse_arrays_ism_inc,axis=0),axis=0))
                    hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk), 
                                              weights=my_railgun.cells_Masse_ISM_exc)
                    my_sky_grid.masse_arrays_ism_exc[:, :,:] += list(hist)
                    my_mv_spectrum.masse_arrays_ism_exc += list(np.sum(np.sum(my_sky_grid.masse_arrays_ism_exc,axis=0),axis=0))
                if(my_variables.floor_jet): 
                    hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk), 
                                              weights=my_railgun.cells_Masse_jet_inc)
                    my_sky_grid.masse_arrays_jet_inc[:,:,:] += list(hist)
                    my_mv_spectrum.masse_arrays_jet_inc += list(np.sum(np.sum(my_sky_grid.masse_arrays_jet_inc,axis=0),axis=0))
                    hist,edges = np.histogramdd(x_y_z_to_hist,
                                              bins=[my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz], 
                                              range=(extenti, extentj, extentk), 
                                              weights=my_railgun.cells_Masse_jet_exc)
                    my_sky_grid.masse_arrays_jet_exc[:, :,:] += list(hist)
                    my_mv_spectrum.masse_arrays_jet_exc += list(np.sum(np.sum(my_sky_grid.masse_arrays_jet_exc,axis=0),axis=0))
            
            #3) échantillonner ces masses projetées sur la grille d'échantillonnage du ciel
           

        del my_railgun.cells_Masse,my_railgun.cells_Masse_ISM, my_railgun.cells_Masse_jet, my_railgun.cells_Volume
        if(my_variables.frame_the_box): 
            del my_railgun.cells_Masse_inc,my_railgun.cells_Masse_exc
            if(my_variables.frame_ism):
                del my_railgun.cells_Masse_box_ism_inc,my_railgun.cells_Masse_box_ism_exc
            if(my_variables.frame_jet):
                del my_railgun.cells_Masse_box_jet_inc,my_railgun.cells_Masse_box_jet_exc                
        if(my_variables.floor_ism): 
            del my_railgun.cells_Masse_ISM_inc,my_railgun.cells_Masse_ISM_exc
        if(my_variables.floor_jet): 
            del my_railgun.cells_Masse_jet_inc,my_railgun.cells_Masse_jet_exc
            
        print("\nMinimal mass found : {0} {1}".format(
            np.str(np.amin(my_sky_grid.masse_arrays)), str(my_variables.unit_masse)))
        print("Maximal mass found : {0} {1}".format(
            np.str(np.amax(my_sky_grid.masse_arrays)), str(my_variables.unit_masse)))
            
        #4) purger les variables de positions projetées et de masses projetées pour ne pas
        #   surcharger la mémoire
        
        os.write(1,b"Progression (in %): Finishing... 100 % \n")
    
        return my_variables,my_railgun,my_sky_grid,my_mv_spectrum