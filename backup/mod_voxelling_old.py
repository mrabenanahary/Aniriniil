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
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                if(i==nx): my_sky_grid.x[i, j, k] = X[nx-1] + 0.5 * my_sky_grid.dx
                else: my_sky_grid.x[i, j, k] = X[i] - 0.5 * my_sky_grid.dx
                if(j==ny): my_sky_grid.y[i, j, k] = Y[ny-1] + 0.5 * my_sky_grid.dy
                else: my_sky_grid.y[i, j, k] = Y[j] - 0.5 * my_sky_grid.dy
                if(k==nz): my_sky_grid.z[i, j, k] = Z[nz-1] + 0.5 * my_sky_grid.dz
                else: my_sky_grid.z[i, j, k] = Z[k] - 0.5 * my_sky_grid.dz
    return my_sky_grid


def sample_sky_plane_grid(my_variables,my_bazooka,my_sky_grid):
    
    if(my_variables.voxelling_scheme=='brute'):
        temp_x_y = np.copy(my_bazooka.geo_position_Gamma_Delta_alpha)
        temp_masse = np.copy(my_bazooka.geo_position_cellsPoints_masse)
        temp_masse_ism = np.copy(my_bazooka.geo_position_cellsPoints_masse_ism)
        temp_masse_jet = np.copy(my_bazooka.geo_position_cellsPoints_masse_jet)
        if(my_variables.frame_the_box): 
            temp_masse_box_inc = np.copy(my_bazooka.geo_position_cellsPoints_masse_box_inc)
            temp_masse_box_exc = np.copy(my_bazooka.geo_position_cellsPoints_masse_box_exc)
        if(my_variables.voxelling_mode == 'Voxel'):
           temp_vobs = np.copy(my_bazooka.geo_vobs)
        test_k = 0
        
        my_sky_grid.redefine_nx_ny_nz_arrays(my_variables,my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz)
        
        for k in range(my_sky_grid.nx):
            a = np.abs(temp_x_y[:, my_variables.alpha_]-my_sky_grid.X[k]) <= my_variables.fac_dx*my_sky_grid.dx
            if(True in a):
                anti_a = ~a
        
                b = temp_x_y[a]
                masse_b = temp_masse[a]
                masse_ism_b = temp_masse_ism[a]
                masse_jet_b = temp_masse_jet[a]
                if(my_variables.frame_the_box): 
                    masse_box_inc_b = temp_masse_box_inc[a]
                    masse_box_exc_b = temp_masse_box_exc[a]
        
                temp_x_y = temp_x_y[anti_a]
                temp_masse = temp_masse[anti_a]
                temp_masse_ism = temp_masse_ism[anti_a]
                temp_masse_jet = temp_masse_jet[anti_a]
                if(my_variables.frame_the_box): 
                    temp_masse_box_inc = temp_masse_box_inc[anti_a]
                    temp_masse_box_exc = temp_masse_box_exc[anti_a]
        
                if(my_variables.voxelling_mode == 'Voxel'):
                    b_vobs = temp_vobs[a]
                    temp_vobs = temp_vobs[anti_a]
                for j in range(my_sky_grid.ny):
                    c = np.abs(b[:, my_variables.Delta_]-my_sky_grid.Y[j]) <= my_variables.fac_dy*my_sky_grid.dy
                    if(True in c):
                        anti_c = ~c
                        # for i in range(len(masse_b[c])):
                        #    masse_arrays[k,j,0] = masse_arrays[k,j,0] + (masse_b[c])[i]
                        if(my_variables.voxelling_mode == 'Masse'):
                            my_sky_grid.masse_arrays[k, j, 0] = np.sum(masse_b[c])
                            my_sky_grid.masse_arrays_ism[k, j, 0] = np.sum(masse_ism_b[c])
                            my_sky_grid.masse_arrays_jet[k, j, 0] = np.sum(masse_jet_b[c])
                            if(my_variables.frame_the_box): 
                                my_sky_grid.masse_arrays_box_inc[k, j, 0] = np.sum(masse_box_inc_b[c])
                                my_sky_grid.masse_arrays_box_exc[k, j, 0] = np.sum(masse_box_exc_b[c])
                        elif(my_variables.voxelling_mode == 'Voxel'):
                            masse_c = masse_b[c]
                            masse_ism_c = masse_ism_b[c]
                            masse_jet_c = masse_jet_b[c]
                            if(my_variables.frame_the_box):  
                                masse_box_inc_c = masse_box_inc_b[c]
                                masse_box_exc_c = masse_box_exc_b[c]
                            c_vobs = b_vobs[c]
                            for h in range(my_sky_grid.nz):
                                d = np.abs(c_vobs-my_sky_grid.Z[h]) <= my_variables.fac_dz*my_sky_grid.dz
                                if(True in d):
                                    anti_d = ~d
                                    my_sky_grid.masse_arrays[k, j, h] = np.sum(masse_c[d])
                                    my_sky_grid.masse_arrays_ism[k, j, h] = np.sum(masse_ism_c[d])
                                    my_sky_grid.masse_arrays_jet[k, j, h] = np.sum(masse_jet_c[d])
                                    if(my_variables.frame_the_box): 
                                        my_sky_grid.masse_arrays_box_inc[k, j, h] = np.sum(masse_box_inc_c[d])
                                        my_sky_grid.masse_arrays_box_exc[k, j, h] = np.sum(masse_box_exc_c[d])
                                    c_vobs = c_vobs[anti_d]
                                    masse_c = masse_c[anti_d]
                                    masse_ism_c = masse_ism_c[anti_d]
                                    masse_jet_c = masse_jet_c[anti_d]
                                    if(my_variables.frame_the_box): 
                                        masse_box_inc_c = masse_box_inc_c[anti_d]
                                        masse_box_exc_c = masse_box_exc_c[anti_d]
                                    del d, anti_d
                                else:
                                    del d
                                #if(voxelling_mode == 'Voxel'):
                                if(((k*my_sky_grid.ny+j)*my_sky_grid.nz+h)/(my_sky_grid.nx*my_sky_grid.ny*my_sky_grid.nz)-test_k > my_variables.each_percent_frac):
                                    #print("k*ny+j= ",k*ny+j)
                                    print("Progression (en %): {0:.2f} ".format(
                                        100*((k*my_sky_grid.ny+j)*my_sky_grid.nz+h)/(my_sky_grid.nx*my_sky_grid.ny*my_sky_grid.nz)))
                                    test_k = ((k*my_sky_grid.ny+j)*my_sky_grid.nz+h)/(my_sky_grid.nx*my_sky_grid.ny*my_sky_grid.nz)
                            b_vobs = b_vobs[anti_c]
                            del masse_c, c_vobs, masse_ism_c, masse_jet_c
                            if(my_variables.frame_the_box): 
                                del masse_box_inc_c,masse_box_exc_c
                        b = b[anti_c]
                        masse_b = masse_b[anti_c]
                        masse_ism_b = masse_ism_b[anti_c]
                        masse_jet_b = masse_jet_b[anti_c]
                        if(my_variables.frame_the_box): 
                            masse_box_inc_b = masse_box_inc_b[anti_c]
                            masse_box_exc_b = masse_box_exc_b[anti_c]
                        # if(len(temp_x_y)==0): break
                        del c, anti_c
                    else:
                        del c
                    if((k*my_sky_grid.ny+j)/(my_sky_grid.nx*my_sky_grid.ny)-test_k > my_variables.each_percent_frac):
                        #print("k*ny+j= ",k*ny+j)
                        print("Progression (en %): {0:.2f} ".format(
                            100*(k*my_sky_grid.ny+j)/(my_sky_grid.nx*my_sky_grid.ny)))
                        #sys.stdout.write('\033[2K\033[1G')
                        test_k = (k*my_sky_grid.ny+j)/(my_sky_grid.nx*my_sky_grid.ny)
                        #if(voxelling_mode == 'Masse'):
                del b, masse_b, masse_ism_b, masse_jet_b
                if(my_variables.frame_the_box): 
                    del masse_box_inc_b,masse_box_exc_b
            else:
                del a
            
            if(k/my_sky_grid.nx-test_k > my_variables.each_percent_frac):
                #print("k*ny+j= ",k*ny+j)
                print("Progression (en %): {0:.2f} ".format(
                    100*k/my_sky_grid.nx))
                #sys.stdout.write('\033[2K\033[1G')
                test_k = k/my_sky_grid.nx
        print("Progression (in %): Finishing... 100 % ")
           
        print("Remaining cells at the end of voxelling : ", np.shape(temp_masse))
        del temp_masse, temp_masse_ism, temp_masse_jet
        if(my_variables.frame_the_box): del temp_masse_box_inc, temp_masse_box_exc
    else:
        raise NotImplementedError("No voxelling_scheme {0} implemented yet! Possibilities are : {1}".format(my_variables.voxelling_scheme,my_variables.voxelling_scheme_list[my_variables.voxelling_method]))
    return my_sky_grid


""" Railgun voxelling functions"""
  
#def voxelling_sky_plane_meshing_bazooka(my_sky_grid,my_railgun,my_variables):
    # maximum bounds
#    return my_sky_grid
