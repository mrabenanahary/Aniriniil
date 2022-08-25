#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""
global globals_has_been_called



import sys
# sys.path.append('./ParaView-5.8.0-RC1-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages')
#from paraview.simple import *
#from math import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import vtk as v
#from vtk import *
from vtk.util.numpy_support import vtk_to_numpy  # thats what you need
from cycler import cycler
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import csv
import os
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
import matplotlib.text as mtext
from math import ceil
from astropy import units as u
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

#importation des packages et modules durhampy personnalisées crées
#from src.postProcessing.plotting import *
#from src.postProcessing.calculus import *

#global global_variables

from mod_global_parameters import * #"""DONE!"""
from mod_read_params import * #"""DONE!"""
from mod_exceptions import *  #"""DONE!"""
from mod_bazooka import *  #"""DONE!"""
from mod_maths import *  #"""DONE!"""
from mod_voxelling import *  #"""DONE!"""
from mod_output import *  #"""DONE!"""

#from src.physics.convert import *

#variables globales :

#global raie_dict, H2O_raie_dict, OI_raie_dict, pH2O_exceptions


"====================PREAMBULE======================"
"""Paramètres Matplotlib"""
mpl.rcdefaults()

Latex=True
font = {'family' : 'serif',
        'size'   : 17,
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


"""============================================================================="""


print("""0) PARAMETRES""")
print("Reading parameters...")

#a) Create a local Variable and set it to default values
my_variables = Variable() 

#b) Udpate my_variables variables with the parameters from input file:
#       - sys.argv contain the vtk input and output file path
#       - filename contains the parameters input  file
my_variables = read_Variable_file(sys.argv,my_variables,print_imported_file=False)

print("All parameters are...")
for item in Variable.list_of_variables : print(item+' = '+str(getattr(my_variables,item)))
print("End of reading parameters...")


"""============================================================================="""


print("""I) Read and get touch on .vtu file data and points in python """\
+"Reading input...")

CellsData, numberOfCells, ndimProblem, numberOfVertixes, numberOfVertixesUnderDimension = \
    extract_vtk_CellsData(my_variables)

#CellsData = extract_vtk_CellsData(my_variables,return_data_props=False)



print("Done...")
print("""II) 2D-part computations""")
print("Doing step 1 over 2...")

my_variables.velocityVector['v1'],\
my_variables.velocityVector['v2'],\
my_variables.velocityVector['v3'],\
my_variables.rhoTracers['ISM'],\
my_variables.rhoTracers['jet'] = \
    (read.extract(CellsData, 'v1', attribute_mode=my_variables.geo_2D_attribute_mode),
     read.extract(CellsData, 'v2', attribute_mode=my_variables.geo_2D_attribute_mode),
     read.extract(CellsData, 'v3', attribute_mode=my_variables.geo_2D_attribute_mode),
     read.extract(CellsData, 'tracer_ism', attribute_mode=my_variables.geo_2D_attribute_mode),
     read.extract(CellsData, 'tracer_jet', attribute_mode=my_variables.geo_2D_attribute_mode))

#deal with very small unphysical negative values of mass due to AMRVAC numerical method
my_variables.rhoTracers['ISM'] = np.maximum(my_variables.rhoTracers['ISM'],0)
my_variables.rhoTracers['jet'] = np.maximum(my_variables.rhoTracers['jet'],0)

print([np.shape(my_variables.velocityVector[el]) for el in ['v1','v2','v3']])
print([np.shape(my_variables.rhoTracers[el]) for el in ['ISM','jet']])

is_mode_implemented(my_variables.voxelling_method,Variable.voxelling_method_list,mode_name='voxelling_method')

"""============================================================================="""
"""The original bazooka-brute approach implemented for geoCompiler"""
if(my_variables.voxelling_method=='bazooka'):

    is_mode_implemented(my_variables.voxelling_scheme,Variable.voxelling_scheme_list[my_variables.voxelling_method],mode_name='voxelling_scheme')
    if(my_variables.voxelling_scheme=='brute'):
        if(my_variables.voxelling_mode=='Voxel'):
            print("Doing step 2 over 2 with voxelling_method : \'{0}\'...".format(my_variables.voxelling_method))
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
            my_bazooka.set_attributes(CellsData,my_variables) #don't forget to update the bazooka attributes according to the vtk input file
            #my_bazooka.output = None
            #Bazooka.output = None
            #del my_bazooka.output, Bazooka.output
            
            #construct the 2D bullets for the 3D bazooka loading
            if(my_variables.frame_the_box):
                my_bazooka, points2D,f_box_2D_inc,f_box_2D_exc = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
            else:
                my_bazooka, points2D = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
            
            my_variables.geo_rho_2D = read.extract(
                CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
            
            print("Done within 2 steps...")
            print("""III) 3D-part computations""")
        
            #load the bazooka for 3D grid computation: here every 3D domain's points are 
            #loaded before they are shot on the sky plane grid
            # Advantage : at each step, every points are treated at the same time
            # Drawback : instancing the 
            if(my_variables.frame_the_box) : my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D,f_box_2D_inc=f_box_2D_inc,f_box_2D_exc=f_box_2D_exc)
            else : my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D)
            
            """DEBUG SCRIPT TO TEST IF 3D COMPUTATIONS HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_3D_computations.py").read())
            """END OF THE DEBUG SCRIPT"""
            
        
            """============================================================================="""
        
            print("""IV) CC points projection on the sky plane""")
            
            
            my_bazooka = full_shoot_with_the_bazooka(my_bazooka,my_variables)
            
            """DEBUG SCRIPT TO TEST IF CC PROJECTIONS ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_CC_projections_sky.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            """============================================================================="""
        
            print("""V) Voxelling """)
            print("""Meshing the sky plane...""")
            
            my_sky_grid = Sky_Grid()
             #don't forget to update the SKy_Grid attributes according to the variables current settings
            my_sky_grid = voxelling_sky_plane_meshing_bazooka(my_sky_grid,my_bazooka,my_variables)
            my_sky_grid.set_attributes(my_variables)
            
            """DEBUG SCRIPT TO TEST IF VOXELLING GRID MESHING ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_voxelling_mesh.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            print("""Sampling projected data on the sky plane grid...""")
            my_sky_grid = sample_sky_plane_grid(my_variables,my_bazooka,my_sky_grid)
            
            """DEBUG SCRIPT TO TEST IF VOXELLING GRID SAMPLING ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_voxelling_sample.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            code_output_to_vtk(my_variables,my_sky_grid)
            
        elif(my_variables.voxelling_mode=='Flow_rates'):
            print("Doing step 2 over 2 with voxelling_method : \'{0}\'...".format(my_variables.voxelling_method))
            my_transverse_layer = Transverse_Layer()
            my_transverse_layer.numberOfCells,\
            my_transverse_layer.ndimProblem,\
            my_transverse_layer.numbeOfVertixes,\
            my_transverse_layer.numberOfVertixesUnderDimension = \
                (numberOfCells, 
                 ndimProblem, 
                 numberOfVertixes, 
                 numberOfVertixesUnderDimension)
            #free some unnecessarily used memory
            my_transverse_layer.set_attributes(CellsData,my_variables) #don't forget to update the bazooka attributes according to the vtk input file
            #my_transverse_layer.output = None
            #Bazooka.output = None
            #del my_transverse_layer.output, Bazooka.output
            
            #construct the 2D bullets for the 3D bazooka loading
            
            print("Done within 2 steps...")
            print("""III) Reference layering""")
            
            
            my_variables.geo_rho_2D = read.extract(
                CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
            my_transverse_layer, points2D, r2 = layer_cells_surfaces(CellsData,my_variables,my_transverse_layer)
            
            print("""IV) Momentum fluxes computing""")
            
            my_transverse_layer.fluxes_surface_center,\
            my_transverse_layer.fluxes_surface_inclination,\
            my_transverse_layer.fluxes_surface_r_z,\
            my_transverse_layer.fluxes_surface_radius_z_inclination,\
            my_transverse_layer.fluxes = \
            compute_momentum_flux_through_circle(my_variables,points2D,r2,
                                                 my_transverse_layer,
                                                 center=my_variables.flux_surface_center,
                                                 r_z_inclination=my_variables.flux_surface_r_z_inclination)
            
            print("""V) Output""")
            code_output_impulsion(my_variables,my_transverse_layer)
            
        else: raise NotImplementedError("No avalaible voxelling mode selected!")
    else: raise NotImplementedError("No avalaible voxelling scheme selected!")
  
    
    """<<< Succesfully checked until this point!"""
    """============================================================================="""
    
    
    """The railgun progressive approach implemented for geoCompiler"""
    """<!> DRAWBACK :  MUCH SLOWER THAN BAZOOKA !"""
    """BIG ADVANTAGE ! MUCH LESS MEMORY-CONSUMING !   """
    
elif(my_variables.voxelling_method=='railgun'): 
    is_mode_implemented(my_variables.voxelling_scheme,Variable.voxelling_scheme_list[my_variables.voxelling_method],mode_name='voxelling_scheme')
    if(my_variables.voxelling_scheme=='progressive'):
        my_railgun = Railgun()
        my_railgun.numberOfCells,\
        my_railgun.ndimProblem,\
        my_railgun.numbeOfVertixes,\
        my_railgun.numberOfVertixesUnderDimension = \
            (numberOfCells, 
             ndimProblem, 
             numberOfVertixes, 
             numberOfVertixesUnderDimension)
        #free some unnecessarily used memory
        my_railgun.set_attributes(CellsData,my_variables) #don't forget to update the bazooka attributes according to the vtk input file
        my_railgun.output = None
        Bazooka.output = None
        del my_railgun.output, Bazooka.output
        
        #calculer les positions 2D-cylindriques, constantes dans le domaine 2D
        CellsNodes = np.empty((my_railgun.numberOfCells,my_railgun.ndimProblem*2,))
        {CellsData.GetCellBounds(i,CellsNodes[i]) for i in range(0,my_railgun.numberOfCells)}
        for idim in range(my_railgun.ndimProblem):
            my_railgun.nodes_r_z_phi[:,:,idim]=np.tile(CellsNodes[:,2*idim:2*idim+2],my_railgun.ndimProblem-1)
            if(idim==my_variables.r_): r_offset = np.amin(my_railgun.nodes_r_z_phi[:,:,idim])
        #r_offset = np.amin(np.tile(CellsNodes[:,2*my_variables.r_:2*my_variables.r_+2],my_railgun.ndimProblem-1))
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
        cells = 0.5*(r_z_phi_max+r_z_phi_min)
        my_railgun.bound_radius = max(np.abs(np.max(r_z_phi_max[:,my_variables.r_])),np.abs(np.min(r_z_phi_max[:,my_variables.r_])))
        my_railgun.bound_z = max(np.abs(np.max(r_z_phi_max[:,my_variables.z_cyl_])),np.abs(np.min(r_z_phi_max[:,my_variables.z_cyl_])))
        
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
            my_railgun.f_box_2D_inc = np.where(((cells[:,my_variables.x_frame]>=my_variables.x_frame_min)&(cells[:,my_variables.x_frame]<=my_variables.x_frame_max))&((cells[:,my_variables.y_frame]>=my_variables.y_frame_min)&(cells[:,my_variables.y_frame]<=my_variables.y_frame_max)),True,False)
            my_railgun.f_box_2D_exc = ~my_railgun.f_box_2D_inc
            my_railgun.f_box_2D_inc = np.where(my_railgun.f_box_2D_inc,1.0,0.0)
            my_railgun.f_box_2D_exc = np.where(my_railgun.f_box_2D_exc,1.0,0.0)
        
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
        del CellsNodes
        
        #fractions nécessaires pour le calcul des 
        my_railgun.cells_f_ISM = my_variables.rhoTracers['ISM']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])
        my_railgun.cells_f_jet = my_variables.rhoTracers['jet']/(my_variables.rhoTracers['ISM']+my_variables.rhoTracers['jet'])

        del my_variables.rhoTracers['ISM'],my_variables.rhoTracers['jet']
        
        #ne pas oublier d'extraire le champ de densite
        my_variables.geo_rho_2D = read.extract(
            CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
        
        # construire les limites de la grille du plan du ciel
        # - abs(alpha) inferior or equal to Rmax
        # - Delta inforior or equal to sqrt(Rmax^2+Hmax^2) and superior or
        #   equal to 
        v_x_max = my_variables.bound_offset_frac*np.max(my_railgun.bound_radius)
        v_x_min = -v_x_max
        
        if(my_variables.geo_angle_i>=0.0):
            v_y_max = np.sqrt(my_railgun.bound_radius**2+my_railgun.bound_z**2) *\
                np.cos((np.pi/2.0)-my_variables.geo_angle_i-np.arccos(1.0/(np.sqrt(1.0+(my_railgun.bound_radius/my_railgun.bound_z)**2))))
            v_y_min = my_railgun.bound_radius * np.sin(np.pi/2.0-my_variables.geo_angle_i)
            if(v_y_min>0.0): v_y_min = - v_y_min
            
        else:
            v_y_min = - np.sqrt(my_railgun.bound_radius**2+my_railgun.bound_z**2) *\
                np.cos((np.pi/2.0)-my_variables.geo_angle_i-np.arccos(1.0/(np.sqrt(1.0+(my_railgun.bound_radius/my_railgun.bound_z)**2))))
            v_y_max = my_railgun.bound_radius * np.sin(np.pi/2.0-my_variables.geo_angle_i)
            if(v_y_max>0.0): v_y_max = - v_y_max
        v_y_max = my_variables.bound_offset_frac*v_y_max
        v_y_min = my_variables.bound_offset_frac*v_y_min
            
        #instance the sky_grid
        my_sky_grid = Sky_Grid()
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
        for k in range(nz+1):
            for j in range(ny+1):
                for i in range(nx+1):
                    if(i==nx): my_sky_grid.x[i, j, k] = X[nx-1] + 0.5 * my_sky_grid.dx
                    else: my_sky_grid.x[i, j, k] = X[i] - 0.5 * my_sky_grid.dx
                    if(j==ny): my_sky_grid.y[i, j, k] = Y[ny-1] + 0.5 * my_sky_grid.dy
                    else: my_sky_grid.y[i, j, k] = Y[j] - 0.5 * my_sky_grid.dy
                    if(k==nz): my_sky_grid.z[i, j, k] = Z[nz-1] + 0.5 * my_sky_grid.dz
                    else: my_sky_grid.z[i, j, k] = Z[k] - 0.5 * my_sky_grid.dz
        
        #itérer l'opération sur l'ensemble des tranches de phi : d'où l'appelation railgun
        for n in range(0, my_variables.nphi):
            
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
            
            #2) projeter sur ces projections les masses correspondantes
            my_railgun.cells_Masse = my_variables.geo_rho_2D * my_railgun.cells_Volume #raw mass
            my_railgun.cells_Masse_ISM = my_railgun.cells_Masse * my_railgun.cells_f_ISM
            my_railgun.cells_Masse_jet = my_railgun.cells_Masse * my_railgun.cells_f_jet
            if(my_variables.frame_the_box):
                my_railgun.cells_Masse_inc = my_railgun.cells_Masse * my_railgun.f_box_2D_inc
                my_railgun.cells_Masse_exc = my_railgun.cells_Masse * my_railgun.f_box_2D_exc
            

            
            
            # convert, if needed, the mass from gram to unit_mass :
            if(my_variables.unit_masse != u.g):
                my_railgun.cells_Masse = (
                    my_railgun.cells_Masse*u.g).to_value(my_variables.unit_masse)
                my_railgun.cells_Masse_ISM = (
                    my_railgun.cells_Masse_ISM*u.g).to_value(my_variables.unit_masse)
                my_railgun.cells_Masse_jet = (
                    my_railgun.cells_Masse_jet*u.g).to_value(my_variables.unit_masse)
                if(my_variables.frame_the_box):
                    my_railgun.cells_Masse_inc = (
                        my_railgun.cells_Masse_inc*u.g).to_value(my_variables.unit_masse)
                    my_railgun.cells_Masse_exc = (
                        my_railgun.cells_Masse_exc*u.g).to_value(my_variables.unit_masse)


            print("\nMinimal mass found : {0} {1}".format(
                np.str(np.amin(my_railgun.cells_Masse)), str(my_variables.unit_masse)))
            print("Maximal mass found : {0} {1}".format(
                np.str(np.amax(my_railgun.cells_Masse)), str(my_variables.unit_masse)))
            
            #3) échantillonner ces masses projetées sur la grille d'échantillonnage du ciel
            
            if(my_variables.voxelling_scheme=='progressive'):
                temp_x_y = np.copy(my_railgun.cells_Delta_Alpha)
                temp_masse = np.copy(my_railgun.cells_Masse)
                temp_masse_ism = np.copy(my_railgun.cells_Masse_ISM)
                temp_masse_jet = np.copy(my_railgun.cells_Masse_jet)
                if(my_variables.frame_the_box): 
                    temp_masse_box_inc = np.copy(my_railgun.cells_Masse_inc)
                    temp_masse_box_exc = np.copy(my_railgun.cells_Masse_exc)
                if(my_variables.voxelling_mode == 'Voxel'):
                   temp_vobs = np.copy(my_railgun.geo_vobs)
                #test_k = 0
                
                my_sky_grid.redefine_nx_ny_nz_arrays(my_variables,my_sky_grid.nx,my_sky_grid.ny,my_sky_grid.nz)
                
                for k in range(my_sky_grid.nx):
                    a = np.abs(temp_x_y[:, my_variables.alpha_-1]-my_sky_grid.X[k]) <= my_variables.fac_dx*my_sky_grid.dx
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
                            c = np.abs(b[:, my_variables.Delta_-1]-my_sky_grid.Y[j]) <= my_variables.fac_dy*my_sky_grid.dy
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

                        del b, masse_b, masse_ism_b, masse_jet_b
                        if(my_variables.frame_the_box): 
                            del masse_box_inc_b,masse_box_exc_b
                    else:
                        del a

                
                del temp_masse_ism, temp_masse_jet
                if(my_variables.frame_the_box): del temp_masse_box_inc, temp_masse_box_exc
            else:
                raise NotImplementedError("No voxelling_scheme {0} implemented yet! Possibilities are : {1}".format(my_variables.voxelling_scheme,my_variables.voxelling_scheme_list[my_variables.voxelling_method]))
            
            print("Progression (en %): {0:.2f} ".format(
                            100*n/my_variables.nphi))
            
            #4) purger les variables de positions projetées et de masses projetées pour ne pas
            #   surcharger la mémoire
            pass #print(n)
        print("Progression (in %): Finishing... 100 % ")
                   
        print("Remaining cells at the end of voxelling : ", np.shape(temp_masse))
        del temp_masse
    else: raise NotImplementedError("No avalaible voxelling scheme selected!")
    code_output_to_vtk(my_variables,my_sky_grid) 
    """============================================================================="""
else: raise NotImplementedError("No avalaible voxelling method selected!")