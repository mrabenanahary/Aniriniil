#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""
from __future__ import print_function
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
#os.environ["PATH"] += os.pathsep + '/obs/mrabenanahary/usr/local/texlive/bin/x86_64-linux'
#print(os.getenv("PATH"))
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
import shutil

#importation des packages et modules durhampy personnalisées crées
#from src.postProcessing.plotting import *
#from src.postProcessing.calculus import *

#global global_variables

from mod_global_parameters import * #"""DONE!"""
from mod_read_params import * #"""DONE!"""
from mod_exceptions import *  #"""DONE!"""
from mod_bazooka import *  #"""DONE!"""
from mod_railgun import *  #"""DONE!"""
from mod_maths import *  #"""DONE!"""
from mod_voxelling import *  #"""DONE!"""
from mod_output import *  #"""DONE!"""
from mod_plot import *


#from src.physics.convert import *

#variables globales :

#global raie_dict, H2O_raie_dict, OI_raie_dict, pH2O_exceptions


"====================PREAMBULE======================"
"""Paramètres Matplotlib"""
mpl.rcdefaults()




"""============================================================================="""


os.write(1,b"0) PARAMETRES\n")
os.write(1,b"Reading parameters...\n")

#a) Create a local Variable and set it to default values
my_variables = Variable() 

#b) Udpate my_variables variables with the parameters from input file:
#       - sys.argv contain the vtk input and output file path
#       - filename contains the parameters input  file
my_variables = read_Variable_file(sys.argv,my_variables,print_imported_file=False)


os.system('mkdir {0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')))
#os.system('cp {0} {1}/.'.format(my_variables.input_par,((my_variables.output_file_vtk).replace('.','')).replace('/','')))
shutil.copyfile('{0}'.format(my_variables.input_par), '{0}/{1}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/',''),(my_variables.input_par).split('/')[-1]))
#my_variables.set_units()

os.write(1,b"All parameters are...\n")
for item in Variable.list_of_variables : print(item+' = '+str(getattr(my_variables,item)))
os.write(1,b"End of reading parameters...\n")


"""============================================================================="""


os.write(1,b"I) Read and get touch on .vtu file data and points in python\n")
os.write(1,b"Reading input...\n")

CellsData, numberOfCells, ndimProblem, numberOfVertixes, numberOfVertixesUnderDimension = \
    extract_vtk_CellsData(my_variables)

#CellsData = extract_vtk_CellsData(my_variables,return_data_props=False)

"""#======Plotting when needed ! ======"""
if(my_variables.plotting):
    plot_maps(my_variables, CellsData, numberOfCells, ndimProblem, numberOfVertixes, numberOfVertixesUnderDimension)



os.write(1,b"Done...\n")

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

my_variables.rhoTracers['ISM'] = my_variables.rho_ism_fac * my_variables.rhoTracers['ISM']
my_variables.rhoTracers['jet'] = my_variables.rho_jet_fac * my_variables.rhoTracers['jet']


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
            os.write(1,b"II) 2D-part computations\n")
            print("Doing step 1 over 1 with voxelling_method : \'{0}\'...".format(my_variables.voxelling_method))
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
                if((my_variables.frame_ism==True)&(my_variables.frame_jet==False)):
                   my_bazooka, points2D,f_box_2D_inc,f_box_2D_exc,\
                       f_box_2D_ism_inc,f_box_2D_ism_exc = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
                if((my_variables.frame_ism==False)&(my_variables.frame_jet==True)):
                   my_bazooka, points2D,f_box_2D_inc,f_box_2D_exc,\
                       f_box_2D_jet_inc,f_box_2D_jet_exc = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
                if((my_variables.frame_ism==True)&(my_variables.frame_jet==True)):
                   my_bazooka, points2D,f_box_2D_inc,f_box_2D_exc,\
                       f_box_2D_ism_inc,f_box_2D_ism_exc,\
                       f_box_2D_jet_inc,f_box_2D_jet_exc = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
                if((my_variables.frame_ism==False)&(my_variables.frame_jet==False)):                
                   my_bazooka, points2D,f_box_2D_inc,f_box_2D_exc = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
            else:
                my_bazooka, points2D = construct_3D_bullets_for_bazooka(CellsData,my_variables,my_bazooka)
            
            my_variables.geo_rho_2D = read.extract(
                CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
            
            os.write(1,b"Done within 1 step...\n")
            os.write(1,b"III) 3D-part computations\n")
        
            #load the bazooka for 3D grid computation: here every 3D domain's points are 
            #loaded before they are shot on the sky plane grid
            # Advantage : at each step, every points are treated at the same time
            # Drawback : instancing the 
            if(my_variables.frame_the_box) : 
                if((my_variables.frame_ism==True)&(my_variables.frame_jet==False)):
                   my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D,f_box_2D_inc=f_box_2D_inc,f_box_2D_exc=f_box_2D_exc,
                                                 f_box_2D_ism_inc=f_box_2D_ism_inc,f_box_2D_ism_exc=f_box_2D_ism_exc)                   
                if((my_variables.frame_ism==False)&(my_variables.frame_jet==True)):
                   my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D,f_box_2D_inc=f_box_2D_inc,f_box_2D_exc=f_box_2D_exc,
                                                 f_box_2D_jet_inc=f_box_2D_jet_inc,f_box_2D_jet_exc=f_box_2D_jet_exc)                   
                if((my_variables.frame_ism==True)&(my_variables.frame_jet==True)):
                   my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D,f_box_2D_inc=f_box_2D_inc,f_box_2D_exc=f_box_2D_exc,
                                                 f_box_2D_ism_inc=f_box_2D_ism_inc,f_box_2D_ism_exc=f_box_2D_ism_exc,
                                                 f_box_2D_jet_inc=f_box_2D_jet_inc,f_box_2D_jet_exc=f_box_2D_jet_exc)
                if((my_variables.frame_ism==False)&(my_variables.frame_jet==False)) :               
                    my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D,f_box_2D_inc=f_box_2D_inc,f_box_2D_exc=f_box_2D_exc)
            else : my_bazooka = load_the_bazooka(my_bazooka,my_variables,points2D)
            
            """DEBUG SCRIPT TO TEST IF 3D COMPUTATIONS HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_3D_computations.py").read())
            """END OF THE DEBUG SCRIPT"""
            
        
            """============================================================================="""
        
            os.write(1,b"IV) CC points projection on the sky plane\n")
            
            
            my_bazooka = full_shoot_with_the_bazooka(my_bazooka,my_variables)
            
            """DEBUG SCRIPT TO TEST IF CC PROJECTIONS ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_CC_projections_sky.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            """============================================================================="""
        
            os.write(1,b"V) Voxelling \n")
            os.write(1,b"Meshing the sky plane...\n")
            
            my_sky_grid = Sky_Grid()
             #don't forget to update the SKy_Grid attributes according to the variables current settings
            my_sky_grid = voxelling_sky_plane_meshing_bazooka(my_sky_grid,my_bazooka,my_variables)
            my_sky_grid.set_attributes(my_variables)
            
            """DEBUG SCRIPT TO TEST IF VOXELLING GRID MESHING ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_voxelling_mesh.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            os.write(1,b"Sampling projected data on the sky plane grid...\n")
            my_sky_grid = sample_sky_plane_grid_bazooka(my_variables,my_bazooka,my_sky_grid)
            
            """DEBUG SCRIPT TO TEST IF VOXELLING GRID SAMPLING ON THE SKY HAS BEEN DONE AS IN THE V1"""
            #exec(open("./test_voxelling_sample.py").read())
            """END OF THE DEBUG SCRIPT"""
            
            code_output_to_vtk(my_variables,my_sky_grid)
            #os.system('mv {0}.vts {1}/.'.format(my_variables.output_file_vtk,((my_variables.output_file_vtk).replace('.','')).replace('/','')))
            if os.path.exists('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.vts'.format(my_variables.output_file_vtk))):
                os.remove('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.vts'.format(my_variables.output_file_vtk)))
            shutil.move('{0}.vts'.format(my_variables.output_file_vtk),'{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')))
            
        elif(my_variables.voxelling_mode=='Flow_rates'):
            print("Doing step 1 over 1 with voxelling_method : \'{0}\'...".format(my_variables.voxelling_method))
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
            
            os.write(1,b"Done within 1 step...\n")
            os.write(1,b"III) Reference layering\n")
            
            
            my_variables.geo_rho_2D = read.extract(
                CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
            my_transverse_layer, points2D, r2 = layer_cells_surfaces(CellsData,my_variables,my_transverse_layer)
            
            os.write(1,b"IV) Momentum fluxes computing\n")
            
            my_transverse_layer.fluxes_surface_center,\
            my_transverse_layer.fluxes_surface_inclination,\
            my_transverse_layer.fluxes_surface_r_z,\
            my_transverse_layer.fluxes_surface_radius_z_inclination,\
            my_transverse_layer.fluxes = \
            compute_momentum_flux_through_circle(my_variables,points2D,r2,
                                                 my_transverse_layer,
                                                 center=my_variables.flux_surface_center,
                                                 r_z_inclination=my_variables.flux_surface_r_z_inclination)
            
            os.write(1,b"V) Output\n")
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
        
        os.write(1,b"II) Compute 2D-grid nodes positions and volumes\n")
        #calculer les positions 2D-cylindriques, constantes dans le domaine 2D
        
        my_variables,my_railgun,cells = compute_2D_grid_nodes_positions_volumes(CellsData,my_variables,my_railgun)
        
        #ne pas oublier d'extraire le champ de densite
        my_variables.geo_rho_2D = read.extract(
            CellsData, 'rho', attribute_mode=my_variables.geo_2D_attribute_mode)
        
        # construire les limites de la grille du plan du ciel
        # - abs(alpha) inferior or equal to Rmax
        # - Delta inforior or equal to sqrt(Rmax^2+Hmax^2) and superior or
        #   equal to 
        
        os.write(1,b"Done...\n")        
        os.write(1,b"III) Building the empty voxel mesh\n")

        #instance the sky_grid
        my_sky_grid = Sky_Grid()
    
        my_variables,my_railgun,my_sky_grid = building_empty_voxel_mesh(my_variables,my_railgun,my_sky_grid)

        os.write(1,b"Done...\n")

        os.write(1,b"IV) Voxelling\n")
    
        #instance the sky_grid
        my_mv_spectrum = MV_spectrum()
    
        my_variables,my_railgun,my_sky_grid,my_mv_spectrum = sample_sky_plane_grid_railgun(my_variables,my_railgun,my_sky_grid,my_mv_spectrum,cells)
        
   
        
        os.write(1,b"V) Output\n")
        code_output_to_vtk(my_variables,my_sky_grid) 
        if(my_variables.output_mv_spectrum==True): code_output_to_vtk_mv(my_variables,my_sky_grid,my_mv_spectrum) 
        #os.system('mv {0}.vts {1}/.'.format(my_variables.output_file_vtk,((my_variables.output_file_vtk).replace('.','')).replace('/','')))
        if os.path.exists('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.vts'.format(my_variables.output_file_vtk))):
            os.remove('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.vts'.format(my_variables.output_file_vtk)))
        shutil.move('{0}.vts'.format(my_variables.output_file_vtk),'{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')))
        if os.path.exists('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.csv'.format(my_variables.output_file_vtk_mv))):
            os.remove('{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')+'/'+'{0}.csv'.format(my_variables.output_file_vtk_mv)))
        shutil.move('{0}.csv'.format(my_variables.output_file_vtk_mv),'{0}'.format(((my_variables.output_file_vtk).replace('.','')).replace('/','')))
        os.write(1,b"Done... Exiting program...\n")
    else: raise NotImplementedError("No avalaible voxelling scheme selected!")
    
    
    """============================================================================="""
else: raise NotImplementedError("No avalaible voxelling method selected!")
