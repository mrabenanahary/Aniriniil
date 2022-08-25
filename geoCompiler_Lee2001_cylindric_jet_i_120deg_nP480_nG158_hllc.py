#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

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


# variables globales :

#global raie_dict, H2O_raie_dict, OI_raie_dict, pH2O_exceptions


"====================PREAMBULE======================"
"""Paramètres Matplotlib"""
mpl.rcdefaults()

Latex = True
font = {'family': 'serif',
        'size': 17,
        'weight': 'extra bold'}
plt.rc('text', usetex=Latex)
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage[squaren,Gray]{SIunits}',
    r'\usepackage{amsmath}',
    r'\usepackage[version=4]{mhchem}']
plt.rcParams['axes.linewidth'] = 2


plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2

plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 2


py_title = \
    """0) PARAMETRES"""


print(py_title)

message = "Reading parameters..."
print(message)

# Conventions (certaines bizzares) des astronomes
not_my_data = set(dir())
geo_angle_i = np.pi/2.0#4.0*np.pi/6.0#15.0*np.pi/18.0 #120°
unit_angle_i = u.rad
geo_angle_i = (geo_angle_i * unit_angle_i).to_value(u.rad)
#print(geo_angle_i)

unit_masse = u.Msun
unit_volume = u.au**3

voxelling_mode = 'Voxel'  # Masse,Voxel
voxelling_mode_list = ['Masse', 'Voxel']
try:
    assert voxelling_mode in voxelling_mode_list
except AssertionError:
    s = ' '
    for el in voxelling_mode_list:
        s += str(el) + ' ||'
    raise AssertionError(
        'Voxelling mode \'{0}\' isn\'t in the list. Possible Modes = {1}'.format(voxelling_mode, s))

x_max = 1.2*8.166666666666666E16
y_max = 1.2*1.2576666666666666E17

epsilon_tol = 1e-6

# geoCompiler-friendly

input_file_vtk = "./Jet_CC_Shang06_0296.vtu"
vtk_Output_Name = "./Cylindric_jet_nG237_nP60_296_ans_Lee2001"
geo_2D_attribute_mode = 'cell'
if(geo_2D_attribute_mode != 'point') and (geo_2D_attribute_mode != 'cell'):
    raise ValueError('in file {0}'.format(str(os.path.basename(sys.argv[0]))) +
                     ' : attribute_mode=geo_2D_attribute_mode is neither a \" point \" or a \" cell \"')
geo_angle_i = geo_angle_i
nphi = 60#5120#512#32#720


if(nphi%2==0): 
    print("You choosed nphi={0} as an even integer.".format(nphi)) 
    print("Symmetry in the data cube after voxelling is obtained for i=+/- pi/2 !")
else:
    print("You choosed nphi={0} as an odd integer.".format(nphi))
    print("Symmetry in the data cube after voxelling may therefore not be obtained for i=+/- pi/2 !")
    print("Beware !\n")

# indexes for coordinates frames:
r_ = 0
z_cyl_ = 1
phi_ = 2

x_ = 0
y_ = 1
z_cart_ = 2

Gamma_ = 0
Delta_ = 1
alpha_ = 2

# voxeling parameters
# ===========================================
voxels_nGamma_cells = 237 #158  # in the biggest 2D-sky plane dimension
# ==> to avoid undersampling (which produces evenly
# space artefact lines where the projected more than one point overlap
# inside the same sampling cell along the dimension) :
#       >>>>>>  voxels_nGamma_cells  must be >= the cell numbers along the biggest
# dimension from the .vtu input

# =========================================
frame_the_box = True
frame_the_box2 = False

#include_frame_box = True
x_frame = r_
y_frame = z_cyl_
x_frame_min = 0.0#2.5e15#0.0
x_frame_max = 1.8e16#3.5e15#1.2e17
y_frame_min = 0.0
y_frame_max = 1.49e17

x_frame2 = r_
y_frame2 = z_cyl_
x_frame_min2 = 0.0#0.0
x_frame_max2 = 3.5e15#2.5e15#1.2e17
y_frame_min2 = 5.0e16
y_frame_max2 = 1.0e17#1.0e16


voxels_dz = 1.0#3.0  # km/s
unit_voxels_dz = u.km/u.s
voxels_z_max = 150.0 #km/s
voxels_z_min = -70.0 #km/s
unit_output_v = u.km/u.s
unit_voxels_z_max = u.km/u.s
unit_voxels_z_min = u.km/u.s
voxels_dz = (voxels_dz*unit_voxels_dz).to_value(u.cm/u.s)
voxels_z_max = (voxels_z_max*unit_voxels_z_max).to_value(u.cm/u.s)
voxels_z_min = (voxels_z_min*unit_voxels_z_min).to_value(u.cm/u.s)

fac_dx = 0.5  # default:0.5
fac_dy = 0.5  # default:0.5
fac_dz = 0.5  # default:0.5

# voxelling progress display for each percent ?
each_percent = 1.0  # in %
each_percent_frac = each_percent/100
number_of_voxels = voxels_nGamma_cells*voxels_nGamma_cells

my_data = set(dir()) - not_my_data
sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())


message = ''.ljust(10,"=")
message += "List of user-defined parameters :"
print(message)
for name in sorted_my_data :
    if(name!='not_my_data'):
        myvalue = eval(name)
        print(name, ' = ', myvalue,"\n")
message = ''.ljust(10,"=")+'End of user-defined parameters list'
print(message)

py_title = \
    """I) Read and get touch on .vtu file data and points in python """


print(py_title)

message = "Reading input..."
print(message)


reader = v.vtkXMLUnstructuredGridReader()

reader.SetFileName(input_file_vtk)
reader.Update()
output = reader.GetOutput()

#potential = vtk_to_numpy(output.GetPointData().GetArray("rho"))
# print(potential)

#nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
points3D_number = output.GetNumberOfCells()

#nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)

points2DData = reader.GetOutput()
#points2D = nodes_nummpy_array

ndim_problem = int(np.shape(points2DData.GetBounds())[0]/2)
number_of_vertixes_problem = 2**ndim_problem
number_of_vertixes_under_dimension = 2**(ndim_problem-1)

del output
message = "Done..."
print(message)



py_title =\
    """II) 2D-part computations"""



print(py_title)


Nstep = 2
step_code = 1
message = "Doing step {0} over {1}...".format(step_code,Nstep)
print(message)



#rr = read.extract(points2DData, 'r_cell', attribute_mode=geo_2D_attribute_mode)
#zz = read.extract(points2DData, 'z_cell', attribute_mode=geo_2D_attribute_mode)
#phiphi = np.full(np.shape(rr), 0)

vr = read.extract(points2DData, 'v1', attribute_mode=geo_2D_attribute_mode)
vz = read.extract(points2DData, 'v2', attribute_mode=geo_2D_attribute_mode)
vphi = read.extract(points2DData, 'v3', attribute_mode=geo_2D_attribute_mode)
rhoISM = read.extract(points2DData, 'tracer_ism', attribute_mode=geo_2D_attribute_mode)
rhoJet = read.extract(points2DData, 'tracer_jet', attribute_mode=geo_2D_attribute_mode)


rhoISM = np.maximum(rhoISM,0.0)
rhoJet = np.maximum(rhoJet,0.0)
#points2D = np.array([[rr[ii], zz[ii], phiphi[ii]] for ii in range(0, len(rr))])
vitesse2D = np.array([[vr[ii], vz[ii], vphi[ii]] for ii in range(0, len(vr))])

del vr, vz, vphi


number_elem_2D = points2DData.GetNumberOfCells()
"""
    cellsPointsId_c = np.full((number_elem_2D,4),0)
    #cellsPointsId_b = np.full((number_elem_2D,4),0)
    cellsPointsCoord_r_z_phi_2D_frame_c = np.empty((number_elem_2D,4,3))
    #cellsPointsCoord_r_z_phi_2D_frame_b = np.empty((number_elem_2D,4,3))
    cellsPoints_drr_dz_dphi_2D_frame_c = np.empty((number_elem_2D,3))
"""
cellsPoints_volume_2D_frame = np.empty((number_elem_2D,))
#cellsPoints_volume_2D_f_ISM = np.empty((number_elem_2D,))
cellsPointsId = np.full((number_elem_2D, 
                         number_of_vertixes_under_dimension), 0)
#cellsPointsId_b_b = np.full((number_elem_2D,4),0)
cellsPointsCoord_r_z_phi_2D_frame = np.empty((number_elem_2D, 
                                              number_of_vertixes_under_dimension, 
                                              ndim_problem))
#cellsPointsCoord_r_z_phi_2D_frame = np.empty((number_elem_2D,4,3))
cellsPoints_drr_dz_dphi_2D_frame = np.empty((number_elem_2D, ndim_problem))
#cellsPoints_volume_2D_frame_b = np.empty((number_elem_2D,))

# initiate iterated variables:
# geo_angle_phi = 0 #rad
geo_dangle_phi = 2*np.pi/nphi

"""#méthode 1 plus lente mais élémentaire

cellIds = v.vtkIdList()
for cellIndex in range(number_elem_2D):
    points2DData.GetCellPoints(cellIndex, cellIds) 
    
    #for i in range(0, cellIds.GetNumberOfIds()): # for every points of the given cell
    #    cellsPointsId[cellIndex,i]=cellIds.GetId(i)
    #    cellsPointsCoord_r_z_phi_2D_frame[cellIndex,i]=points2DData.GetPoint(cellsPointsId_c[cellIndex,i]) # get coordinates of the given point of the given cell, type: class 'tuple'
    
    cellsPointsId_c[cellIndex]=[cellIds.GetId(i) for i in range(0,cellIds.GetNumberOfIds())]
    cellsPointsCoord_r_z_phi_2D_frame_c[cellIndex]=[points2DData.GetPoint(cellsPointsId_c[cellIndex,i]) for i in range(0,cellIds.GetNumberOfIds())]# get coordinates of the given point of the given cell, type: class 'tuple'
    cellsPointsCoord_r_z_phi_2D_frame_c[cellIndex,:,phi_]=cellsPointsCoord_r_z_phi_2D_frame_c[cellIndex,:,phi_]-\
        0.5*geo_dangle_phi
    r_z_phi_max_c = np.amax(cellsPointsCoord_r_z_phi_2D_frame_c[cellIndex],axis=0)
    r_z_phi_min_c = np.amin(cellsPointsCoord_r_z_phi_2D_frame_c[cellIndex],axis=0)
    cellsPoints_drr_dz_dphi_2D_frame_c[cellIndex,r_]=np.abs(r_z_phi_max_c[r_]*r_z_phi_max_c[r_]-\
                                                    r_z_phi_min_c[r_]*r_z_phi_min_c[r_])
    cellsPoints_drr_dz_dphi_2D_frame_c[cellIndex,z_cyl_]=np.abs(r_z_phi_max_c[z_cyl_]-\
                                                    r_z_phi_min_c[z_cyl_])
cellsPoints_drr_dz_dphi_2D_frame_c[:,phi_]=geo_dangle_phi
"""
step_code = 2
message = "Doing step {0} over {1}...".format(step_code,Nstep)
print(message)

# méthode 2 plus compliquée mais plus rapide
#cellIndex_array = [i for i in range(number_elem_2D)]

temp_CellsPoint = np.empty((number_elem_2D,ndim_problem*2,))
{points2DData.GetCellBounds(i,temp_CellsPoint[i]) for i in range(0,number_elem_2D)}
for idim in range(ndim_problem):
    cellsPointsCoord_r_z_phi_2D_frame[:,:,idim]=np.tile(temp_CellsPoint[:,2*idim:2*idim+2],ndim_problem-1)
r_offset = np.amin(np.tile(temp_CellsPoint[:,2*r_:2*r_+2],ndim_problem-1))
cellsPointsCoord_r_z_phi_2D_frame[:,:,phi_]=cellsPointsCoord_r_z_phi_2D_frame[:, :, phi_]-0.5*geo_dangle_phi
if(r_offset<0): cellsPointsCoord_r_z_phi_2D_frame[:,:,r_]=cellsPointsCoord_r_z_phi_2D_frame[:, :, r_]-r_offset

r_z_phi_max = np.array([list(np.amax(
    cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, number_elem_2D)])#+r_offset
r_z_phi_min = np.array([list(np.amin(
    cellsPointsCoord_r_z_phi_2D_frame[i], axis=0)) for i in range(0, number_elem_2D)])#+r_offset
# cellsPoints_drr_dz_dphi_2D_frame[cellIndex,r_]=np.abs(r_z_phi_max[r_]*r_z_phi_max[r_]-\
#                                                    _z_phi_min[r_]*r_z_phi_min[r_])
cellsPoints_drr_dz_dphi_2D_frame[:, r_] = np.abs(r_z_phi_max[:, r_]*r_z_phi_max[:, r_] -
                                                 r_z_phi_min[:, r_]*r_z_phi_min[:, r_])
cellsPoints_drr_dz_dphi_2D_frame[:, z_cyl_] = np.abs(r_z_phi_max[:, z_cyl_] -
                                                     r_z_phi_min[:, z_cyl_])
cellsPoints_drr_dz_dphi_2D_frame[:, phi_] = geo_dangle_phi

points2D = 0.5*(r_z_phi_max+r_z_phi_min)
if(frame_the_box):
    #if(include_frame_box):
    f_box_2D_inc = np.where(((points2D[:,x_frame]>=x_frame_min)&(points2D[:,x_frame]<=x_frame_max))&((points2D[:,y_frame]>=y_frame_min)&(points2D[:,y_frame]<=y_frame_max)),True,False)
    #else:
    f_box_2D_exc = ~f_box_2D_inc#np.where(((points2D[:,x_frame]>=x_frame_min)&(points2D[:,x_frame]<=x_frame_max))&((points2D[:,y_frame]>=y_frame_min)&(points2D[:,y_frame]<=y_frame_max)),False,True)
if(frame_the_box2):
    f_box_2D_inc = f_box_2D_inc | np.where(((points2D[:,x_frame2]>=x_frame_min2)&(points2D[:,x_frame2]<=x_frame_max2))&((points2D[:,y_frame2]>=y_frame_min2)&(points2D[:,y_frame2]<=y_frame_max2)),True,False)
    f_box_2D_exc = ~f_box_2D_inc#f_box_2D_exc | np.where(((points2D[:,x_frame2]>=x_frame_min)&(points2D[:,x_frame2]<=x_frame_max))&((points2D[:,y_frame2]>=y_frame_min)&(points2D[:,y_frame2]<=y_frame_max)),False,True)


f_box_2D_inc = np.where(f_box_2D_inc,1.0,0.0)
f_box_2D_exc = np.where(f_box_2D_exc,1.0,0.0)

points2D[:,phi_] = points2D[:,phi_] + geo_dangle_phi

temp_dangle_t = cellsPoints_drr_dz_dphi_2D_frame[:, phi_]
temp_drr = cellsPoints_drr_dz_dphi_2D_frame[:, r_]  # r2**2-r1**2
temp_dz = cellsPoints_drr_dz_dphi_2D_frame[:, z_cyl_]
cellsPoints_volume_2D_frame = 0.5 * temp_drr * temp_dangle_t * \
    temp_dz  # 0.5 * (r2**2-r1**2) * (phi2-phi1) * (z2 - z1)
del temp_dangle_t, temp_drr, temp_dz, temp_CellsPoint


"""#tests
MyMin = np.min
MyMax = np.max
a=cellsPoints_drr_dz_dphi_2D_frame[:,r_]
b=cellsPoints_drr_dz_dphi_2D_frame_c[:,r_]
print('min=','-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=cellsPoints_drr_dz_dphi_2D_frame[:,z_cyl_]
b=cellsPoints_drr_dz_dphi_2D_frame_c[:,z_cyl_]
print('min=','-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=cellsPoints_drr_dz_dphi_2D_frame[:,phi_]
b=cellsPoints_drr_dz_dphi_2D_frame_c[:,phi_]
print('min=','-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=cellsPointsCoord_r_z_phi_2D_frame[:,:,r_]
b=cellsPointsCoord_r_z_phi_2D_frame_c[:,:,r_]
print('min=','-->',MyMin(a-b,axis=1),'-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=1),'-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=cellsPointsCoord_r_z_phi_2D_frame[:,:,z_cyl_]
b=cellsPointsCoord_r_z_phi_2D_frame_c[:,:,z_cyl_]
print('min=','-->',MyMin(a-b,axis=1),'-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=1),'-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=cellsPointsCoord_r_z_phi_2D_frame[:,:,phi_]
b=cellsPointsCoord_r_z_phi_2D_frame_c[:,:,phi_]
print('min=','-->',MyMin(a-b,axis=1),'-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=','-->',MyMax(a-b,axis=1),'-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
"""
"""
print('>>> max<<<< = ', \
          np.max(cellsPointsCoord_r_z_phi_2D_frame_b-cellsPointsCoord_r_z_phi_2D_frame),\
          '>>> min<<<< = ', \
          np.min(cellsPointsCoord_r_z_phi_2D_frame_b-cellsPointsCoord_r_z_phi_2D_frame))"""


message = "Done with {0} steps...".format(Nstep)
print(message)

"""Reproduce the 2D plane in a 3D axisymmetric grid by rotation-filling
This place is supposed to be placed at phi=0, i.e. plane (O,x,z)"""


py_title = \
    """III) 3D-part computations"""


print(py_title)


number_elem_3D = number_elem_2D * nphi
# variables to contain each 3D variables
geo_position_x_y_z = np.empty((number_elem_3D, ndim_problem))
geo_position_r_z_phi = np.empty((number_elem_3D, ndim_problem))
geo_position_f_ism = np.empty((number_elem_3D,))
geo_position_f_jet = np.empty((number_elem_3D,))
if(frame_the_box):
    geo_position_f_box_inc = np.empty((number_elem_3D,))
    geo_position_f_box_exc = np.empty((number_elem_3D,))
    #geo_position_masse_box_inc = np.empty((number_elem_3D,))
    #geo_position_masse_box_exc = np.empty((number_elem_3D,))
    geo_position_cellsPoints_masse_box_inc = np.empty((number_elem_3D,))
    geo_position_cellsPoints_masse_box_exc = np.empty((number_elem_3D,))
geo_position_cellsPoints_volume = np.empty((number_elem_3D,))
geo_position_cellsPoints_masse = np.empty((number_elem_3D,))
geo_position_cellsPoints_masse_ism = np.empty((number_elem_3D,))
geo_position_cellsPoints_masse_jet = np.empty((number_elem_3D,))

geo_cells_r_z_phi_v = np.empty((number_elem_3D, ndim_problem))
geo_position_3D_cell_index_in_2D_frame = np.full((number_elem_3D,), 0)
geo_index_first_cell_in_2D_frame_iphi = np.full((nphi,), 0)
geo_rho = np.empty((number_elem_3D,))
# index to iterate over 3D-geo_* array's 1st dimension (= ID n° of the 3D-cell)
geo_3D_index = 0

# get 2D rho
geo_rho_2D = read.extract(
    points2DData, 'rho', attribute_mode=geo_2D_attribute_mode)

# iterate over each of the nphi phi-layers of the 3D domain
for n in range(0, nphi):  # phi \in [0,2\pi[

    # instance n
    geo_index_first_cell_in_2D_frame_iphi[n] = n*number_elem_2D
    # instance phi:
    geo_angle_phi = n*geo_dangle_phi

    geo_angle_phi_array = np.full((number_elem_2D), geo_angle_phi)
    del geo_angle_phi

    # arrays components for which we assign the same values in this loop level

    geo_position_3D_cell_index_in_2D_frame[geo_index_first_cell_in_2D_frame_iphi[n]:
                                           geo_index_first_cell_in_2D_frame_iphi[n] +
                                           number_elem_2D] = np.arange(0, number_elem_2D)

    geo_position_r_z_phi[geo_index_first_cell_in_2D_frame_iphi[n]:
                         geo_index_first_cell_in_2D_frame_iphi[n] +
                         number_elem_2D, r_] = points2D[:, r_]  # r in 3D
    geo_position_r_z_phi[geo_index_first_cell_in_2D_frame_iphi[n]:
                         geo_index_first_cell_in_2D_frame_iphi[n] +
                         number_elem_2D, z_cyl_] = points2D[:, z_cyl_]  # z in 3D
    geo_position_r_z_phi[geo_index_first_cell_in_2D_frame_iphi[n]:
                         geo_index_first_cell_in_2D_frame_iphi[n] +
                         number_elem_2D, phi_] = geo_angle_phi_array  # phi in 3D
    geo_position_cellsPoints_volume[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D] = cellsPoints_volume_2D_frame 
    geo_position_f_ism[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D] = rhoISM/(rhoISM+rhoJet)
    geo_position_f_jet[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D] = rhoJet/(rhoISM+rhoJet)
    if(frame_the_box):
        geo_position_f_box_inc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D] = f_box_2D_inc
        geo_position_f_box_exc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D] = f_box_2D_exc
    # radial velocity
    geo_cells_r_z_phi_v[geo_index_first_cell_in_2D_frame_iphi[n]:
                        geo_index_first_cell_in_2D_frame_iphi[n] +
                        number_elem_2D, r_] = vitesse2D[:, r_]  # r in 3D
    geo_cells_r_z_phi_v[geo_index_first_cell_in_2D_frame_iphi[n]:
                        geo_index_first_cell_in_2D_frame_iphi[n] +
                        number_elem_2D, z_cyl_] = vitesse2D[:, z_cyl_]  # z in 3D
    geo_cells_r_z_phi_v[geo_index_first_cell_in_2D_frame_iphi[n]:
                        geo_index_first_cell_in_2D_frame_iphi[n] +
                        number_elem_2D, phi_] = vitesse2D[:, phi_]  # z in 3D
    del geo_angle_phi_array
    # arrays components for which we use intermediate arrays
    geo_r_array = np.array(geo_position_r_z_phi[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                geo_index_first_cell_in_2D_frame_iphi[n] +
                                                number_elem_2D, r_])
    geo_phi_array = np.array(geo_position_r_z_phi[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                  geo_index_first_cell_in_2D_frame_iphi[n] +
                                                  number_elem_2D, phi_])
    # 3D cartesian coordinates
    geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:
                       geo_index_first_cell_in_2D_frame_iphi[n] +
                       number_elem_2D, x_] = list(geo_r_array*np.cos(geo_phi_array))  # x in 3D
    geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:
                       geo_index_first_cell_in_2D_frame_iphi[n] +
                       number_elem_2D, y_] = list(geo_r_array*np.sin(geo_phi_array))  # y in 3D
    geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:
                       geo_index_first_cell_in_2D_frame_iphi[n] +
                       number_elem_2D, z_cart_] = points2D[:, z_cyl_]  # z in 3D
    geo_rho[geo_index_first_cell_in_2D_frame_iphi[n]:
            geo_index_first_cell_in_2D_frame_iphi[n] +
            number_elem_2D] = geo_rho_2D
    geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:
                                   geo_index_first_cell_in_2D_frame_iphi[n] +
                                   number_elem_2D] = geo_rho[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                             geo_index_first_cell_in_2D_frame_iphi[n] +
                                                             number_elem_2D] * \
        geo_position_cellsPoints_volume[geo_index_first_cell_in_2D_frame_iphi[n]:
                                        geo_index_first_cell_in_2D_frame_iphi[n] +
                                        number_elem_2D] #raw mass
    geo_position_cellsPoints_masse_ism[geo_index_first_cell_in_2D_frame_iphi[n]:
                                   geo_index_first_cell_in_2D_frame_iphi[n] +
                                   number_elem_2D] = geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                             geo_index_first_cell_in_2D_frame_iphi[n] +
                                                             number_elem_2D] *\
    geo_position_f_ism[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D]
    geo_position_cellsPoints_masse_jet[geo_index_first_cell_in_2D_frame_iphi[n]:
                                   geo_index_first_cell_in_2D_frame_iphi[n] +
                                   number_elem_2D] = geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                             geo_index_first_cell_in_2D_frame_iphi[n] +
                                                             number_elem_2D] *\
    geo_position_f_jet[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D]
    if(frame_the_box):
        geo_position_cellsPoints_masse_box_inc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                   geo_index_first_cell_in_2D_frame_iphi[n] +
                                   number_elem_2D] = geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                             geo_index_first_cell_in_2D_frame_iphi[n] +
                                                             number_elem_2D] *\
    geo_position_f_box_inc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D]
        geo_position_cellsPoints_masse_box_exc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                   geo_index_first_cell_in_2D_frame_iphi[n] +
                                   number_elem_2D] = geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:
                                                             geo_index_first_cell_in_2D_frame_iphi[n] +
                                                             number_elem_2D] *\
    geo_position_f_box_exc[geo_index_first_cell_in_2D_frame_iphi[n]:
                                    geo_index_first_cell_in_2D_frame_iphi[n] +
                                    number_elem_2D]
                    
    message="Progression : {0:.2f} % --> point n° : {1} ".format(100*geo_3D_index/number_elem_3D,geo_3D_index)
    print(message)
    geo_3D_index += number_elem_2D
    
print(''.ljust(len(message),' '))
print("Progression : Finished.... 100 %")
print("number of 3D points = {0}".format(len(geo_rho)))


# convert, if needed, the mass from gram to unit_mass :
if(unit_masse != u.g):
    geo_position_cellsPoints_masse = (
        geo_position_cellsPoints_masse*u.g).to_value(unit_masse)
    geo_position_cellsPoints_masse_ism = (
        geo_position_cellsPoints_masse_ism*u.g).to_value(unit_masse)
    geo_position_cellsPoints_masse_jet = (
        geo_position_cellsPoints_masse_jet*u.g).to_value(unit_masse)
    if(frame_the_box):
        geo_position_cellsPoints_masse_box_inc = (
            geo_position_cellsPoints_masse_box_inc*u.g).to_value(unit_masse)
        geo_position_cellsPoints_masse_box_exc = (
            geo_position_cellsPoints_masse_box_exc*u.g).to_value(unit_masse)
# if(unit_masse!=u.g or unit_volume!=u.cm**3):
#    geo_rho=(geo_rho*u.g/(u.cm**3)).to_value(unit_masse/(u.au**3))


print("\nMinimal mass found : {0} {1}".format(
    np.str(np.amin(geo_position_cellsPoints_masse)), str(unit_masse)))
print("Maximal mass found : {0} {1}".format(
    np.str(np.amax(geo_position_cellsPoints_masse)), str(unit_masse)))

"""
points = v.vtkPoints()
quad = v.vtkQuad()
cells = v.vtkCellArray()

numPoints = ds.south_north.size*ds.west_east.size
   
print('Write points \n')
for i,j in product(ds.south_north.values,ds.west_east.values):
        points.InsertNextPoint(ds.lat.isel(south_north=i,west_east=j), ds.lon.isel(south_north=i,west_east=j), ds.HGT.sel(south_north=i,west_east=j).values/6370000.0)

print('Write cells \n')
for idx in range(points.GetNumberOfPoints()-ds.west_east.size):
    if (idx%ds.west_east.size != 0):
        quad.GetPointIds().SetId(0,idx)
        quad.GetPointIds().SetId(1,idx+1)
        quad.GetPointIds().SetId(2,idx+ds.west_east.size+1)
        quad.GetPointIds().SetId(3,idx+ds.west_east.size)
        cells.InsertNextCell(quad)

print('Create unstructured grid \n') 
grid = v.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(v.VTK_QUAD, cells)

writer = v.vtkXMLUnstructuredGridWriter()
writer.SetFileName('cosipy.vtu')
writer.SetInputData(grid)
writer.Write() 
"""

"""
geo_position_x_y_z_b = np.empty((number_elem_3D,3))
geo_position_r_z_phi_b = np.empty((number_elem_3D,3))
geo_position_3D_cell_index_in_2D_frame_b = np.full((number_elem_3D,),0)
geo_index_first_cell_in_2D_frame_iphi_b = np.full((nphi,),0)
geo_rho_b = np.empty((number_elem_3D,))
geo_3D_index_b = 0 #index to iterate over 3D-geo_* array's 1st dimension (= ID n° of the 3D-cell)

#initiate iterated variables:
#geo_angle_phi_b = 0 #rad
geo_dangle_phi_b = geo_dangle_phi

#get 2D rho
geo_rho_2D_b = read.extract(points2DData,'rho',attribute_mode=geo_2D_attribute_mode)

#iterate over each of the nphi phi-layers of the 3D domain
for n in range(0,nphi): #phi \in [0,2\pi[
    
    #instance n
    geo_index_first_cell_in_2D_frame_iphi_b[n]=n*number_elem_2D
    #instance phi:
    geo_angle_phi_b = n*geo_dangle_phi_b
    
    #iterate over each 2D phi-layer/plane's pixel
    for i,pixels in enumerate(points2D):
        #--> i-th pixel's (r,z) coordinates
        geo_position_3D_cell_index_in_2D_frame_b[geo_3D_index_b]=i
        
        geo_position_r_z_phi_b[geo_3D_index_b,r_],geo_r_b = pixels[r_],pixels[r_] #r in 3D
        geo_position_r_z_phi_b[geo_3D_index_b,z_cyl_],geo_z_b = pixels[z_cyl_],pixels[z_cyl_] #z in 3D
        geo_position_r_z_phi_b[geo_3D_index_b,phi_] = geo_angle_phi_b #phi in 3D
        #trivial case phi = 0 :
        if n==0:
            geo_phi_b = 0.0 #to forcefully ensure that i-th pixel's phi = 0
            
            #signal if something is wrong with the input .vtu file for which it must have phi=0
            try:
                assert np.abs(geo_phi_b-points2D[i,phi_])<epsilon_tol
            except AssertionError:
                raise AssertionError('Something fishy about the .vtu : the 2D plane\'s phi is != 0 : phi(pixel n° {0}) = {1}'.format(np.str(i),np.str(points2D[i,2])))    
            
            #--> i-th pixel's (x,y,z) co-frame coordinates (trivial case phi=0)
            geo_position_x_y_z_b[geo_3D_index_b,x_]=geo_r_b #x in 3D
            geo_position_x_y_z_b[geo_3D_index_b,y_]=0.0   #y in 3D
            geo_position_x_y_z_b[geo_3D_index_b,z_cart_]=geo_z_b #z in 3D
            
        #else, for any other non-trivial phi-plane
        else:
            #get the i-th pixel's phi coordinate
            geo_phi_b = geo_angle_phi_b 
            geo_position_x_y_z_b[geo_3D_index_b,x_]=geo_r_b*np.cos(geo_phi_b) #x in 3D
            geo_position_x_y_z_b[geo_3D_index_b,y_]=geo_r_b*np.sin(geo_phi_b)   #y in 3D
            geo_position_x_y_z_b[geo_3D_index_b,z_cart_]=geo_z_b #z in 3D
        
        
        #get rho
        #print(geo_position_3D_cell_index_in_2D_frame[geo_3D_index])
        geo_rho_b[geo_3D_index_b] = geo_rho_2D_b[geo_position_3D_cell_index_in_2D_frame_b[geo_3D_index_b]] 
        geo_3D_index_b+=1
    print(geo_3D_index_b)

#tests
MyMin = np.min
MyMax = np.max
a=geo_rho_b
b=geo_rho
print('min=',MyMin(a-b),' et ','max=',MyMax(a-b))        
a=geo_position_x_y_z_b
b=geo_position_x_y_z
print('min=',MyMin(a-b,axis=1),'-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=',MyMax(a-b,axis=1),'-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=geo_position_r_z_phi_b
b=geo_position_r_z_phi
print('min=',MyMin(a-b,axis=1),'-->',MyMin(a-b,axis=0),'-->',MyMin(a-b),' et ','max=',MyMax(a-b,axis=1),'-->',MyMax(a-b,axis=0),'-->',MyMax(a-b))
a=geo_index_first_cell_in_2D_frame_iphi_b
b=geo_index_first_cell_in_2D_frame_iphi
print('min=',MyMin(a-b),' et ','max=',MyMax(a-b))
a=geo_position_3D_cell_index_in_2D_frame_b
b=geo_position_3D_cell_index_in_2D_frame
print('min=',MyMin(a-b),' et ','max=',MyMax(a-b))
"""

#del points2DData



py_title = \
    """IV) CC points projection on the sky plane"""


print(py_title)

geo_position_Gamma_Delta_alpha = np.empty((number_elem_3D, 3))

"""
#méthode 1 plus lente mais élémentaire
for i,cellCoords in enumerate(geo_position_r_z_phi):
    geo_position_Gamma_Delta_alpha[i,Gamma_] = cellCoords[r_]*np.sin(geo_angle_i)*\
        np.cos(cellCoords[phi_])+cellCoords[z_cyl_]*np.cos(geo_angle_i) #r in 3D
    geo_position_Gamma_Delta_alpha[i,Delta_] = -cellCoords[r_]*np.cos(geo_angle_i)*\
        np.cos(cellCoords[phi_])+cellCoords[z_cyl_]*np.sin(geo_angle_i) # z in 3D
    geo_position_Gamma_Delta_alpha[i,alpha_] = cellCoords[r_]*np.sin(cellCoords[phi_])#phi in 3D
"""

# méthode 2 plus compliquée mais plus rapide
geo_position_Gamma_Delta_alpha[:, Gamma_] = geo_position_r_z_phi[:, r_]*np.sin(geo_angle_i) *\
    np.cos(geo_position_r_z_phi[:, phi_]) +\
    geo_position_r_z_phi[:, z_cyl_]*np.cos(geo_angle_i)  # r in 3D
geo_position_Gamma_Delta_alpha[:, Delta_] = -geo_position_r_z_phi[:, r_]*np.cos(geo_angle_i) *\
    np.cos(geo_position_r_z_phi[:, phi_]) + \
    geo_position_r_z_phi[:, z_cyl_]*np.sin(geo_angle_i)  # z in 3D
geo_position_Gamma_Delta_alpha[:, alpha_] = geo_position_r_z_phi[:,
                                                                 r_]*np.sin(geo_position_r_z_phi[:, phi_])  # phi in 3D
if(voxelling_mode == 'Voxel'):
    geo_vobs = np.empty((number_elem_3D,))
    geo_vobs = geo_cells_r_z_phi_v[:, phi_] *\
        np.sin(geo_angle_i)*np.sin(geo_position_r_z_phi[:, phi_]) -\
        -geo_cells_r_z_phi_v[:, r_]*np.sin(geo_angle_i)*np.cos(geo_position_r_z_phi[:, phi_]) -\
        geo_cells_r_z_phi_v[:, z_cyl_]*np.cos(geo_angle_i)


py_title = \
    """V) Voxelling """


print(py_title)

# maximum bounds
v_x_max = np.max(geo_position_Gamma_Delta_alpha[:, alpha_])
v_y_max = np.max(geo_position_Gamma_Delta_alpha[:, Delta_])

v_x_min = np.min(geo_position_Gamma_Delta_alpha[:, alpha_])
v_y_min = np.min(geo_position_Gamma_Delta_alpha[:, Delta_])

# Dimensions
#nx, ny, nz = voxels_nGamma_cells, voxels_nGamma_cells, 1
if(v_x_max > v_y_max):
    nx = int(voxels_nGamma_cells)
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
    ny = int(voxels_nGamma_cells)
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
if(voxelling_mode == 'Masse'):
    nz = 1
    lz = 1
    dz = lz/nz
    Z = np.arange(0, lz+0.1*dz, dz, dtype='float64')
elif(voxelling_mode == 'Voxel'):
    v_max = voxels_z_max #max(np.abs(np.amin(geo_vobs)), np.abs(np.amax(geo_vobs)))
    v_min = voxels_z_min
    nz = int(ceil(v_max/voxels_dz))
    #v_max = nz*voxels_dz
    #dz = voxels_dz
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
x = np.zeros((nx + 1, ny + 1, nz + 1))
y = np.zeros((nx + 1, ny + 1, nz + 1))
z = np.zeros((nx + 1, ny + 1, nz + 1))
# We add some random fluctuation to make the grid more interesting
for k in range(nz+1):
    for j in range(ny+1):
        for i in range(nx+1):
            if(i==nx): x[i, j, k] = X[nx-1] + 0.5 * dx
            else: x[i, j, k] = X[i] - 0.5 * dx
            if(j==ny): y[i, j, k] = Y[ny-1] + 0.5 * dy
            else: y[i, j, k] = Y[j] - 0.5 * dy
            if(k==nz): z[i, j, k] = Z[nz-1] + 0.5 * dz
            else: z[i, j, k] = Z[k] - 0.5 * dz



# Variables
"""
masse_arrays = np.zeros((nx,ny,nz))
test_k = 0
for k in range(0,nx):
    a=list(np.abs(geo_position_Gamma_Delta_alpha[:,alpha_]-X[k])<0.5*dx)
    for j in range(0,ny):
        b=list(np.abs(geo_position_Gamma_Delta_alpha[:,Delta_]-Y[j])<0.5*dx)
        c = [a[i] and b[i] for i in range(number_elem_3D)]
        masse_arrays[k,j,0] = np.sum(geo_position_cellsPoints_masse[c])
        if((k*ny+j)//int(each_percent_frac*nx*ny)>test_k): 
            #print("k*ny+j= ",k*ny+j)
            print("Progression (en %): {0}".format(str(100*(k*ny+j)/(nx*ny))))
            test_k+=1
        del a,b,c
 
array2D_vtk = masse_arrays.reshape( (nx, ny, nz)) 
#temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1)) 
gridToVTK("./Masse_sky_plane", x, y, z, cellData = {"array2D_vtk" : array2D_vtk})
"""


masse_arrays = np.zeros((nx, ny, nz), dtype=float)
masse_arrays_ism = np.zeros((nx, ny, nz), dtype=float)
masse_arrays_jet = np.zeros((nx, ny, nz), dtype=float)
if(frame_the_box): 
    masse_arrays_box_inc = np.zeros((nx, ny, nz), dtype=float)
    masse_arrays_box_exc = np.zeros((nx, ny, nz), dtype=float)
# np.copy(geo_rho)#np.copy(geo_position_cellsPoints_masse)#
#variable_to_voxel = np.copy(geo_position_cellsPoints_masse)
temp_x_y = np.copy(geo_position_Gamma_Delta_alpha)
temp_masse = np.copy(geo_position_cellsPoints_masse)
temp_masse_ism = np.copy(geo_position_cellsPoints_masse_ism)
temp_masse_jet = np.copy(geo_position_cellsPoints_masse_jet)
if(frame_the_box): 
    temp_masse_box_inc = np.copy(geo_position_cellsPoints_masse_box_inc)
    temp_masse_box_exc = np.copy(geo_position_cellsPoints_masse_box_exc)
if(voxelling_mode == 'Voxel'):
    temp_vobs = np.copy(geo_vobs)
test_k = 0
#del variable_to_voxel

"""
for k in range(nx):
    a=np.abs(geo_position_Gamma_Delta_alpha[:,alpha_]-X[k])<0.5*dx
    for j in range(0,ny):
        b=np.abs(geo_position_Gamma_Delta_alpha[:,Delta_]-Y[j])<0.5*dx
        c = a & b
        masse_arrays[k,j,0] = np.sum(geo_position_cellsPoints_masse[c])
        if((k*ny+j)//int(each_percent_frac*nx*ny)>test_k): 
            #print("k*ny+j= ",k*ny+j)
            print("Progression (en %): {0}".format(str(100*(k*ny+j)/(nx*ny))))
            test_k+=1
        #del b,c
    #del a
array2D_vtk = masse_arrays.reshape( (nx, ny, nz)) 
#temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1)) 
gridToVTK("./Masse_sky_plane", x, y, z, cellData = {"array2D_vtk" : array2D_vtk})
"""

for k in range(nx):
    a = np.abs(temp_x_y[:, alpha_]-X[k]) <= fac_dx*dx
    if(True in a):
        anti_a = ~a

        b = temp_x_y[a]
        masse_b = temp_masse[a]
        masse_ism_b = temp_masse_ism[a]
        masse_jet_b = temp_masse_jet[a]
        if(frame_the_box): 
            masse_box_inc_b = temp_masse_box_inc[a]
            masse_box_exc_b = temp_masse_box_exc[a]

        temp_x_y = temp_x_y[anti_a]
        temp_masse = temp_masse[anti_a]
        temp_masse_ism = temp_masse_ism[anti_a]
        temp_masse_jet = temp_masse_jet[anti_a]
        if(frame_the_box): 
            temp_masse_box_inc = temp_masse_box_inc[anti_a]
            temp_masse_box_exc = temp_masse_box_exc[anti_a]

        if(voxelling_mode == 'Voxel'):
            b_vobs = temp_vobs[a]
            temp_vobs = temp_vobs[anti_a]
        for j in range(ny):
            c = np.abs(b[:, Delta_]-Y[j]) <= fac_dy*dy
            if(True in c):
                anti_c = ~c
                # for i in range(len(masse_b[c])):
                #    masse_arrays[k,j,0] = masse_arrays[k,j,0] + (masse_b[c])[i]
                if(voxelling_mode == 'Masse'):
                    masse_arrays[k, j, 0] = np.sum(masse_b[c])
                    masse_arrays_ism[k, j, 0] = np.sum(masse_ism_b[c])
                    masse_arrays_jet[k, j, 0] = np.sum(masse_jet_b[c])
                    if(frame_the_box): 
                        masse_arrays_box_inc[k, j, 0] = np.sum(masse_box_inc_b[c])
                        masse_arrays_box_exc[k, j, 0] = np.sum(masse_box_exc_b[c])
                elif(voxelling_mode == 'Voxel'):
                    masse_c = masse_b[c]
                    masse_ism_c = masse_ism_b[c]
                    masse_jet_c = masse_jet_b[c]
                    if(frame_the_box):  
                        masse_box_inc_c = masse_box_inc_b[c]
                        masse_box_exc_c = masse_box_exc_b[c]
                    c_vobs = b_vobs[c]
                    for h in range(nz):
                        d = np.abs(c_vobs-Z[h]) <= fac_dz*dz
                        if(True in d):
                            anti_d = ~d
                            masse_arrays[k, j, h] = np.sum(masse_c[d])
                            masse_arrays_ism[k, j, h] = np.sum(masse_ism_c[d])
                            masse_arrays_jet[k, j, h] = np.sum(masse_jet_c[d])
                            if(frame_the_box): 
                                masse_arrays_box_inc[k, j, h] = np.sum(masse_box_inc_c[d])
                                masse_arrays_box_exc[k, j, h] = np.sum(masse_box_exc_c[d])
                            c_vobs = c_vobs[anti_d]
                            masse_c = masse_c[anti_d]
                            masse_ism_c = masse_ism_c[anti_d]
                            masse_jet_c = masse_jet_c[anti_d]
                            if(frame_the_box): 
                                masse_box_inc_c = masse_box_inc_c[anti_d]
                                masse_box_exc_c = masse_box_exc_c[anti_d]
                            del d, anti_d
                        else:
                            del d
                        #if(voxelling_mode == 'Voxel'):
                        if(((k*ny+j)*nz+h)/(nx*ny*nz)-test_k > each_percent_frac):
                            #print("k*ny+j= ",k*ny+j)
                            print("Progression (en %): {0:.2f} ".format(
                                100*((k*ny+j)*nz+h)/(nx*ny*nz)))
                            test_k = ((k*ny+j)*nz+h)/(nx*ny*nz)
                    b_vobs = b_vobs[anti_c]
                    del masse_c, c_vobs, masse_ism_c, masse_jet_c
                    if(frame_the_box): 
                        del masse_box_inc_c,masse_box_exc_c
                b = b[anti_c]
                masse_b = masse_b[anti_c]
                masse_ism_b = masse_ism_b[anti_c]
                masse_jet_b = masse_jet_b[anti_c]
                if(frame_the_box): 
                    masse_box_inc_b = masse_box_inc_b[anti_c]
                    masse_box_exc_b = masse_box_exc_b[anti_c]
                # if(len(temp_x_y)==0): break
                del c, anti_c
            else:
                del c
            if((k*ny+j)/(nx*ny)-test_k > each_percent_frac):
                #print("k*ny+j= ",k*ny+j)
                print("Progression (en %): {0:.2f} ".format(
                    100*(k*ny+j)/(nx*ny)))
                #sys.stdout.write('\033[2K\033[1G')
                test_k = (k*ny+j)/(nx*ny)
                #if(voxelling_mode == 'Masse'):
        del b, masse_b, masse_ism_b, masse_jet_b
        if(frame_the_box): 
            del masse_box_inc_b,masse_box_exc_b
    else:
        del a
    
    if(k/nx-test_k > each_percent_frac):
        #print("k*ny+j= ",k*ny+j)
        print("Progression (en %): {0:.2f} ".format(
            100*k/n))
        #sys.stdout.write('\033[2K\033[1G')
        test_k = k/nx
print("Progression (in %): Finishing... 100 % ")
   
print("Remaining cells at the end of voxelling : ", np.shape(temp_masse))
del temp_masse, temp_masse_ism, temp_masse_jet
if(frame_the_box): del temp_masse_box_inc, temp_masse_box_exc


if(voxelling_mode == 'Voxel'):
    del temp_vobs 

        # if(len(temp_x_y)==0): break
array2D_vtk = masse_arrays.reshape((nx, ny, nz))
array2D_vtk_ism = masse_arrays_ism.reshape((nx, ny, nz))
array2D_vtk_jet = masse_arrays_jet.reshape((nx, ny, nz))
if(frame_the_box):
    array2D_vtk_box_inc = masse_arrays_box_inc.reshape((nx, ny, nz))
    array2D_vtk_box_exc = masse_arrays_box_exc.reshape((nx, ny, nz))
#temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1))
if(not frame_the_box):
    gridToVTK(vtk_Output_Name, x, y, (z *
              u.cm/u.s).to_value(unit_output_v), cellData={"masse totale": array2D_vtk,
                                                           "masse ism": array2D_vtk_ism,
                                                           "masse jet": array2D_vtk_jet})
else:
    #if(include_frame_box):
    name_box_inc = "masse box (included)"
    #else:
    name_box_exc = "masse box (excluded)"
    gridToVTK(vtk_Output_Name, x, y, (z *
              u.cm/u.s).to_value(unit_output_v), cellData={"masse totale": array2D_vtk,
                                                           "masse ism": array2D_vtk_ism,
                                                           "masse jet": array2D_vtk_jet,
                                                           name_box_inc: array2D_vtk_box_inc,
                                                           name_box_exc: array2D_vtk_box_exc})
# test plot :

#py_title = \
#    """VI) Plotting """


#print(py_title)

"""for n in range(0, 1):
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111)
plt.tight_layout()

# variables to plot
# n=0
"""
# arrays components for which we assign the same values in this loop level

"""    
x = geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D,x_] #x in 3D


"""
"""
y = geo_position_Gamma_Delta_alpha[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D,Delta_] #alpha in 3D
x = -geo_position_Gamma_Delta_alpha[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D,alpha_] #alpha in 3D
"""

# y = geo_position_Gamma_Delta_alpha[:,Delta_] #alpha in 3D
# x = geo_position_Gamma_Delta_alpha[:,alpha_] #alpha in 3D

"""
y = geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D,y_] #y in 3D
"""
"""
y = geo_position_x_y_z[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D,z_cart_]
"""
"""
array2D = geo_rho[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D]
"""

#array2D = geo_rho

#array2D = (geo_position_cellsPoints_masse*u.g).to_value(u.Msun)

"""

array2D = geo_position_cellsPoints_volume[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D]

"""
"""
array2D = (geo_position_cellsPoints_masse[geo_index_first_cell_in_2D_frame_iphi[n]:\
                     geo_index_first_cell_in_2D_frame_iphi[n]+\
                         number_elem_2D]*u.g).to_value(u.Msun) #in Msun
"""
"""
x=points2D[:,0]
y=points2D[:,1]
"""

"""
ax.set_xlim(-x_max,x_max)
ax.set_ylim(-y_max,y_max)

ax.grid(True)
#ax.imshow(ints2D,array2D)
ax.set_title(r'2D domain (cell-centered pixels) n={0} over {1}'.format(np.str(n),np.str(nphi-1)))
ax.set_xlabel(r'$\Gamma$ (cm)')
ax.set_ylabel(r'$\alpha$ (cm)')
im = ax.scatter(x,y,s=0.2,c=array2D, marker='s',cmap=cm.jet,norm=mpl.colors.LogNorm())
#ax.set_aspect('equal', 'box')

#fig.colorbar(im, ax=ax, label=r'$\rho~(g.cm^{-3})$')
fig.colorbar(im, ax=ax, label=r'$\rho~(g.cm^{-3})$')
print("Printing the JPEG :")
plt.savefig('test_plot_vtu_n={0} over {1}.jpg'.format(np.str(n),np.str(nphi-1)))
print("Printing the PDF :")
plt.savefig('test_plot_vtu_n={0} over {1}.pdf'.format(np.str(n),np.str(nphi-1)))
#plt.savefig('test_plot_vtu_s=4.jpg')
print("Plotting on Matplotlib :")
plt.show()
plt.close('all')
print('n=',n)
"""

del message