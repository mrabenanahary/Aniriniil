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
geo_angle_i = 4.0*np.pi/6.0#15.0*np.pi/18.0 #120°
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

input_file_vtk = "./input/Jet_cylindric_pulser_Lee2001/hllc/old/Jet_CC_Shang06_0610.vtu"
vtk_Output_Name = "./hllc/DE/Lee_2001_us_300yrs_unstratified/Cylindric_jet_nG237_nP60_296_ans"
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
x_frame_max = 6.66667e16#3.5e15#1.2e17
y_frame_min = 0.0
y_frame_max = 1.45e17

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


rhoISM = np.maximum(rhoISM,0)
rhoJet = np.maximum(rhoJet,0)


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

points2D[:,phi_] = points2D[:,phi_] + 0.5*geo_dangle_phi

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