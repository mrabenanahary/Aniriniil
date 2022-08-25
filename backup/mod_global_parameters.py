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
= Module contenant tous les noms de variables globales utilisées dans les scripts
= qui font appel à ce module
= 
=
=================================================================
"""



from cycler import cycler
import numpy as np
import vtk as v
import csv
import os
from math import ceil
from astropy import units as u
import time
import struct
from mod_exceptions import *


class Variable:
    #below are the list of global parameters AND their default value:
    not_my_data = set(dir()) #add your global variables after this line
    
    """==================================================================="""
    """Parameters default values"""
    
    dummy_test = {}
    geo_angle_i = 4.0*np.pi/6.0
    test_angle_i_value(geo_angle_i,varname='geo_angle_i',vmin=-np.pi,vmax=np.pi)
    unit_angle_i = u.rad
    geo_angle_i = (geo_angle_i * unit_angle_i).to_value(u.rad)
    unit_masse = u.Msun
    unit_volume = u.au**3

    voxelling_mode = 'Voxel'  # Masse,Voxel
    voxelling_mode_list = ['Voxel','Flow_rates']
    
    voxelling_method = 'bazooka'
    voxelling_method_list = ['bazooka','railgun']
    voxelling_scheme = 'brute'
    voxelling_scheme_list = {'bazooka':['brute'],'railgun':['progressive']}
    
    
    test_voxelling_mode(voxelling_mode,voxelling_mode_list,voxelling_parameter_name="voxelling_mode")
    test_voxelling_mode(voxelling_method,voxelling_method_list,voxelling_parameter_name="voxelling_method")
    test_voxelling_mode(voxelling_scheme,voxelling_scheme_list[voxelling_method],voxelling_parameter_name="voxelling_scheme")
    
    x_max = 1.2*8.166666666666666E16
    y_max = 1.2*1.2576666666666666E17
    epsilon_tol = 1e-6
    
    input_par = "./test_input.par"
    input_file_vtk = "./input/Jet_cylindric_pulser_Lee2001/hllc/old/Jet_CC_Shang06_0610.vtu"#"./Jet_Shang06_B4_0090_LEVEL5.vtu"#"./input/Jet_cylindric_pulser_Lee2001/hllc/old/Jet_CC_Shang06_0610.vtu"
    output_file_vtk = "./hllc/DE/Lee_2001_us_300yrs_unstratified/Cylindric_jet_nG237_nP60_296_ans"#"./Cylindric_jet_nG237_nP60_10000_ans_i120deg"#"./hllc/DE/Lee_2001_us_300yrs_unstratified/Cylindric_jet_nG237_nP60_296_ans"

    #impulsion mode
    output_impulsion_ascii = "./data.dat"
    unit_momentum_flux_output = u.Msun*u.km/(u.s*u.yr)
    unit_mass_flux_output = u.Msun/u.yr
    flux_surface_center=[[0,2.6875e16,0]]
    flux_surface_r_z_inclination=[[1.4e16,0.0]]    
    
    geo_2D_attribute_mode = 'cell'
    possible_geo_2D_attribute_mode = ['cell']
    
    test_attribute_mode(input_file_vtk,geo_2D_attribute_mode)
    
    nphi = 60#5120#512#32#720
    
    if(nphi%2==0): 
        nphi_choice_msg = "You choosed nphi={0} as an even integer.".format(nphi)
        nphi_choice_msg = nphi_choice_msg + "\n" + "Symmetry in the data cube after voxelling is obtained for i=+/- pi/2 !"
    else:
        nphi_choice_msg = "You choosed nphi={0} as an odd integer.".format(nphi)
        nphi_choice_msg = nphi_choice_msg + "\n" +"Symmetry in the data cube after voxelling may therefore not be obtained for i=+/- pi/2 !"
        nphi_choice_msg = nphi_choice_msg + "\n" +"Beware !\n"
    
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
    
    #include_frame_box = True
    x_frame = r_
    y_frame = z_cyl_
    x_frame_min = 0.0#2.5e15#0.0
    x_frame_max = 6.66667e16#3.5e15#1.2e17
    y_frame_min = 0.0
    y_frame_max = 1.45e17
    
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
    
    floor_ism = True
    threshold_f_ism = 1.0e-1
    
    floor_jet = True
    threshold_f_jet = 1.0e-1
    
    fac_dx = 0.5  # default:0.5
    fac_dy = 0.5  # default:0.5
    fac_dz = 0.5  # default:0.5
    
    # voxelling progress display for each percent ?
    each_percent = 1.0  # in %
    each_percent_frac = each_percent/100
    number_of_voxels = voxels_nGamma_cells*voxels_nGamma_cells
    
    #vtk computation related variables
    velocityVector = {'v1':0,'v2':0,'v3':0}
    rhoTracers = {'ISM':[0,0],'jet':[0,0]}
    
    #for 3D treatment global variables
    geo_rho_2D = None

    """End of parameters default values"""
    """==================================================================="""
    
    my_data = set(dir()) - not_my_data
    sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())
    
    #avoid putting list of variable between not_my_data and my_data, since
    #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
    #only for the Variable class attribute
    
    list_of_variables = []
    list_of_variables_print = ''.ljust(10,"=")
    list_of_variables_print += "List of user-defined parameters :\n"
    #print(list_of_variables_print)
    for name in sorted_my_data :
        if(name!='not_my_data'):
            myvalue = eval(name)
            list_of_variables.append(name)
            list_of_variables_print = list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
    list_of_variables_print = list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    #print(list_of_variables_print)
    
    def __init__(self):
        for varname in Variable.list_of_variables:
            setattr(self,varname,getattr(Variable,varname))
            
    def __repr__(self):
        o='{'
        for i,varname in enumerate(Variable.list_of_variables):
            o=o+varname+' : '+str(getattr(self,varname))
            if(i!=np.size(Variable.list_of_variables)-1): o+=' , '
        o+='}'
        return o
    
    def __str__(self):
        o=r""
        for i,varname in enumerate(Variable.list_of_variables):
            o=o+varname+r" : "+str(getattr(self,varname))
            if(i!=np.size(Variable.list_of_variables)-1): o+="\n"
        return o        
    
class Bazooka:
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(Variable.input_file_vtk)
    reader.Update()
    output = reader.GetOutput()
    
    
    #below are the list of global parameters AND their default value:
    not_my_data = set(dir()) #add your global variables after this line
    
    
    """==================================================================="""
    """Parameters default values"""
    
    numberOfCells = output.GetNumberOfCells()
    ndimProblem = int(np.shape(output.GetBounds())[0]/2)
    numbeOfVertixes = 2**ndimProblem
    numberOfVertixesUnderDimension = 2**(ndimProblem-1)
    
    cellsPoints_volume_2D_frame = np.empty((numberOfCells,))
    cellsPointsId = np.full((numberOfCells, 
                             numberOfVertixesUnderDimension), 0)
    cellsPointsCoord_r_z_phi_2D_frame = np.empty((numberOfCells, 
                                                  numberOfVertixesUnderDimension, 
                                                  ndimProblem))
    cellsPoints_drr_dz_dphi_2D_frame = np.empty((numberOfCells, ndimProblem))
    geo_dangle_phi = 2*np.pi/Variable.nphi
    
    
    #3D-variables needed for bazooka method
    numberOf3Delem = numberOfCells * Variable.nphi
    # variables to contain each 3D variables
    geo_index_first_cell_in_2D_frame_iphi = 0
    geo_position_x_y_z = np.empty((numberOf3Delem, ndimProblem))
    geo_position_r_z_phi = np.empty((numberOf3Delem, ndimProblem))
    geo_position_f_ism = np.empty((numberOf3Delem,))
    geo_position_f_jet = np.empty((numberOf3Delem,))

    geo_position_cellsPoints_volume = np.empty((numberOf3Delem,))
    geo_position_cellsPoints_masse = np.empty((numberOf3Delem,))
    geo_position_cellsPoints_masse_ism = np.empty((numberOf3Delem,))
    geo_position_cellsPoints_masse_jet = np.empty((numberOf3Delem,))

    if(Variable.floor_ism):
        geo_position_f_ism_inc = np.empty((numberOf3Delem,))
        geo_position_f_ism_exc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_ism_inc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_ism_exc = np.empty((numberOf3Delem,))
    if(Variable.floor_jet):
        geo_position_f_jet_inc = np.empty((numberOf3Delem,))
        geo_position_f_jet_exc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_jet_inc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_jet_exc = np.empty((numberOf3Delem,))
    if(Variable.frame_the_box):
        geo_position_f_box_inc = np.empty((numberOf3Delem,))
        geo_position_f_box_exc = np.empty((numberOf3Delem,))
        #geo_position_masse_box_inc = np.empty((numberOf3Delem,))
        #geo_position_masse_box_exc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_box_inc = np.empty((numberOf3Delem,))
        geo_position_cellsPoints_masse_box_exc = np.empty((numberOf3Delem,))        
    
    geo_cells_r_z_phi_v = np.empty((numberOf3Delem, ndimProblem))
    geo_position_3D_cell_index_in_2D_frame = np.full((numberOf3Delem,), 0)
    geo_index_first_cell_in_2D_frame_iphi = np.full((Variable.nphi,), 0)
    geo_rho = np.empty((numberOf3Delem,))
    # index to iterate over 3D-geo_* array's 1st dimension (= ID n° of the 3D-cell)
    geo_3D_index = 0
    
    #2D variables for the voxelling
    geo_position_Gamma_Delta_alpha = np.empty((numberOf3Delem, 3))
    if(Variable.voxelling_mode == 'Voxel'):
        geo_vobs = np.empty((numberOf3Delem,))
        
    """End of parameters default values"""
    """==================================================================="""
    
    my_data = set(dir()) - not_my_data
    sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())
    

    #avoid putting list of variable between not_my_data and my_data, since
    #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
    #only for the Bazooka class attribute
    
    
    list_of_variables = []
    list_of_variables_print = ''.ljust(10,"=")
    list_of_variables_print += "List of user-defined parameters :\n"
    #print(list_of_variables_print)
    
    for name in sorted_my_data :
        if(name!='not_my_data'):
            myvalue = eval(name)
            list_of_variables.append(name)
            list_of_variables_print = list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
    list_of_variables_print = list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    #print(list_of_variables_print)
    
    
    def __init__(self):
        for varname in Bazooka.list_of_variables:
            setattr(self,varname,getattr(Bazooka,varname))
            
    def __repr__(self):
        o='{'
        for i,varname in enumerate(Bazooka.list_of_variables):
            o=o+varname+' : '+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+=' , '
        o+='}'
        return o
    
    def __str__(self):
        o=r""
        for i,varname in enumerate(Bazooka.list_of_variables):
            o=o+varname+r" : "+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+="\n"
        return o
    
    def set_attributes(self,output,my_variables):   
        #below are the list of the parameters AND their updated value:
        self.not_my_data = set(dir()) #add your global updated variables after this line
        
        """==================================================================="""
        """Parameters updated values"""
        
        self.numberOfCells = output.GetNumberOfCells()
        self.ndimProblem = int(np.shape(output.GetBounds())[0]/2)
        self.numbeOfVertixes = 2**self.ndimProblem
        self.numberOfVertixesUnderDimension = 2**(self.ndimProblem-1)
        
        self.cellsPoints_volume_2D_frame = np.empty((self.numberOfCells,))
        self.cellsPointsId = np.full((self.numberOfCells, 
                                 self.numberOfVertixesUnderDimension), 0)
        self.cellsPointsCoord_r_z_phi_2D_frame = np.empty((self.numberOfCells, 
                                                      self.numberOfVertixesUnderDimension, 
                                                      self.ndimProblem))
        self.cellsPoints_drr_dz_dphi_2D_frame = np.empty((self.numberOfCells, self.ndimProblem))
        self.geo_dangle_phi = 2*np.pi/my_variables.nphi
        
        
        #3D-variables needed for bazooka method
        self.numberOf3Delem = self.numberOfCells * my_variables.nphi
        # variables to contain each 3D variables
        self.geo_index_first_cell_in_2D_frame_iphi = 0
        self.geo_position_x_y_z = np.empty((self.numberOf3Delem, self.ndimProblem))
        self.geo_position_r_z_phi = np.empty((self.numberOf3Delem, self.ndimProblem))
        self.geo_position_f_ism = np.empty((self.numberOf3Delem,))
        self.geo_position_f_jet = np.empty((self.numberOf3Delem,))

        self.geo_position_cellsPoints_volume = np.empty((self.numberOf3Delem,))
        self.geo_position_cellsPoints_masse = np.empty((self.numberOf3Delem,))
        self.geo_position_cellsPoints_masse_ism = np.empty((self.numberOf3Delem,))
        self.geo_position_cellsPoints_masse_jet = np.empty((self.numberOf3Delem,))
        
        if(my_variables.floor_ism):
            self.geo_position_f_ism_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_f_ism_exc = np.empty((self.numberOf3Delem,))
            self.geo_position_cellsPoints_masse_ism_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_cellsPoints_masse_ism_exc = np.empty((self.numberOf3Delem,))
        else: del self.geo_position_f_ism_inc, self.geo_position_f_ism_exc,self.geo_position_cellsPoints_masse_ism_inc,self.geo_position_cellsPoints_masse_ism_exc
        
        if(my_variables.floor_jet):
            self.geo_position_f_jet_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_f_jet_exc = np.empty((self.numberOf3Delem,))
            self.geo_position_cellsPoints_masse_jet_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_cellsPoints_masse_jet_exc = np.empty((self.numberOf3Delem,))
        else: del self.geo_position_f_jet_inc, self.geo_position_f_jet_exc,self.geo_position_cellsPoints_masse_jet_inc,self.geo_position_cellsPoints_masse_jet_exc
        
        if(my_variables.frame_the_box):
            self.geo_position_f_box_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_f_box_exc = np.empty((self.numberOf3Delem,))
            #geo_position_masse_box_inc = np.empty((numberOf3Delem,))
            #geo_position_masse_box_exc = np.empty((numberOf3Delem,))
            self.geo_position_cellsPoints_masse_box_inc = np.empty((self.numberOf3Delem,))
            self.geo_position_cellsPoints_masse_box_exc = np.empty((self.numberOf3Delem,))
        else: del self.geo_position_f_box_inc, self.geo_position_f_box_exc,self.geo_position_cellsPoints_masse_box_inc,self.geo_position_cellsPoints_masse_box_exc        
        
        self.geo_cells_r_z_phi_v = np.empty((self.numberOf3Delem, self.ndimProblem))
        self.geo_position_3D_cell_index_in_2D_frame = np.full((self.numberOf3Delem,), 0)
        self.geo_index_first_cell_in_2D_frame_iphi = np.full((my_variables.nphi,), 0)
        self.geo_rho = np.empty((self.numberOf3Delem,))
        # index to iterate over 3D-geo_* array's 1st dimension (= ID n° of the 3D-cell)
        self.geo_3D_index = 0
        
        #2D variables for the voxelling
        self.geo_position_Gamma_Delta_alpha = np.empty((self.numberOf3Delem, 3))
        if(my_variables.voxelling_mode == 'Voxel'):
            self.geo_vobs = np.empty((self.numberOf3Delem,))
            
        """End of parameters updated values"""
        """==================================================================="""
    
        self.my_data = set(dir()) - self.not_my_data
        self.sorted_my_data = sorted(list(self.my_data), key=lambda x:x.lower())
        
    
        #avoid putting list of variable between not_my_data and my_data, since
        #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
        #only for the Bazooka class attribute
        
        self.list_of_variables = []
        self.list_of_variables_print = ''.ljust(10,"=")
        self.list_of_variables_print += "List of user-defined parameters :\n"
        #print(list_of_variables_print)
        for name in self.sorted_my_data :
            if(name!='not_my_data'):
                myvalue = eval(name)
                self.list_of_variables.append(name)
                self.list_of_variables_print = self.list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
        self.list_of_variables_print = self.list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'

class Railgun:
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(Variable.input_file_vtk)
    reader.Update()
    output = reader.GetOutput()
    
    
    #below are the list of global parameters AND their default value:
    not_my_data = set(dir()) #add your global variables after this line
    
    
    """==================================================================="""
    """Parameters default values"""

    numberOfCells = output.GetNumberOfCells()
    ndimProblem = int(np.shape(output.GetBounds())[0]/2)
    numbeOfVertixes = 2**ndimProblem
    numberOfVertixesUnderDimension = 2**(ndimProblem-1)
    
    nodes_r_z_phi = np.empty((numberOfCells, 
                             numberOfVertixesUnderDimension, 
                             ndimProblem))
    cells_drr_dz_dphi = np.empty((numberOfCells, ndimProblem))
    geo_dangle_phi = 2*np.pi/Variable.nphi
    
    cells_Volume = None
    cells_Masse = None
    cells_Masse_ISM = None
    cells_Masse_jet = None
    cells_Masse_inc = None
    cells_Masse_exc = None
    
    geo_vobs = None
    
    numberOf3Delem = numberOfCells * Variable.nphi
    
    cells_Delta_Alpha = np.empty((numberOfCells, ndimProblem-1))
    
    cells_f_ISM = None
    cells_f_jet = None
    f_box_2D_inc = None
    f_box_2D_exc = None
    
    bound_radius = 0
    bound_z = 0 
        
    """End of parameters default values"""
    """==================================================================="""
    
    my_data = set(dir()) - not_my_data
    sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())
    

    #avoid putting list of variable between not_my_data and my_data, since
    #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
    #only for the Railgun class attribute
    
    
    list_of_variables = []
    list_of_variables_print = ''.ljust(10,"=")
    list_of_variables_print += "List of user-defined parameters :\n"
    #print(list_of_variables_print)
    
    for name in sorted_my_data :
        if(name!='not_my_data'):
            myvalue = eval(name)
            list_of_variables.append(name)
            list_of_variables_print = list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
    list_of_variables_print = list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    #print(list_of_variables_print)
    
    
    def __init__(self):
        for varname in Railgun.list_of_variables:
            setattr(self,varname,getattr(Railgun,varname))
            
    def __repr__(self):
        o='{'
        for i,varname in enumerate(Railgun.list_of_variables):
            o=o+varname+' : '+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+=' , '
        o+='}'
        return o
    
    def __str__(self):
        o=r""
        for i,varname in enumerate(Railgun.list_of_variables):
            o=o+varname+r" : "+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+="\n"
        return o
    
    def set_attributes(self,output,my_variables):   
        #below are the list of the parameters AND their updated value:
        self.not_my_data = set(dir()) #add your global updated variables after this line
        
        """==================================================================="""
        """Parameters updated values"""
        
        self.numberOfCells = output.GetNumberOfCells()
        self.ndimProblem = int(np.shape(output.GetBounds())[0]/2)
        self.numbeOfVertixes = 2**self.ndimProblem
        self.numberOfVertixesUnderDimension = 2**(self.ndimProblem-1)
        
        self.nodes_r_z_phi = np.empty((self.numberOfCells, 
                                       self.numberOfVertixesUnderDimension, 
                                       self.ndimProblem))
        self.cells_drr_dz_dphi = np.empty((self.numberOfCells, self.ndimProblem))
        self.geo_dangle_phi = 2*np.pi/my_variables.nphi
        
        
        
        
        
        
        self.numberOf3Delem = self.numberOfCells * my_variables.nphi
    
        self.cells_Delta_Alpha = np.empty((self.numberOfCells, self.ndimProblem-1))
        
            
        """End of parameters updated values"""
        """==================================================================="""
    
        self.my_data = set(dir()) - self.not_my_data
        self.sorted_my_data = sorted(list(self.my_data), key=lambda x:x.lower())
        
    
        #avoid putting list of variable between not_my_data and my_data, since
        #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
        #only for the Railgun class attribute
        
        self.list_of_variables = []
        self.list_of_variables_print = ''.ljust(10,"=")
        self.list_of_variables_print += "List of user-defined parameters :\n"
        #print(list_of_variables_print)
        for name in self.sorted_my_data :
            if(name!='not_my_data'):
                myvalue = eval(name)
                self.list_of_variables.append(name)
                self.list_of_variables_print = self.list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
        self.list_of_variables_print = self.list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'

        
 
    

    
class Sky_Grid:
     
    #below are the list of global parameters AND their default value:
    not_my_data = set(dir()) #add your global variables after this line

    """==================================================================="""
    """Parameters default values"""
    
    #2D sky plane sampling grid nodes
    x = 0
    y = 0
    z = 0
    dx = 0
    dy = 0
    dz = 0
    
    X = 0
    Y = 0
    Z = 0
    
    nx=2
    ny=2
    nz=2
    
    masse_arrays = np.zeros((nx, ny, nz), dtype=float)
    masse_arrays_ism = np.zeros((nx, ny, nz), dtype=float)
    masse_arrays_jet = np.zeros((nx, ny, nz), dtype=float)
    if(Variable.frame_the_box): 
        masse_arrays_box_inc = np.zeros((nx, ny, nz), dtype=float)
        masse_arrays_box_exc = np.zeros((nx, ny, nz), dtype=float)
    if(Variable.floor_ism):
        masse_arrays_ism_inc = np.zeros((nx, ny, nz), dtype=float)
        masse_arrays_ism_exc = np.zeros((nx, ny, nz), dtype=float)
    if(Variable.floor_jet):
        masse_arrays_jet_inc = np.zeros((nx, ny, nz), dtype=float)
        masse_arrays_jet_exc = np.zeros((nx, ny, nz), dtype=float)
    
    
    """End of parameters default values"""
    """==================================================================="""
    
    my_data = set(dir()) - not_my_data
    sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())
    
    list_of_variables = []
    list_of_variables_print = ''.ljust(10,"=")
    list_of_variables_print += "List of user-defined parameters :\n"
    
    for name in sorted_my_data :
        if(name!='not_my_data'):
            myvalue = eval(name)
            list_of_variables.append(name)
            list_of_variables_print = list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
    list_of_variables_print = list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    
    
    def __init__(self):
        for varname in Sky_Grid.list_of_variables:
            setattr(self,varname,getattr(Sky_Grid,varname))
            
    def __repr__(self):
        o='{'
        for i,varname in enumerate(Sky_Grid.list_of_variables):
            o=o+varname+' : '+str(getattr(self,varname))
            if(i!=np.size(Sky_Grid.list_of_variables)-1): o+=' , '
        o+='}'
        return o
    
    def __str__(self):
        o=r""
        for i,varname in enumerate(Sky_Grid.list_of_variables):
            o=o+varname+r" : "+str(getattr(self,varname))
            if(i!=np.size(Sky_Grid.list_of_variables)-1): o+="\n"
        return o
    
    def set_attributes(self,my_variables):
        self.not_my_data = set(dir()) #add your updated variables after this line 
        #2D sky plane sampling grid nodes
 
        """==================================================================="""
        """Parameters updated values"""

        self.masse_arrays = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        self.masse_arrays_ism = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        self.masse_arrays_jet = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        if(my_variables.frame_the_box): 
            self.masse_arrays_box_inc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
            self.masse_arrays_box_exc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        if(my_variables.floor_ism):
            self.masse_arrays_ism_inc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
            self.masse_arrays_ism_exc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        if(my_variables.floor_jet):
            self.masse_arrays_jet_inc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
            self.masse_arrays_jet_exc = np.zeros((self.nx, self.ny, self.nz), dtype=float)
        

        """End of parameters updated values"""
        """==================================================================="""
        
        self.my_data = set(dir()) - self.not_my_data
        self.sorted_my_data = sorted(list(self.my_data), key=lambda x:x.lower())
        
        self.list_of_variables = []
        self.list_of_variables_print = ''.ljust(10,"=")
        self.list_of_variables_print += "List of user-defined parameters :\n"
        
        for name in self.sorted_my_data :
            if(name!='not_my_data'):
                myvalue = eval(name)
                self.list_of_variables.append(name)
                self.list_of_variables_print = self.list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
        self.list_of_variables_print = self.list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    
    def redefine_nx_ny_nz_arrays(self,my_variables,nx,ny,nz):
        self.masse_arrays = np.zeros((nx, ny, nz), dtype=float)
        self.masse_arrays_ism = np.zeros((nx, ny, nz), dtype=float)
        self.masse_arrays_jet = np.zeros((nx, ny, nz), dtype=float)
        if(my_variables.frame_the_box): 
            self.masse_arrays_box_inc = np.zeros((nx, ny, nz), dtype=float)
            self.masse_arrays_box_exc = np.zeros((nx, ny, nz), dtype=float)
        
    
class Transverse_Layer:
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(Variable.input_file_vtk)
    reader.Update()
    output = reader.GetOutput()
    
    
    #below are the list of global parameters AND their default value:
    not_my_data = set(dir()) #add your global variables after this line
    
    
    """==================================================================="""
    """Parameters default values"""
    
    numberOfCells = output.GetNumberOfCells()
    ndimProblem = int(np.shape(output.GetBounds())[0]/2)
    numbeOfVertixes = 2**ndimProblem
    numberOfVertixesUnderDimension = 2**(ndimProblem-1)
    

    cells_Id = np.full((numberOfCells, 
                             numberOfVertixesUnderDimension), 0)
    cells_r_z_phi = np.empty((numberOfCells, 
                                                  numberOfVertixesUnderDimension, 
                                                  ndimProblem))
    cells_drr = np.empty((numberOfCells,))
    cells_surface = np.empty((numberOfCells,))
    
    fluxes = []
    fluxes_surface_center = []
    fluxes_surface_inclination =[] 
    fluxes_surface_r_z = []
    fluxes_surface_radius_z_inclination = []
        
    """End of parameters default values"""
    """==================================================================="""
    
    my_data = set(dir()) - not_my_data
    sorted_my_data = sorted(list(my_data), key=lambda x:x.lower())
    

    #avoid putting list of variable between not_my_data and my_data, since
    #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
    #only for the Bazooka class attribute
    
    
    list_of_variables = []
    list_of_variables_print = ''.ljust(10,"=")
    list_of_variables_print += "List of user-defined parameters :\n"
    #print(list_of_variables_print)
    
    for name in sorted_my_data :
        if(name!='not_my_data'):
            myvalue = eval(name)
            list_of_variables.append(name)
            list_of_variables_print = list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
    list_of_variables_print = list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'
    #print(list_of_variables_print)
    
    
    def __init__(self):
        for varname in Transverse_Layer.list_of_variables:
            setattr(self,varname,getattr(Transverse_Layer,varname))
            
    def __repr__(self):
        o='{'
        for i,varname in enumerate(Transverse_Layer().list_of_variables):
            o=o+varname+' : '+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+=' , '
        o+='}'
        return o
    
    def __str__(self):
        o=r""
        for i,varname in enumerate(Transverse_Layer.list_of_variables):
            o=o+varname+r" : "+str(getattr(self,varname))
            if(i!=np.size(self.list_of_variables)-1): o+="\n"
        return o
    
    def set_attributes(self,output,my_variables):   
        #below are the list of the parameters AND their updated value:
        self.not_my_data = set(dir()) #add your global updated variables after this line
        
        """==================================================================="""
        """Parameters updated values"""
        
        self.numberOfCells = output.GetNumberOfCells()
        self.ndimProblem = int(np.shape(output.GetBounds())[0]/2)
        self.numbeOfVertixes = 2**self.ndimProblem
        self.numberOfVertixesUnderDimension = 2**(self.ndimProblem-1)
        

        self.cells_Id = np.full((self.numberOfCells, 
                                 self.numberOfVertixesUnderDimension), 0)
        self.cells_r_z_phi = np.empty((self.numberOfCells, 
                                                      self.numberOfVertixesUnderDimension, 
                                                      self.ndimProblem))
        self.cells_drr = np.empty((self.numberOfCells,))
        self.cells_surface = np.empty((self.numberOfCells,))
        
            
        """End of parameters updated values"""
        """==================================================================="""
    
        self.my_data = set(dir()) - self.not_my_data
        self.sorted_my_data = sorted(list(self.my_data), key=lambda x:x.lower())
        
    
        #avoid putting list of variable between not_my_data and my_data, since
        #it may cause infinite loop. Maybe! In any case, it isn't useful since it is
        #only for the Bazooka class attribute
        
        self.list_of_variables = []
        self.list_of_variables_print = ''.ljust(10,"=")
        self.list_of_variables_print += "List of user-defined parameters :\n"
        #print(list_of_variables_print)
        for name in self.sorted_my_data :
            if(name!='not_my_data'):
                myvalue = eval(name)
                self.list_of_variables.append(name)
                self.list_of_variables_print = self.list_of_variables_print + name + ' = ' + str(myvalue) + "\n"
        self.list_of_variables_print = self.list_of_variables_print + ''.ljust(10,"=")+'End of user-defined parameters list'    