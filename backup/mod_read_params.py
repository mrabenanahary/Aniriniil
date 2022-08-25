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
    




def read_file_param(filename,
                    return_params_dict=True,
                    a_variable = Variable(),
                    print_imported_file = False):
    #global global_variables
    test_if_isfile(filename)
    with open(filename) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        params_dict =  yaml.full_load(file)
        list_of_unvalid_parameters = {}
        all_parameters_valid = True
        if(print_imported_file) : 
            print('\n')
            print("Non-defaulted imported parameters are :")
            for item,doc in params_dict.items():
                print(item+" : "+str(doc))
            print('\n')
            
        for item,doc in params_dict.items():
            if(item in Variable.list_of_variables):
                setattr(a_variable,item,doc)
            else:
                all_parameters_valid = False
                list_of_unvalid_parameters[item]=doc
                
        if all_parameters_valid==False : raise AttributeError('Invalid unexisting parameter entered for a Variable: {0}'.format(list_of_unvalid_parameters)) 
        is_mode_implemented(a_variable.geo_2D_attribute_mode,Variable.possible_geo_2D_attribute_mode,mode_name='geo_attribute_mode')
        is_mode_implemented(a_variable.voxelling_mode,Variable.voxelling_mode_list,mode_name='voxelling_mode')
        is_mode_implemented(a_variable.voxelling_method,Variable.voxelling_method_list,mode_name='voxelling_method')
        is_mode_implemented(a_variable.voxelling_scheme,Variable.voxelling_scheme_list[a_variable.voxelling_method],mode_name='voxelling_scheme')
        
        test_angle_i_value(a_variable.geo_angle_i,varname='geo_angle_i',vmin=-np.pi,vmax=np.pi,beginning_msg='From/in the parameters input file {0} : '.format(filename))
    
        #print(">>>>>\n", global_variables)
        if return_params_dict==True: return a_variable,params_dict
        else : return a_variable
        
def read_Variable_file(sysargv,my_variables,print_imported_file=False):
    if len(sysargv) >= 2:
        test_if_isfile(sysargv[1])
        my_variables.input_par = str(sysargv[1])
    else:
        test_if_there_is_directory(Variable.input_par)
        my_variables.input_par = str(Variable.input_par)

    varfilename = my_variables.input_par                
    #read the parameters file
    my_variables = read_file_param(varfilename,
                        return_params_dict=False,a_variable=my_variables,
                        print_imported_file=print_imported_file)
    
    if len(sysargv) >= 3:
        test_if_isfile(sysargv[2])
        my_variables.input_file_vtk = str(sysargv[2])    
    
    if len(sysargv) >= 4:
        test_if_isfile(sysargv[3])
        my_variables.output_file_vtk = str(sysargv[3])
    
    return my_variables
    
def extract_vtk_CellsData(my_variables,return_data_props=True):
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(my_variables.input_file_vtk)
    reader.Update()
    output = reader.GetOutput()
    
    numberOfCells = output.GetNumberOfCells()
    CellsData = reader.GetOutput()
    
    ndimProblem = int(np.shape(CellsData.GetBounds())[0]/2)
    numbeOfVertixes = 2**ndimProblem
    numberOfVertixesUnderDimension = 2**(ndimProblem-1)
    
    if(return_data_props): return CellsData, numberOfCells, ndimProblem, numbeOfVertixes, numberOfVertixesUnderDimension
    else : return CellsData
    
def transfert_local_var_to_global_var(local_var):
    return local_var
            
def transfert_global_var_to_local_var(global_var):
    return global_var
"""
if(not __name__ == "__main__"):
    #global global_variables 
    global_variables.dummy_test['mod_read_params']=True
    global_variables.fac_dx = None
else:
    print("From mod_read_params.py :", Variable.nphi_choice_msg)
    a,b = read_file_param('test_input.par')
    c = Variable()
    global_variables = transfert_local_var_to_global_var(c)
    global_variables.fac_dy = False
    d= Variable()
    global_variables = transfert_local_var_to_global_var(d)
    e= Variable()
    e.fac_dy = None
    e = transfert_global_var_to_local_var(global_variables)
    print("global var before :")
    print(global_variables)
    e.fac_dz = False
    e = transfert_global_var_to_local_var(global_variables)
    print("global var after :")
    print(global_variables)
    
"""
    
    
    #for names in Variable.list_of_variables:
    #    print(names,' : ',getattr(a,names))
    #    print(names,' : ',getattr(Variable,names))