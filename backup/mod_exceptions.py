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
#global globals_has_been_called


from cycler import cycler
import numpy as np
import csv
import os
from math import ceil
from astropy import units as u
import time
import struct


def test_voxelling_mode(voxelling_mode,voxelling_mode_list,voxelling_parameter_name='Voxelling mode'):
    try:
        assert voxelling_mode in voxelling_mode_list
    except AssertionError:
        s = ' '
        for el in voxelling_mode_list:
            s += str(el) + ' ||'
        raise AssertionError(
            ' {2} \'{0}\' isn\'t in the list. Possible Modes = {1}'.format(voxelling_mode, s,voxelling_parameter_name))
        
def test_attribute_mode(filename,attribute_mode):
    if(attribute_mode != 'point') and (attribute_mode != 'cell'):
        raise ValueError('In input parfile {0}'.format(str(filename)) +
                         ' : attribute_mode=geo_2D_attribute_mode is neither a \" point \" or a \" cell \"')
        
def test_if_isfile(filename):
    if(not os.path.isfile(filename)):
                raise FileNotFoundError("Parameters file {0} not found !".format(filename))
                
def test_if_there_is_directory(filename):
    modified_filename = filename.split(sep='/')
    modified_filename = '/'.join(modified_filename[:-1])
    try:
        assert os.path.isdir(modified_filename)
    except  NotADirectoryError:
        raise NotADirectoryError('File {0} can\'t be opened. {1} doesn\'t exist or is\'n t a directory'.format(filename,modified_filename))
    
def is_mode_implemented(string_mode,possible_mode,mode_name=''):
    if(string_mode not in possible_mode):
            raise NotImplementedError("{0} \'{1}\' isn\'t implemented yet !\n Current possible choices are {2}".format(mode_name,string_mode,possible_mode))
            
            
def test_angle_i_value(angle_i,varname='angle_i',vmin=-np.pi,vmax=np.pi, beginning_msg=''):
    try:
        assert ((angle_i>vmin)&(angle_i<=vmax))
    except  ValueError:
        raise ValueError('{4}The angle {0} hasn\' correct value : it must be between {1} and {2} and here it values {3}'.format(varname,vmin,vmax,angle_i,beginning_msg))