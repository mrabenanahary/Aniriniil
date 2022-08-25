#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:40:51 2019

@author: mialy94
"""
import os
originalDir = os.getcwd()
os.chdir("../../")
from durhampy.physics.convert import *                                               #import d'éléments du main package durhampy
os.chdir(originalDir)

defDocString = lambda classe, nom_de_methode, docstring: exec("{0:}.{1:}.__doc__ = docstring".format(classe,nom_de_methode),locals(),globals())



#def helpDocString(dictionnaire):
    
    