#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

exec(open("./voxelling_mesh.py").read())


to_3D_compar_found_errors = True

v1_list = []
v1_list.append('x')
v1_list.append('y')
v1_list.append('z')





v2_list = []
v2_list.append('my_sky_grid.x')
v2_list.append('my_sky_grid.y')
v2_list.append('my_sky_grid.z')


for i,el in enumerate(v1_list):
    truth = False in (getattr(my_sky_grid,v2_list[i].split(sep='.')[-1])==eval(el))
    if(truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", not truth)

v1_list = []
v1_list.append('dx')
v1_list.append('dy')
v1_list.append('dz')
v2_list = []
v2_list.append('my_sky_grid.dx')
v2_list.append('my_sky_grid.dy')
v2_list.append('my_sky_grid.dz')

for i,el in enumerate(v1_list):
    truth = (getattr(my_sky_grid,v2_list[i].split(sep='.')[-1])==eval(el))
    if(not truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", truth)

if(not to_3D_compar_found_errors):
    print("The comparaison between hardcoded v1 and harmonised v2 code for voxelling grid meshing on the sky show erros !!")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Please debug it !  <<<<<<<<<<<<<<<<<<<<<<<<<")
else:
    print("The comparaison between hardcoded v1 and harmonised v2 code for voxelling grid meshing on the sky show no errors !!")
    print(" Everything is fine ! You may continue your way....")
    
