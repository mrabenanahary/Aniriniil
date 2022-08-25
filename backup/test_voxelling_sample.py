#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

exec(open("./voxelling_sample.py").read())


to_3D_compar_found_errors = True



v1_list = []
v1_list.append('masse_arrays')
v1_list.append('masse_arrays_ism')
v1_list.append('masse_arrays_jet')
if(frame_the_box): 
    v1_list.append('masse_arrays_box_inc')
    v1_list.append('masse_arrays_box_exc')




v2_list = []
v2_list.append('my_sky_grid.masse_arrays')
v2_list.append('my_sky_grid.masse_arrays_ism')
v2_list.append('my_sky_grid.masse_arrays_jet')
if(frame_the_box): 
    v2_list.append('my_sky_grid.masse_arrays_box_inc')
    v2_list.append('my_sky_grid.masse_arrays_box_exc')


for i,el in enumerate(v1_list):
    truth = False in (getattr(my_sky_grid,v2_list[i].split(sep='.')[-1])==eval(el))
    if(truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", not truth)


if(not to_3D_compar_found_errors):
    print("The comparaison between hardcoded v1 and harmonised v2 code for voxelling grid meshing on the sky show erros !!")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Please debug it !  <<<<<<<<<<<<<<<<<<<<<<<<<")
else:
    print("The comparaison between hardcoded v1 and harmonised v2 code for voxelling grid meshing on the sky show no errors !!")
    print(" Everything is fine ! You may continue your way....")
    
