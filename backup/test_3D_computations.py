#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

exec(open("./3D_computations.py").read())



to_3D_compar_found_errors = True

v1_list = []
v1_list.append('geo_index_first_cell_in_2D_frame_iphi')
v1_list.append('geo_position_3D_cell_index_in_2D_frame')
v1_list.append('geo_position_r_z_phi')
v1_list.append('geo_position_cellsPoints_volume')
v1_list.append('geo_position_f_ism')
v1_list.append('geo_position_f_jet')
v1_list.append('cellsPoints_drr_dz_dphi_2D_frame')
if(my_variables.frame_the_box):
    v1_list.append('geo_position_f_box_inc')
    v1_list.append('geo_position_f_box_exc')
v1_list.append('geo_cells_r_z_phi_v')
v1_list.append('geo_position_x_y_z')
v1_list.append('geo_rho')
v1_list.append('geo_position_cellsPoints_masse')
v1_list.append('geo_position_cellsPoints_masse_ism')
v1_list.append('geo_position_cellsPoints_masse_jet')
if(my_variables.frame_the_box):
    v1_list.append('geo_position_cellsPoints_masse_box_inc')
    v1_list.append('geo_position_cellsPoints_masse_box_exc')




v2_list = []
v2_list.append('my_bazooka.geo_index_first_cell_in_2D_frame_iphi')
v2_list.append('my_bazooka.geo_position_3D_cell_index_in_2D_frame')
v2_list.append('my_bazooka.geo_position_r_z_phi')
v2_list.append('my_bazooka.geo_position_cellsPoints_volume')
v2_list.append('my_bazooka.geo_position_f_ism')
v2_list.append('my_bazooka.geo_position_f_jet')
v2_list.append('my_bazooka.cellsPoints_drr_dz_dphi_2D_frame')
if(my_variables.frame_the_box):
    v2_list.append('my_bazooka.geo_position_f_box_inc')
    v2_list.append('my_bazooka.geo_position_f_box_exc')
v2_list.append('my_bazooka.geo_cells_r_z_phi_v')
v2_list.append('my_bazooka.geo_position_x_y_z')
v2_list.append('my_bazooka.geo_rho')
v2_list.append('my_bazooka.geo_position_cellsPoints_masse')
v2_list.append('my_bazooka.geo_position_cellsPoints_masse_ism')
v2_list.append('my_bazooka.geo_position_cellsPoints_masse_jet')
if(my_variables.frame_the_box):
    v2_list.append('my_bazooka.geo_position_cellsPoints_masse_box_inc')
    v2_list.append('my_bazooka.geo_position_cellsPoints_masse_box_exc')


for i,el in enumerate(v1_list):
    truth = False in (getattr(my_bazooka,v2_list[i].split(sep='.')[-1])==eval(el))
    if(truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", not truth)

v1_list = []    
v1_list.append('geo_3D_index')
v2_list = []
v2_list.append('my_bazooka.geo_3D_index')

for i,el in enumerate(v1_list):
    truth = (getattr(my_bazooka,v2_list[i].split(sep='.')[-1])==eval(el))
    if(not truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", truth)
    print(el, " : ", truth)
    
if(not to_3D_compar_found_errors):
    print("The comparaison between hardcoded v1 and harmonised v2 code for 3D computations show erros !!")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Please debug it !  <<<<<<<<<<<<<<<<<<<<<<<<<")
else:
    print("The comparaison between hardcoded v1 and harmonised v2 code for 3D computations show no errors !!")
    print(" Everything is fine ! You may continue your way....")
    
