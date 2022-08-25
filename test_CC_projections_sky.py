#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

exec(open("./CC_projections_sky.py").read())


to_3D_compar_found_errors = True

v1_list = []
v1_list.append('geo_position_Gamma_Delta_alpha')
if(my_variables.voxelling_mode == 'Voxel'):
    v1_list.append('geo_vobs')





v2_list = []
v2_list.append('my_bazooka.geo_position_Gamma_Delta_alpha')
if(my_variables.voxelling_mode == 'Voxel'):
    v2_list.append('my_bazooka.geo_vobs')



for i,el in enumerate(v1_list):
    truth = False in (getattr(my_bazooka,v2_list[i].split(sep='.')[-1])==eval(el))
    if(truth):
        to_3D_compar_found_errors = to_3D_compar_found_errors & False    
    print(el," : ", not truth)
    
if(not to_3D_compar_found_errors):
    print("The comparaison between hardcoded v1 and harmonised v2 code for CC projections on the sky plane show erros !!")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Please debug it !  <<<<<<<<<<<<<<<<<<<<<<<<<")
else:
    print("The comparaison between hardcoded v1 and harmonised v2 code for CC projections on the sky plane show no errors !!")
    print(" Everything is fine ! You may continue your way....")
    
