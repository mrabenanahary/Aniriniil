#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 23:24:13 2019

@author: mialy94
=================================================================
= Auteur : Mialy RABENANAHARY
= 
= units.py
= 
= Module contenant toutes les constantes globales (k_B,c,h,etc....)
= dans les différents systemes d'unités disponibles
= 
=
=================================================================
"""

from astropy import units as u
from astropy import constants as astropy_const

global c_lumiere_cgs, sample_of_dimensionless_number, hPlanck_cgs, k_B_cgs

c_lumiere_cgs = astropy_const.c.to_value(u.cm/u.s)#2.99792458e10 #cm/s
hPlanck =  astropy_const.h.to(u.erg*u.s) #cgs
hPlanck_cgs = hPlanck.to_value(u.erg*u.s)
k_B = astropy_const.k_B
k_B_cgs = k_B.to_value(u.erg/u.K)

#definitions d'unites
sample_of_dimensionless_number = 1 * u.m/u.micrometer
dixMoins18Watt = u.def_unit(r'\times10^{-18} W',(1e-18)*u.watt)
KelvinDEnergy = u.def_unit(r'K',k_B.to_value(u.Joule/u.K) * u.Joule)