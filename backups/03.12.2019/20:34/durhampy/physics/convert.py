#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 00:06:00 2019

@author: mialy94

=================================================================
= Auteur : Mialy RABENANAHARY
= 
= convert.py
= 
= Module contenant toutes les fonctions qui permettent de convertir les grandeurs
= physiques en d'autres 
=
=================================================================
"""

from astropy import units as u
from astropy import constants as const

global c_lumiere_cgs

c_lumiere_cgs = 2.99792458e10 #cm/s

def WavelengthTONu(wavelength,
                   unit_of_wavelength = u.meter,
                   unit_of_nu = u.GHz,
                   return_quantities = False):

    c = c_lumiere_cgs * u.cm/u.s
    
    lbd = wavelength *unit_of_wavelength 
    nu = c.to(u.m/u.s) / lbd.to(u.m)
    if(return_quantities):
        output = nu.to(unit_of_nu)
        try:
            assert type(output)==type(1*u.m)
        except AssertionError:
                raise AssertionError('Erreur de conversion vers une Quantité Astropy avec unités!! La quantité rétournée n\'en est pas une.')
    else:
        output = (nu.to(unit_of_nu)).value
        try:
            assert type(output)==(type((c.to(u.m/u.s).value)) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantité Astropy avec unités vers une quantité sans unité (e.g. float) !!! ')
    #print('Continue?')
    return output
    

"""def BrillanceSurfaceTOIntegratedIntensity(brillance_de_surface,
                                          unit_of_input = 'cgs',
                                          unit_of_output = 'cgs', 
                                          wavelength = ):
"""

a = WavelengthTONu(wavelength=145,unit_of_wavelength=u.micrometer,unit_of_nu=u.GHz)
print(a)