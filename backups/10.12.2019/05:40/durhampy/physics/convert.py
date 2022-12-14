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
from astropy import constants as astropy_const
import os
originalDir = os.getcwd()
os.chdir("../../")
from durhampy.data.units import *
os.chdir(originalDir)
import inspect


#liste des convertisseurs disponible dans la fonction globale de conversion physique
global acceptable_inputs, num_to_acceptable_inputs_map

acceptable_inputs =     ["frequency-to-wavelength",
                         "frequency-to-energy",
                         "wavelength-to-frequency",
                         "wavelength-to-energy",
                         "energy-to-frequency",
                         "energy-to-wavelength",
                         "surfaceBrightness-to-integratedIntensity",
                         "surfaceFlux-to-integratedIntensity",
                         "integratedIntensity-to-surfaceFlux",
                         "surfaceBrightness-to-surfaceFlux",
                         "surfaceFlux-to-surfaceBrightness",
                         "integratedIntensity-to-surfaceBrightness"]

#pour utiliser les switch case sous Python
acceptable_inputs_to_num_map = {} 
i = 0
for entrees in acceptable_inputs:
    acceptable_inputs_to_num_map[entrees] = i
    i=i+1

num_to_acceptable_inputs_map = {}

for entrees in acceptable_inputs_to_num_map.keys():
    num_to_acceptable_inputs_map[acceptable_inputs_to_num_map[entrees]] = entrees
    
try:
    bool_map = True
    for entrees in acceptable_inputs:
        bool_map = bool_map and (entrees==num_to_acceptable_inputs_map[acceptable_inputs_to_num_map[entrees]])
    assert bool_map
except AssertionError:
    raise AssertionError("Erreur dans le mapping des conversions physiques possibles >>> {0:} >>>> {1:}".format(acceptable_inputs_to_num_map,num_to_acceptable_inputs_map))

"""=======================Auto-coh??rence unit??s-grandeurs physiques========================"""

def raiseInvalidUnit(quantity,correct_unit,nom_de_fonction):
    try:
        assert ((quantity/(1*correct_unit)).decompose()).unit==(sample_of_dimensionless_number.decompose()).unit
    except AssertionError:
        raise AssertionError("Erreur d'unit?? renseign??e dans la fonction {2:}(): la quantit?? astropy {0:} a ??t?? entr??e alors que la grandeur physique \n correcte attendue correspond ?? celle en unit?? de {1:}".format(repr(quantity),repr(correct_unit),nom_de_fonction))






"""=======================Conversion entre longueur d'onde et fr??quence========================"""
def WavelengthTONu(wavelength,
                   unit_of_wavelength = u.meter,
                   unit_of_nu = u.GHz,
                   return_quantities = False):
    
    """
    Fonction qui convertit une longueur d'onde l dans une unit?? donn??e en 
    fr??quence nu en une autre unit?? donn??e, via la formule :
        nu = c/l
        
    Utilisation:
            
        nu = WavelengthTONu(wavelength,**kwargs):                
                
    Param??tres d'entr??e:
    ===================
        Param??tres obligatoires:
        ------------------------
        
        wavelength : r??el indiquant la valeur de la longueur d'onde ?? convertir exprim??e 
            en l'unit?? indiqu??e par unit_of_wavelength
        
        Param??tres facultatifs, **kwargs:
        ---------------------------------
            
        unit_of_wavelength: unit?? de la longueur d'onde renseign??e via le param??tre 
            wavelength. Type : astropy.units.core.IrreducibleUnit. 
            Valeur par d??faut : u.m (pour le m??tre S.I.)

        unit_of_nu: unit?? de la fr??quence issue de la conversion et retourn??e
            par cette fonction. Type : astropy.units.core.IrreducibleUnit. 
            Valeur par d??faut : u.GHz (pour le gigahertz d??riv?? du S.I.)
            
        return_quantities : bool??en indiquant si la fonction doit retourner
            une valeur de fr??quence de type astropy.units.quantity.Quantity (si return_quantities==True)
            dans l'unit?? de unit_of_nu ou simplement la valeur sous type float (si return_quantities==False) dans l'unit??
            de unit_of_nu. Par d??faut False.
            
    Param??tres de sortie:
    ===================
    
    nu : valeur en unit?? de unit_of_nu de la fr??quence convertie ?? partir de la valeur de longueur d'onde
        renseign??e dans wavelength en unit?? de unit_of_wavelength.
        Si return_quantities==True : nu est de type astropy.units.quantity.Quantity
        Si return_quantities==False : nu est de type float
                        
        """
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(wavelength*unit_of_wavelength,u.m,"WavelengthTONu")
    raiseInvalidUnit(1*unit_of_nu,u.Hz,"WavelengthTONu")
    
    c = c_lumiere_cgs * u.cm/u.s
    
    lbd = wavelength *unit_of_wavelength 
    nu = c.to(u.m/u.s) / lbd.to(u.m)
    if(return_quantities):
        output = nu.to(unit_of_nu)
        try:
            assert type(output)==type(1*u.m)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output = (nu.to_value(unit_of_nu))
        try:
            assert type(output)==(type((c.to_value(u.m/u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output

def NuTOWavelength(nu,
                   unit_of_nu = u.GHz,
                   unit_of_wavelength = u.meter,
                   return_quantities = False):
    """
    Fonction qui convertit une fr??quence nu dans une unit?? donn??e en 
    longueur d'onde l en une autre unit?? donn??e, via la formule :
        l = c/nu
        
    Utilisation:
            
        wavelength = NuTOWavelength(nu,**kwargs):                
                
    Param??tres d'entr??e:
    ===================
        Param??tres obligatoires:
        ------------------------
        
        nu : r??el indiquant la valeur de la fr??quence ?? convertir exprim??e 
            en l'unit?? indiqu??e par unit_of_nu
        
        Param??tres facultatifs, **kwargs:
        ---------------------------------

        unit_of_nu: unit?? de la fr??quence renseign??e via le param??tre nu.
            Type : astropy.units.core.IrreducibleUnit. 
            Valeur par d??faut : u.GHz (pour le gigahertz d??riv?? du S.I.)
            
        unit_of_wavelength: unit?? de la longueur d'onde issue de la conversion et retourn??e
            par cette fonction. Type : astropy.units.core.IrreducibleUnit. 
            Valeur par d??faut : u.m (pour le m??tre S.I.)
            
        return_quantities : bool??en indiquant si la fonction doit retourner
            une valeur de longueur d'onde de type astropy.units.quantity.Quantity (si return_quantities==True)
            dans l'unit?? de unit_of_wavelength ou simplement la valeur sous type float (si return_quantities==False) dans l'unit??
            de unit_of_wavelength. Par d??faut False.
            
    Param??tres de sortie:
    ===================
    
    wavelength : valeur en unit?? de unit_of_wavelength de la longueur d'onde convertie ?? partir de la valeur de fr??quence
        renseign??e dans nu en unit?? de unit_of_nu.
        Si return_quantities==True : wavelength est de type astropy.units.quantity.Quantity
        Si return_quantities==False : wavelength est de type float
                        
        """
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(1*unit_of_wavelength,u.m,"NuTOWavelength")
    raiseInvalidUnit(nu*unit_of_nu,u.Hz,"NuTOWavelength")
    
    c = c_lumiere_cgs * u.cm/u.s
    
    frequency = nu *unit_of_nu
    wavelength = c.to(u.m/u.s) / frequency.to(u.Hz)
    if(return_quantities):
        output = wavelength.to(unit_of_wavelength)
        try:
            assert type(output)==type(1*u.Hz)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output = (wavelength.to_value(unit_of_wavelength))
        try:
            assert type(output)==(type((c.to_value(u.m/u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output

"""=======================(END) Conversion entre longueur d'onde et fr??quence (END)========================"""



"""=======================Conversion entre fr??quence et ??nergie (en unit??s de Kelvin)========================"""

def NuTOEnergie(   nu,
                   unit_of_nu = u.GHz,
                   unit_of_energy = KelvinDEnergy,
                   return_quantities = False):
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(1*unit_of_nu,u.GHz,"NuTOEnergie")
    raiseInvalidUnit(1*unit_of_energy,KelvinDEnergy,"NuTOEnergie")
    
    h = hPlanck_cgs * (u.erg*u.s) # in cgs units system above in the 50th first lines of this scripts
    f = nu * unit_of_nu
    
    Eup = h.to(u.Joule*u.s) * f.to(u.Hz) #in SI units system for convenience of manipulation
    Eup = Eup.to(u.Joule) # in SI for convenience
    kB = k_B.to(u.Joule/u.K) # in SI for convenience
    
    Eup_in_Kelvin = (Eup.to_value(u.Joule) / kB.to_value(u.Joule/u.K)) # in units of Kelvin of energy
    
    #print(Eup, Eup_in_Kelvin)
    #print(Eup.to(u.Joule), Eup.to(u.eV))
    
    if(return_quantities):
        output = (Eup_in_Kelvin*KelvinDEnergy).to(unit_of_energy)
        try:
            assert type(output)==type(1*u.Hz)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output = (Eup_in_Kelvin*KelvinDEnergy).to_value(unit_of_energy)
        try:
            assert type(output)==(type((h.to_value(u.Joule * u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output    

def EnergieTONu(   energy,
                   unit_of_energy = KelvinDEnergy,
                   unit_of_nu = u.GHz,
                   return_quantities = False):
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(1*unit_of_energy,KelvinDEnergy,"EnergieTONu")
    raiseInvalidUnit(1*unit_of_nu,u.GHz,"EnergieTONu")
    
    h = hPlanck_cgs * (u.erg*u.s) # in cgs units system above in the 50th first lines of this scripts
    energyu = energy * unit_of_energy
    
    #energyu = h.to(u.Joule*u.s) * f.to(u.Hz) #in SI units system for convenience of manipulation
    nu = energyu.to(u.Joule) / h.to(u.Joule*u.s) # en Hz
    nu_GHz = nu.to(u.GHz)
       
    #print(Eup, Eup_in_Kelvin)
    #print(Eup.to(u.Joule), Eup.to(u.eV))
    
    if(return_quantities):
        output = nu_GHz.to(unit_of_nu)
        try:
            assert type(output)==type(1*u.Hz)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output =  nu_GHz.to_value(unit_of_nu)
        try:
            assert type(output)==(type((h.to_value(u.Joule * u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output    
    
"""=======================(END) Conversion entre fr??quence et ??nergie (en unit??s de Kelvin) (END)========================"""






"""=======================Conversion entre longueur d'onde et ??nergie (en unit??s de Kelvin)========================"""

def WavelengthTOEnergie(   wavelength,
                   unit_of_wavelength = u.micrometer,
                   unit_of_energy = KelvinDEnergy,
                   return_quantities = False):
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(1*unit_of_wavelength,u.micrometer,"WavelengthTOEnergie")
    raiseInvalidUnit(1*unit_of_energy,KelvinDEnergy,"WavelengthTOEnergie")
    
    h = hPlanck_cgs * (u.erg*u.s) # in cgs units system above in the 50th first lines of this scripts
    c = c_lumiere_cgs * (u.cm/u.s)
    l = wavelength * unit_of_wavelength
    
    Eup = h.to(u.Joule*u.s) * c.to(u.m/u.s) / l.to(u.m) #in SI units system for convenience of manipulation
    Eup = Eup.to(u.Joule) # in SI for convenience
    kB = k_B.to(u.Joule/u.K) # in SI for convenience
    
    Eup_in_Kelvin = (Eup.to_value(u.Joule) / kB.to_value(u.Joule/u.K)) # in units of Kelvin of energy
    
    #print(Eup, Eup_in_Kelvin)
    #print(Eup.to(u.Joule), Eup.to(u.eV))
    
    if(return_quantities):
        output = (Eup_in_Kelvin*KelvinDEnergy).to(unit_of_energy)
        try:
            assert type(output)==type(1*u.Hz)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output = (Eup_in_Kelvin*KelvinDEnergy).to_value(unit_of_energy)
        try:
            assert type(output)==(type((h.to_value(u.Joule * u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output    

def EnergieTOWavelength(   energy,
                   unit_of_energy = KelvinDEnergy,
                   unit_of_wavelength = u.micrometer,
                   return_quantities = False):
    #v??rification d'entr??e des bonnes unit??s physiques    
    raiseInvalidUnit(1*unit_of_energy,KelvinDEnergy,"EnergieTOWavelength")
    raiseInvalidUnit(1*unit_of_wavelength,u.micrometer,"EnergieTOWavelength")
    
    h = hPlanck_cgs * (u.erg*u.s) # in cgs units system above in the 50th first lines of this scripts
    c = c_lumiere_cgs * (u.cm/u.s)
    energyu = energy * unit_of_energy
    
    #energyu = h.to(u.Joule*u.s) * f.to(u.Hz) #in SI units system for convenience of manipulation
    wavelength = h.to(u.Joule*u.s) * c.to(u.m/u.s) / energyu.to(u.Joule)#energyu.to(u.Joule) / h.to(u.Joule*u.s) # en Hz
    wavelength_micrometer = wavelength.to(u.micrometer)
       
    #print(energyu, energyu_in_Kelvin)
    #print(energyu.to(u.Joule), energyu.to(u.eV))
    
    if(return_quantities):
        output = wavelength_micrometer.to(unit_of_wavelength)
        try:
            assert type(output)==type(1*u.Hz)
        except AssertionError:
            raise AssertionError('Erreur de conversion vers une Quantit?? Astropy avec unit??s!! La quantit?? r??tourn??e n\'en est pas une.')
    else:
        output =  wavelength_micrometer.to_value(unit_of_wavelength)
        try:
            assert type(output)==(type((h.to_value(u.Joule * u.s))) or float)
        except AssertionError:
            raise AssertionError('Erreur de conversion d\'une Quantit?? Astropy avec unit??s vers une quantit?? sans unit?? (e.g. float) !!! ')
    #print('Continue?')
    return output    
    
"""=======================(END) Conversion entre longueur d'onde et ??nergie (en unit??s de Kelvin) (END)========================"""


"""=================================Conversion entre brillance de surface, intensit?? int??gr??e, et flux surfacique ========================"""    
""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Fonctions de conversion en intensit?? int??gr??e<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

def SurfaceBrightnessTOIntegratedIntensity(surface_brightness,
                                          wavelength,
                                          unit_of_surface_brightness,
                                          unit_of_wavelength,
                                          unit_of_integrated_intensity,
                                          return_quantities = False):
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """
    nom_de_cette_fonction = "SurfaceBrightnessTOIntegratedIntensity"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    raiseInvalidUnit(surface_brightness*unit_of_surface_brightness, u.erg/u.sr/u.cm**2/u.s, nom_de_cette_fonction)
    raiseInvalidUnit(wavelength*unit_of_wavelength, u.m, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_integrated_intensity, u.K*u.km/u.s, nom_de_cette_fonction)
    
    I_nu_d_nu_cgs = (surface_brightness * unit_of_surface_brightness)      
    #print(I_nu_d_nu_cgs)
    I_nu_d_nu_cgs = I_nu_d_nu_cgs.to(u.erg/u.sr/u.cm**2/u.s)                    #brillance de surface en erg/s/cm^2/sr
    #print(I_nu_d_nu_cgs)
    lbd = (wavelength * unit_of_wavelength)                                     
    lbd = lbd.to(u.micrometer)                                                  #longueur d'onde en microm??tre
    
    coeff_CGS = (3.62146E-2) * u.K * u.sr/u.erg                                 #voir references : bloc_notes.pdf
    
    #print(I_nu_d_nu_cgs.to_value(u.erg/u.sr/u.cm**2/u.s))
    #print( ((lbd.to_value(u.micrometer))**3))
    rTdV = coeff_CGS.to_value(u.K * u.sr/u.erg) * I_nu_d_nu_cgs.to_value(u.erg/u.sr/u.cm**2/u.s) * ((lbd.to_value(u.micrometer))**3)
    if(return_quantities==True) : rTdV = rTdV * u.K * u.km / u.s
    return rTdV

def SurfaceFluxTOIntegratedIntensity(     surface_flux,
                                          wavelength,
                                          theta_s,
                                          unit_of_surface_flux,
                                          unit_of_wavelength,
                                          unit_of_theta_s,
                                          unit_of_integrated_intensity, 
                                          return_quantities = False):
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """
    nom_de_cette_fonction = "SurfaceFluxTOIntegratedIntensity"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    raiseInvalidUnit(surface_flux*unit_of_surface_flux, dixMoins18Watt/u.m**2, nom_de_cette_fonction)
    raiseInvalidUnit(wavelength*unit_of_wavelength, u.m, nom_de_cette_fonction)
    raiseInvalidUnit(theta_s*unit_of_theta_s, u.arcsecond, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_integrated_intensity, u.K*u.km/u.s, nom_de_cette_fonction)
    
    Flux = (surface_flux * unit_of_surface_flux)
    Flux = Flux.to(dixMoins18Watt/u.m**2)                                       #flux en 10^{-18} W/m^2
    lbd = (wavelength * unit_of_wavelength) 
    lbd = lbd.to(u.micrometer)                                                  #longueur d'onde en microm??tre
    thetaS = (theta_s * unit_of_theta_s) 
    thetaS = thetaS.to(u.arcsecond)                                             #dimension angulaire en arcsecondes
    
    
    coeff_CGS = (1.54076e-6) * u.K * u.arcsecond**2 / u.erg                     #voir references : bloc_notes.pdf
    
    rTdV = coeff_CGS.to_value(u.K * u.arcsecond**2 / u.erg) * Flux.to_value(dixMoins18Watt/u.m**2) \
                                                            * ((lbd.to_value(u.micrometer))**3) \
                                                            / ((thetaS.to_value(u.arcsecond))**2)
    if(return_quantities==True) : rTdV = rTdV * u.K * u.km / u.s
    return rTdV

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Fonctions de conversion en flux surfacique<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

def IntegratedIntensityTOSurfaceFlux(     integrated_intensity,
                                          theta_s,
                                          wavelength,
                                          unit_of_integrated_intensity,
                                          unit_of_theta_s,
                                          unit_of_wavelength,
                                          unit_of_surface_flux, 
                                          return_quantities = False):
    
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """
    nom_de_cette_fonction = "IntegratedIntensityTOSurfaceFlux"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    raiseInvalidUnit(integrated_intensity*unit_of_integrated_intensity, u.K*u.km/u.s, nom_de_cette_fonction)
    raiseInvalidUnit(theta_s*unit_of_theta_s, u.arcsecond, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_surface_flux, dixMoins18Watt/u.m**2, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_wavelength, u.micrometer, nom_de_cette_fonction)
    
    TdV_Kkmps = (integrated_intensity * unit_of_integrated_intensity)      
    TdV_Kkmps = TdV_Kkmps.to(u.K*u.km/u.s)                                          #intensit?? int??gr??e en K.km/s
    thetaS = (theta_s * unit_of_theta_s)                                     
    thetaS = thetaS.to(u.arcsecond)                                             #dimension angulaire en arcsecondes
    lbd = (wavelength * unit_of_wavelength) 
    lbd = lbd.to(u.micrometer)                                                  #longueur d'onde en microm??tre
    
    coeff_CGS = (6.49028e5) * u.erg/u.arcsecond**2/u.K                             #voir references : bloc_notes.pdf
    

    rSurfaceFlux = coeff_CGS.to_value(u.erg/u.arcsecond**2/u.K) * TdV_Kkmps.to_value(u.K*u.km/u.s) * (((thetaS.to_value(u.arcsecond))**2)/((lbd.to_value(u.micrometer))**3))
    #print((((u.erg/u.arcsecond**2/u.K*u.K*u.km/u.s*(u.arcsecond**2)/((u.micrometer)**3))).decompose()).to(u.watt/u.m**2))
    if(return_quantities==True) : rSurfaceFlux = rSurfaceFlux * dixMoins18Watt/u.m**2
    return rSurfaceFlux


def SurfaceBrightnessTOSurfaceFlux(       surface_brightness,
                                          theta_s,
                                          unit_of_surface_brightness,
                                          unit_of_theta_s,
                                          unit_of_surface_flux, 
                                          return_quantities = False):
    
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """
    nom_de_cette_fonction = "SurfaceBrightnessTOSurfaceFlux"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    #print(unit_of_surface_brightness)
    #print(surface_brightness*unit_of_surface_brightness)
    raiseInvalidUnit(surface_brightness*unit_of_surface_brightness, u.erg/u.sr/u.cm**2/u.s, nom_de_cette_fonction)
    raiseInvalidUnit(theta_s*unit_of_theta_s, u.arcsecond, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_surface_flux, dixMoins18Watt/u.m**2, nom_de_cette_fonction)
    
    brillanceDeSurface = (surface_brightness * unit_of_surface_brightness)      
    brillanceDeSurface = brillanceDeSurface.to(u.erg/u.sr/u.cm**2/u.s)                                          #intensit?? int??gr??e en K.km/s
    thetaS = (theta_s * unit_of_theta_s)                                     
    thetaS = thetaS.to(u.arcsecond)                                             #dimension angulaire en arcsecondes                                       #longueur d'onde en microm??tre
    
    coeff_CGS = (2.35044e4) * u.sr/u.arcsecond**2                               #voir references : bloc_notes.pdf
    

    rSurfaceFlux = coeff_CGS.to_value(u.sr/u.arcsecond**2) * brillanceDeSurface.to_value(u.erg/u.sr/u.cm**2/u.s) * ((thetaS.to_value(u.arcsecond))**2)
    #print((((u.sr/u.arcsecond**2*u.erg/u.sr/u.cm**2/u.s)*u.arcsecond**2).decompose()).to(u.watt/u.m**2))
    if(return_quantities==True) : rSurfaceFlux = rSurfaceFlux * dixMoins18Watt/u.m**2
    return rSurfaceFlux

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Fonctions de conversion en brillance de surface<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

def SurfaceFluxTOSurfaceBrightness(       surface_flux,
                                          theta_s,
                                          unit_of_surface_flux,
                                          unit_of_theta_s,
                                          unit_of_brightness_surface,
                                          return_quantities = False):
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """    
    nom_de_cette_fonction = "SurfaceFluxTOSurfaceBrightness"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    raiseInvalidUnit(surface_flux*unit_of_surface_flux, dixMoins18Watt/u.m**2, nom_de_cette_fonction)
    raiseInvalidUnit(theta_s*unit_of_theta_s, u.arcsecond, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_brightness_surface, u.erg/u.sr/u.cm**2/u.s, nom_de_cette_fonction)
    
    Flux = (surface_flux * unit_of_surface_flux)
    Flux = Flux.to(dixMoins18Watt/u.m**2)                                       #flux en 10^{-18} W/m^2
    thetaS = (theta_s * unit_of_theta_s) 
    thetaS = thetaS.to(u.arcsecond)                                             #dimension angulaire en arcsecondes
    
    
    coeff_CGS = (4.25452e-5) * u.arcsecond**2 / u.sr                     #voir references : bloc_notes.pdf
    
    rBrillanceSurface = coeff_CGS.to_value(u.arcsecond**2 / u.sr) * Flux.to_value(dixMoins18Watt/u.m**2) \
                                                            / ((thetaS.to_value(u.arcsecond))**2)
    if(return_quantities==True) : rBrillanceSurface = rBrillanceSurface * u.erg/u.sr/u.cm**2/u.s
    return rBrillanceSurface

def IntegratedIntensityTOSurfaceBrightness(integrated_intensity,
                                          wavelength,
                                          unit_of_integrated_intensity,
                                          unit_of_wavelength,
                                          unit_of_brightness_surface,
                                          return_quantities = False):
    
    """ Fonction valid??e comme convertissant correctement les unit??s par Mialy le 04/12/2019 """
    nom_de_cette_fonction = "IntegratedIntensityTOSurfaceBrightness"
    
    #on commence par v??rifier que l'user entre des unit??s valables
    raiseInvalidUnit(integrated_intensity*unit_of_integrated_intensity, u.K*u.km/u.s, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_brightness_surface, u.erg/u.sr/u.cm**2/u.s, nom_de_cette_fonction)
    raiseInvalidUnit(1*unit_of_wavelength, u.micrometer, nom_de_cette_fonction)
    
    TdV_Kkmps = (integrated_intensity * unit_of_integrated_intensity)      
    TdV_Kkmps = TdV_Kkmps.to(u.K*u.km/u.s)                                          #intensit?? int??gr??e en K.km/s                                           #dimension angulaire en arcsecondes
    lbd = (wavelength * unit_of_wavelength) 
    lbd = lbd.to(u.micrometer)                                                  #longueur d'onde en microm??tre
    
    coeff_CGS =  27.6131 * u.erg/u.sr/u.K                             #voir references : bloc_notes.pdf
    

    rBrillanceSurface = coeff_CGS.to_value(u.erg/u.sr/u.K) * TdV_Kkmps.to_value(u.K*u.km/u.s) / ((lbd.to_value(u.micrometer))**3)
    if(return_quantities==True) : rBrillanceSurface = rBrillanceSurface * u.erg/u.sr/u.cm**2/u.s
    return rBrillanceSurface

"""=======================(END) Conversion entre brillance de surface, intensit?? int??gr??e, et flux surfacique (END)========================"""






















"""=================================Fonction globale effectuant automatiquement la conversion entre ??nergie, fr??quence et longueur d'onde ========================"""    



def convertisseurPhysique(
                      type_de_conversion,
                      energie = 500 ,                                           #param??tres pour la conversion entre longueur d'onde, fr??quence et ??negie (et les param??tres de conversion entre brillance de surface, flux surfacique et intensit?? int??gr??e)
                      unite_d_energie = KelvinDEnergy,
                      longueur_d_onde = 145,
                      unite_de_longueur_d_onde = u.micrometer,
                      frequence = 756,
                      unite_de_frequence = u.GHz,
                      brillance_de_surface = 4e-5,                              #param??tres pour la conversion entre brillance de surface, flux surfacique et intensit?? int??gr??e + les param??tres pr??c??dents
                      unite_de_brillance_de_surface = u.erg/u.cm**2/u.s/u.sr,
                      flux_de_surface = 595,
                      unite_de_flux_surfacique = dixMoins18Watt/u.m**2,
                      intensite_integree = 77.32,
                      unite_d_intensite_integree = u.K * u.km / u.s,
                      theta_source = 10,
                      unite_de_theta_source = u.arcsecond,
                      retourner_qtes_physiques=False):
    
    #verifier la validite des unites entrees directement dans les fonctions pr??c??dentes appel??es ci-dessous

    
    #test de validit?? de l'input entr??
    try : 
        assert type_de_conversion in acceptable_inputs
    except AssertionError:
        raise ValueError("L'input renseign?? dans l'argument type_de_conversion n'est pas dans la liste des conversions disponibles: {0:}".format([el for el in acceptable_inputs]))
    
    #modifier les switch cases en cas d'ajout ou de retrait de types de conversions
    
    i_input = acceptable_inputs_to_num_map[type_de_conversion]
    
    if(  i_input == acceptable_inputs_to_num_map["frequency-to-wavelength"]):
        out_put = NuTOWavelength(                        nu=frequence, 
                                                         unit_of_nu=unite_de_frequence,
                                                         unit_of_wavelength=unite_de_longueur_d_onde,
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["frequency-to-energy"]):
        out_put = NuTOEnergie(                           nu=frequence,
                                                         unit_of_nu=unite_de_frequence,
                                                         unit_of_energy=unite_d_energie,
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["wavelength-to-frequency"]):
        out_put = WavelengthTONu(                        wavelength=longueur_d_onde,
                                                         unit_of_wavelength=unite_de_longueur_d_onde,
                                                         unit_of_nu=unite_de_frequence,
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["wavelength-to-energy"]):
        out_put = WavelengthTOEnergie(                   wavelength=longueur_d_onde,
                                                         unit_of_wavelength=unite_de_longueur_d_onde, 
                                                         unit_of_energy=unite_d_energie, 
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["energy-to-frequency"]):
        out_put = EnergieTONu(                           energy=energie,
                                                         unit_of_energy=unite_d_energie,
                                                         unit_of_nu=unite_de_frequence,
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["energy-to-wavelength"]):
        out_put = EnergieTOWavelength(                   energy=energie,
                                                         unit_of_energy=unite_d_energie,
                                                         unit_of_wavelength=unite_de_longueur_d_onde, 
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["surfaceBrightness-to-integratedIntensity"]):
        out_put = SurfaceBrightnessTOIntegratedIntensity(surface_brightness=brillance_de_surface,
                                                         wavelength=longueur_d_onde, 
                                                         unit_of_surface_brightness = unite_de_brillance_de_surface,
                                                         unit_of_wavelength = unite_de_longueur_d_onde,
                                                         unit_of_integrated_intensity = unite_d_intensite_integree, 
                                                         return_quantities = retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["surfaceFlux-to-integratedIntensity"]):
        out_put = SurfaceFluxTOIntegratedIntensity(      surface_flux=flux_de_surface,
                                                         wavelength=longueur_d_onde, 
                                                         theta_s=theta_source,
                                                         unit_of_surface_flux = unite_de_flux_surfacique,
                                                         unit_of_wavelength = unite_de_longueur_d_onde,
                                                         unit_of_theta_s= unite_de_theta_source,
                                                         unit_of_integrated_intensity = unite_d_intensite_integree, 
                                                         return_quantities = retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["integratedIntensity-to-surfaceFlux"]):
        out_put = IntegratedIntensityTOSurfaceFlux(      integrated_intensity=intensite_integree,
                                                         theta_s=theta_source,
                                                         wavelength=longueur_d_onde,
                                                         unit_of_integrated_intensity = unite_d_intensite_integree,
                                                         unit_of_theta_s = unite_de_theta_source,
                                                         unit_of_wavelength=unite_de_longueur_d_onde,
                                                         unit_of_surface_flux = unite_de_flux_surfacique, 
                                                         return_quantities = retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["surfaceBrightness-to-surfaceFlux"]):
        out_put = SurfaceBrightnessTOSurfaceFlux(        surface_brightness=brillance_de_surface,
                                                         theta_s=theta_source,
                                                         unit_of_surface_brightness = unite_de_brillance_de_surface,
                                                         unit_of_theta_s = unite_de_theta_source,
                                                         unit_of_surface_flux = unite_de_flux_surfacique, 
                                                         return_quantities = retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["surfaceFlux-to-surfaceBrightness"]):
        out_put = SurfaceFluxTOSurfaceBrightness(        surface_flux=flux_de_surface,
                                                         theta_s=theta_source,
                                                         unit_of_surface_flux=unite_de_flux_surfacique,
                                                         unit_of_theta_s=unite_de_theta_source,
                                                         unit_of_brightness_surface=unite_de_brillance_de_surface,
                                                         return_quantities=retourner_qtes_physiques)
        
    elif(i_input == acceptable_inputs_to_num_map["integratedIntensity-to-surfaceBrightness"]):
        out_put = IntegratedIntensityTOSurfaceBrightness(integrated_intensity=intensite_integree,
                                                         wavelength=longueur_d_onde,
                                                         unit_of_integrated_intensity=unite_d_intensite_integree,
                                                         unit_of_wavelength=unite_de_longueur_d_onde,
                                                         unit_of_brightness_surface=unite_de_brillance_de_surface,
                                                         return_quantities=retourner_qtes_physiques)
    else:
        raise("Probl??me insoluble, aucun des cas de conversions physiques disponibles n'a ??t?? appel?? !!") #message fatal
    return out_put






















"""=====================================test du module convert.py=============================================="""
if __name__ == "__main__":
    
    a = WavelengthTONu(wavelength=145,unit_of_wavelength=u.micrometer,unit_of_nu=u.GHz)
    print("a=",a)
    
    
    b = NuTOWavelength(nu=2067.5342,unit_of_wavelength=u.micrometer,unit_of_nu=u.GHz,return_quantities=True)
    print("b=",b)
    
    #raiseInvalidUnit(1*u.erg/u.sr/u.cm**2/u.s,u.watt/u.s,"WavelengthTONu")
    
    retQtt = True
    
    
    c = SurfaceBrightnessTOIntegratedIntensity(1,2,unit_of_surface_brightness=u.erg/u.sr/u.cm**2/u.s,unit_of_wavelength=u.micrometer,return_quantities=retQtt, unit_of_integrated_intensity = u.K * u.km/u.s)
    print("c=",c)
    
    d = SurfaceBrightnessTOIntegratedIntensity(1,2,unit_of_surface_brightness=u.watt/u.sr/u.m**2,unit_of_wavelength=u.m, return_quantities=retQtt, unit_of_integrated_intensity = u.K * u.km/u.s)
    print("d=",d)
    print((3.62146e-2)*((1e7)/(1e4))*((2e6)**3))
    
    
    e =  SurfaceFluxTOIntegratedIntensity(1e-6,145,10,unit_of_surface_flux = dixMoins18Watt/u.m**2, unit_of_wavelength = u.micrometer, unit_of_theta_s = u.arcsecond, unit_of_integrated_intensity = u.K*u.km/u.s)
    print("e=",e)
    
    f = SurfaceFluxTOIntegratedIntensity(1e-6,145,10,unit_of_surface_flux = u.erg/u.s/u.cm**2,unit_of_wavelength=u.millimeter,unit_of_theta_s=u.radian,return_quantities=True, unit_of_integrated_intensity = u.K*u.km/u.s)
    print("f=",f)
    print(1.54076*1e-6*1e-6*((1e-7*1e18)/1e-4)*1e-9*((145*1e6)**3) / ((10*206265)**2))
    
    g = IntegratedIntensityTOSurfaceFlux(1,1,1,                                          unit_of_integrated_intensity = u.K*u.km/u.s,
                                          unit_of_theta_s = u.arcsecond,
                                          unit_of_wavelength = u.micrometer,
                                          unit_of_surface_flux = dixMoins18Watt/u.m**2)
    print('g=',g)
    
    h = IntegratedIntensityTOSurfaceFlux(1e7,2,12,unit_of_integrated_intensity = u.K*u.m/u.s,unit_of_theta_s = u.radian,unit_of_wavelength = u.mm,return_quantities=True,unit_of_surface_flux = dixMoins18Watt/u.m**2)
    print("h=",h)
    print(6.49028e5*1e7*1e-3*(((2*206265)**2)/((12*1e-3*1e6)**3)), "<=> h_unit / h_a_la_main = ", h/(6.49028e5*1e7*1e-3*(((2*206265)**2)/((12*1e-3*1e6)**3))))
    
    i = SurfaceBrightnessTOSurfaceFlux(1,1,unit_of_surface_brightness = u.erg/u.sr/u.cm**2/u.s,
                                          unit_of_theta_s = u.arcsecond,
                                          unit_of_surface_flux = dixMoins18Watt/u.m**2)
    print('i=',i)
    
    j = SurfaceBrightnessTOSurfaceFlux(1e7,12,unit_of_surface_brightness = u.watt/u.sr/u.m**2,unit_of_theta_s = u.radian,unit_of_surface_flux=dixMoins18Watt/u.m**2, return_quantities=True)
    print("j=",j)
    print((2.35044e4)*1e7*(1e7/1e4)*(12*206265)**2)
    
    k = SurfaceFluxTOSurfaceBrightness(1,1,unit_of_surface_flux = dixMoins18Watt/u.m**2,
                                          unit_of_theta_s = u.arcsecond,
                                          unit_of_brightness_surface = u.erg/u.sr/u.cm**2/u.s)
    print('k=',k)
    
    l = SurfaceFluxTOSurfaceBrightness(1e7,12,unit_of_surface_flux = u.erg/u.s/u.cm**2,unit_of_theta_s = u.radian,return_quantities=True,unit_of_brightness_surface = u.erg/u.sr/u.cm**2/u.s)
    print("l=",l)
    print(4.25452e-5*1e7*((1e-7*1e18)/(1e-4))*1/((12*206265)**2))
    
    m = IntegratedIntensityTOSurfaceBrightness(1,1,unit_of_integrated_intensity = u.K*u.km/u.s,
                                          unit_of_wavelength = u.micrometer,
                                          unit_of_brightness_surface = u.erg/u.sr/u.cm**2/u.s)
    print('m=',m)
    
    n = IntegratedIntensityTOSurfaceBrightness(1e7,12,unit_of_integrated_intensity = u.K*u.m/u.s,unit_of_wavelength = u.millimeter,return_quantities=True,unit_of_brightness_surface = u.erg/u.sr/u.cm**2/u.s)
    print("n=",n)
    print(27.6131*1e7*(1e-3)*1/((12*1e-3*1e6)**3))
    
    aa = NuTOEnergie(10418.241360755546,return_quantities=False)
    bb = NuTOEnergie(10418.241360755546,return_quantities=True)
    print(aa)
    print(bb)
    
    cc = EnergieTONu(500,return_quantities=False)
    dd = EnergieTONu(500,return_quantities=True)
    print(cc)
    print(dd)
    
    
    ee = NuTOEnergie(10418.241360755546e9,unit_of_nu = u.Hz,return_quantities=False)
    ff = NuTOEnergie(10418.241360755546e9,unit_of_nu = u.Hz,return_quantities=True)
    print(ee)
    print(ff)
    
    gg = EnergieTONu(0.043086404716694296,return_quantities=False,unit_of_energy=u.eV)
    hh = EnergieTONu(0.043086404716694296,return_quantities=True,unit_of_energy=u.eV)
    print(gg)
    print(hh)
    
    del aa,bb,cc,dd,ee,ff,gg,hh
    
    aa = WavelengthTOEnergie(28.77583,return_quantities=False)
    bb = WavelengthTOEnergie(28.77583,return_quantities=True)
    print(aa)
    print(bb)
    
    cc = EnergieTOWavelength(500,return_quantities=False)
    dd = EnergieTOWavelength(500,return_quantities=True)
    print(cc)
    print(dd)
    
    
    ee = WavelengthTOEnergie(28.77583e-6,unit_of_wavelength = u.m,return_quantities=False)
    ff = WavelengthTOEnergie(28.77583e-6,unit_of_wavelength = u.m,return_quantities=True)
    print(ee)
    print(ff)
    
    gg = EnergieTOWavelength(0.043086404716694296,return_quantities=False,unit_of_energy=u.eV)
    hh = EnergieTOWavelength(0.043086404716694296,return_quantities=True,unit_of_energy=u.eV)
    print(gg)
    print(hh)
    
    #test de la conversion globale:
    
    """frequency-to-wavelength"""
    Type = "frequency-to-wavelength"
    print(convertisseurPhysique(type_de_conversion=Type, frequence=10418.241360755546e9,unite_de_frequence=u.Hz,unite_de_longueur_d_onde=u.m, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, frequence=10418.241360755546e9,unite_de_frequence=u.Hz,unite_de_longueur_d_onde=u.m, retourner_qtes_physiques=True),
          NuTOWavelength(nu=10418.241360755546e9,unit_of_wavelength=u.m,unit_of_nu=u.Hz,return_quantities=True))
    
    """frequency-to-energy"""
    Type = "frequency-to-energy"
    print(convertisseurPhysique(type_de_conversion=Type, frequence=10418.241360755546e9,unite_de_frequence=u.Hz,unite_d_energie=u.erg, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, frequence=10418.241360755546e9,unite_de_frequence=u.Hz,unite_d_energie=u.erg, retourner_qtes_physiques=True),
          NuTOEnergie(nu=10418.241360755546e9,unit_of_energy=u.erg,unit_of_nu=u.Hz,return_quantities=True))
    
    """wavelength-to-frequency"""
    Type = "wavelength-to-frequency"
    print(convertisseurPhysique(type_de_conversion=Type, longueur_d_onde=145e-6, unite_de_longueur_d_onde=u.m, unite_de_frequence=u.Hz, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, longueur_d_onde=145e-6, unite_de_longueur_d_onde=u.m, unite_de_frequence=u.Hz, retourner_qtes_physiques=True),
          WavelengthTONu(wavelength=145e-6,unit_of_wavelength=u.m,unit_of_nu=u.Hz, return_quantities=True))
    
    """wavelength-to-energy"""
    Type = "wavelength-to-energy"
    print(convertisseurPhysique(type_de_conversion=Type, longueur_d_onde=145e-6, unite_de_longueur_d_onde=u.m, unite_d_energie=u.erg, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, longueur_d_onde=145e-6, unite_de_longueur_d_onde=u.m, unite_d_energie=u.erg, retourner_qtes_physiques=True),
          WavelengthTOEnergie(wavelength=145e-6,unit_of_wavelength=u.m,unit_of_energy=u.erg, return_quantities=True))
    
    """energy-to-frequency"""
    Type = "energy-to-frequency"
    print(convertisseurPhysique(type_de_conversion=Type, energie=0.043086404716694296, unite_d_energie=u.eV, unite_de_frequence=u.Hz, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, energie=0.043086404716694296, unite_d_energie=u.eV, unite_de_frequence=u.Hz, retourner_qtes_physiques=True),
          EnergieTONu(0.043086404716694296,return_quantities=True,unit_of_energy=u.eV,unit_of_nu=u.Hz))
    
    """energy-to-wavelength"""
    Type = "energy-to-wavelength"
    print(convertisseurPhysique(type_de_conversion=Type, energie=0.043086404716694296, unite_d_energie=u.eV, unite_de_longueur_d_onde=u.m, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, energie=0.043086404716694296, unite_d_energie=u.eV, unite_de_longueur_d_onde=u.m, retourner_qtes_physiques=True),
          EnergieTOWavelength(0.043086404716694296,return_quantities=True,unit_of_energy=u.eV,unit_of_wavelength=u.m))


    """Les fonctions :
        
        - NuTOWavelength, NuTOEnergie, WavelengthTONu, WavelengthTOEnergie, EnergieTONu, EnergieTOWavelength, 
        - convertisseurPhysique,
        
    ont ??t?? v??rifi??es comme valides et fonctionnant correctement par Mialy le 09/12/2019 ?? 03h00"""
    
    
    """surfaceBrightness-to-integratedIntensity"""
    Type = "surfaceBrightness-to-integratedIntensity"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, brillance_de_surface=1, unite_de_brillance_de_surface=u.watt/u.m**2/u.sr, longueur_d_onde=2, unite_de_longueur_d_onde=u.m, unite_d_intensite_integree=u.K * u.km/ u.s, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, brillance_de_surface=1, unite_de_brillance_de_surface=u.watt/u.m**2/u.sr, longueur_d_onde=2, unite_de_longueur_d_onde=u.m, unite_d_intensite_integree=u.K * u.km/ u.s, retourner_qtes_physiques=True),
          (3.62146e-2)*((1e7)/(1e4))*((2e6)**3))
    """surfaceFlux-to-integratedIntensity"""
    Type = "surfaceFlux-to-integratedIntensity"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, flux_de_surface=1e-6, unite_de_flux_surfacique=u.erg/u.s/u.cm**2, longueur_d_onde=145, unite_de_longueur_d_onde=u.millimeter, theta_source=10, unite_de_theta_source=u.radian, unite_d_intensite_integree=u.K * u.km/ u.s, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, flux_de_surface=1e-6, unite_de_flux_surfacique=u.erg/u.s/u.cm**2, longueur_d_onde=145, unite_de_longueur_d_onde=u.millimeter, theta_source=10, unite_de_theta_source=u.radian, unite_d_intensite_integree=u.K * u.km/ u.s, retourner_qtes_physiques=True),
          1.54076*1e-6*1e-6*((1e-7*1e18)/1e-4)*1e-9*((145*1e6)**3) / ((10*206265)**2))
    
    """integratedIntensity-to-surfaceFlux"""
    Type = "integratedIntensity-to-surfaceFlux"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, intensite_integree=1e7, unite_d_intensite_integree=u.K * u.m/u.s, theta_source=2, unite_de_theta_source=u.radian, longueur_d_onde=12, unite_de_longueur_d_onde=u.mm, unite_de_flux_surfacique=dixMoins18Watt/u.m**2, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, intensite_integree=1e7, unite_d_intensite_integree=u.K * u.m/u.s, theta_source=2, unite_de_theta_source=u.radian, longueur_d_onde=12, unite_de_longueur_d_onde=u.mm, unite_de_flux_surfacique=dixMoins18Watt/u.m**2, retourner_qtes_physiques=True),
          6.49028e5*1e7*1e-3*(((2*206265)**2)/((12*1e-3*1e6)**3)), "<=> h_unit / h_a_la_main = ", h/(6.49028e5*1e7*1e-3*(((2*206265)**2)/((12*1e-3*1e6)**3))))
    """surfaceBrightness-to-surfaceFlux"""
    Type = "surfaceBrightness-to-surfaceFlux"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, brillance_de_surface=1e7, unite_de_brillance_de_surface=u.watt/u.m**2/u.sr, theta_source=12, unite_de_theta_source=u.radian, unite_de_flux_surfacique=dixMoins18Watt/u.m**2, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, brillance_de_surface=1e7, unite_de_brillance_de_surface=u.watt/u.m**2/u.sr, theta_source=12, unite_de_theta_source=u.radian, unite_de_flux_surfacique=dixMoins18Watt/u.m**2, retourner_qtes_physiques=True),
          (2.35044e4)*1e7*(1e7/1e4)*(12*206265)**2)
    """surfaceFlux-to-surfaceBrightness"""
    Type = "surfaceFlux-to-surfaceBrightness"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, flux_de_surface=1e7, unite_de_flux_surfacique=u.erg/u.s/u.cm**2, theta_source=12, unite_de_theta_source=u.radian, unite_de_brillance_de_surface=u.erg/u.s/u.sr/u.cm**2, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, flux_de_surface=1e7, unite_de_flux_surfacique=u.erg/u.s/u.cm**2, theta_source=12, unite_de_theta_source=u.radian, unite_de_brillance_de_surface=u.erg/u.s/u.sr/u.cm**2, retourner_qtes_physiques=True),
          4.25452e-5*1e7*((1e-7*1e18)/(1e-4))*1/((12*206265)**2))
    """"integratedIntensity-to-surfaceBrightness"""
    Type = "integratedIntensity-to-surfaceBrightness"
    print('\n',
          convertisseurPhysique(type_de_conversion=Type, intensite_integree=1e7, unite_d_intensite_integree=u.K * u.m/u.s, longueur_d_onde=12, unite_de_longueur_d_onde=u.mm, unite_de_brillance_de_surface=u.erg/u.s/u.sr/u.cm**2, retourner_qtes_physiques=False),
          convertisseurPhysique(type_de_conversion=Type, intensite_integree=1e7, unite_d_intensite_integree=u.K * u.m/u.s, longueur_d_onde=12, unite_de_longueur_d_onde=u.mm, unite_de_brillance_de_surface=u.erg/u.s/u.sr/u.cm**2, retourner_qtes_physiques=True),
          27.6131*1e7*(1e-3)*1/((12*1e-3*1e6)**3))