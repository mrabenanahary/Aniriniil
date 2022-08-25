#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import sys
import scipy as scipy
import matplotlib.pyplot as plt
from matplotlib import pyplot
import matplotlib as mpl
from math import pi
from scipy.special import gammainc
from scipy.special import gamma
from scipy.interpolate import griddata

#############################
# constants in cgs
#############################

c = 3.0e+10

m_e = 9.1e-28

m_p = 1.67e-24

sigma_T = 6.65e-25

e = 4.803e-10

#############################
# constants for g/h emission function
#############################

c_1 = 0.78

c_2 = 0.25

c_3 = 2.175


#######################################################################################
# GET SYNCHROTRON EMISSION & ABSORPTION
#######################################################################################

def synchrotron(
    data,
    filename,
    fileid,
    theta,
    npixel,
    nu,
    alpha,
    frac_L,
    ):

    # ########################
    # constants
    # ########################

    Rnorm = 3.085677e+17

    L_tot = 1e+46

    lfac_in = 10.

    eta_rho = 1e-3

    gamma_ad = 5. / 3.

    epsilon = 0.1

    c_e = 1e+3

    energy_ratio = 0.1

    index = 2. * alpha + 1.

    # ########################
    # extracting points
    # ########################

    points = data.get('points')

    x = points[:, 0]

    y = points[:, 1]

    z = points[:, 2]

    # ########################
    # extracting tracers
    # ########################

    tracer_1 = data.get('flrho1')

    tracer_2 = data.get('flrho2')

    tracer_3 = data.get('flrho3')

    # ########################
    # extracting magnetic field
    # ########################

    if filename == 'toro' or filename == 'polo' or filename == 'helico':

        magnetic_field_x = data.get('b1')

        magnetic_field_y = data.get('b2')

        magnetic_field_z = data.get('b3')

        magnetic_field = np.sqrt(magnetic_field_x ** 2.
                                 + magnetic_field_y ** 2.
                                 + magnetic_field_z ** 2.)
    else:

        magnetic_field = data.get('bturb')

    # ########################
    # extracting velocities
    # ########################

    velocity_x = data.get('v1')

    velocity_y = data.get('v2')

    velocity_z = data.get('v3')

    velocity = np.sqrt(velocity_x ** 2 + velocity_y ** 2. + velocity_z
                       ** 2.)

    # ########################
    # extracting pressure
    # ########################

    pressure = data.get('p')

    # ########################
    # extracting density
    # ########################

    rho = data.get('rho')

    # ########################
    # extracting lorentz factor
    # ########################

    lfac = data.get('lfac')

    # ########################
    # calculating entalpy
    # ########################

    h = 1. + gamma_ad / (gamma_ad - 1.) * (pressure / (rho * c ** 2.))

    # ########################
    # getting normalization
    # ########################

    velocity_norm = c

    rho_norm = L_tot * frac_L / (Rnorm ** 2. * c ** 3. * pi
                                 * np.sqrt(lfac_in ** 2. - 1.)
                                 * (lfac_in * h - 1.)) * eta_rho ** -1.

    density_part_norm = rho_norm / m_p

    pressure_norm = m_p * density_part_norm * c ** 2.

    magnetic_field_norm = np.sqrt(8. * pi * pressure_norm)

    del h

    # ########################
    # applying normalization
    # ########################

    rho_phys = rho_norm * rho

    nelem = len(rho_phys)

    density_part_phys = density_part_norm * rho

    pressure_phys = pressure_norm * pressure

    magnetic_field_phys = magnetic_field_norm * magnetic_field

    velocity_phys = velocity_norm * velocity

    velocity_x_phys = velocity_norm * velocity_x

    velocity_y_phys = velocity_norm * velocity_y

    velocity_z_phys = velocity_norm * velocity_z

    del rho_norm
    del rho
    del density_part_norm
    del pressure_norm
    del pressure
    del magnetic_field_norm
    del magnetic_field
    del velocity_norm
    del velocity_x
    del velocity_y
    del velocity_z

    # ########################
    # neglects the AM
    # ########################

    tracer_index = np.where((tracer_1 <= 0.) & (tracer_2 < 1e-5)
                            & (tracer_3 < 1.))

    pressure_phys[tracer_index] = 0.

    energy_intern = pressure_phys / (gamma_ad - 1.)

    ei_index = np.where(energy_intern < 0.)

    energy_intern[ei_index] = 0.

    del tracer_1
    del tracer_2
    del tracer_3

    # ########################
    # calculating magnetic energy
    # ########################

    magnetic_energy = magnetic_field_phys ** 2. / (8. * pi)

    # ########################
    # calculating larmor frequency
    # ########################

    nu_L = e * magnetic_field_phys / (2. * pi * m_e * c)

    del magnetic_field_phys

    # ########################
    # table initialization
    # ########################

    emission = np.zeros(nelem)

    absorption = np.zeros(nelem)

    lfac_min = np.zeros(nelem)

    lfac_max = np.zeros(nelem)

    t_min = np.zeros(nelem)

    t_max = np.zeros(nelem)

    K = np.zeros(nelem)

    # ########################
    # calculating lfac/t min and lfac/t absorption_max
    # ########################

    index_em = np.where((nu_L > 0.) & (pressure_phys > 0.))

    lfac_min[index_em] = energy_ratio * energy_intern[index_em] \
        / (epsilon * density_part_phys[index_em]) * ((index - 2.)
            / (index - 1.)) * ((1. - c_e ** (1. - index)) / (1. - c_e
                               ** (2. - index)))

    lfac_min[np.where(lfac_min < 1.)] = 1.

    lfac_max[index_em] = c_e * lfac_min[index_em]

    t_min[index_em] = nu / (3. * lfac_min[index_em] ** 2.
                            * nu_L[index_em])

    t_max[index_em] = nu / (3. * lfac_max[index_em] ** 2.
                            * nu_L[index_em])

    # ########################
    # calculating normalization factor
    # ########################

    K[index_em] = (energy_ratio * energy_intern[index_em] * (index
                   - 2.) / (1. - c_e ** (2. - index))) ** (index - 1.) \
        * ((1. - c_e ** (1. - index)) / (epsilon
           * density_part_phys[index_em] * (index - 1.))) ** (index
            - 2.)

    # ########################
    # calculating emission & absorption
    # ########################

    emission[index_em] = 9. * sigma_T * c * magnetic_energy[index_em] \
        * c_1 / (24. * pi ** 2. * nu_L[index_em]) * np.sqrt(nu
            / nu_L[index_em]) * (g_emission(nu, index, K[index_em],
                                 nu_L[index_em], t_max[index_em])
                                 - g_emission(nu, index, K[index_em],
                                 nu_L[index_em], t_min[index_em]))

    absorption[index_em] = 3. * np.sqrt(3.) * sigma_T * c \
        * magnetic_energy[index_em] * c_1 / (16. * pi ** 2. * m_e * nu
            ** 2. * nu_L[index_em]) * (h_absorption(nu, index,
            K[index_em], nu_L[index_em], t_max[index_em])
            - h_absorption(nu, index, K[index_em], nu_L[index_em],
            t_min[index_em]))

    # ########################
    # del old variables
    # ########################

    del lfac_min
    del lfac_max
    del t_min
    del t_max
    del K
    del pressure_phys
    del nu_L
    del energy_intern
    del magnetic_energy
    del density_part_phys

    # ########################
    # applying doppler factor (simple way)
    # ########################

    print('emission (max) =', emission.max(), 'cgs')
    print('absorption (max) =', absorption.max(), 'cgs')

    before_boost = emission.max()

    beta = np.sqrt(velocity_x_phys ** 2. + velocity_y_phys ** 2
                   + velocity_z_phys ** 2.) / c

    lfac[np.where(lfac < 1.)] = 1.

    doppler = 1. / (lfac - np.sqrt(lfac ** 2. - 1.) * np.cos((90.
                    - theta) * pi / 180.))

    emission = emission * doppler ** 2.

    absorption = absorption * doppler ** -1.

    after_boost = emission.max()

    print('Boost ratio on emission =', after_boost / before_boost)

    print('emission (max) =', emission.max(), 'cgs')

    print('absorption (max) =', absorption.max(), 'cgs')

    del beta
    del lfac
    del velocity_x_phys
    del velocity_y_phys
    del velocity_z_phys
    del doppler

    # ########################
    # return
    # #######################

    return {'emission': emission, 'absorption': absorption,
            'points': data.get('points')}


def shotgun(
    inputdata,
    npixel,
    redshift,
    index_start,
    index_end,
    delta_file,
    surface_ref,
    number_file,
    file_num,
    previous_I,
    phi,
    delta,
    theta,
    ):

    # ########################
    # extracting
    # ########################

    points     = inputdata.get('points'    )

    emission   = inputdata.get('emission'  )

    absorption = inputdata.get('absorption')

    # ########################
    # constants
    # #######################

    parsec = 3.086e+17

    x0     = [0., 0., 0.]

    costheta = np.cos(theta * pi / 180.)

    sintheta = np.sin(theta * pi / 180.)

    # ########################
    # make interpolation grid
    # ########################

    x = points[:, 0]

    y = points[:, 1]

    z = points[:, 2]

    L = [max(x), max(y)*costheta + 2.*max(z)*sintheta, max(z)*costheta]
    
    print('Ly max:', L[0])
    print('Lx max:', L[1])
    print('Lz max:', L[2])

    grid = make_grid(
        phi,
        theta,
        npixel,
        L,
        delta=delta,
        x0=x0,
        )

    # ########################
    # calculating typical distance / surface
    # ########################
    
    print('et y :', np.var(y))
    print('et x :', np.var(x))
    print('et z :', np.var(z))

    delta_los = np.var(y) ** 0.5 * parsec * L[1] / npixel[1]

    print('delta_los =', delta_los, 'cm')

    surface_los = np.var(x) ** 0.5 * np.var(z) ** 0.5 * parsec ** 2. \
        * (L[0] * L[2]) / (npixel[0] * npixel[2])

    print('surface_los :', surface_los, 'cm^2')

    H_0 = 70.

    distance_lum = 2. * c * (redshift + 1. - np.sqrt(redshift + 1.)) \
        / (H_0 * 3.24076e-20)

    print('distance lumineuse =', distance_lum, 'cm')

    # ########################
    # emission interpolation
    # ########################

    emission = griddata((x, y, z), emission, (grid[0], grid[1],
                        grid[2]), method='nearest')

    emission = np.reshape(emission, (npixel[0], npixel[1], npixel[2]))

    print('emission interpolation DONE')
    
    plt.imshow(emission.sum(axis = 1), aspect = 'auto')
    plt.show()

    # ########################
    # absorption interpolation
    # ########################

    absorption = griddata((x, y, z), absorption, (grid[0], grid[1],
                          grid[2]), method='nearest')

    absorption = np.reshape(absorption, (npixel[0], npixel[1],
                            npixel[2]))

    print('absorption interpolation DONE')

    # ########################
    # opacity
    # ########################

    tau = absorption * delta_los

    # ########################
    # integrating between surface
    # ########################

    intensity = np.zeros((npixel[0], npixel[2]))

    #if number_file != file_num:

    #    intensity = previous_I

    print ('reference surface : number =', surface_ref)

    surface_start = index_start + int(delta_file * surface_ref)

    surface_end = index_start + int(delta_file * (1. + surface_ref))

    print ('index start is', index_start)

    print ('index end is', index_end)

    print ('surface start is', surface_start)

    print ('surface end is', surface_end)

    total_intensity = 0.

    for i in range(npixel[0]):
        for j in range(npixel[2]):
            for k in range(surface_start, surface_end):

                if absorption[i, k, j] > 0.:

                        intensity[i, j] = intensity[i, j] * np.exp(-tau[i,k, j]) + emission[i, k, j] / absorption[i,k, j] * (1. - np.exp(-tau[i, k, j]))
                
                    #intensity[i,j] = intensity[i,j] + emission[i,k,j] * delta_los

    # ########################
    # del and return
    # ########################

    del emission
    del absorption
    del tau

    plt.imshow(intensity, aspect = 'auto', cmap = 'jet')
    plt.colorbar()
    plt.show()

    return {'intensity': intensity, 'surface_los': surface_los,
            'distance_lum': distance_lum}


#######################################################################################
# grid function
#######################################################################################


def make_grid(phi,theta,npixel,L,delta=0,x0=[0,0,0]):
    
    xi,yi,zi=np.mgrid[0:npixel[0],0:npixel[1],0:npixel[2]]
    
    
    xi=xi/(npixel[0]-1.)
    yi=yi/(npixel[1]-1.)
    zi=zi/(npixel[2]-1.)

    
    xi = xi*L[0] - L[0]/2.
    yi = yi*L[1] - L[1]/2.
    zi = zi*L[2] + 1.#+ L[2]/2. 
    
    
    xi = xi.flatten()
    yi = yi.flatten()
    zi = zi.flatten()
    
    nelem = len(xi)

# Since the standard phi counts from x (and not -y as we do by rotating around x first)
# add 90 to the internal phi.  
    phii = phi + 90.
# Now rotate the points around theta and then phi:

    cosphi = np.cos(phii*np.pi/180.)
    sinphi = np.sin(phii*np.pi/180.)
    costheta = np.cos(theta*np.pi/180.)
    sintheta = np.sin(theta*np.pi/180.)
    
    xid = xi
    yid = np.empty(nelem)
    zid = np.empty(nelem)
# Rotate by theta around the x-axis:
    #yid = costheta*yi - sintheta*zi
    #zid = sintheta*yi + costheta*zi
    yid  = costheta*yi - sintheta*zi
    zid  = sintheta*yi + costheta*zi

    xidd = np.empty(nelem)
    yidd = np.empty(nelem)
    zidd = zid
# Rotate by phi around the z-axis:
    xidd =  cosphi*xid - sinphi*yid
    yidd =  sinphi*xid + cosphi*yid

    return [xidd+x0[0],yidd+x0[1],zidd+x0[2]]

#######################################################################################
# emission functions
#######################################################################################

def g_emission(
    nu,
    index,
    K,
    nu_L,
    t,
    ):

    g_function = K * (nu / (3. * nu_L)) ** (-index / 2.) * t ** (c_2
            + (index - 1.) / 2.) * (c_3 * t) ** ((1. - index) / 2.
            - c_2) * gamma(c_2 + (index - 1.) / 2.) * (1.
            - gammainc(c_2 + (index - 1.) / 2., c_3 * t))

    return g_function


def h_absorption(
    nu,
    index,
    K,
    nu_L,
    t,
    ):

    h_function = K * (index + 2.) * 3. ** (index / 2.) * (nu / nu_L) \
        ** (-index / 2.) * t ** (c_2 + index / 2.) * (c_3 * t) ** (-c_2
            - index / 2.) * gamma(c_2 + index / 2.) * (1.
            - gammainc(c_2 + index / 2., c_3 * t))

    return h_function











