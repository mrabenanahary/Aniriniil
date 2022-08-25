#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import syncmap
import read
from math import pi
import rotation_data
import matplotlib.pyplot as plt
from pprint import pprint
import numpy as np
import time
start_time = time.time()

#############################
# angle l.o.s
#############################

theta = 45.

phi = 0

nphi = 10

delta = 0

costheta = np.cos(theta * pi / 180.)

sintheta = np.sin(theta * pi / 180.)

#############################
# size and resolution of the box
#############################

npixel = [125, int(125*costheta + 500*sintheta),500]

# npixel = [120, 600, 200]

x0 = [0, 0, 0]

#############################
# constants
#############################

alpha      = 0.6

redshift   = 0.00428

frac_L     = 1e-3

nu         = 1e+9

erg_to_mJy = 1e+26

R_jet      = 3.086e+17

#############################
# select snapshot
#############################

fileid = 0

filename = 'hydro'

fileroot = 'hydro_open_2E_'

root = '../../PHD_1/CODE/Code_2020_09/PP/DATA/'

#root = ''

filename_id = root + 'jet2D_srmhd_' + fileroot

#############################
# LTD file
#############################

delta_t_obs = 10.  

gamma_mean  = (10. * 1. + 3 * 2.)/3.

print('gamma mean in jet :', gamma_mean)

number_max  = 1

print(number_max, 'point in light curve')

I_boxed     = np.zeros((npixel[0], npixel[2]))

index_start = 0

index_end   = npixel[1]

if theta == 0. : 

    delta_y     = index_end 

    number_file = index_end / delta_y

else : 

    costheta    = np.cos(theta * pi / 180.)
    
   # delta_y     = (delta_t_obs / gamma_mean) * (costheta / ( 1. - costheta)) 

    delta_y = index_end 

    number_file =  int(1. / delta_y)

    if number_file == 0. : 

        number_file = 1
    
    delta_y     = int(index_end / number_file)

#number_file = int((index_end - index_start) / delta_y)

print ('delta_y =', delta_y)

print ('one point =', number_file, 'files !')

#############################
# start integrating
#############################

for number in range(number_max):

    index_surface = 0

    file_num = number

    number_file = number + 1

    # ########################
    # calculating one point
    # ########################

    for i in range(number, number_file):

        print('opening and rotation of the file')

        data_3D = rotation_data.rot2D(filename, i, nphi,
                filename_id, type='vtu')

        print('calculating emission, absorption and intensity')

        emission = syncmap.synchrotron(
            data_3D,
            filename,
            fileid,
            theta,
            npixel,
            nu,
            alpha,
            frac_L,
            )

        intensity = syncmap.shotgun(
            emission,
            npixel,
            redshift,
            index_start,
            index_end,
            delta_y,
            index_surface,
            number_file,
            file_num,
            I_boxed,
            phi,
            delta,
            theta,
            )

        I_boxed = intensity.get('intensity')

        # ####################
        # export constants
        # ####################

        surface_los = intensity.get('surface_los')

        distance_lum = intensity.get('distance_lum')

        # ####################
        # looking for the next file
        # ####################

        index_surface = index_surface + 1

    # ########################
    # sum up
    # ########################

    total_intensity = 0.

    flux = surface_los / distance_lum ** 2. * (1. + redshift) * I_boxed

    total_intensity = erg_to_mJy * flux.sum()

    print('total intensity =', total_intensity, 'mJy')
    
    np.savetxt(filename + '_' + str(number) + '.txt', flux)

    # ########################
    # next point
    # ########################

    number_file = number_file + 1

print('--- %s seconds ---' % (time.time() - start_time))
