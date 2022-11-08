# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define an initial configuration.
Grains are a combination of square and disk.
We have 2 temporary classes about grains and contact."""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from multiprocessing import Pool
from functools import partial
import Create_LG_IC_Square
import Create_LG_IC_Disk

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    #Create an initial configuration with tempo grains

    #simulation with squares and disks
    if dict_geometry['N_grain_disk']*dict_geometry['N_grain_square'] != 0:
        simulation_report.write('A sample composed of disks and square is not available for the moment !')
        raise ValueError('A sample composed of disks and square is not available for the moment !')

    #simulation with squares
    elif dict_geometry['N_grain_disk'] == 0:
        Create_LG_IC_Square.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #simulation with disk
    elif dict_geometry['N_grain_square'] == 0:
        Create_LG_IC_Disk.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample, simulation_report):
    #Convert an initial configuration with tempo grains to current configuration with real grain

    #simulation with squares and disks
    if dict_geometry['N_grain_disk']*dict_geometry['N_grain_square'] != 0:
        simulation_report.write('A sample composed of disks and square is not available for the moment !')
        raise ValueError('A sample composed of disks and square is not available for the moment !')

    #simulation with squares
    elif dict_geometry['N_grain_disk'] == 0:
        Create_LG_IC_Square.From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample)

    #simulation with disk
    elif dict_geometry['N_grain_square'] == 0:
        Create_LG_IC_Disk.From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample)
