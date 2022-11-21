# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define an initial configuration.
Grains are a combination of square and disk."""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import Create_LG_IC_Disk
import Create_LG_IC_Square
import Create_LG_IC_Mix

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    '''Create an initial configuration with tempo grains'''

    #simulation with disk
    if dict_geometry['type'] == 'AllDisk':
        Create_LG_IC_Disk.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #simulation with squares
    elif dict_geometry['type'] == 'AllSquare':
        Create_LG_IC_Square.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #simulation with squares and disks
    elif dict_geometry['type'] == 'Mix':
        Create_LG_IC_Mix.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample, simulation_report):
    '''Convert an initial configuration with tempo grains to current configuration with real grain'''

    #simulation with disk
    if dict_geometry['type'] == 'AllDisk':
        Create_LG_IC_Disk.From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample)

    #simulation with squares
    elif dict_geometry['type'] == 'AllSquare':
        Create_LG_IC_Square.From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample)

    #simulation with squares and disks
    elif dict_geometry['type'] == 'Mix':
        Create_LG_IC_Mix.From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample)
