# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file to restart a simulation after a crash.
There is a save at the end of each PFDEM iteration.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import os
import shutil
from datetime import datetime
from pathlib import Path
import pickle

#Own function and class
from Write_txt import Write_txt
from Create_i_AC import Create_i_AC_local
import Owntools
import Grain
import Contact
import Contact_gw
import Report
import User
import main

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

name_to_load = 'LG_frac_15_run_2_save_dicts_before_pf'

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open(name_to_load,'rb')
dict_save = pickle.load(toload)
toload.close()
dict_algorithm = dict_save['algorithm']
dict_geometry = dict_save['geometry']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitations = dict_save['sollicitations']
dict_tracker = dict_save['tracker']
simulation_report = dict_save['report']

#-------------------------------------------------------------------------------
#plan simulation
#-------------------------------------------------------------------------------

simulation_report.write('\nA crash occurs...\n\n')
dict_algorithm['Debug_DEM'] = True
dict_algorithm['i_print_plot'] = 50


#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

if name_to_load[-10:] =='_before_pf':
    main.main_iteration_from_pf(dict_algorithm, dict_geometry, dict_material, dict_sollicitations, dict_sample, dict_tracker, simulation_report)

while not User.Criteria_StopSimulation(dict_algorithm):

    main.main_iteration_until_pf(dict_algorithm, dict_geometry, dict_material, dict_sollicitations, dict_sample, dict_tracker, simulation_report)
    main.main_iteration_from_pf(dict_algorithm, dict_geometry, dict_material, dict_sollicitations, dict_sample, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_simulation(dict_algorithm, dict_tracker, simulation_report)
