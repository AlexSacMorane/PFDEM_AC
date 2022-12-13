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
from Create_LG_IC import LG_tempo, From_LG_tempo_to_usable
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

name_to_load = 'frac_5_run_1_save_dicts'
name_report = 'Debug/Report_after_crash'

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open(name_to_load,'rb')
dict_save = pickle.load(toload,encoding = 'bytes')
toload.close()
dict_algorithm = dict_save['algorithm']
dict_geometry = dict_save['geometry']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitations = dict_save['sollicitations']
dict_tracker = dict_save['tracker']

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

#create a simulation report
simulation_report = Report.Report(name_report,datetime.now())

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

while not User.Criteria_StopSimulation(dict_algorithm):

    main.main_iteration(dict_algorithm, dict_geometry, dict_material, dict_sollicitations, dict_sample, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_simulation(dict_algorithm, dict_tracker, simulation_report)
