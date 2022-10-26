# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
import os
import shutil
from datetime import datetime
from pathlib import Path
import random
from multiprocessing import Pool
from functools import partial
import time

#Own function and class
from Write_txt import Write_txt
from PFtoDEM_Multi import PFtoDEM_Multi
from Create_i_AC import Create_i_AC
from Create_LG_IC import LG_tempo, From_LG_tempo_to_usable_f
import Owntools
import Grain
import Etai
import Contact
import Contact_gw
import Report
import User

#-------------------------------------------------------------------------------
#Plan simulation
#-------------------------------------------------------------------------------

if Path('Debug').exists():
    shutil.rmtree('Debug')
if Path('Input').exists():
    shutil.rmtree('Input')
if Path('Output').exists():
    shutil.rmtree('Output')
if Path('Data').exists():
    shutil.rmtree('Data')

os.mkdir('Debug')
os.mkdir('Debug/DEM_ite')
os.mkdir('Debug/DEM_ite/Init')
os.mkdir('Input')
os.mkdir('Output')
os.mkdir('Data')

#-------------------------------------------------------------------------------
# tic
#-------------------------------------------------------------------------------

simulation_report = Report.Report('Debug/Report',datetime.now())

#-------------------------------------------------------------------------------
#Initial conditions
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

if dict_algorithm['Debug'] or dict_algorithm['Debug_DEM'] :
    simulation_report.write('This simulation can be debugged\n')
if dict_algorithm['SaveData'] :
    if not Path('../'+dict_algorithm['main_folder_name']).exists():
        os.mkdir('../'+dict_algorithm['main_folder_name'])
    simulation_report.write('This simulation is saved\n')
    i_run = 1
    folderpath = Path('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['template_simulation_name']+str(i_run))
    while folderpath.exists():
        i_run = i_run + 1
        folderpath = Path('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['template_simulation_name']+str(i_run))

    #add element in dict
    dict_algorithm['name_folder'] = dict_algorithm['template_simulation_name']+str(i_run)

if dict_algorithm['SaveData'] or dict_algorithm['Debug'] or dict_algorithm['Debug_DEM']:
    simulation_report.write('\n')

#Creation of the grain (without PF)
simulation_report.write_and_print('Creation of the grains\n','\nCREATION OF THE GRAINS\n')
LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#Spatial discretisation
User.Add_SpatialDiscretisation(dict_geometry,dict_sample)

#PF parameters
User.Add_WidthInt_DoubleWellBarrier(dict_material, dict_sample)

#Dissolution energy
User.Add_DissolutionEnergy(dict_algorithm,dict_geometry,dict_material,dict_sollicitations)

simulation_report.tac_tempo(datetime.now(),'Initialisation')

#-------------------------------------------------------------------------------
#Creation of polygonal particles
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

# Creation of the real list of grains
From_LG_tempo_to_usable_f(dict_algorithm, dict_ic, dict_geometry, dict_material, dict_sample, simulation_report)

#creation pf the dissolution .txt
Owntools.Write_e_dissolution_txt(dict_sample,dict_sollicitations)

#Study the grain and saving the grain surface
Owntools.Study_Polygonal_and_init_f(dict_algorithm,dict_geometry,dict_sample,dict_sollicitations,simulation_report)


S_grains = 0
for grain in dict_sample['L_g']:
    S_grains = S_grains + grain.surface
simulation_report.write('Total Surface '+str(round(S_grains,0))+' µm2\n')

simulation_report.tac_tempo(datetime.now(),'Creation and study of polygonal particles')

#-------------------------------------------------------------------------------
#Distribution etai
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

Etai.etai_distribution_dissolution(dict_algorithm,dict_sample, dict_sollicitations, simulation_report)

#create eta for undissolved
L_etai_undissolved = []
for etai in range(0,len(dict_sample['L_ig_etai_undissolved'])):
    L_etai_undissolved.append(Etai.Etai(etai,dict_sample['L_ig_etai_undissolved'][etai]))
#create eta for dissolved
L_etai_dissolved = []
for etai in range(0,len(dict_sample['L_ig_etai_dissolved'])):
    L_etai_dissolved.append(Etai.Etai(len(dict_sample['L_ig_etai_undissolved'])+etai,dict_sample['L_ig_etai_dissolved'][etai]))

#Add elements in dict
dict_sample['L_etai_undissolved'] = L_etai_undissolved
dict_sample['L_etai_dissolved'] = L_etai_dissolved

#Saving the grain surface
S_grains_dissolvable = 0
for grain in dict_sample['L_g']:
    if grain.dissolved :
        S_grains_dissolvable = S_grains_dissolvable + grain.surface
simulation_report.write('Total Surface dissolvable '+str(round(S_grains_dissolvable,0))+' µm2\n')

#Plot etai distribution
if dict_algorithm['Debug'] :
    Owntools.Debug_etai_f(dict_sample)

simulation_report.tac_tempo(datetime.now(),'Etai distribution')

#-------------------------------------------------------------------------------
#Main
#-------------------------------------------------------------------------------

#Tracker and create an new dict
dict_tracker = {
    't_L' : [0],
    'S_grains_L' : [S_grains],
    'S_grains_dissolvable_L' : [S_grains_dissolvable],
    'S_dissolved_L' : [0],
    'S_dissolved_perc_L' : [0],
    'S_dissolved_perc_dissolvable_L' : [0],
    'n_grains_L' : [len(dict_sample['L_g'])],
    'k0_xmin_L' : [],
    'k0_xmax_L' : []
}

# Preparation and add elements in dicts
dict_algorithm['i_PF'] = 0
dict_sample['L_contact_gw'] = []
dict_sample['L_ij_contact_gw'] = []
dict_sample['id_contact_gw'] = 0
dict_sample['L_contact'] = []
dict_sample['L_ij_contact'] = []
dict_sample['id_contact'] = 0

if dict_algorithm['Debug'] :
    Owntools.Debug_f2(dict_algorithm,dict_sample)

while not User.Criteria_StopSimulation(dict_algorithm):
      # update element in dict
      dict_algorithm['i_PF'] = dict_algorithm['i_PF'] + 1

      simulation_report.write_and_print('\nIteration '+str(dict_algorithm['i_PF'])+' / '+str(dict_algorithm['n_t_PFDEM'])+'\n','\nITERATION PF '+str(dict_algorithm['i_PF'])+' / '+str(dict_algorithm['n_t_PFDEM'])+'\n')

      #prepare iteration
      os.mkdir('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF']))
      os.mkdir('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/txt')
      os.mkdir('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/png')

      if dict_algorithm['MovePF_selector'] == 'Interpolation':
          # Saving to compute a rigid body motion
          L_center_g = []
          for grain in dict_sample['L_g']:
              L_center_g.append(grain.center.copy())

      # Compute kinetic energy criteria
      Ecin_stop = 0
      for grain in dict_sample['L_g']:
          Ecin_stop = Ecin_stop + 0.5*grain.m*(dict_algorithm['Ecin_ratio']*grain.r_mean/dict_algorithm['dt_DEM'])**2/len(dict_sample['L_g'])

      #Update element in dict
      dict_algorithm['Ecin_stop'] = Ecin_stop

      #Trackers and add element in dict
      dict_tracker['Ecin'] = []
      dict_tracker['Force_applied'] = []
      dict_tracker['k0_xmin'] = []
      dict_tracker['k0_xmax'] = []
      dict_tracker['y_box_max'] = [dict_sample['y_box_max']]
      dict_tracker['Force_on_upper_wall'] = []

      dict_algorithm['i_DEM'] = - 1
      DEM_loop_statut = True
      simulation_report.write('\n')
      simulation_report.tic_tempo(datetime.now())

      #-----------------------------------------------------------------------------
      # DEM iteration
      #-----------------------------------------------------------------------------

      while DEM_loop_statut :
          # update element in dict
          dict_algorithm['i_DEM'] = dict_algorithm['i_DEM'] + 1

          for grain in dict_sample['L_g']:
              grain.init_f_control(dict_sollicitations)

          # Detection of contacts between grains
          if dict_algorithm['i_DEM'] % dict_algorithm['i_update_neighborhoods']  == 0:
              Contact.Update_Neighborhoods_f(dict_algorithm,dict_sample)
          Contact.Grains_Polyhedral_contact_Neighborhoods_f(dict_algorithm,dict_material,dict_sample)

          # Detection of contacts between grain and walls
          if dict_algorithm['i_DEM'] % dict_algorithm['i_update_neighborhoods']  == 0:
              Contact_gw.Update_wall_Neighbouroods(dict_algorithm, dict_sample)
          Contact_gw.Grains_Polyhedral_Wall_contact_Neighbourood(dict_material,dict_sample)

          #Compute contact interactions (g-g and g-w)
          for contact in dict_sample['L_contact']:
              if dict_algorithm['Spring_type'] == 'Ponctual':
                  contact.DEM_2grains_Polyhedral_normal()
                  contact.DEM_2grains_Polyhedral_tangential(dict_algorithm['dt_DEM'])
              elif dict_algorithm['Spring_type'] == 'Surface':
                  contact.DEM_2grains_Polyhedral_normal_surface()
                  contact.DEM_2grains_Polyhedral_tangential_surface(dict_algorithm['dt_DEM'])
              else :
                  simulation_report.write('Spring type not available !')
                  raise ValueError('Spring type not available !')
          for contact in dict_sample['L_contact_gw'] :
              if dict_algorithm['Spring_type'] == 'Ponctual':
                  contact.DEM_gw_Polyhedral_normal()
                  contact.DEM_gw_Polyhedral_tangential(dict_algorithm['dt_DEM'])
              else : #Surface must be coded for contact gw
                  simulation_report.write('Spring type not available !')
                  raise ValueError('Spring type not available !')

          #Move particles and trackers
          #Semi implicit euler scheme
          Ecin = 0
          Force_applied = 0
          for grain in dict_sample['L_g']:
              a_i = grain.f/grain.m
              v_i = grain.v + a_i*dict_algorithm['dt_DEM']
              dw_i = grain.mz/grain.inertia
              w_i = grain.w + dw_i*dict_algorithm['dt_DEM']
              grain.update_geometry_kinetic(v_i,a_i,w_i,dict_algorithm['dt_DEM']) #Move grains
              Ecin = Ecin + 0.5*grain.m*np.linalg.norm(grain.v)**2/len(dict_sample['L_g'])
              Force_applied = Force_applied + np.linalg.norm(grain.f)/len(dict_sample['L_g'])

          #Control the y_max to verify vertical confinement
          Owntools.Control_y_max_NR(dict_sample,dict_sollicitations)
          Owntools.Compute_k0(dict_sample,dict_sollicitations)

          #trackers
          dict_tracker['Ecin'].append(Ecin)
          dict_tracker['Force_applied'].append(Force_applied)
          dict_tracker['y_box_max'].append(dict_sample['y_box_max'])
          dict_tracker['Force_on_upper_wall'].append(dict_sollicitations['Force_on_upper_wall'])
          dict_tracker['k0_xmin'].append(dict_sample['k0_xmin'])
          dict_tracker['k0_xmax'].append(dict_sample['k0_xmax'])

          if dict_algorithm['i_DEM'] %dict_algorithm['i_print_plot'] == 0:
              print('\nPF '+str(dict_algorithm['i_PF'])+' -> i_DEM '+str(dict_algorithm['i_DEM']+1)+' / '+str(dict_algorithm['i_DEM_stop']+1)+' (max)')
              print('Ecin',int(Ecin),'/',int(dict_algorithm['Ecin_stop']),'('+str(int(100*Ecin/dict_algorithm['Ecin_stop'])),' %)')
              print('F_confinement',int(dict_sollicitations['Force_on_upper_wall']),'/',int(dict_sollicitations['Vertical_Confinement_Force']),'('+str(int(100*dict_sollicitations['Force_on_upper_wall']/dict_sollicitations['Vertical_Confinement_Force'])),' %)')

              Owntools.save_DEM_tempo(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker)

              if dict_algorithm['Debug_DEM'] :
                Owntools.Debug_DEM_f(dict_algorithm, dict_sample)
                Write_txt(dict_algorithm,dict_sample)

          #-----------------------------------------------------------------------------
          # Stop conditions
          #-----------------------------------------------------------------------------

          if dict_algorithm['i_DEM'] >= dict_algorithm['i_DEM_stop'] :
              DEM_loop_statut = False
              print("DEM loop stopped by too many iterations.")
              simulation_report.write('/!\ End of DEM steps with '+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+'/!\ \n')
          if Ecin < dict_algorithm['Ecin_stop'] and dict_algorithm['i_DEM'] > dict_algorithm['n_window_stop'] and (dict_sollicitations['Vertical_Confinement_Force']*0.95<dict_sollicitations['Force_on_upper_wall'] and dict_sollicitations['Force_on_upper_wall']<dict_sollicitations['Vertical_Confinement_Force']*1.05):
              k0_xmin_window = dict_tracker['k0_xmin'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              k0_xmax_window = dict_tracker['k0_xmax'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              y_box_max_window = dict_tracker['y_box_max'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              if max(k0_xmin_window) - min(k0_xmin_window) < dict_algorithm['dk0_stop'] and max(k0_xmax_window) - min(k0_xmax_window) < dict_algorithm['dk0_stop'] and max(y_box_max_window) - min(y_box_max_window) < dict_algorithm['dy_box_max_stop']:
                  DEM_loop_statut = False
                  print("DEM loop stopped by steady state reached.")
                  simulation_report.write("DEM loop stopped by steady state reached with "+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+"\n")

      #-----------------------------------------------------------------------------
      # Debugging at the end of DEM step
      #-----------------------------------------------------------------------------

      if dict_algorithm['Debug'] :
        Owntools.Debug_f(dict_algorithm,dict_sample)
        Owntools.Debug_Trackers(dict_algorithm,dict_sollicitations,dict_tracker)
        Write_txt(dict_algorithm,dict_sample)
        Owntools.Plot_chain_force(dict_algorithm['i_PF'],dict_algorithm['i_DEM'])

        Owntools.save_DEM_final(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker)

      #-----------------------------------------------------------------------------
      # Compute Vertical and horizontal sollicitations to compute k0
      #-----------------------------------------------------------------------------

      Owntools.Control_y_max_NR(dict_sample,dict_sollicitations)
      Owntools.Compute_k0(dict_sample,dict_sollicitations)

      simulation_report.write('k0_xmin : '+str(round(dict_sample['k0_xmin'],2))+' / k0_xmax : '+str(round(dict_sample['k0_xmax'],2))+'\n')

      #Update element in dict
      dict_tracker['k0_xmin_L'].append(dict_sample['k0_xmin'])
      dict_tracker['k0_xmax_L'].append(dict_sample['k0_xmax'])

      simulation_report.tac_tempo(datetime.now(),'DEM loop '+str(dict_algorithm['i_PF']))

      #-----------------------------------------------------------------------------
      # Move PF
      #-----------------------------------------------------------------------------

      simulation_report.tic_tempo(datetime.now())

      if dict_algorithm['MovePF_selector'] == 'DeconstructRebuild':
          Grain.DEMtoPF_Decons_rebuild_f(dict_algorithm,dict_material,dict_sample)

          for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
              etai.update_etai_M(dict_sample['L_g'])
              etai.Write_txt_Decons_rebuild(dict_algorithm,dict_sample)
      #does not work for multiple grain by etai
      #elif MovePF_selector == 'Interpolation':
      #  for id_grain in range(len(dict_sample['L_g'])):
      #      rbm = dict_sample['L_g'][id_grain].center - L_center_g[id_grain]
      #      grain = dict_sample['L_g'][id_grain]
      #      grain.DEMtoPF_Interpolation(x_L,y_L,i_PF,rbm)
      else :
          simulation_report.write('Method to move phase field not available !')
          raise ValueError('Method to move phase field not available !')

      simulation_report.tac_tempo(datetime.now(),'Move phase-field iteration '+str(dict_algorithm['i_PF']))

      #-----------------------------------------------------------------------------
      # PF Simulation
      #-----------------------------------------------------------------------------

      simulation_report.tic_tempo(datetime.now())

      #Create the .i for moose
      Create_i_AC(dict_algorithm, dict_material, dict_sample)

      #Running moose
      os.system('mpiexec -n '+str(dict_algorithm['np_proc'])+' ~/projects/moose/modules/combined/combined-opt -i PF_'+str(dict_algorithm['i_PF'])+'.i')

      #sorting files
      j_str = Owntools.Sort_Files(dict_algorithm)

      simulation_report.tac_tempo(datetime.now(),'PF '+str(dict_algorithm['i_PF']))

      #-----------------------------------------------------------------------------
      # Conversion PF to DEM
      #-----------------------------------------------------------------------------

      simulation_report.tic_tempo(datetime.now())

      #Convert PF data into DEM data
      FileToRead = 'Output/PF_'+str(dict_algorithm['i_PF'])+'/PF_'+str(dict_algorithm['i_PF'])+'_other_'+j_str
      PFtoDEM_Multi(FileToRead,dict_algorithm,dict_sample)

      #Geometric study
      Owntools.Study_Polygonal_and_init_f(dict_algorithm,dict_geometry,dict_sample,dict_sollicitations,simulation_report)
      S_grains = 0
      S_grains_dissolvable = 0
      for grain in dict_sample['L_g']:
          S_grains = S_grains + grain.surface
          if grain.dissolved :
              S_grains_dissolvable = S_grains_dissolvable + grain.surface
      simulation_report.write('Total Surface '+str(int(S_grains))+' µm2\n')
      simulation_report.write('Total Surface dissolvable '+str(int(S_grains_dissolvable))+' µm2\n')

      # Tracker
      dict_tracker['t_L'].append(dict_tracker['t_L'][-1] + dict_algorithm['dt_PF']*dict_algorithm['n_t_PF'])
      dict_tracker['S_grains_L'].append(S_grains)
      dict_tracker['S_dissolved_L'].append(dict_tracker['S_grains_L'][0]-S_grains)
      dict_tracker['S_dissolved_perc_L'].append((dict_tracker['S_grains_L'][0]-S_grains)/(dict_tracker['S_grains_L'][0])*100)
      dict_tracker['n_grains_L'].append(len(dict_sample['L_g']))
      dict_tracker['S_grains_dissolvable_L'].append(S_grains_dissolvable)
      dict_tracker['S_dissolved_perc_dissolvable_L'].append((dict_tracker['S_grains_dissolvable_L'][0]-S_grains_dissolvable)/(dict_tracker['S_grains_dissolvable_L'][0])*100)

      simulation_report.tac_tempo(datetime.now(),'From PF to DEM iteration '+str(dict_algorithm['i_PF']))

      #-----------------------------------------------------------------------------
      # Reinitialisation of contact for the next step
      #-----------------------------------------------------------------------------

      for contact in dict_sample['L_contact']:
          contact.init_contact(dict_sample['L_g'])
      for contact in dict_sample['L_contact_gw']:
          contact.init_contact_gw(dict_sample['L_g'])

      #-----------------------------------------------------------------------------
      # Print Grains configuration
      #-----------------------------------------------------------------------------

      if dict_algorithm['Debug'] :
        Owntools.Debug_f2(dict_algorithm,dict_sample)

      #-----------------------------------------------------------------------------
      # Save tempo
      #-----------------------------------------------------------------------------

      if dict_algorithm['SaveData'] :
        Owntools.save_tempo(dict_algorithm,dict_tracker)

#-------------------------------------------------------------------------------
# toc
#-------------------------------------------------------------------------------

simulation_report.end(datetime.now())

#-------------------------------------------------------------------------------
# Debugging and Output
#-------------------------------------------------------------------------------

if dict_algorithm['Debug'] :

    #Making movies
    Owntools.make_mp4(int(2*dict_algorithm['i_PF']))

    #Trackers
    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['t_L'],dict_tracker['S_grains_L'])
    plt.title('Evolution of the grains surface')
    fig.savefig('Debug/Evolution_Surface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_L'][:-1],dict_tracker['k0_xmin_L'],label='k0 with xmin')
    plt.plot(dict_tracker['S_dissolved_L'][:-1],dict_tracker['k0_xmax_L'],label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Grains surface dissolved (µm2)')
    fig.savefig('Debug/Evolution_k0_with_TotalSurface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_perc_L'][:-1],dict_tracker['k0_xmin_L'],label='k0 with xmin')
    plt.plot(dict_tracker['S_dissolved_perc_L'][:-1],dict_tracker['k0_xmax_L'],label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Percentage of grains surface dissolved (%)')
    fig.savefig('Debug/Evolution_k0_with_percentage_dissolved.png')
    plt.close(1)

#-------------------------------------------------------------------------------
#Saving data
#-------------------------------------------------------------------------------

if dict_algorithm['SaveData'] :

    name_actual_folder = os.path.dirname(os.path.realpath(__file__))
    shutil.copytree(name_actual_folder, '../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder'])

    Owntools.save_final(dict_algorithm,dict_tracker)
