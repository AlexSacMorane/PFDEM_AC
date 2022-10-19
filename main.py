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

#Own function and class
from Write_txt import Write_txt
from PFtoDEM_Multi import PFtoDEM_Multi
from Create_i_AC import Create_i_AC
from Create_LG_IC import LG_tempo, From_LG_tempo_to_usable, Grain_Tempo
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
if dict_algorithm['SaveData'] or dict_algorithm['Debug'] or dict_algorithm['Debug_DEM']:
    simulation_report.write('\n')

#Creation of the grain (without PF)
simulation_report.write_and_print('Creation of the grains\n','\nCREATION OF THE GRAINS\n')
LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#Spatial discretisation
User.Add_SpatialDiscretisation(dict_geometry,dict_sample)

# PF parameters
User.Add_WidthInt_DoubleWellBarrier(dict_material, dict_sample)

simulation_report.tac_tempo(datetime.now(),'Initialisation')

#-------------------------------------------------------------------------------
#Creation of polygonal particles
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

# Creation of the real list of grains
From_LG_tempo_to_usable(dict_ic, dict_material, dict_sample)

#creation pf the dissolution .txt
Owntools.Write_e_dissolution_txt(dict_sample,dict_sollicitations)

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#load data needed
L_g = dict_sample['L_g']
gravity = dict_sollicitations['gravity']
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

#Saving the grain surface
S_grains = 0
for grain in L_g:
    grain.Geometricstudy(dict_geometry,dict_sample,simulation_report)
    grain.init_f_control(dict_sollicitations)
    S_grains = S_grains + grain.surface
simulation_report.write('Total Surface '+str(round(S_grains,0))+' µm2\n')

simulation_report.tac_tempo(datetime.now(),'Creation and study of polygonal particles')

#-------------------------------------------------------------------------------
#Distribution etai
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

Etai.etai_distribution_dissolution(dict_algorithm,dict_sample, dict_sollicitations, simulation_report)

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
#load data needed
L_ig_etai_undissolved = dict_sample['L_ig_etai_undissolved']
L_ig_etai_dissolved = dict_sample['L_ig_etai_dissolved']
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

#create eta for undissolved
L_etai_undissolved = []
for etai in range(0,len(L_ig_etai_undissolved)):
    L_etai_undissolved.append(Etai.Etai(etai,L_ig_etai_undissolved[etai]))
#create eta for dissolved
L_etai_dissolved = []
for etai in range(0,len(L_ig_etai_dissolved)):
    L_etai_dissolved.append(Etai.Etai(len(L_ig_etai_undissolved)+etai,L_ig_etai_dissolved[etai]))

#Add elements in dict
L_etai_undissolved = dict_sample['L_ig_etai_undissolved']
L_etai_dissolved = dict_sample['L_ig_etai_dissolved']

if dict_algorithm['Debug'] :
    Owntools.Debug_etai_f(dict_sample)

simulation_report.tac_tempo(datetime.now(),'Etai distribution')

#-------------------------------------------------------------------------------
#Main
#-------------------------------------------------------------------------------

#preparation
L_contact_gw = []
L_ij_contact_gw = []
id_contact_gw = 0
L_contact = []
L_ij_contact = []
id_contact = 0
i_PF = 0

# Tracker
t_L = [0]
S_grains_L = [S_grains]
S_dissolved_L = [0]
S_dissolved_perc_L = [0]
n_grains_L = [len(L_g)]
k0_xmin_L = []
k0_xmax_L = []

#create an new dict
dict_tracker = {
    't_L' : t_L,
    'S_grains_L' : S_grains_L,
    'S_dissolved_L' : S_dissolved_L,
    'S_dissolved_perc_L' : S_dissolved_perc_L,
    'n_grains_L' : n_grains_L,
    'k0_xmin_L' : k0_xmin_L,
    'k0_xmax_L' : k0_xmax_L
}

# add elements in dicts
dict_algorithm['i_PF'] = i_PF
dict_sample['L_contact_gw'] = L_contact_gw
dict_sample['L_ij_contact_gw'] = L_ij_contact_gw
dict_sample['id_contact_gw'] = id_contact_gw
dict_sample['L_contact'] = L_contact
dict_sample['L_ij_contact'] = L_ij_contact
dict_sample['id_contact'] = id_contact

if dict_algorithm['Debug'] :
    Owntools.Debug_f2(dict_algorithm,dict_sample)

while not User.Criteria_StopSimulation(dict_algorithm):
      i_PF = i_PF + 1

      # update element in dict
      dict_algorithm['i_PF'] = i_PF

      simulation_report.write_and_print('\nIteration '+str(i_PF)+' / '+str(dict_algorithm['n_t_PF'])+'\n','\nITERATION PF '+str(i_PF)+' / '+str(dict_algorithm['n_t_PF'])+'\n')

      #prepare iteration
      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF))
      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF)+'/txt')
      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF)+'/png')

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

      #Trackers
      Ecin_tracker = []
      Force_applied_tracker = []
      k0_xmin_tracker = []
      k0_xmax_tracker = []
      y_box_max_tracker = [dict_sample['y_box_max']]
      F_on_ymax_tracker = []

      #add element in dict
      dict_tracker['Ecin'] = Ecin_tracker
      dict_tracker['Force_applied'] = Force_applied_tracker
      dict_tracker['k0_xmin'] = k0_xmin_tracker
      dict_tracker['k0_xmax'] = k0_xmin_tracker
      dict_tracker['y_box_max'] = y_box_max_tracker
      dict_tracker['Force_on_upper_wall'] = F_on_ymax_tracker

      i_DEM = - 1
      DEM_loop_statut = True
      simulation_report.write('\n')
      simulation_report.tic_tempo(datetime.now())

      #-----------------------------------------------------------------------------
      # DEM iteration
      #-----------------------------------------------------------------------------

      while DEM_loop_statut :
          i_DEM = i_DEM + 1

          # update element in dict
          dict_algorithm['i_DEM'] = i_DEM

          for grain in dict_sample['L_g']:
              grain.init_f_control(dict_sollicitations)

          # Detection of contacts between grains
          if dict_algorithm['i_DEM'] % dict_algorithm['i_update_neighborhoods']  == 0:
              Contact.Update_Neighbouroods(dict_algorithm,dict_sample)
          Contact.Grains_Polyhedral_contact_Neighbouroods(dict_material,dict_sample)
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

              Owntools.save_tempo(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker)

              if dict_algorithm['Debug_DEM'] :
                Owntools.Debug_DEM_f(dict_algorithm, dict_sample)
                Write_txt(dict_algorithm,dict_sample)

          #-----------------------------------------------------------------------------
          # Stop conditions
          #-----------------------------------------------------------------------------

          if dict_algorithm['i_DEM'] >= dict_algorithm['i_DEM_stop'] :
              DEM_loop_statut = False
              print("DEM loop stopped by too many iterations.")
              simulation_report.tac_tempo(datetime.now(),'DEM loop '+str(dict_algorithm['i_PF']))
              simulation_report.write('/!\ End of DEM steps with '+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+'/!\ \n')
          if Ecin < dict_algorithm['Ecin_stop'] and dict_algorithm['i_DEM'] > dict_algorithm['n_window_stop'] and (dict_sollicitations['Vertical_Confinement_Force']*0.95<dict_sollicitations['Force_on_upper_wall'] and dict_sollicitations['Force_on_upper_wall']<dict_sollicitations['Vertical_Confinement_Force']*1.05):
              k0_xmin_window = dict_tracker['k0_xmin'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              k0_xmax_window = dict_tracker['k0_xmax'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              y_box_max_window = dict_tracker['y_box_max'][dict_algorithm['i_DEM']+1-dict_algorithm['n_window_stop']:dict_algorithm['i_DEM']+1]
              if max(k0_xmin_window) - min(k0_xmin_window) < dict_algorithm['dk0_stop'] and max(k0_xmax_window) - min(k0_xmax_window) < dict_algorithm['dk0_stop'] and max(y_box_max_window) - min(y_box_max_window) < dict_algorithm['dy_box_max_stop']:
                  DEM_loop_statut = False
                  print("DEM loop stopped by steady state reached.")
                  simulation_report.tac_tempo(datetime.now(),'DEM loop '+str(dict_algorithm['i_PF']))
                  simulation_report.write("DEM loop stopped by steady state reached with "+str(dict_algorithm['i_DEM']+1)+' iterations / '+str(dict_algorithm['i_DEM_stop']+1)+"\n")

      #-----------------------------------------------------------------------------
      # Debugging at the end of DEM step
      #-----------------------------------------------------------------------------

      if dict_algorithm['Debug'] :
        Owntools.Debug_f(dict_algorithm,dict_sample)
        Owntools.Debug_Trackers(dict_algorithm,dict_sollicitations,dict_tracker)
        Write_txt(dict_algorithm,dict_sample)
        Owntools.Plot_chain_force(dict_algorithm['i_PF'],dict_algorithm['i_DEM'])

        Owntools.save_final(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker)

      #-----------------------------------------------------------------------------
      # Compute Vertical and horizontal sollicitations to compute k0
      #-----------------------------------------------------------------------------

      Owntools.Control_y_max_NR(dict_sample,dict_sollicitations)
      Owntools.Compute_k0(dict_sample,dict_sollicitations)

      simulation_report.write('k0_xmin : '+str(round(dict_sample['k0_xmin'],2))+' / k0_xmax : '+str(round(dict_sample['k0_xmax'],2))+'\n')

      #Update element in dict
      dict_tracker['k0_xmin_L'].append(dict_sample['k0_xmin'])
      dict_tracker['k0_xmax_L'].append(dict_sample['k0_xmax'])

      #-----------------------------------------------------------------------------
      # Computing rigid body motion
      #-----------------------------------------------------------------------------

      L_rbm = []
      for id_grain in range(len(dict_sample['L_g'])):
          rbm = dict_sample['L_g'][id_grain].center - L_center_g[id_grain]
          L_rbm.append(rbm)

      #-----------------------------------------------------------------------------
      # Move PF
      #-----------------------------------------------------------------------------

      Owntools.Stop_Debug(simulation_report)

      if MovePF_selector == 'DeconstructRebuild':
          for grain in L_g:
              grain.DEMtoPF_Decons_rebuild(w,x_L,y_L,i_PF)
          for etai in L_etai_undissolved+L_etai_dissolved:
              etai.update_etai_M(L_g)
              etai.Write_txt_Decons_rebuild(i_PF,x_L,y_L)
      #does not work for multiple grain by etai
      #elif MovePF_selector == 'Interpolation':
      #  for id_grain in range(len(L_g)):
      #      grain = L_g[id_grain]
      #      rbm = L_rbm[id_grain]
      #      grain.DEMtoPF_Interpolation(x_L,y_L,i_PF,rbm)
      else :
          simulation_report.write('Method to move phase field not available !')
          raise ValueError('Method to move phase field not available !')

      #-----------------------------------------------------------------------------
      # PF Simulation
      #-----------------------------------------------------------------------------

      #Create the .i for moose
      Create_i_AC(i_PF,dt_PF,x_L,y_L,M_pf,kc_pf,n_t_PF_moose,double_well_height,L_etai_undissolved,L_etai_dissolved)

      simulation_report.tic_tempo(datetime.now())
      #Running moose
      os.system('mpiexec -n '+str(np_proc)+' ~/projects/moose/modules/combined/combined-opt -i PF_'+str(i_PF)+'.i')
      simulation_report.tac_tempo(datetime.now(),'PF '+str(i_PF))

      #sorting files
      j_str = Owntools.Sort_Files(np_proc,n_t_PF_moose,i_PF,L_g)

      #Convert PF data into DEM data
      simulation_report.tic_tempo(datetime.now())
      FileToRead = 'Output/PF_'+str(i_PF)+'/PF_'+str(i_PF)+'_other_'+j_str
      PFtoDEM_Multi(FileToRead,x_L,y_L,L_g,L_etai_undissolved+L_etai_dissolved,np_proc)

      #Geometric study
      S_grains = 0
      for grain in L_g:
          grain.Geometricstudy(self,dict_geometry,dict_sample,simulation_report)
          grain.init_f_control(dict_sollicitations)
          S_grains = S_grains + grain.surface

      simulation_report.tac_tempo(datetime.now(),'From PF to DEM '+str(i_PF))

      # Tracker
      t_L.append(t_L[-1] + dt_PF*n_t_PF_moose)
      S_grains_L.append(S_grains)
      S_dissolved_L.append(S_grains_L[0]-S_grains)
      S_dissolved_perc_L.append((S_grains_L[0]-S_grains)/(S_grains_L[0])*100)
      n_grains_L.append(len(L_g))

      simulation_report.write('Total Surface '+str(int(S_grains))+' µm2\n')

      #-----------------------------------------------------------------------------
      # Reinitialisation of contact for the next step
      #-----------------------------------------------------------------------------

      for contact in L_contact:
          contact.init_contact(L_g)
      for contact in L_contact_gw:
          contact.init_contact_gw(L_g)

      #-----------------------------------------------------------------------------
      # Print Grains configuration
      #-----------------------------------------------------------------------------

      if dict_algorithm['Debug'] :
        Owntools.Debug_f2(dict_algorithm,dict_sample)

      #-----------------------------------------------------------------------------
      # Save tempo
      #-----------------------------------------------------------------------------

      if SaveData :
        outfile = open('../Data_MG_Box_AC_M/'+name_folder+'_save_tempo','wb')
        dict = {}
        dict['k0_xmin_L'] = k0_xmin_L
        dict['k0_xmax_L'] = k0_xmax_L
        dict['S_dissolved_L'] = S_dissolved_L
        dict['S_dissolved_perc_L'] = S_dissolved_perc_L
        dict['S_grains_L'] = S_grains_L
        pickle.dump(dict,outfile)
        outfile.close()

#-------------------------------------------------------------------------------
# toc
#-------------------------------------------------------------------------------

simulation_report.end(datetime.now())

#-------------------------------------------------------------------------------
# Debugging
#-------------------------------------------------------------------------------

if Debug :

    #Making movies
    Owntools.make_mp4(int(2*i_PF))

    #Trackers
    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(t_L,S_grains_L)
    plt.title('Evolution of the grains surface')
    fig.savefig('Debug/Evolution_Surface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(S_dissolved_L[:-1],k0_xmin_L,label='k0 with xmin')
    plt.plot(S_dissolved_L[:-1],k0_xmax_L,label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Grains surface dissolved (µm2)')
    fig.savefig('Debug/Evolution_k0_with_TotalSurface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(S_dissolved_perc_L[:-1],k0_xmin_L,label='k0 with xmin')
    plt.plot(S_dissolved_perc_L[:-1],k0_xmax_L,label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Percentage of grains surface dissolved (%)')
    fig.savefig('Debug/Evolution_k0_with_percentage_dissolved.png')
    plt.close(1)

#-------------------------------------------------------------------------------
#Saving data
#-------------------------------------------------------------------------------

if SaveData :

    shutil.copytree('../Simu_MG_Box_AC_M', '../Data_MG_Box_AC_M/'+name_folder)

    os.remove('../Data_MG_Box_AC_M/'+name_folder+'_save_tempo')

    outfile = open('../Data_MG_Box_AC_M/'+name_folder+'_save','wb')
    dict = {}
    dict['k0_xmin_L'] = k0_xmin_L
    dict['k0_xmax_L'] = k0_xmax_L
    dict['S_dissolved_L'] = S_dissolved_L
    dict['S_dissolved_perc_L'] = S_dissolved_perc_L
    dict['S_grains_L'] = S_grains_L
    pickle.dump(dict,outfile)
    outfile.close()
