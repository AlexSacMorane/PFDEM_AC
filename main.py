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
#Global parameters
#-------------------------------------------------------------------------------

simulation_report.tic_tempo(datetime.now())

dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

#Creation of the grain (without PF)
print('\nCREATION OF THE GRAINS\n')
LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#Spatial discretisation
User.Add_SpatialDiscretisation(dict_geometry,dict_sample)

# PF parameters
User.Add_WidthInt_DoubleWellBarrier(dict_material, dict_sample)

Owntools.Stop_Debug(simulation_report)

# Creation of the real list of grains
L_g = From_LG_tempo_to_usable(L_g_tempo,x_L,y_L,w)

#creation pf the dissolution .txt
Owntools.Write_e_dissolution_txt(Dissolution_Energy,x_L,y_L)

#Criteria to stop simulation (PF and DEM)
def Criteria_StopSimulation(i_PF,n_t_PF):
    Criteria_Verified = False
    if i_PF >= n_t_PF:
        Criteria_Verified = True
    return Criteria_Verified

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
    name_folder = dict_algorithm['template_simulation_name']+str(i_run)
if dict_algorithm['SaveData'] or dict_algorithm['Debug'] or dict_algorithm['Debug_DEM']:
    simulation_report.write('\n')

Owntools.Stop_Debug(simulation_report)

#-------------------------------------------------------------------------------
#Initial condition
#-------------------------------------------------------------------------------

#Saving the grain surface before pressure solution
S_grains = 0
for grain in L_g:
    grain.Geometricstudy(x_L,y_L,20,simulation_report)
    grain.init_f_control(gravity)
    S_grains = S_grains + grain.surface

simulation_report.tac_tempo(datetime.now(),'Creation of grains and initialisation')

#Distribution etai
simulation_report.tic_tempo(datetime.now())
L_ig_etai_undissolved, L_ig_etai_dissolved = Etai.etai_distribution_dissolution(L_g,frac_dissolved)
n_dissolved = 0
for L_ig in L_ig_etai_dissolved:
    n_dissolved = n_dissolved + len(L_ig)
n_total = 0
for L_ig in L_ig_etai_undissolved + L_ig_etai_dissolved:
    n_total = n_total + len(L_ig)

simulation_report.write_and_print('Real fraction dissolved : '+str(int(100*n_dissolved/n_total))+'\n','Real fraction dissolved : '+str(int(100*n_dissolved/n_total)))
simulation_report.write('Number of eta for undissolved grain : '+str(len(L_ig_etai_undissolved))+'\n'+\
                        'Number of eta for dissolved grain : '+str(len(L_ig_etai_dissolved))+'\n')

#create eta for undissolved
L_etai_undissolved = []
for etai in range(0,len(L_ig_etai_undissolved)):
    L_etai_undissolved.append(Etai.Etai(etai,L_ig_etai_undissolved[etai]))
#create eta for dissolved
L_etai_dissolved = []
for etai in range(0,len(L_ig_etai_dissolved)):
    L_etai_dissolved.append(Etai.Etai(len(L_ig_etai_undissolved)+etai,L_ig_etai_dissolved[etai]))
if Debug :
    Owntools.Debug_etai_f(L_g,x_box_min,x_box_max,y_box_min,y_box_max,x_L,y_L)

simulation_report.tac_tempo(datetime.now(),'Etai distribution')

#Detection of initial contact
L_contact_gw = []
L_ij_contact_gw = []
id_contact_gw = 0
L_contact = []
L_ij_contact = []
id_contact = 0
n_contact_deleted = 0

# Tracker
t_L = [0]
S_grains_L = [S_grains]
S_dissolved_L = [0]
S_dissolved_perc_L = [0]
n_grains_L = [len(L_g)]
k0_xmin_L = []
k0_xmax_L = []

simulation_report.write('Total Surface '+str(round(S_grains,0))+' µm2\n')

#-------------------------------------------------------------------------------
#Main
#-------------------------------------------------------------------------------

i_PF = 0
if Debug :
    Owntools.Debug_f2(L_g,i_PF,x_box_min,x_box_max,y_box_min,y_box_max,x_L,y_L)

while not Criteria_StopSimulation(i_PF,n_t_PF):

      i_PF = i_PF + 1
      simulation_report.write('\nIteration '+str(i_PF)+' / '+str(n_t_PF)+'\n')
      print('\nITERATION PF '+str(i_PF)+' / '+str(n_t_PF)+'\n')

      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF))
      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF)+'/txt')
      os.mkdir('Debug/DEM_ite/PF_'+str(i_PF)+'/png')

      # Saving to compute a rigid body motion
      L_center_g = []
      for grain in L_g:
          L_center_g.append(grain.center.copy())

      # Compute kinetic energy criteria
      Ecin_stop = 0
      for grain in L_g:
          Ecin_stop = Ecin_stop + 0.5*grain.m*(Ecin_ratio*grain.r_mean/dt_DEM)**2/len(L_g)

      #-----------------------------------------------------------------------------
      # DEM iteration
      #-----------------------------------------------------------------------------

      i_DEM = - 1

      #Trackers
      Ecin_tracker = []
      Force_applied_tracker = []
      k0_xmin_tracker = []
      k0_xmax_tracker = []
      y_box_max_tracker = [y_box_max]
      F_on_ymax_tracker = []

      simulation_report.write('\n')
      DEM_loop_statut = True
      simulation_report.tic_tempo(datetime.now())

      while DEM_loop_statut :
          i_DEM = i_DEM + 1

          for grain in L_g:
              grain.init_f_control(gravity)

          # Detection of contacts between grains
          #L_contact, L_ij_contact, id_contact = Contact.Grains_Polyhedral_contact(L_g,L_ij_contact,L_contact,id_contact,mu_friction,coeff_restitution)
          if i_DEM % i_update_neighbouroods  == 0:
              Contact.Update_Neighbouroods(L_g,1.5)
          L_contact, L_ij_contact, id_contact = Contact.Grains_Polyhedral_contact_Neighbouroods(L_g,L_ij_contact,L_contact,id_contact,mu_friction,coeff_restitution)
          # Detection of contacts between grain and walls
          #L_contact_gw, L_ij_contact_gw, id_contact_gw = Contact_gw.Grains_Polyhedral_Wall_contact(L_g,L_contact_gw,L_ij_contact_gw,id_contact_gw,x_box_min,x_box_max,y_box_min,y_box_max,0,coeff_restitution)
          if i_DEM % i_update_neighbouroods  == 0:
              wall_neighborhood = Contact_gw.Update_wall_Neighbouroods(L_g,1.5,x_box_min,x_box_max,y_box_min,y_box_max)
          L_contact_gw, L_ij_contact_gw, id_contact_gw = Contact_gw.Grains_Polyhedral_Wall_contact_Neighbourood(wall_neighborhood,L_contact_gw,L_ij_contact_gw,id_contact_gw,x_box_min,x_box_max,y_box_min,y_box_max,0,coeff_restitution)

          #Compute contact interactions (g-g and g-w)
          for contact in L_contact:
              if Spring_type == 'Ponctual':
                  contact.DEM_2grains_Polyhedral_normal()
                  contact.DEM_2grains_Polyhedral_tangential(dt_DEM)
              elif Spring_type == 'Surface':
                  contact.DEM_2grains_Polyhedral_normal_surface()
                  contact.DEM_2grains_Polyhedral_tangential_surface(dt_DEM)
              else :
                  simulation_report.write('Spring type not available !')
                  raise ValueError('Spring type not available !')
          for contact in L_contact_gw :
              if Spring_type == 'Ponctual':
                  contact.DEM_gw_Polyhedral_normal()
                  contact.DEM_gw_Polyhedral_tangential(dt_DEM)
              else : #Surface must be coded for contact gw
                  simulation_report.write('Spring type not available !')
                  raise ValueError('Spring type not available !')

          #Move particles and trackers
          #Semi implicit euler scheme
          Ecin = 0
          Force_applied = 0
          for grain in L_g:
              a_i = grain.f/grain.m
              v_i = grain.v + a_i*dt_DEM
              dw_i = grain.mz/grain.inertia
              w_i = grain.w + dw_i*dt_DEM
              grain.update_geometry_kinetic(v_i,a_i,w_i,dt_DEM) #Move grains
              Ecin = Ecin + 0.5*grain.m*np.linalg.norm(grain.v)**2/len(L_g)
              Force_applied = Force_applied + np.linalg.norm(grain.f)/len(L_g)
          Ecin_tracker.append(Ecin)
          Force_applied_tracker.append(Force_applied)

          #Control the y_max to verify vertical confinement
          y_box_max, F_on_ymax = Owntools.Control_y_max_NR(y_box_max,Vertical_Confinement_Force,L_contact_gw,L_g)
          y_box_max_tracker.append(y_box_max)
          F_on_ymax_tracker.append(F_on_ymax)
          F_on_xmin, F_on_xmax = Owntools.Compute_horizontal_sollicitations(L_contact_gw)
          k0_xmin, k0_xmax = Owntools.Compute_k0(F_on_xmin, F_on_xmax, F_on_ymax, x_box_min, x_box_max, y_box_min, y_box_max)
          k0_xmin_tracker.append(k0_xmin)
          k0_xmax_tracker.append(k0_xmax)

          if i_DEM %i_print_plot == 0:
              print('\nPF '+str(i_PF)+' -> i_DEM '+str(i_DEM+1)+' / '+str(i_DEM_stop+1)+' (max)')
              print('Ecin',int(Ecin),'/',int(Ecin_stop),'('+str(int(100*Ecin/Ecin_stop)),' %)')
              print('F_confinement',int(F_on_ymax),'/',int(Vertical_Confinement_Force),'('+str(int(100*F_on_ymax/Vertical_Confinement_Force)),' %)')

              outfile = open('Debug/DEM_ite/PF_'+str(i_PF)+'/save_tempo','wb')
              dict = {}
              dict['E_cin'] = Ecin_tracker
              dict['E_cin_stop'] = Ecin_stop
              dict['Force'] = Force_applied_tracker
              dict['n_window_stop'] = n_window_stop
              dict['k0_xmin_tracker'] = k0_xmin_tracker
              dict['k0_xmax_tracker'] = k0_xmax_tracker
              dict['dk0_stop'] = dk0_stop
              dict['y_box_max'] = y_box_max_tracker[:-1]
              dict['dy_box_max_stop'] = dy_box_max_stop
              dict['F_on_ymax'] = F_on_ymax_tracker
              dict['L_g'] = L_g
              dict['L_contact'] = L_contact
              dict['L_ij_contact'] = L_ij_contact
              dict['L_contact_gw'] = L_contact_gw
              dict['L_ij_contact_gw'] = L_ij_contact_gw
              pickle.dump(dict,outfile)
              outfile.close()

              if Debug_DEM :
                Owntools.Debug_DEM_f(L_g,L_contact,L_contact_gw,i_PF,i_DEM,x_box_min,x_box_max,y_box_min,y_box_max,x_L,y_L)
                Write_txt(L_g,L_contact,L_contact_gw,i_PF,i_DEM)

          #-----------------------------------------------------------------------------
          # Stop conditions
          #-----------------------------------------------------------------------------

          if i_DEM >= i_DEM_stop :
              DEM_loop_statut = False
              print("DEM loop stopped by too many iterations.")
              simulation_report.tac_tempo(datetime.now(),'DEM loop '+str(i_PF))
              simulation_report.write('/!\ End of DEM steps with '+str(i_DEM+1)+' iterations / '+str(i_DEM_stop+1)+'/!\ \n')
          if Ecin < Ecin_stop and i_DEM > n_window_stop and (Vertical_Confinement_Force*0.95<F_on_ymax and F_on_ymax<Vertical_Confinement_Force*1.05):
              k0_xmin_window = k0_xmin_tracker[i_DEM+1-n_window_stop:i_DEM+1]
              k0_xmax_window = k0_xmax_tracker[i_DEM+1-n_window_stop:i_DEM+1]
              y_box_max_window = y_box_max_tracker[i_DEM+1-n_window_stop:i_DEM+1]
              if max(k0_xmin_window) - min(k0_xmin_window) < dk0_stop and max(k0_xmax_window) - min(k0_xmax_window) < dk0_stop and max(y_box_max_window) - min(y_box_max_window) < dy_box_max_stop:
                  DEM_loop_statut = False
                  print("DEM loop stopped by steady state reached.")
                  simulation_report.tac_tempo(datetime.now(),'DEM loop '+str(i_PF))
                  simulation_report.write("DEM loop stopped by steady state reached with "+str(i_DEM+1)+' iterations / '+str(i_DEM_stop+1)+"\n")

      #-----------------------------------------------------------------------------
      # Computing rigid body motion
      #-----------------------------------------------------------------------------

      L_rbm = []
      for id_grain in range(len(L_g)):
          rbm = L_g[id_grain].center - L_center_g[id_grain]
          L_rbm.append(rbm)

      #-----------------------------------------------------------------------------
      # Debugging at the end of DEM step
      #-----------------------------------------------------------------------------

      if Debug :
        Owntools.Debug_f(L_g,i_PF,x_box_min,x_box_max,y_box_min,y_box_max,x_L,y_L)
        Owntools.Debug_Trackers(Ecin_tracker,Ecin_stop,Force_applied_tracker,k0_xmin_tracker,k0_xmax_tracker,y_box_max_tracker,F_on_ymax_tracker,Vertical_Confinement_Force,i_PF)
        Write_txt(L_g,L_contact,L_contact_gw,i_PF,i_DEM)
        Owntools.Plot_chain_force(i_PF,i_DEM)

        os.remove('Debug/DEM_ite/PF_'+str(i_PF)+'/save_tempo')
        outfile = open('Debug/DEM_ite/PF_'+str(i_PF)+'/save_final','wb')
        dict = {}
        dict['E_cin'] = Ecin_tracker
        dict['E_cin_stop'] = Ecin_stop
        dict['Force'] = Force_applied_tracker
        dict['n_window_stop'] = n_window_stop
        dict['k0_xmin_tracker'] = k0_xmin_tracker
        dict['k0_xmax_tracker'] = k0_xmax_tracker
        dict['dk0_stop'] = dk0_stop
        dict['y_box_max'] = y_box_max_tracker[:-1]
        dict['dy_box_max_stop'] = dy_box_max_stop
        dict['F_on_ymax'] = F_on_ymax_tracker
        dict['L_g'] = L_g
        dict['L_contact'] = L_contact
        dict['L_ij_contact'] = L_ij_contact
        dict['L_contact_gw'] = L_contact_gw
        dict['L_ij_contact_gw'] = L_ij_contact_gw
        pickle.dump(dict,outfile)
        outfile.close()

      #-----------------------------------------------------------------------------
      # Compute Vertical and horizontal sollicitations to compute k0
      #-----------------------------------------------------------------------------

      F_on_xmin, F_on_xmax = Owntools.Compute_horizontal_sollicitations(L_contact_gw)
      k0_xmin, k0_xmax = Owntools.Compute_k0(F_on_xmin, F_on_xmax, F_on_ymax, x_box_min, x_box_max, y_box_min, y_box_max)
      k0_xmin_L.append(k0_xmin)
      k0_xmax_L.append(k0_xmax)

      simulation_report.write('F_xmin : '+str(int(F_on_xmin))+' / F_xmax : '+str(int(F_on_xmax))+'\n')
      simulation_report.write('F_ymax : '+str(int(F_on_ymax))+' / ymax : '+str(round(y_max,2))+'\n')
      simulation_report.write('k0_xmin : '+str(round(k0_xmin,2))+' / k0_xmax : '+str(round(k0_xmax,2))+'\n')

      #-----------------------------------------------------------------------------
      # Move PF
      #-----------------------------------------------------------------------------

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
          grain.Geometricstudy(x_L,y_L,20,simulation_report)
          grain.init_f_control(gravity)
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

      if Debug :
        Owntools.Debug_f2(L_g,i_PF,x_box_min,x_box_max,y_box_min,y_box_max,x_L,y_L)

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
