# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    #this function is called in main.py to have all the parameters needed in the simulation

    #---------------------------------------------------------------------------
    #Geometric parameters

    N_grain = 51 #number of grains
    R_mean = 350 #µm radius
    L_R = [1.1*R_mean, 1*R_mean, 0.9*R_mean] #from larger to smaller
    L_percentage_R = [1/3, 1/3, 1/3] #distribution of the different radius
    grain_discretisation = 20 #approximatively the number of vertices for one grain

    #write dict
    dict_geometry = {
    'N_grain' : N_grain,
    'R_mean' : R_mean,
    'L_R' : L_R,
    'L_percentage_R' : L_percentage_R,
    'grain_discretisation' : grain_discretisation
    }

    #---------------------------------------------------------------------------
    #Material parameters

    Y = 70*(10**9)*(10**6)*(10**(-12)) #Young Modulus µN/µm2
    nu = 0.3 #Poisson's ratio
    rho = 2500*10**(-6*3) #density kg/µm3
    rho_surf = 4/3*rho*R_mean #kg/µm2
    mu_friction_gg = 0.5 #grain-grain
    mu_friction_gw = 0 #grain-wall
    coeff_restitution = 0.2 #1 is perfect elastic
    # PF parameters
    M_pf = 1 # mobility
    kc_pf = 3 #graident coefficient

    #write dict
    dict_material = {
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'rho_surf' : rho_surf,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution,
    'M_pf' : M_pf,
    'kc_pf' : kc_pf
    }

    #---------------------------------------------------------------------------
    #Sample definition

    #Box définition
    x_box_min = 0 #µm
    x_box_max = 2*R_mean*math.sqrt(N_grain/0.6) #µm 0.6 from Santamarina, 2014 to avoid boundaries effect
    y_box_min = 0 #µm

    #write dict
    dict_sample = {
    'x_box_min' : x_box_min,
    'x_box_max' : x_box_max,
    'y_box_min' : y_box_min
    }

    #---------------------------------------------------------------------------
    #External sollicitations

    Vertical_Confinement_Pressure = 500*10**5 #Pa
    Vertical_Confinement_Force = Vertical_Confinement_Pressure*(x_box_max-x_box_min)*(2*R_mean)*10**(-6) #µN
    gravity = 0 #µm/s2
    Dissolution_Energy = 0.0025 #e-12 J
    frac_dissolved = 0.15 #Percentage of grain dissolved

    #write dict
    dict_sollicitations = {
    'Vertical_Confinement_Pressure' : Vertical_Confinement_Pressure,
    'Vertical_Confinement_Force' : Vertical_Confinement_Force,
    'gravity' : gravity,
    'Dissolution_Energy' : Dissolution_Energy,
    'frac_dissolved' : frac_dissolved
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    #Phase field
    dt_PF = 0.075 #s time step during MOOSE simulation
    n_t_PF = 5 #number of iterations PF-DEM
    factor_distribution_etai = 1.5 #margin to distribute etai
    MovePF_selector = 'DeconstructRebuild' #Move PF

    #DEM parameters
    dt_DEM_crit = math.pi*min(L_R)/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/5 #s time step during DEM simulation
    factor_neighborhood = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 100 #the frequency of the update of the neighborhood of the grains and the walls
    Spring_type = 'Ponctual' #Kind of contact
    #Stop criteria of the DEM
    i_DEM_stop = 3000 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0002
    n_window_stop = 50
    dk0_stop = 0.05
    dy_box_max_stop = 0.5

    #PF-DEM
    n_t_PF = 3 #number of cycle DEM-PF

    #Number of processor
    np_proc = 4

    #Debugging
    Debug = True #plot configuration before and after DEM simulation
    Debug_DEM = False #plot configuration inside DEM
    i_print_plot = 50 #frenquency of the print and plot in DEM step
    SaveData = False #save data
    main_folder_name = 'Data_MG_Box_AC_M' #where data are saved
    template_simulation_name = 'Run_' #template of the simulation name

    #write dict
    dict_algorithm = {
    'dt_PF' : dt_PF,
    'n_t_PF' : n_t_PF,
    'dt_DEM_crit' : dt_DEM_crit,
    'dt_DEM' : dt_DEM,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window_stop' : n_window_stop,
    'dk0_stop' : dk0_stop,
    'dy_box_max_stop' : dy_box_max_stop,
    'n_t_PF' : n_t_PF,
    'MovePF_selector' : MovePF_selector,
    'Spring_type' : Spring_type,
    'np_proc' : np_proc,
    'Debug' : Debug,
    'Debug_DEM' : Debug_DEM,
    'SaveData' : SaveData,
    'main_folder_name' : main_folder_name,
    'template_simulation_name' : template_simulation_name,
    'i_print_plot' : i_print_plot,
    'factor_neighborhood' : factor_neighborhood,
    'factor_distribution_etai' : factor_distribution_etai
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    n_generation = 2 #number of grains generation /!\ Work only for 2 /!\
    factor_ymax_box = 1.5 #margin to generate grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 3000 #stop criteria for DEM during IC
    i_update_neighborhoods_gen = 5 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 100 #the frequency of the update of the neighborhood of the grains and the walls during IC combination
    i_print_plot_IC = 100 #frenquency of the print and plot for IC
    dt_DEM_IC = 2*dt_DEM_crit/5 #s time step during IC
    Ecin_ratio_IC = 0.0005
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood

    #write dict
    dict_ic = {
    'n_generation' : n_generation,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'factor_ymax_box' : factor_ymax_box,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'N_test_max' : N_test_max
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations

#-------------------------------------------------------------------------------

def Add_SpatialDiscretisation(dict_geometry,dict_sample):

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    R_mean = dict_geometry['R_mean']
    N_grain = dict_geometry['N_grain']
    x_box_min = dict_sample['x_box_min']
    x_box_max = dict_sample['x_box_max']
    y_box_min = dict_sample['y_box_min']
    y_box_max = dict_sample['y_box_max']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

    #Spatial discretisation
    x_min = x_box_min-R_mean*0.25 #µm
    x_max = x_box_max+R_mean*0.25 #µm
    n_x =  int(math.sqrt(N_grain/0.6)*30) #30 node on a mean grain diameter, 0.6 from Santamarian 2014 to avoid boundaries effect
    x_L = list(np.linspace(x_min,x_max,n_x))
    y_min = y_box_min-R_mean*0.25 #µm
    y_max = y_box_max+R_mean*0.25 #µm
    n_y = int(n_x*0.6)
    y_L = list(np.linspace(y_min,y_max,n_y))

    #add elements in dict
    dict_sample['x_min'] = x_min
    dict_sample['x_max'] = x_max
    dict_sample['n_x'] = n_x
    dict_sample['x_L'] = x_L
    dict_sample['y_min'] = y_min
    dict_sample['y_max'] = y_max
    dict_sample['n_y'] = n_y
    dict_sample['y_L'] = y_L

#-------------------------------------------------------------------------------

def Add_WidthInt_DoubleWellBarrier(dict_material, dict_sample):
    #need Spatial dicretisation to run

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    x_L = dict_sample['x_L']
    y_L = dict_sample['y_L']
    kc_pf = dict_material['kc_pf']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    w = math.sqrt((x_L[4]-x_L[0])**2+(y_L[4]-y_L[0])**2)
    double_well_height = 20*kc_pf/w/w #double well height J/µm2
    #the factor 20 was found by try and error

    #add elements in dict
    dict_material['w'] = w
    dict_material['double_well_height'] = double_well_height

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    #Criteria to stop simulation (PF and DEM)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    i_PF = dict_algorithm['i_PF']
    n_t_PF = dict_algorithm['n_t_PF']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

    Criteria_Verified = False
    if i_PF >= n_t_PF:
        Criteria_Verified = True
    return Criteria_Verified
