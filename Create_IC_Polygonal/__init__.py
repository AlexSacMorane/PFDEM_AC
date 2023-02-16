# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation for initial conditions with polygonal grains.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#Own
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gw_ic
from Create_IC_Polygonal.Grain_ic_polygonal import Grain_Tempo_Polygonal
from Create_IC_Polygonal.Contact_gg_ic_polygonal import Contact_Tempo_Polygonal, Update_Neighborhoods, Grains_contact_Neighborhoods
from Create_IC_Polygonal.Contact_gw_ic_polygonal import Contact_gw_Tempo_Polygonal, Update_wall_Neighborhoods, Grains_Polyhedral_Wall_contact_Neighborhood
import Grain

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def Discretize_Grains(dict_ic, n_border):
    """
    Convert sphere particle to polygonal particule.

        Input :
            an ic dictionnary (a dict)
            a number of vertices (an int)
        Output :
            an ic discrete grain (a dict)
    """
    dict_ic_discrete = dict_ic.copy()
    for i_grain in range(len(dict_ic['L_g_tempo'])):
        dict_ic_discrete['L_g_tempo'][i_grain] = Grain_Tempo_Polygonal(dict_ic['L_g_tempo'][i_grain], n_border)

    return dict_ic_discrete

#-------------------------------------------------------------------------------

def DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    i_DEM_0 = dict_ic['i_DEM_IC']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    k0_tracker = []
    Fv_tracker = []
    for grain in dict_ic['L_g_tempo']:
        Force_stop = Force_stop + 0.5*grain.mass*dict_sollicitation['gravity']
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_ic['Ecin_ratio_IC']*grain.radius/dict_ic['dt_DEM_IC'])**2
        grain.v = np.array([0, 0])

    while DEM_loop_statut :

        dict_ic['i_DEM_IC'] = dict_ic['i_DEM_IC'] + 1

        #Contact detection
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            Update_Neighborhoods(dict_ic)
        Grains_contact_Neighborhoods(dict_ic,dict_material)

        # Detection of contacts between grain and walls
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            wall_neighborhood = Update_wall_Neighborhoods(dict_ic['L_g_tempo'],dict_ic['factor_neighborhood_IC'],dict_sample['x_box_min'],dict_sample['x_box_max'],dict_sample['y_box_min'],dict_sample['y_box_max'])
        Grains_Polyhedral_Wall_contact_Neighborhood(wall_neighborhood,dict_sample['x_box_min'],dict_sample['x_box_max'],dict_sample['y_box_min'],dict_sample['y_box_max'], dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in  dict_ic['L_contact']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Move grains
        for grain in dict_ic['L_g_tempo']:
            grain.euler_semi_implicite(dict_ic['dt_DEM_IC'],10*dict_ic['Ecin_ratio_IC'])

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_ic['L_g_tempo'])):
            if dict_ic['L_g_tempo'][id_grain].center[0] < dict_sample['x_box_min'] :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[0] > dict_sample['x_box_max'] :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] < dict_sample['y_box_min'] :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[1] > dict_sample['y_box_max'] :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
            dict_ic['L_g_tempo'].pop(id_grain)

        #Control the y_max to have the pressure target
        dict_sample['y_box_max'], Fv = Control_y_max_NR(dict_sample['y_box_max'],dict_sollicitation['Vertical_Confinement_Force'],dict_ic['L_contact_gw'],dict_ic['L_g_tempo'])

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(dict_sample['y_box_max'])
        Fv_tracker.append(Fv)
        F_xmin = 0
        F_xmax = 0
        for contact in dict_ic['L_contact_gw']:
            if contact.nature == 'gwx_min' :
                F_xmin = F_xmin + contact.Fwg_n
            if contact.nature == 'gwx_max' :
                F_xmax = F_xmax + contact.Fwg_n
        if Fv != 0:
            k0 = np.mean([F_xmin/(dict_sample['y_box_max'] - dict_sample['y_box_min'])*(dict_sample['x_box_max'] - dict_sample['x_box_min'])/Fv,
                          F_xmax/(dict_sample['y_box_max'] - dict_sample['y_box_min'])*(dict_sample['x_box_max'] - dict_sample['x_box_min'])/Fv])
        else :
            k0 = 0
        k0_tracker.append(k0)

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] ==0:
            if dict_sollicitation['gravity'] > 0 :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            else :
                print('i_DEM',dict_ic['i_DEM_IC'],'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            if dict_ic['Debug_DEM'] :
                Plot_Config_Loaded(dict_ic['L_g_tempo'],dict_sample['x_box_min'],dict_sample['x_box_max'],dict_sample['y_box_min'],dict_sample['y_box_max'],dict_ic['i_DEM_IC'])

        #Check stop conditions for DEM
        if dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC'] + i_DEM_0:
             DEM_loop_statut = False
        if dict_sollicitation['gravity'] > 0:
            if Ecin < Ecin_stop and F < Force_stop and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                  DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC']*0.1 + i_DEM_0 and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #update dict
    dict_ic['Ecin_tracker'] = Ecin_tracker
    dict_ic['Ymax_tracker'] = Ymax_tracker
    dict_ic['k0_tracker'] = k0_tracker
    dict_ic['Fv_tracker'] = Fv_tracker

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    """
    Compute total kinetic energy.

        Input :
            a list of temporary grains (a list)
        Output :
            the total kinetic energy (a float)
    """
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    """
    Compute total force applied on grains in the sample.

        Input :
            a list of temporary grains (a list)
        Output :
            the total force applied (a float)
    """
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a coordinate of the upper wall (a float)
            a confinement value (a float)
            a list of contact grain - wall (a list)
            a list of temporary grain (a list)
        Output :
            the coordinate of the upper wall (a float)
            a force applied on the upper wall before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        y_max = y_max + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        y_max = Reset_y_max(L_g,Force_target)

    for contact in L_contact_gw:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = y_max

    return y_max, F

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_ymax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    """
    The upper wall is located as a single contact verify the target value.

        Input :
            the list of temporary grains (a list)
            the confinement force (a float)
        Output :
            the upper wall position (a float)
    """
    y_max = None
    id_grain_max = None
    for id_grain in range(len(L_g)):
        grain = L_g[id_grain]
        y_max_grain = max(grain.l_border_y)

        if y_max != None and y_max_grain > y_max:
            y_max = y_max_grain
            id_grain_max = id_grain
        elif y_max == None:
            y_max = y_max_grain
            id_grain_max = id_grain

    factor = 5
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].radius)
    y_max = y_max - (Force/k)**(2/3)

    return y_max

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(L_g,x_min,x_max,y_min,y_max,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/DEM_ite/Init/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Config_Loaded_End(L_g,x_min,x_max,y_min,y_max):
    """
    Plot loaded configuration at the end of the initial configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/DEM_ite/ConfigLoaded.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample):
    """
    Create a real grain from a temporary grain.

    The phase variable is built. The distance between the point of the mesh and the particle center determines the value of the variable.
    A cosine profile is applied inside the interface.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with the list of real grains
    """
    L_g = []
    for grain_tempo in dict_ic['L_g_tempo']:

        dict_ic_to_real = {
        'Id' : grain_tempo.id,
        'Dissolved' : grain_tempo.dissolved,
        'Y' : grain_tempo.y,
        'Nu' : grain_tempo.nu,
        'Rho_surf' : grain_tempo.rho_surf,
        'Center' : grain_tempo.center,
        'L_border' : grain_tempo.l_border,
        'L_border_x' : grain_tempo.l_border_x,
        'L_border_y' : grain_tempo.l_border_y,
        'L_r' : grain_tempo.l_r,
        'L_theta_r' : grain_tempo.l_theta_r,
        'R_min' : min(grain_tempo.l_r),
        'R_max' : max(grain_tempo.l_r),
        'R_mean' : np.mean(grain_tempo.l_r),
        'Surface' : grain_tempo.surface,
        'Mass' : grain_tempo.mass,
        'Inertia' : grain_tempo.inertia
        }
        #create real grain
        L_g.append(Grain.Grain(dict_ic_to_real))

    #Add element in dict
    dict_sample['L_g'] = L_g

#-------------------------------------------------------------------------------
