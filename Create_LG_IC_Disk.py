# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define an initial configuration.
Grains are considered as circular to be faster.
We have 2 temporary classes about grains and contact."""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
import matplotlib.pyplot as plt
import Grain

#-------------------------------------------------------------------------------
#Classes definition
#-------------------------------------------------------------------------------

class Grain_Tempo:

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Radius, dict_material):
    '''
    defining the grain

    each grain is described by a id (an integer class)
                               a center (a array class [X,Y])
                               a radius (a float)
                               material properties (a dict)
    '''
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    Rho_surf = dict_material['rho_surf']
    Y = dict_material['Y']
    Nu = dict_material['nu']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    self.id = ID
    self.center = Center
    self.radius = Radius
    self.r_max = Radius
    self.y = Y
    self.nu = Nu
    self.g = Y /2/(1+Nu) #shear modulus
    self.rho_surf = Rho_surf
    self.mass = math.pi*Radius**2*Rho_surf
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.plot_preparation()
    self.neighbourood = []

#-------------------------------------------------------------------------------

  def add_F(self,F):
      '''add a force (an array [Fx,Fy]) to the grain'''
      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      '''
      initialize the force applied to the grain
      a gravity of g is applied
      '''
      self.fx = 0
      self.fy = -g*self.mass

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM):
    '''move the grain following a semi implicit euler scheme'''
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    self.center = self.center + self.v*dt_DEM
    for i in range(len(self.l_border_x)):
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*dt_DEM
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*dt_DEM

#-------------------------------------------------------------------------------

  def plot_preparation(self):
    '''prepare the plot of a grain'''
    N_p_border = 180 #number of nodes
    L_border_x = []
    L_border_y = []
    for i in range(N_p_border):
       theta = 2*math.pi*i/N_p_border
       L_border_x.append(self.center[0]+self.radius*math.cos(theta))
       L_border_y.append(self.center[1]+self.radius*math.sin(theta))
    L_border_x.append(L_border_x[0])
    L_border_y.append(L_border_y[0])
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y

#-------------------------------------------------------------------------------

class Contact_Tempo:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, Mu, Coeff_Restitution, Nature, G2 = None, Limit = None):
    '''defining the contact

    the contact can be grain-grain or grain-wall
    each contact is described by a id (an integer class)
                                 a first grain (a grain_tempo class)
                                 a friction coefficient (a float)
                                 a restitution coefficient (a float) : 1 all restitued / 0 any restitued
                                 a nature of the contact (a string)
                                 a second grain (a grain_tempo class), not used in the case of grain-wall
                                 a coordinate of the wall (a float), not used in the case of the grain-grain
    '''
    self.id = ID
    self.nature = Nature
    self.tangential_old_statut = False
    self.ft = 0
    self.overlap_tangential = 0
    self.mu = Mu
    self.coeff_restitution = Coeff_Restitution
    if Nature == 'gg':
        self.g1 = G1
        self.g2 = G2
        Y_eq = 1/((1-G1.nu*G1.nu)/G1.y+(1-G2.nu*G2.nu)/G2.y)
        R_eq = 1/(1/G1.radius+1/G2.radius)
        k = 4/3*Y_eq*math.sqrt(R_eq) #unlinear spring
        self.k = k
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        self.eta = eta
    else :
        self.g = G1
        self.limit = Limit
        factor = 5 #increase the stiffness of the wall
        k = factor*4/3*G1.y/(1-G1.nu*G1.nu)*math.sqrt(G1.radius) #unlinear spring
        self.k = k
        self.kt = 0
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        self.eta = eta

#-------------------------------------------------------------------------------

  def normal(self):
    '''compute the normal reaction of the contact

    2 cases : grain-grain and grain-wall (one "elif" by wall)'''
    if self.nature == 'gg':
        #unlinear stiffness
        n12 = (self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center)
        self.n12 = n12
        overlap = self.g1.radius + self.g2.radius - np.linalg.norm(self.g1.center - self.g2.center)
        self.overlap = overlap
        F12_n = self.k*overlap**(3/2)
        F12 = F12_n*n12
        self.F12_n = F12_n
        self.g1.add_F(-F12)
        self.g2.add_F( F12)
        #damping
        F12_damp = -np.dot(self.g2.v - self.g1.v,n12)*self.eta*n12
        self.g1.add_F(-F12_damp)
        self.g2.add_F( F12_damp)

    elif self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        overlap = self.limit-(self.g.center[1]-self.g.radius)
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg)
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.g.add_F(Fwg_damp)

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        overlap = (self.g.center[1]+self.g.radius)-self.limit
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg)
        #damping

    elif self.nature == 'gwx_min':
        #unlinear stiffness
        nwg = np.array([1,0])
        self.nwg = nwg
        overlap = self.limit-(self.g.center[0]-self.g.radius)
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg)
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.g.add_F(Fwg_damp)

    elif self.nature == 'gwx_max':
        #linear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        overlap = self.g.center[0]+self.g.radius-self.limit
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg)
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.g.add_F(Fwg_damp)

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
    '''compute the tangential reaction of the contact

    2 cases : grain-grain and grain-wall (one "elif" by wall)'''
    if self.nature == 'gg':
        if self.mu > 0 :
            #unlinear stiffness
            G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
            R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
            kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap))
            kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F12_n),0))
            self.kt = kt

            t12 = np.array([-self.n12[1], self.n12[0]])
            self.t12 = t12
            if self.tangential_old_statut:
                #if a reaction has been already computed
                #need to project the tangential reaction on the new tangential plane
                self.ft = self.ft*np.dot(self.t12_old,self.t12)
            else:
                self.tangential_old_statut = True
            Delta_Us = np.dot(self.g1.v-self.g2.v,self.t12) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            self.t12_old = self.t12
            if abs(self.ft) > abs(self.mu*self.F12_n) or kt == 0: #Coulomb criteria
                self.ft = self.mu * abs(self.F12_n) * np.sign(self.ft)
            F12 = -self.ft*t12
            self.g1.add_F(-F12)
            self.g2.add_F( F12)
            #damping
            gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
            mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
            eta = 2 * gamma * math.sqrt(mass_eq*kt)
            F12_damp = -np.dot(self.g2.v - self.g1.v,t12)*eta/2*t12
            self.g1.add_F(-F12_damp)
            self.g2.add_F( F12_damp)
        else :
            self.kt = 0
            self.t12 = np.array([-self.n12[1], self.n12[0]])

    elif self.nature == 'gwy_min':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([-1, 0])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.add_F(Fwg)
        else :
            twg = np.array([-1, 0])
            self.twg = twg

    elif self.nature == 'gwy_max':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([1, 0])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.add_F(Fwg)
        else :
            twg = np.array([1, 0])
            self.twg = twg

    elif self.nature == 'gwx_min':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([0, 1])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.add_F(Fwg)
        else :
            twg = np.array([0, 1])
            self.twg = twg

    elif self.nature == 'gwx_max':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([0, -1])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.add_F(Fwg)
        else :
            twg = np.array([0, 1])
            self.twg = twg

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    '''create an initial condition with disk grain'''
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    n_generation = dict_ic['n_generation']
    if n_generation != 2:
        simulation_report.write('n_generation must be equal to 2 !')
        raise ValueError('n_generation must be equal to 2 !')
    factor = dict_ic['factor_ymax_box']
    N_grain = dict_geometry['N_grain_disk']
    L_radius = dict_geometry['L_R']
    L_percentage_radius = dict_geometry['L_percentage_R']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #define the y_max for the grains generation
    radius_mean = 0
    for i in range(len(L_radius)):
        radius_mean = radius_mean + L_radius[i]*L_percentage_radius[i]
    dy_creation = N_grain / n_generation * factor*(2*radius_mean)**2/(x_max-x_min)

    #plan the grains generation
    L_n_grain_radius_try_one = []
    L_n_grain_radius = []
    L_n_grain_radius_done = []
    for percentage in L_percentage_radius:
        L_n_grain_radius_try_one.append(int(N_grain*percentage/n_generation))
        L_n_grain_radius.append(int(N_grain*percentage))
        L_n_grain_radius_done.append(0)

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    i_DEM = 0
    L_L_g_tempo = []

    #---------------------------------------------------------------------------

    print('First generation of grains')
    L_g_tempo = []

    #add elements in dicts
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_ic['L_L_g_tempo'] = L_L_g_tempo
    dict_ic['i_DEM_IC'] = i_DEM
    dict_ic['L_n_grain_radius_try_one'] = L_n_grain_radius_try_one
    dict_ic['L_n_grain_radius'] = L_n_grain_radius
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done
    dict_sample['y_box_min_ic'] = dict_sample['y_box_min']
    dict_sample['dy_creation'] = dy_creation

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 1, simulation_report)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in L_g_tempo:
        if grain.center[1]+grain.radius > y_max:
            y_max = grain.center[1]+grain.radius

    #add element in dict
    dict_sample['y_box_max'] = y_max

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #---------------------------------------------------------------------------

    print('Second generation of grains')
    L_g_tempo = []

    #update elements un dict
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_sample['y_box_min_ic'] = dict_sample['y_box_max']

    Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, 2, simulation_report)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = dict_sample['y_box_min_ic']
    for grain in L_g_tempo:
        if grain.center[1]+grain.radius > y_max:
            y_max = grain.center[1]+grain.radius

    #update element in dict
    dict_sample['y_box_max'] = y_max

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, False, simulation_report)

    #---------------------------------------------------------------------------

    print('Combine generations of grains')

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_L_g_tempo = dict_ic['L_L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #combine all smaller sample
    L_g = []
    for L_g_tempo in L_L_g_tempo:
        for g_tempo in L_g_tempo:
            L_g.append(g_tempo)

    #update element in dict
    dict_ic['L_g_tempo'] = L_g

    DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, True, simulation_report)

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    simulation_report.write_and_print(str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n','\n'+str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n')

    return L_g_tempo, y_max

#-------------------------------------------------------------------------------

def DEM_loading(dict_ic, dict_material, dict_sample, dict_sollicitations, multi_generation, simulation_report):
    '''loading the granular system'''
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    dt_DEM = dict_ic['dt_DEM_IC']
    i_DEM_stop = dict_ic['i_DEM_stop_IC']
    i_DEM = dict_ic['i_DEM_IC']
    Ecin_ratio_IC = dict_ic['Ecin_ratio_IC']
    i_print_plot_IC = dict_ic['i_print_plot_IC']
    factor_neighborhood_IC = dict_ic['factor_neighborhood_IC']
    if multi_generation :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_com']
        y_min = dict_sample['y_box_min']
    else :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_gen']
        y_min = dict_sample['y_box_min_ic']
    mu_gg = 0
    mu_gw = 0
    e_gg = dict_material['coeff_restitution']
    e_gw = dict_material['coeff_restitution']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_max = dict_sample['y_box_max']
    Forcev_target = dict_sollicitations['Vertical_Confinement_Force']
    gravity = dict_sollicitations['gravity']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


    i_DEM_0 = i_DEM
    DEM_loop_statut = True

    #Initialisation
    L_contact_gg = []
    L_contact_ij = []
    L_contact_gw = []
    L_contact_gw_ij = []
    id_contact = 0

    #trackers and stop conditions
    if gravity > 0:
        Force_stop = 0
        for grain in L_g_tempo:
            Force_stop = Force_stop + 0.8*grain.mass*gravity

    Force_tracker = []
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    for grain in L_g_tempo:
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(Ecin_ratio_IC*grain.radius/dt_DEM)**2

    while DEM_loop_statut :

        i_DEM = i_DEM + 1

        #Contact detection
        if (i_DEM-i_DEM_0-1) % i_update_neighborhoods  == 0:
            Update_Neighborhoods(dict_ic)
        L_contact_gg, L_contact_ij, id_contact = Grains_Disk_contact_Neighborhoods(L_g_tempo,L_contact_gg,L_contact_ij,id_contact,mu_gg,e_gg)
        # Detection of contacts between grain and walls
        if (i_DEM-i_DEM_0-1) % i_update_neighborhoods  == 0:
            wall_neighborhood = Update_wall_Neighborhoods(L_g_tempo,factor_neighborhood_IC,x_min,x_max,y_min,y_max)
        L_contact_gw, L_contact_gw_ij, id_contact = Grains_Disk_Wall_contact_Neighborhood(wall_neighborhood,L_contact_gw,L_contact_gw_ij,id_contact,x_min,x_max,y_min,y_max,mu_gw,e_gw)

        #Sollicitation computation
        for grain in L_g_tempo:
             grain.init_F_control(gravity)
        for contact in  L_contact_gg+L_contact_gw:
            contact.normal()
            contact.tangential(dt_DEM)

        #Move grains
        for grain in L_g_tempo :
            grain.euler_semi_implicite(dt_DEM)

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(L_g_tempo)):
            if L_g_tempo[id_grain].center[0] < x_min :
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[0] > x_max :
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[1] < y_min :
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[1] > y_max :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(L_g_tempo[id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(L_g_tempo[id_grain].id)+' has been deleted because it is out of the box')
            L_g_tempo.pop(id_grain)

        #Control the y_max to have the pressure target
        y_max, Fv = Control_y_max_NR(y_max,Forcev_target,L_contact_gw,L_g_tempo)

        #Tracker
        F = F_total(L_g_tempo)
        Ecin = E_cin_total(L_g_tempo)
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(y_max)

        if i_DEM % i_print_plot_IC == 0:
            if gravity > 0 :
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            else:
                print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            if dict_ic['Debug_DEM']:
                Plot_Config_Loaded(L_g_tempo,x_min,x_max,y_min,y_max,i_DEM)

        #Check stop conditions for DEM
        if i_DEM >= i_DEM_stop + i_DEM_0:
             DEM_loop_statut = False
        if gravity > 0 :
            if Ecin < Ecin_stop and F < Force_stop and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                  DEM_loop_statut = False
        else :
            if Ecin < Ecin_stop and i_DEM >= i_DEM_stop*0.1 + i_DEM_0 and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
                  DEM_loop_statut = False
        if L_g_tempo == []:
            DEM_loop_statut = False

    #trackers
    Plot_Trackers(Force_tracker, Ecin_tracker, Ymax_tracker, i_DEM)

    #Update dict
    dict_ic['L_g_tempo'] = L_g_tempo
    L_L_g_tempo = dict_ic['L_L_g_tempo']
    L_L_g_tempo.append(L_g_tempo.copy())
    dict_ic['L_L_g_tempo'] = L_L_g_tempo
    dict_sample['y_box_max'] = y_max
    dict_ic['i_DEM_IC'] = i_DEM

#-------------------------------------------------------------------------------

def Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, id_generation, simulation_report):
    '''generate the grains

    a position is tried, then we verify this new grain does not overlap with previously created ones
    '''
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g_tempo = dict_ic['L_g_tempo']
    N_test_max = dict_ic['N_test_max']
    if id_generation == 1:
        L_n_grain_radius = dict_ic['L_n_grain_radius_try_one']
    else :
        L_n_grain_radius = dict_ic['L_n_grain_radius']
    L_n_grain_radius_done = dict_ic['L_n_grain_radius_done']
    L_radius = dict_geometry['L_R']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_min = dict_sample['y_box_min_ic']
    dy_creation = dict_sample['dy_creation']
    y_max = y_min + dy_creation
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #Parameters for the method
    n_not_created = 0

    for i in range(len(L_radius)):
        radius = L_radius[i]
        n_grain = L_n_grain_radius[i]
        n_grain_done = L_n_grain_radius_done[i]
        last_id_grain_created = np.sum(L_n_grain_radius_done)
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*radius,x_max-1.1*radius),random.uniform(y_min+1.1*radius,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,radius,dict_material)
                grain_created = True
                for grain in L_g_tempo:
                    if Intersection(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_radius_done[i] = L_n_grain_radius_done[i] + 1

    #Update dict
    dict_ic['L_g_tempo'] = L_g_tempo
    dict_ic['L_n_grain_radius_done'] = L_n_grain_radius_done

#-------------------------------------------------------------------------------

def Intersection(g1,g2):
    '''verify if there is an overlap between two grains or not'''
    d_12 = np.linalg.norm(g1.center - g2.center)
    return d_12 < g1.radius + g2.radius

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    '''compute total kinetic energy'''
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    '''compute total force applied on grain'''
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    '''
    Control the upper wall to apply force

    a Newton-Raphson method is applied
    '''
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
    '''
    compute the function f to control the upper wall

    This function represents the difference between the force applied and the target value
    '''
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    '''compute the derivative function df to control the upper wall'''
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    '''the upper wall is moved as a single contact verify the target value'''
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

def Update_Neighborhoods(dict_ic):
    '''
    determine a neighborhood for each grain. This function is called every x time step

    grain contact is determined by Grains_Polyhedral_contact_Neighborhoods()
    notice that if there is a potential contact between grain_i and grain_j
    grain_i is not in the neighborhood of grain_j
    whereas grain_j is in the neighborhood of grain_i
    with i_grain < j_grain
    '''
    #factor determines the size of the neighborhood window
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g = dict_ic['L_g_tempo']
    factor = dict_ic['factor_neighborhood_IC']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    for i_grain in range(len(L_g)-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(L_g)):
            if np.linalg.norm(L_g[i_grain].center-L_g[j_grain].center) < factor*(L_g[i_grain].radius+L_g[j_grain].radius):
                neighborhood.append(L_g[j_grain])
        L_g[i_grain].neighbourood = neighborhood

    #Update dict
    dict_ic['L_g_tempo'] = L_g

#-------------------------------------------------------------------------------

def Grains_Disk_contact_Neighborhoods(L_g,L_contact,L_ij_contact,id_contact,mu_gg,e_gg):
    '''
    detect contact between a grain and grains from its neighborhood

    the neighborhood is updated with Update_Neighborhoods()
    '''
    for i_grain in range(len(L_g)-1) :
        grain_i = L_g[i_grain]
        for neighbour in L_g[i_grain].neighbourood:
            grain_j = neighbour
            j_grain = neighbour.id
            if Intersection(grain_i,grain_j):
                if (i_grain,j_grain) not in L_ij_contact:  #contact not detected previously
                   #creation of contact
                   L_ij_contact.append((i_grain,j_grain))
                   L_contact.append(Contact_Tempo(id_contact, grain_i, mu_gg, e_gg, 'gg', grain_j))
                   id_contact = id_contact + 1

            else :
                if (i_grain,j_grain) in L_ij_contact : #contact detected previously is not anymore
                       L_contact.pop(L_ij_contact.index((i_grain,j_grain)))
                       L_ij_contact.remove((i_grain,j_grain))

    return L_contact, L_ij_contact, id_contact

#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(L_g,factor,x_min,x_max,y_min,y_max):
    '''
    determine a neighborhoods for wall. This function is called every x time step

    grain_wall contact is determined by Grains_Polyhedral_Wall_contact_Neighborhood
    '''
    #factor determines the size of the neighborhood window
    wall_neighborhood = []
    for grain in L_g:

        p_x_min = min(grain.l_border_x)
        p_x_max = max(grain.l_border_x)
        p_y_min = min(grain.l_border_y)
        p_y_max = max(grain.l_border_y)

        #grain-wall x_min
        if abs(p_x_min-x_min) < factor*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-x_max) < factor*grain.radius :
            wall_neighburhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-y_min) < factor*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-y_max) < factor*grain.radius :
            wall_neighborhood.append(grain)

    return wall_neighborhood

#-------------------------------------------------------------------------------

def Grains_Disk_Wall_contact_Neighborhood(wall_neighborhood,L_contact_gw,L_contact_gw_ij,id_contact,x_min,x_max,y_min,y_max,mu_gw,e_gw):
  '''
  detect contact grain in the neighborhood of the wall and  the wall

  the neighborhood is updated with Update_wall_Neighborhoods()
  we realize iterations on the grain list and compare with the coordinate of the different walls
  '''
  for grain in wall_neighborhood:

      # contact grain-wall x_min
      if grain.center[0] < x_min + grain.radius and (grain.id,-1) not in L_contact_gw_ij:
          L_contact_gw.append(Contact_Tempo(id_contact, grain, mu_gw, e_gw, 'gwx_min', None, x_min))
          id_contact = id_contact + 1
          L_contact_gw_ij.append((grain.id,-1))
      elif grain.center[0] > x_min + grain.radius and (grain.id,-1) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      # contact grain-wall x_max
      if grain.center[0] > x_max - grain.radius and (grain.id,-2) not in L_contact_gw_ij:
          L_contact_gw.append(Contact_Tempo(id_contact, grain, mu_gw, e_gw, 'gwx_max', None, x_max))
          id_contact = id_contact + 1
          L_contact_gw_ij.append((grain.id,-2))
      elif grain.center[0] < x_max - grain.radius and (grain.id,-2) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      # contact grain-wall y_min
      if grain.center[1] < y_min + grain.radius and (grain.id,-3) not in L_contact_gw_ij:
          L_contact_gw.append(Contact_Tempo(id_contact, grain, mu_gw, e_gw, 'gwy_min', None, y_min))
          id_contact = id_contact + 1
          L_contact_gw_ij.append((grain.id,-3))
      elif grain.center[1] > y_min + grain.radius and (grain.id,-3) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      # contact grain-wall y_max
      if grain.center[1] > y_max - grain.radius and (grain.id,-4) not in L_contact_gw_ij:
          L_contact_gw.append(Contact_Tempo(id_contact, grain, mu_gw, e_gw, 'gwy_max', None, y_max))
          id_contact = id_contact + 1
          L_contact_gw_ij.append((grain.id,-4))
      elif grain.center[1] < y_max - grain.radius and (grain.id,-4) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-4))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)

  return L_contact_gw, L_contact_gw_ij, id_contact

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(L_g,x_min,x_max,y_min,y_max,i):
    '''plot loading configuration'''
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        grain.plot_preparation()
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
    '''plot final loading configuration'''
    plt.figure(1,figsize=(16,9))
    L_x = []
    L_y = []
    L_u = []
    L_v = []
    for grain in L_g:
        grain.plot_preparation()
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
        plt.plot(grain.center[0],grain.center[1],'xk')
        L_x.append(grain.center[0])
        L_y.append(grain.center[1])
        L_u.append(grain.fx)
        L_v.append(grain.fy)
    plt.plot([x_min,x_min,x_max,x_max,x_min],[y_max,y_min,y_min,y_max,y_max],'k')
    plt.axis('equal')
    plt.savefig('Debug/ConfigLoaded.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Trackers(Force_tracker, Ecin_tracker, Ymax_tracker, i_DEM):
    '''plot trackers'''
    plt.figure(1,figsize=(16,9))
    plt.subplot(321)
    plt.plot(Force_tracker)
    plt.vlines(int(2/3*len(Force_tracker)),min(Force_tracker),max(Force_tracker))
    plt.title('Force total')
    plt.subplot(322)
    plt.plot(Force_tracker[int(2/3*len(Force_tracker)):])
    plt.title('End force total')
    plt.subplot(323)
    plt.plot(Ecin_tracker)
    plt.vlines(int(2/3*len(Ecin_tracker)),min(Ecin_tracker),max(Ecin_tracker))
    plt.title('Kinetic energy')
    plt.subplot(324)
    plt.plot(Ecin_tracker[int(2/3*len(Ecin_tracker)):])
    plt.title('End Kinetic energy')
    plt.subplot(325)
    plt.plot(Ymax_tracker)
    plt.vlines(int(2/3*len(Ymax_tracker)),min(Ymax_tracker),max(Ymax_tracker))
    plt.title('Upper wall position')
    plt.subplot(326)
    plt.plot(Ymax_tracker[int(2/3*len(Ymax_tracker)):])
    plt.title('End upper wall position')
    plt.savefig('Debug/DEM_ite/Init/Trackers_'+str(i_DEM)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_geometry, dict_material, dict_sample):
    '''from a tempo configuration (circular grains), an initial configuration (polygonal grains) is generated'''
    L_g = []
    for grain_tempo in dict_ic['L_g_tempo']:

        L_border = []
        L_border_x = []
        L_border_y = []
        L_r = []
        L_theta_r = []
        for i_border in range(dict_geometry['grain_discretisation']):
            theta = 2*math.pi*i_border/dict_geometry['grain_discretisation']
            p = np.array(grain_tempo.center)+grain_tempo.radius*np.array([math.cos(theta),math.sin(theta)])
            L_border.append(p)
            L_border_x.append(p[0])
            L_border_y.append(p[1])
            L_r.append(grain_tempo.radius)
            L_theta_r = [theta]
        L_border.append(L_border[0])
        L_border_x.append(L_border_x[0])
        L_border_y.append(L_border_y[0])

        dict_ic_to_real = {
        'Id' : grain_tempo.id,
        'Y' : grain_tempo.y,
        'Nu' : grain_tempo.nu,
        'Rho_surf' : grain_tempo.rho_surf,
        'Center' : grain_tempo.center,
        'L_border' : L_border,
        'L_border_x' : L_border_x,
        'L_border_y' : L_border_y,
        'L_r' : L_r,
        'L_theta_r' : L_theta_r,
        'R_min' : grain_tempo.radius,
        'R_max' : grain_tempo.radius,
        'R_mean' : grain_tempo.radius,
        'Surface' : math.pi*grain_tempo.radius**2,
        'Mass' : math.pi*grain_tempo.radius**2*grain_tempo.rho_surf,
        'Inertia' : math.pi*grain_tempo.radius**2*grain_tempo.rho_surf*grain_tempo.radius**2
        }
        #create real grain
        L_g.append(Grain.Grain(dict_ic_to_real))

    #Add element in dict
    dict_sample['L_g'] = L_g

#-------------------------------------------------------------------------------

def IC(x_L,y_L,w,grain_tempo):
  '''create initial phase field, assuming a circular grain.'''
  #init the phase field
  etai_M_IC = np.zeros((len(y_L),len(x_L)))

  #extract a part focused on the grain
  x_extract_min = grain_tempo.center[0] - grain_tempo.radius - w
  x_extract_max = grain_tempo.center[0] + grain_tempo.radius + w
  y_extract_min = grain_tempo.center[1] - grain_tempo.radius - w
  y_extract_max = grain_tempo.center[1] + grain_tempo.radius + w

  #look for this part inside the global mesh
  #create search list
  x_L_search_min = abs(np.array(x_L)-x_extract_min)
  x_L_search_max = abs(np.array(x_L)-x_extract_max)
  y_L_search_min = abs(np.array(y_L)-y_extract_min)
  y_L_search_max = abs(np.array(y_L)-y_extract_max)

  #get index
  i_x_min = list(x_L_search_min).index(min(x_L_search_min))
  i_x_max = list(x_L_search_max).index(min(x_L_search_max))
  i_y_min = list(y_L_search_min).index(min(y_L_search_min))
  i_y_max = list(y_L_search_max).index(min(y_L_search_max))

  for y in y_L[i_y_min:i_y_max+1]:
    for x in x_L[i_x_min:i_x_max+1]:
      r = np.linalg.norm(np.array([x,y])-grain_tempo.center)
      if r<grain_tempo.radius-w/2:
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 1
      elif r>grain_tempo.radius+w/2:
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 0
      else :
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 0.5*(1 + np.cos(math.pi*(r-grain_tempo.radius+w/2)/w))

  return etai_M_IC

#-------------------------------------------------------------------------------
