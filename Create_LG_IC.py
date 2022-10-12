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
import pickle
import matplotlib.pyplot as plt
import Grain

#-------------------------------------------------------------------------------
#Classes definition
#-------------------------------------------------------------------------------

class Grain_Tempo:

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Radius, Y, Nu, Rho_surf):
    #defining the grain
    #each grain is described by a id (an integer class)
    #                           a center (a array class [X,Y])
    #                           a radius (a float)
    #                           a Young modulus (a float)
    #                           a Poisson's ratio (a float)
    #                           a surface mass (a float)

    self.id = ID
    self.center = Center
    self.radius = Radius
    self.y = Y
    self.nu = Nu
    self.g = Y /2/(1+Nu) #shear modulus
    self.rho_surf = Rho_surf
    self.mass = math.pi*Radius**2*Rho_surf
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.plot_preparation()

#-------------------------------------------------------------------------------

  def add_F(self,F):
      #add a force (an array [Fx,Fy]) to the grain

      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      #initialize the force applied to the grain
      #a gravity of g is applied

      self.fx = 0
      self.fy = -g*self.mass

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM):
    #move the grain following a semi implicit euler scheme

    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    self.center = self.center + self.v*dt_DEM
    for i in range(len(self.l_border_x)):
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*dt_DEM
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*dt_DEM

#-------------------------------------------------------------------------------

  def plot_preparation(self):
    #prepare the plot of a grain

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
    #defining the contact
    #the contact can be grain-grain or grain-wall
    #each contact is described by a id (an integer class)
    #                             a first grain (a grain_tempo class)
    #                             a friction coefficient (a float)
    #                             a restitution coefficient (a float) : 1 all restitued / 0 any restitued
    #                             a nature of the contact (a string)
    #                             a second grain (a grain_tempo class), not used in the case of grain-wall
    #                             a coordinate of the wall (a float), not used in the case of the grain-grain

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
    #compute the normal reaction of the contact

    #2 cases : grain-grain and grain-wall (one "elif" by wall)
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
    #compute the tangential reaction of the contact

    #2 cases : grain-grain and grain-wall (one "elif" by wall)
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

def LG_tempo(x_min,x_max,y_min,N_grain,L_radius,L_percentage_radius,rho_surf,Y,nu,mu_gg,e_gg,mu_gw,e_gw,Force_target,gravity,dt_DEM,simulation_report):
    #create an initial condition

    #number of grains generation
    n_generation = 2 #Work only for 2

    #define the y_max for the grains generation
    factor = 1.5
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

    #Parameters for the DEM
    i_DEM_stop = 10000
    i_DEM = 0

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    simulation_report.write('Creation of the grains\n')
    L_L_g_tempo = []
    y_min_init = y_min #save for rebuild

    print('First generation of grains')
    L_g_tempo = []
    L_g_tempo, L_n_grain_radius_done = Create_grains(L_g_tempo,L_n_grain_radius_try_one,L_n_grain_radius_done,L_radius,x_min,x_max,y_min,y_min+dy_creation,Y,nu,rho_surf,simulation_report)
    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = y_min
    for grain in L_g_tempo:
        if grain.center[1]+grain.radius > y_max:
            y_max = grain.center[1]+grain.radius
    L_g_tempo, y_min, i_DEM = DEM_loading(L_g_tempo, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min, y_max, dt_DEM, Force_target, gravity, i_DEM_stop, i_DEM, simulation_report)
    L_L_g_tempo.append(L_g_tempo.copy())

    print('Second generation of grains')
    L_g_tempo = []
    L_g_tempo, L_n_grain_radius_done = Create_grains(L_g_tempo,L_n_grain_radius,L_n_grain_radius_done,L_radius,x_min,x_max,y_min,y_min+dy_creation,Y,nu,rho_surf,simulation_report)
    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = y_min
    for grain in L_g_tempo:
        if grain.center[1]+grain.radius > y_max:
            y_max = grain.center[1]+grain.radius
    L_g_tempo, y_max, i_DEM = DEM_loading(L_g_tempo, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min, y_max, dt_DEM, Force_target, gravity, i_DEM_stop, i_DEM, simulation_report)
    L_L_g_tempo.append(L_g_tempo.copy())

    print('Combine generations of grains')
    L_g = []
    for L_g_tempo in L_L_g_tempo:
        for g_tempo in L_g_tempo:
            L_g.append(g_tempo)
    L_g_tempo, y_max, i_DEM = DEM_loading(L_g, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min_init, y_max, dt_DEM, Force_target, gravity, i_DEM_stop, i_DEM, simulation_report)

    simulation_report.write_and_print(str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n','\n'+str(len(L_g_tempo))+' / '+str(N_grain)+' grains have been created\n')

    return L_g_tempo, y_max

#-------------------------------------------------------------------------------

def DEM_loading(L_g_tempo, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min, y_max, dt_DEM, Forcev_target, gravity, i_DEM_stop, i_DEM, simulation_report):
    #loading the granular

    i_DEM_0 = i_DEM
    DEM_loop_statut = True

    #Initialisation
    L_contact_gg = []
    L_contact_ij = []
    L_contact_gw = []
    L_contact_gw_ij = []
    id_contact = 0

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Ymax_tracker = []
    Ymax_stop = 0
    for grain in L_g_tempo:
        Force_stop = Force_stop + 0.8*grain.mass*gravity
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(0.00005*grain.radius/dt_DEM)**2

    while DEM_loop_statut :

        i_DEM = i_DEM + 1

        #Contact detection
        for i_grain in range(len(L_g_tempo)-1):
              for j_grain in range(i_grain+1,len(L_g_tempo)):
                  #contact grain-grain
                  if Intersection(L_g_tempo[i_grain],L_g_tempo[j_grain]) and (i_grain,j_grain) not in L_contact_ij:
                      L_contact_gg.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], mu_gg, e_gg, 'gg', L_g_tempo[j_grain]))
                      id_contact = id_contact + 1
                      L_contact_ij.append((i_grain,j_grain))
                  elif not Intersection(L_g_tempo[i_grain],L_g_tempo[j_grain]) and (i_grain,j_grain) in L_contact_ij:
                      i_contact = L_contact_ij.index((i_grain,j_grain))
                      L_contact_gg.pop(i_contact)
                      L_contact_ij.pop(i_contact)
        for i_grain in range(len(L_g_tempo)):
              # contact grain-wall x_min
              if L_g_tempo[i_grain].center[0] < x_min + L_g_tempo[i_grain].radius and (i_grain,-1) not in L_contact_gw_ij:
                  L_contact_gw.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwx_min', None, x_min))
                  L_contact_gw_ij.append((i_grain,-1))
              elif L_g_tempo[i_grain].center[0] > x_min + L_g_tempo[i_grain].radius and (i_grain,-1) in L_contact_gw_ij:
                  i_contact = L_contact_gw_ij.index((i_grain,-1))
                  L_contact_gw.pop(i_contact)
                  L_contact_gw_ij.pop(i_contact)
              # contact grain-wall x_max
              if L_g_tempo[i_grain].center[0] > x_max - L_g_tempo[i_grain].radius and (i_grain,-2) not in L_contact_gw_ij:
                  L_contact_gw.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwx_max', None, x_max))
                  id_contact = id_contact + 1
                  L_contact_gw_ij.append((i_grain,-2))
              elif L_g_tempo[i_grain].center[0] < x_max - L_g_tempo[i_grain].radius and (i_grain,-2) in L_contact_gw_ij:
                  i_contact = L_contact_gw_ij.index((i_grain,-2))
                  L_contact_gw.pop(i_contact)
                  L_contact_gw_ij.pop(i_contact)
              # contact grain-wall y_min
              if L_g_tempo[i_grain].center[1] < y_min + L_g_tempo[i_grain].radius and (i_grain,-3) not in L_contact_gw_ij:
                  L_contact_gw.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwy_min', None, y_min))
                  id_contact = id_contact + 1
                  L_contact_gw_ij.append((i_grain,-3))
              elif L_g_tempo[i_grain].center[1] > y_min + L_g_tempo[i_grain].radius and (i_grain,-3) in L_contact_gw_ij:
                  i_contact = L_contact_gw_ij.index((i_grain,-3))
                  L_contact_gw.pop(i_contact)
                  L_contact_gw_ij.pop(i_contact)
              # contact grain-wall y_max
              if L_g_tempo[i_grain].center[1] > y_max - L_g_tempo[i_grain].radius and (i_grain,-4) not in L_contact_gw_ij:
                  L_contact_gw.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwy_max', None, y_max))
                  id_contact = id_contact + 1
                  L_contact_gw_ij.append((i_grain,-4))
              elif L_g_tempo[i_grain].center[1] < y_max - L_g_tempo[i_grain].radius and (i_grain,-4) in L_contact_gw_ij:
                  i_contact = L_contact_gw_ij.index((i_grain,-4))
                  L_contact_gw.pop(i_contact)
                  L_contact_gw_ij.pop(i_contact)

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
            if L_g_tempo[id_grain].center[0] < x_min - L_g_tempo[id_grain].radius:
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[0] > x_max + L_g_tempo[id_grain].radius:
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[1] < y_min - L_g_tempo[id_grain].radius:
                L_ig_to_delete.append(id_grain)
            elif L_g_tempo[id_grain].center[1] > y_max + L_g_tempo[id_grain].radius:
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write('Grain '+str(L_g_tempo[id_grain].id)+' has been deleted because it is out of the box\n')
            L_g_tempo.pop(id_grain)

        #Control the y_max to have the pressure target
        y_max, Fv = Control_y_max_NR(y_max,Forcev_target,L_contact_gw,L_g_tempo)

        #Tracker
        F = F_total(L_g_tempo)
        Ecin = E_cin_total(L_g_tempo)
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Ymax_tracker.append(y_max)

        if i_DEM%100==0:
            print('i_DEM',i_DEM,'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/Forcev_target),'%')
            Plot_Config_Loaded(L_g_tempo,x_min,x_max,y_min,y_max,i_DEM)

        #Check stop conditions for DEM
        if i_DEM >= i_DEM_stop + i_DEM_0:
             DEM_loop_statut = False
        if Ecin < Ecin_stop and F < Force_stop and (0.95*Forcev_target<Fv and Fv<1.05*Forcev_target):
              DEM_loop_statut = False
        if L_g_tempo == []:
            DEM_loop_statut = False

    return L_g_tempo, y_max, i_DEM

#-------------------------------------------------------------------------------

def Create_grains(L_g_tempo,L_n_grain_radius,L_n_grain_radius_done,L_radius,x_min,x_max,y_min,y_max,Y,nu,rho_surf,simulation_report):
    #generate the grains
    #a position is tried, then we verify this new grain does not overlap with previously created ones

    #Parameters for the method
    N_test_max = 5000
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
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,radius,Y,nu,rho_surf)
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

    return L_g_tempo, L_n_grain_radius_done

#-------------------------------------------------------------------------------

def Intersection(g1,g2):
    #verify is two grains overlap or not

    d_12 = np.linalg.norm(g1.center - g2.center)
    return d_12 < g1.radius + g2.radius

#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    #compute total kinetic energy

    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def F_total(L_g):
    #compute total force applied on grain

    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy])
    return F

#-------------------------------------------------------------------------------

def Control_y_max_NR(y_max,Force_target,L_contact_gw,L_g):
    #Control the upper wall to apply force
    #a Newton-Raphson method is applied

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
    #compute the function f to control the upper wall
    #difference between the force applied and the target value

    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    #compute the derivative function df to control the upper wall

    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g,Force):
    #the upper wall is located as a single contact verify the target value

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
    #plot loading configuration

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
    #plot final loading configuration

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

def From_LG_tempo_to_usable(L_g_tempo,x_L,y_L,w):
    #from a tempo configuration (circular grains), an initial configuration (polygonal grains) is generated

    L_g = []
    for grain_tempo in L_g_tempo:
        ci_M_IC = IC(x_L,y_L,grain_tempo.radius,w,grain_tempo.center)
        L_g.append(Grain.Grain(grain_tempo.id,ci_M_IC,None,np.array([0,0]),np.array([0,0]),grain_tempo.y,grain_tempo.nu,grain_tempo.rho_surf))
    return L_g

#-------------------------------------------------------------------------------

def IC(x_L,y_L,R,w,C):
  #create initial phase field, assuming a circular grain.

  etai_M_IC = np.zeros((len(y_L),len(x_L)))

  for y in y_L:
    for x in x_L:
      r = np.linalg.norm(np.array([x,y])-C)
      if r<R-w/2:
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 1
      elif r>R+w/2:
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 0
      else :
        etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = 0.5*(1 + np.cos(math.pi*(r-R+w/2)/w))

  return etai_M_IC

#-------------------------------------------------------------------------------
