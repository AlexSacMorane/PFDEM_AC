# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between a wall and a grain
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math
import matplotlib.pyplot as plt
import random

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gw:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, Mu, Coeff_Restitution, Nature, Limit, Overlap):
    #defining the contact grain-wall
    #each contact is described by a id (an integer class)
    #                             a grain (a grain class)
    #                             a friction coefficient (a float)
    #                             a restitution coefficient (a float) : 1 = restitution of all energy, 0 = reestituion of any energy
    #                             a nature description (a float) : identify the wall
    #                             a coordinate of the wall (a float) : x or y depending of the wall nature
    #                             an overlap (a float)

    self.id = ID
    self.g = G
    factor = 5 #factor just to increase the stiffness
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.r_mean) #Hertz law
    self.kt = 0
    self.ft = 0
    self.limit = Limit
    self.nature = Nature
    self.mu = Mu
    self.coeff_restitution = Coeff_Restitution
    self.overlap = Overlap
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def update_overlap(self,new_overlap):
    #update the overlap of a contact already created.

    self.overlap = new_overlap

#-------------------------------------------------------------------------------

  def init_contact_gw(self,L_g):
    #initialize the contact with updating the grain,
    #                            putting at 0 the tangential reaction
    #                            saying the boolean at False (new contact grain-w&all)

    self.g = L_g[self.g.id]
    self.ft = 0
    self.tangential_old_statut = False

#-------------------------------------------------------------------------------

  def  DEM_gw_Polyhedral_normal(self):
    #compute the normal reaction of a contact grain-wall
    #Here a pontual spring is considered
    #conditions "if" are defined and same for each wall nature

    if self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.m
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])
        #damping
        Fwg_damp_n = 0
        self.Fwg_damp_n = Fwg_damp_n

    elif self.nature == 'gwx_min':
        #unlinear stiffness
        nwg = np.array([1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.m
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

    elif self.nature == 'gwx_max':
        #unlinear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.m
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.update_f(Fwg_damp[0],Fwg_damp[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------

  def DEM_gw_Polyhedral_tangential(self, dt_DEM):
   #compute the tangential reaction of a contact grain-wall
   #Here a pontual spring is considered
   #conditions "if" are defined and same for each wall nature

   if self.nature == 'gwy_min':
       #unlinear stiffness
       twg = np.array([-1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(min(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

   elif self.nature == 'gwy_max':
       #unlinear stiffness
       twg = np.array([1, 0])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_y.index(max(self.g.l_border_y))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])

   elif self.nature == 'gwx_min':
       #unlinear stiffness
       twg = np.array([0, 1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(min(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

   elif self.nature == 'gwx_max':
       #linear stiffness
       twg = np.array([0, -1])
       self.twg = twg
       r = np.linalg.norm(self.g.l_border[:-1][self.g.l_border_x.index(max(self.g.l_border_x))] - self.g.center) - self.overlap
       Delta_Us = (np.dot(self.g.v,self.twg) - r*self.g.w) * dt_DEM
       self.overlap_tangential = self.overlap_tangential + Delta_Us
       self.ft = self.ft - self.kt*Delta_Us
       if abs(self.ft) > abs(self.mu*self.Fwg_n) :
           self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
       Fwg = self.ft*twg
       self.g.update_f(Fwg[0],Fwg[1],self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------

  def DEM_gw_Polyhedral_normal_surface(self,simulation_report):
   #compute the tangential reaction of a contact grain-wall
   #Here a surface spring is considered
   #conditions "if" are defined and same for each wall nature
   pass

#-------------------------------------------------------------------------------

  def DEM_gw_Polyhedral_tangential_surface(self,simulation_report):
   #compute the tangential reaction of a contact grain-wall
   #Here a surface spring is considered
   #conditions "if" are defined and same for each wall nature
   pass

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Grains_Polyhedral_Wall_contact(L_g,L_contact_gw,L_contact_gw_ij,id_contact_gw,x_min,x_max,y_min,y_max,mu,coeff_restitution):
  #detect contact grain-wall
  #we realize iterations on the grain list and compare with the coordinate of the different walls

  for grain in L_g:

      p_x_min = min(grain.l_border_x)
      p_x_max = max(grain.l_border_x)
      p_y_min = min(grain.l_border_y)
      p_y_max = max(grain.l_border_y)

      #grain-wall x_min
      if p_x_min < x_min and (grain.id,-1) not in L_contact_gw_ij:
          overlap = x_min - p_x_min
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwx_min', x_min, overlap))
          L_contact_gw_ij.append((grain.id,-1))
          id_contact_gw = id_contact_gw + 1
      elif p_x_min < x_min and (grain.id,-1) in L_contact_gw_ij:
          overlap = x_min - p_x_min
          L_contact_gw[L_contact_gw_ij.index((grain.id,-1))].update_overlap(overlap)
      elif p_x_min > x_min and (grain.id,-1) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall x_max
      if p_x_max > x_max and (grain.id,-2) not in L_contact_gw_ij:
          overlap = p_x_max - x_max
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwx_max', x_max, overlap))
          L_contact_gw_ij.append((grain.id,-2))
          id_contact_gw = id_contact_gw + 1
      elif p_x_max > x_max and (grain.id,-2) in L_contact_gw_ij:
          overlap = p_x_max - x_max
          L_contact_gw[L_contact_gw_ij.index((grain.id,-2))].update_overlap(overlap)
      elif p_x_max < x_max and (grain.id,-2) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall y_min
      if p_y_min < y_min and (grain.id,-3) not in L_contact_gw_ij:
          overlap = y_min - p_y_min
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwy_min', y_min, overlap))
          L_contact_gw_ij.append((grain.id,-3))
          id_contact_gw = id_contact_gw + 1
      elif p_y_min < y_min and (grain.id,-3) in L_contact_gw_ij:
          overlap = y_min - p_y_min
          L_contact_gw[L_contact_gw_ij.index((grain.id,-3))].update_overlap(overlap)
      elif p_y_min > y_min and (grain.id,-3) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall y_max
      if p_y_max > y_max and (grain.id,-4) not in L_contact_gw_ij:
          overlap = p_y_max - y_max
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwy_max', y_max, overlap))
          L_contact_gw_ij.append((grain.id,-4))
          id_contact_gw = id_contact_gw + 1
      elif p_y_max > y_max and (grain.id,-4) in L_contact_gw_ij:
          overlap = p_y_max - y_max
          L_contact_gw[L_contact_gw_ij.index((grain.id,-4))].update_overlap(overlap)
      elif p_y_max < y_max and (grain.id,-4) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-4))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)

  return L_contact_gw, L_contact_gw_ij, id_contact_gw

#-------------------------------------------------------------------------------

def Grains_Polyhedral_Wall_contact_Neighbourood(wall_neighbourood,L_contact_gw,L_contact_gw_ij,id_contact_gw,x_min,x_max,y_min,y_max,mu,coeff_restitution):
  #detect contact grain in the neighbourood of the wall and  the wall
  #the neighbourood is updated with Update_wall_Neighbouroods()
  #we realize iterations on the grain list and compare with the coordinate of the different walls

  for grain in wall_neighbourood:

      p_x_min = min(grain.l_border_x)
      p_x_max = max(grain.l_border_x)
      p_y_min = min(grain.l_border_y)
      p_y_max = max(grain.l_border_y)

      #grain-wall x_min
      if p_x_min < x_min and (grain.id,-1) not in L_contact_gw_ij:
          overlap = x_min - p_x_min
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwx_min', x_min, overlap))
          L_contact_gw_ij.append((grain.id,-1))
          id_contact_gw = id_contact_gw + 1
      elif p_x_min < x_min and (grain.id,-1) in L_contact_gw_ij:
          overlap = x_min - p_x_min
          L_contact_gw[L_contact_gw_ij.index((grain.id,-1))].update_overlap(overlap)
      elif p_x_min > x_min and (grain.id,-1) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall x_max
      if p_x_max > x_max and (grain.id,-2) not in L_contact_gw_ij:
          overlap = p_x_max - x_max
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwx_max', x_max, overlap))
          L_contact_gw_ij.append((grain.id,-2))
          id_contact_gw = id_contact_gw + 1
      elif p_x_max > x_max and (grain.id,-2) in L_contact_gw_ij:
          overlap = p_x_max - x_max
          L_contact_gw[L_contact_gw_ij.index((grain.id,-2))].update_overlap(overlap)
      elif p_x_max < x_max and (grain.id,-2) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall y_min
      if p_y_min < y_min and (grain.id,-3) not in L_contact_gw_ij:
          overlap = y_min - p_y_min
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwy_min', y_min, overlap))
          L_contact_gw_ij.append((grain.id,-3))
          id_contact_gw = id_contact_gw + 1
      elif p_y_min < y_min and (grain.id,-3) in L_contact_gw_ij:
          overlap = y_min - p_y_min
          L_contact_gw[L_contact_gw_ij.index((grain.id,-3))].update_overlap(overlap)
      elif p_y_min > y_min and (grain.id,-3) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)
      #grain-wall y_max
      if p_y_max > y_max and (grain.id,-4) not in L_contact_gw_ij:
          overlap = p_y_max - y_max
          L_contact_gw.append(Contact_gw(id_contact_gw, grain, mu, coeff_restitution, 'gwy_max', y_max, overlap))
          L_contact_gw_ij.append((grain.id,-4))
          id_contact_gw = id_contact_gw + 1
      elif p_y_max > y_max and (grain.id,-4) in L_contact_gw_ij:
          overlap = p_y_max - y_max
          L_contact_gw[L_contact_gw_ij.index((grain.id,-4))].update_overlap(overlap)
      elif p_y_max < y_max and (grain.id,-4) in L_contact_gw_ij:
          i_contact = L_contact_gw_ij.index((grain.id,-4))
          L_contact_gw.pop(i_contact)
          L_contact_gw_ij.pop(i_contact)

  return L_contact_gw, L_contact_gw_ij, id_contact_gw

#-------------------------------------------------------------------------------

def Update_wall_Neighbouroods(L_g,factor,x_min,x_max,y_min,y_max):
    #determine a neighbouroods for wall. This function is called every x time step
    #grain_wall contact is determined by Grains_Polyhedral_Wall_contact_Neighbourood
    #factor determines the size of the neighbourood window

    wall_neighbourood = []
    for grain in L_g:

        p_x_min = min(grain.l_border_x)
        p_x_max = max(grain.l_border_x)
        p_y_min = min(grain.l_border_y)
        p_y_max = max(grain.l_border_y)

        #grain-wall x_min
        if abs(p_x_min-x_min) < factor*grain.r_max :
            wall_neighbourood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-x_max) < factor*grain.r_max :
            wall_neighbourood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-y_min) < factor*grain.r_max :
            wall_neighbourood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-y_max) < factor*grain.r_max :
            wall_neighbourood.append(grain)

    return wall_neighbourood
