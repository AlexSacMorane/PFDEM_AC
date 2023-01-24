# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be


The goal of this file is to define a new class.
The new class is about the grain during initial condition generation with perfect disks.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain_Tempo:
  """
  A temporary grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Radius, Dissolvable, Type, dict_material):
    """Defining the grain.

        Input :
            itself (a grain_tempo)
            an id (a int)
            a center coordinate (a 1 x 2 numpy array)
            a radius (a float)
            a material dictionnary (a dict)
            a grain type, disk or square (a float)
        Output :
            Nothing, but a temporary grain is generated (a grain_tempo)
    """
    n_border = 90 #number of vertices
    L_border = []
    L_border_x = []
    L_border_y = []
    #Build the border (for print)
    for i in range(n_border):
        theta = 2*math.pi*i/n_border
        p = np.array(Center) + np.array([Radius*math.cos(theta),Radius*math.sin(theta)])
        L_border.append(p)
        L_border_x.append(p[0])
        L_border_y.append(p[1])
    L_border.append(L_border[0])
    L_border_x.append(L_border_x[0])
    L_border_y.append(L_border_y[0])
    #save
    self.id = ID
    self.radius = Radius
    self.dissolved = Dissolvable
    self.type = Type
    self.theta = 0
    if Dissolvable :
        self.rho_surf = dict_material['rho_surf_dissolvable']
    else :
        self.rho_surf = dict_material['rho_surf_undissolvable']
    self.surface = math.pi*Radius**2
    self.mass = self.rho_surf*self.surface
    self.inertia = self.mass*Radius**2
    self.center = np.array(Center)
    self.l_border = L_border
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.y = dict_material['Y']
    self.nu = dict_material['nu']
    self.g = dict_material['Y']/2/(1+dict_material['nu']) #shear modulus
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.w = 0

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      """
      Add a force to the grain.

        Input :
            itself (a grain_tempo)
            a force applied (a 1 x 2 numpy array)
            a application point (a 1 x 2 numpy array)
        Output :
            Nothing, but attributes are updated (three floats)
      """
      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]
      v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
      v2 = np.array([F[0], F[1], 0])
      self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      """
      Initialize the force applied to the grain.

      A gravity is assumed.

        Input :
            itself (a grain_tempo)
            a gravity value (a float)
        Output :
            Nothing, but attributes concerning the force applied are initialized (three floats)
      """
      self.fx = 0
      self.fy = -g*self.mass
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM,factor):
    """
    Move the grain following a semi implicit euler scheme.

        Input :
            itself (a grain_tempo)
            a time step (a float)
            a factor to limite the displacement (a float)
        Output :
            Nothing, but the grain is moved
    """
    #translation
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    if np.linalg.norm(self.v) > self.radius*factor/dt_DEM: #limitation of the speed
        self.v = self.v * self.radius*factor/dt_DEM/np.linalg.norm(self.v)
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + self.v*dt_DEM
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*dt_DEM
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*dt_DEM
    self.center = self.center + self.v*dt_DEM

    #rotation
    dw_i = self.mz/self.inertia
    self.w = self.w + dw_i*dt_DEM
    self.theta = self.theta + self.w*dt_DEM
