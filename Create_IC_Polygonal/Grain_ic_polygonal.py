# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the polygonal grains during initial condition generation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import random
import numpy as np

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain_Tempo_Polygonal:
  """
  A temporary grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, grain_sphere, n_border):
    """Defining the grain.

        Input :
            a grain tempo (a tempo grain)
            a discretization of the grain (an int)
        Output :
            a grain tempo polygonal (a tempo polygonal grain)
    """
    L_border = []
    L_border_x = []
    L_border_y = []
    L_r = []
    L_theta_r = []
    if grain_sphere.type == 'Disk':
        #Build the border
        for i in range(n_border):
            theta = 2*math.pi*i/n_border
            p = np.array(grain_sphere.center) + np.array([grain_sphere.radius*math.cos(theta),grain_sphere.radius*math.sin(theta)])
            L_border.append(p)
            L_border_x.append(p[0])
            L_border_y.append(p[1])
            L_r.append(grain_sphere.radius)
            L_theta_r.append(theta)
        L_border.append(L_border[0])
        L_border_x.append(L_border_x[0])
        L_border_y.append(L_border_y[0])
        #save
        self.radius = grain_sphere.radius
        self.surface = math.pi*grain_sphere.radius**2
        self.mass = math.pi*grain_sphere.radius**2*grain_sphere.rho_surf
        self.inertia = self.mass*grain_sphere.radius**2
    elif grain_sphere.type == 'Square':
        Lenght = 2*grain_sphere.radius/math.sqrt(2)
        #initial random rotation
        Theta = random.uniform(0,math.pi/2)
        #Build the border
        P_rot = np.array([[math.cos(Theta), -math.sin(Theta)],
                          [math.sin(Theta),  math.cos(Theta)]]) #initial rotation
        #Build the border
        for i in range(n_border):
            angle = 2*math.pi*i/n_border #angle of the vertex
            angle_to_determine_R = angle%(math.pi/2)
            if angle_to_determine_R > math.pi/4:
                angle_to_determine_R = angle_to_determine_R - math.pi/2
            p = np.array([Lenght/2/math.cos(angle_to_determine_R),0]) #compute the radius, considering the grain as a square
            P_rot_p = np.array([[math.cos(angle), -math.sin(angle)],
                                [math.sin(angle),  math.cos(angle)]])
            p = np.dot(P_rot_p,p)
            p = np.dot(P_rot,p)
            p = p + np.array(grain_sphere.center)
            L_border.append(p)
            L_border_x.append(p[0])
            L_border_y.append(p[1])
            L_r.append(Lenght/2/math.cos(angle_to_determine_R))
            angle_r = Theta + angle
            while angle_r >= 2*math.pi:
                angle_r = angle_r - 2*math.pi
            while angle_r < 0:
                angle_r = angle_r + 2*math.pi
            L_theta_r.append(angle_r)
        L_border.append(L_border[0])
        L_border_x.append(L_border_x[0])
        L_border_y.append(L_border_y[0])
        #save
        self.radius = grain_sphere.radius
        self.surface = Lenght**2
        self.mass = Lenght**2*grain_sphere.rho_surf
        self.inertia = 1/6*self.mass*Lenght**2
    #common
    self.id = grain_sphere.id
    self.center = np.array(grain_sphere.center)
    self.dissolved = grain_sphere.dissolved
    self.type = grain_sphere.type
    self.l_border = L_border
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.l_r = L_r
    self.l_theta_r = L_theta_r
    self.r_min = min(L_r)
    self.r_max = max(L_r)
    self.rho_surf = grain_sphere.rho_surf
    self.y = grain_sphere.y
    self.nu = grain_sphere.nu
    self.g = grain_sphere.g #shear modulus
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.theta = 0
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

    for i_theta_r in range(len(self.l_theta_r)) :
        theta_r = self.l_theta_r[i_theta_r]
        theta_r = theta_r + self.w*dt_DEM
        while theta_r >= 2*math.pi:
            theta_r = theta_r - 2*math.pi
        while theta_r < 0 :
            theta_r = theta_r + 2*math.pi
        self.l_theta_r[i_theta_r] = theta_r

    for i in range(len(self.l_border)):
        p = self.l_border[i] - self.center
        Rot_Matrix = np.array([[math.cos(self.w*dt_DEM), -math.sin(self.w*dt_DEM)],
                               [math.sin(self.w*dt_DEM),  math.cos(self.w*dt_DEM)]])
        p = np.dot(Rot_Matrix,p)
        self.l_border[i] = p + self.center
        self.l_border_x[i] = p[0] + self.center[0]
        self.l_border_y[i] = p[1] + self.center[1]
