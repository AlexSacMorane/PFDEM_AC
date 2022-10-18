# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define an initial configuration.
Grains are considered as a square.
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

  def __init__(self, ID, Center, Dimension, Theta, Y, Nu, Rho_surf):
    #defining the grain
    #each grain is described by a id (an integer class)
    #                           a center (a array class [X,Y])
    #                           a dimension (a float)
    #                           an initial rotation (a float between 0 and 2pi)
    #                           a Young modulus (a float)
    #                           a Poisson's ratio (a float)
    #                           a surface mass (a float)

    #Build the border
    n_border = 20 #number of vertices
    L_border = []
    L_border_x = []
    L_border_y = []
    L_r = []
    L_theta_r = []
    P_rot = np.array([[math.cos(Theta), -math.sin(Theta)],
                      [math.sin(Theta),  math.cos(Theta)]]) #initial rotation
    for i in range(n_border):
        angle = 2*math.pi*i/n_border #angle of the vertex
        angle_to_determine_R = angle%(math.pi/2)
        if angle_to_determine_R > math.pi/4:
            angle_to_determine_R = angle_to_determine_R - math.pi/2
        p = np.array([Dimension/2/math.cos(angle_to_determine_R),0]) #compute the radius, considering the grain as a square
        P_rot_p = np.array([[math.cos(angle), -math.sin(angle)],
                            [math.sin(angle),  math.cos(angle)]])
        p = np.dot(P_rot_p,p)
        p = np.dot(P_rot,p)
        p = p + np.array(Center)
        L_border.append(p)
        L_border_x.append(p[0])
        L_border_y.append(p[1])
        L_r.append(Dimension/2/math.cos(angle_to_determine_R))
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
    self.id = ID
    self.center = np.array(Center)
    self.dimension = Dimension
    self.theta = Theta
    self.l_border = L_border
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.l_theta_r = L_theta_r
    self.l_r = L_r
    self.r_min = min(L_r)
    self.r_max = max(L_r)
    self.y = Y
    self.nu = Nu
    self.g = Y /2/(1+Nu) #shear modulus
    self.rho_surf = Rho_surf
    self.mass = Dimension**2*Rho_surf
    self.inertia = 1/6*self.mass*Dimension**2
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.w = 0

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      #add a force (an array [Fx,Fy]) to the grain

      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]
      v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
      v2 = np.array([F[0], F[1], 0])
      self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      #initialize the force applied to the grain
      #a gravity of g is applied

      self.fx = 0
      self.fy = -g*self.mass
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM):
    #move the grain following a semi implicit euler scheme

    #translation
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    if np.linalg.norm(self.v) > self.dimension*0.01/dt_DEM: #limitation of the speed
        self.v = self.v * self.dimension*0.01/dt_DEM/np.linalg.norm(self.v)
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
    #rigib body rotation
    for i in range(len(self.l_border)):
        p = self.l_border[i] - self.center
        Rot_Matrix = np.array([[math.cos(self.w*dt_DEM), -math.sin(self.w*dt_DEM)],
                               [math.sin(self.w*dt_DEM),  math.cos(self.w*dt_DEM)]])
        p = np.dot(Rot_Matrix,p)
        self.l_border[i] = p + self.center
        self.l_border_x[i] = p[0] + self.center[0]
        self.l_border_y[i] = p[1] + self.center[1]

#-------------------------------------------------------------------------------

class Contact_Tempo:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, Mu, Coeff_Restitution):
    #defining the contact grain-grain
    #each contact is described by a id (an integer class)
    #                             two grains (a grain class)
    #                             a friction coefficient (a float)
    #                             a restitution coefficient (a float) : 1 = restitution of all energy, 0 = reestituion of any energy

    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = Mu
    self.coeff_restitution = Coeff_Restitution
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def normal(self):
    #compute the normal reaction of a contact grain-grain
    #Here a pontual spring is considered

    #looking for the nearest nodes
    d_virtual = max(self.g1.r_max,self.g2.r_max)
    ij_min = [0,0]
    d_ij_min = 100*d_virtual #Large
    for i in range(len(self.g1.l_border[:-1])):
        for j in range(len(self.g2.l_border[:-1])):
            d_ij = np.linalg.norm(self.g2.l_border[:-1][j]-self.g1.l_border[:-1][i]+d_virtual*(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]
    self.ij_min = ij_min

    #-----------------------------------------------------------------------------
    #Computing CP
    #-----------------------------------------------------------------------------

    M = (self.g1.l_border[:-1][ij_min[0]]+self.g2.l_border[:-1][ij_min[1]])/2
    # 5 candidates for CP
    N = np.array([self.g1.l_border[:-1][ij_min[0]][0] - self.g2.l_border[:-1][ij_min[1]][0],
                  self.g1.l_border[:-1][ij_min[0]][1] - self.g2.l_border[:-1][ij_min[1]][1]])
    N = N/np.linalg.norm(N)
    PB = np.array([-N[1] ,N[0]])
    PB = PB/np.linalg.norm(PB)

    #candidats from grain 1
    if ij_min[0] <len(self.g1.l_border[:-1]) - 1:
        M1 = self.g1.l_border[:-1][ij_min[0]+1]-self.g1.l_border[:-1][ij_min[0]]
    else :
        M1 = self.g1.l_border[:-1][0]-self.g1.l_border[:-1][ij_min[0]]
    M1 = M1/np.linalg.norm(M1)
    M3 = self.g1.l_border[:-1][ij_min[0]-1]-self.g1.l_border[:-1][ij_min[0]]
    M3 = M3/np.linalg.norm(M3)
    #reorganize the candidats
    if np.dot(M1,PB) < 0:
        Mtempo = M1.copy()
        M1 = M3.copy()
        M3 = Mtempo.copy()

    #candidats from grain 2
    if ij_min[1] <len(self.g2.l_border[:-1]) - 1:
        M2 = self.g2.l_border[:-1][ij_min[1]+1]-self.g2.l_border[:-1][ij_min[1]]
    else :
        M2 = self.g2.l_border[:-1][0]-self.g2.l_border[:-1][ij_min[1]]
    M2 = M2/np.linalg.norm(M2)
    M4 = self.g2.l_border[:-1][ij_min[1]-1]-self.g2.l_border[:-1][ij_min[1]]
    M4 = M4/np.linalg.norm(M4)
    #reorganize the candidats
    if np.dot(M2,PB) < 0:
      Mtempo = M2.copy()
      M2 = M4.copy()
      M4 = Mtempo.copy()

    #compute the different angles
    theta_PB = math.pi/2
    theta_M1 =  math.acos(np.dot(M1,N))
    theta_M2 =  math.acos(np.dot(M2,N))
    theta_M3 = -math.acos(np.dot(M3,N))
    theta_M4 = -math.acos(np.dot(M4,N))

    #find the PC
    if theta_M2 < theta_PB and theta_PB < theta_M1\
       and theta_M3 < -theta_PB and -theta_PB < theta_M4:
       PC = PB
    else:
      L_Mi = [M1,M2,M3,M4]
      L_theta_Mi_PB=[theta_M1-theta_PB, theta_PB-theta_M2, -theta_M3-theta_PB, theta_PB+theta_M4]
      PC = L_Mi[L_theta_Mi_PB.index(min(L_theta_Mi_PB))]

    #-----------------------------------------------------------------------------
    # Compute the normal and tangential planes
    #-----------------------------------------------------------------------------

    PC_normal = np.array([PC[1],-PC[0]])
    if np.dot(PC_normal,(self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center))<0 :
        PC_normal = np.array([-PC[1],PC[0]])
    self.pc_normal = PC_normal #n12
    self.pc_tangential = np.array([-PC_normal[1],PC_normal[0]])

    #-----------------------------------------------------------------------------
    # Compute the overlap
    #-----------------------------------------------------------------------------

    d_b = np.dot(M-self.g2.l_border[:-1][ij_min[1]],PC_normal)
    d_a = np.dot(M-self.g1.l_border[:-1][ij_min[0]],PC_normal)
    overlap = d_b - d_a
    self.overlap_normal = overlap

    if overlap > 0:
    #-----------------------------------------------------------------------------
    # Compute the reaction
    #-----------------------------------------------------------------------------

        #Spring term
        Y_eq = 1/((1-self.g1.nu*self.g1.nu)/self.g1.y+(1-self.g2.nu*self.g2.nu)/self.g2.y)
        R_eq = 1/(1/np.linalg.norm(self.g1.center-self.g1.l_border[self.ij_min[0]])+1/np.linalg.norm(self.g2.center-self.g2.l_border[self.ij_min[1]]))
        k = 4/3*Y_eq*math.sqrt(R_eq)
        F_2_1_n = -k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.Ep_n = 2/5 * k * overlap**(5/2) #-dEp/dx = F_2_1_n
        self.g1.add_F( F_2_1, self.g1.l_border[:-1][ij_min[0]])
        self.g2.add_F(-F_2_1, self.g2.l_border[:-1][ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.add_F( F_2_1_damp, self.g1.l_border[:-1][ij_min[0]])
        self.g2.add_F(-F_2_1_damp, self.g2.l_border[:-1][ij_min[1]])

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0
        self.Ep_n = 0

#-------------------------------------------------------------------------------

  def tangential(self,dt_DEM):
    #compute the tangential reaction of a contact grain-grain
    #Here a pontual spring is considered

    if self.overlap_normal > 0 and self.mu > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/np.linalg.norm(self.g1.center-self.g1.l_border[self.ij_min[0]])+1/np.linalg.norm(self.g2.center-self.g2.l_border[self.ij_min[1]]))
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap_normal))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F_2_1_n),0)) #not linear spring

        r1 = np.linalg.norm(self.g1.l_border[:-1][self.ij_min[0]] - self.g1.center) - self.overlap_normal/2
        r2 = np.linalg.norm(self.g2.l_border[:-1][self.ij_min[1]] - self.g2.center) - self.overlap_normal/2
        Delta_Us = (np.dot(self.g1.v-self.g2.v,self.pc_tangential) + r1*self.g1.w + r2*self.g2.w)*dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - kt*Delta_Us
        self.tangential_old = self.pc_tangential
        if abs(self.ft) > abs(self.mu*self.F_2_1_n) or kt == 0: #Coulomb criteria
            self.ft = self.mu * abs(self.F_2_1_n) * np.sign(self.ft)

        self.g1.add_F( self.ft*self.pc_tangential, self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.add_F(-self.ft*self.pc_tangential, self.g2.l_border[:-1][self.ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.ft_damp = F_2_1_damp_t
        self.g1.add_F( F_2_1_damp, self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.add_F(-F_2_1_damp, self.g2.l_border[:-1][self.ij_min[1]])

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------

class Contact_gw_Tempo:

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
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.r_max) #Hertz law
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

  def  normal(self):
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
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])
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
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

    elif self.nature == 'gwx_max':
        #unlinear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,nwg)*eta
        Fwg_damp = Fwg_damp_n*nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
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
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(min(self.g.l_border_y))])

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
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_y.index(max(self.g.l_border_y))])

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
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(min(self.g.l_border_x))])

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
       self.g.add_F(Fwg, self.g.l_border[self.g.l_border_x.index(max(self.g.l_border_x))])

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def LG_tempo(x_min,x_max,y_min,N_grain,L_dimension,L_percentage_dimension,rho_surf,Y,nu,mu_gg,e_gg,mu_gw,e_gw,Force_target,gravity,dt_DEM,simulation_report):
    #create an initial condition

    #number of grains generation
    n_generation = 2 #Work only for 2

    #define the y_max for the grains generation
    factor = 1.2
    dimension_mean = 0
    for i in range(len(L_dimension)):
        dimension_mean = dimension_mean + L_dimension[i]*L_percentage_dimension[i]
    dy_creation = N_grain/n_generation * factor*(dimension_mean)**2/(x_max-x_min)

    #plan the grains generation
    L_n_grain_dimension_try_one = []
    L_n_grain_dimension = []
    L_n_grain_dimension_done = []
    for percentage in L_percentage_dimension:
        L_n_grain_dimension_try_one.append(int(N_grain*percentage/n_generation))
        L_n_grain_dimension.append(int(N_grain*percentage))
        L_n_grain_dimension_done.append(0)

    #Parameters for the DEM
    i_DEM_stop = 100000
    i_DEM = 0

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    simulation_report.write('Creation of the grains\n')
    L_L_g_tempo = []
    y_min_init = y_min #save for rebuild

    print('First generation of grains')
    L_g_tempo = []
    L_g_tempo, L_n_grain_dimension_done = Create_grains(L_g_tempo,L_n_grain_dimension_try_one,L_n_grain_dimension_done,L_dimension,x_min,x_max,y_min,y_min+dy_creation,Y,nu,rho_surf,simulation_report)
    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = y_min
    for grain in L_g_tempo:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)
    L_g_tempo, y_min, i_DEM = DEM_loading(L_g_tempo, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min, y_max, dt_DEM, Force_target, gravity, i_DEM_stop, i_DEM, simulation_report)
    L_L_g_tempo.append(L_g_tempo.copy())

    print('Second generation of grains')
    L_g_tempo = []
    L_g_tempo, L_n_grain_dimension_done = Create_grains(L_g_tempo,L_n_grain_dimension,L_n_grain_dimension_done,L_dimension,x_min,x_max,y_min,y_min+dy_creation,Y,nu,rho_surf,simulation_report)
    #DEM to find the steady-state configuration after loading
    #find the maximum y (center+radius)
    y_max = y_min
    for grain in L_g_tempo:
        if max(grain.l_border_y) > y_max:
            y_max = max(grain.l_border_y)
    L_g_tempo, y_min, i_DEM = DEM_loading(L_g_tempo, mu_gg, e_gg, mu_gw, e_gw, x_min, x_max, y_min, y_max, dt_DEM, Force_target, gravity, i_DEM_stop, i_DEM, simulation_report)
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
        Force_stop = Force_stop + 0.5*grain.mass*gravity
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(0.0001*grain.r_max/dt_DEM)**2

    while DEM_loop_statut :

        i_DEM = i_DEM + 1

        #Contact detection
        for i_grain in range(len(L_g_tempo)-1):
              for j_grain in range(i_grain+1,len(L_g_tempo)):
                  #contact grain-grain
                  if Grains_Polyhedral_contact_f(L_g_tempo[i_grain],L_g_tempo[j_grain]) and (i_grain,j_grain) not in L_contact_ij:
                      L_contact_gg.append(Contact_Tempo(id_contact, L_g_tempo[i_grain], L_g_tempo[j_grain], mu_gg, e_gg))
                      id_contact = id_contact + 1
                      L_contact_ij.append((i_grain,j_grain))
                  elif not Grains_Polyhedral_contact_f(L_g_tempo[i_grain],L_g_tempo[j_grain]) and (i_grain,j_grain) in L_contact_ij:
                      i_contact = L_contact_ij.index((i_grain,j_grain))
                      L_contact_gg.pop(i_contact)
                      L_contact_ij.pop(i_contact)

        #Contact detection with walls
        for i_grain in range(len(L_g_tempo)):
          #coordinate extremum of a grain
          p_x_min = min(L_g_tempo[i_grain].l_border_x)
          p_x_max = max(L_g_tempo[i_grain].l_border_x)
          p_y_min = min(L_g_tempo[i_grain].l_border_y)
          p_y_max = max(L_g_tempo[i_grain].l_border_y)

          # contact grain-wall x_min
          if p_x_min < x_min and (i_grain,-1) not in L_contact_gw_ij:
              overlap = x_min - p_x_min
              L_contact_gw.append(Contact_gw_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwx_min', x_min, overlap))
              id_contact = id_contact + 1
              L_contact_gw_ij.append((i_grain,-1))
          elif p_x_min < x_min and (i_grain,-1) in L_contact_gw_ij:
              overlap = x_min - p_x_min
              L_contact_gw[L_contact_gw_ij.index((i_grain,-1))].update_overlap(overlap)
          elif p_x_min > x_min and (i_grain,-1) in L_contact_gw_ij:
              i_contact = L_contact_gw_ij.index((i_grain,-1))
              L_contact_gw.pop(i_contact)
              L_contact_gw_ij.pop(i_contact)
          #grain-wall x_max
          if p_x_max > x_max and (i_grain,-2) not in L_contact_gw_ij:
              overlap = p_x_max - x_max
              L_contact_gw.append(Contact_gw_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwx_max', x_max, overlap))
              L_contact_gw_ij.append((i_grain,-2))
              id_contact = id_contact + 1
          elif p_x_max > x_max and (i_grain,-2) in L_contact_gw_ij:
              overlap = p_x_max - x_max
              L_contact_gw[L_contact_gw_ij.index((i_grain,-2))].update_overlap(overlap)
          elif p_x_max < x_max and (i_grain,-2) in L_contact_gw_ij:
              i_contact = L_contact_gw_ij.index((i_grain,-2))
              L_contact_gw.pop(i_contact)
              L_contact_gw_ij.pop(i_contact)
          # contact grain-wall y_min
          if p_y_min < y_min and (i_grain,-3) not in L_contact_gw_ij:
              overlap = y_min - p_y_min
              L_contact_gw.append(Contact_gw_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwy_min', y_min, overlap))
              id_contact = id_contact + 1
              L_contact_gw_ij.append((i_grain,-3))
          elif p_y_min < y_min and (i_grain,-3) in L_contact_gw_ij:
              overlap = y_min - p_y_min
              L_contact_gw[L_contact_gw_ij.index((i_grain,-3))].update_overlap(overlap)
          elif p_y_min > y_min and (i_grain,-3) in L_contact_gw_ij:
              i_contact = L_contact_gw_ij.index((i_grain,-3))
              L_contact_gw.pop(i_contact)
              L_contact_gw_ij.pop(i_contact)
          # contact grain-wall y_max
          if p_y_max > y_max and (i_grain,-4) not in L_contact_gw_ij:
              overlap = p_y_max - y_max
              L_contact_gw.append(Contact_gw_Tempo(id_contact, L_g_tempo[i_grain], mu_gw, e_gw, 'gwy_max', y_max, overlap))
              id_contact = id_contact + 1
              L_contact_gw_ij.append((i_grain,-4))
          elif p_y_max > y_max and (i_grain,-4) in L_contact_gw_ij:
              overlap = p_y_max - y_max
              L_contact_gw[L_contact_gw_ij.index((i_grain,-4))].update_overlap(overlap)
          elif p_y_max < y_max and (i_grain,-4) in L_contact_gw_ij:
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
            #coordinate extremum of a grain
            p_x_min = min(L_g_tempo[i_grain].l_border_x)
            p_x_max = max(L_g_tempo[i_grain].l_border_x)
            p_y_min = min(L_g_tempo[i_grain].l_border_y)
            p_y_max = max(L_g_tempo[i_grain].l_border_y)
            if p_x_max < x_min :
                L_ig_to_delete.append(id_grain)
            elif p_x_min > x_max :
                L_ig_to_delete.append(id_grain)
            elif p_y_max < y_min :
                L_ig_to_delete.append(id_grain)
            elif p_y_min > y_max :
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

        if i_DEM%50==0:
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

def Create_grains(L_g_tempo,L_n_grain_dimension,L_n_grain_dimension_done,L_dimension,x_min,x_max,y_min,y_max,Y,nu,rho_surf,simulation_report):
    #generate the grains
    #a position is tried, then we verify this new grain does not overlap with previously created ones

    #Parameters for the method
    N_test_max = 5000
    n_not_created = 0

    for i in range(len(L_dimension)):
        dimension = L_dimension[i]
        n_grain = L_n_grain_dimension[i]
        n_grain_done = L_n_grain_dimension_done[i]
        last_id_grain_created = np.sum(L_n_grain_dimension_done)
        for id_grain in range(last_id_grain_created, last_id_grain_created + n_grain - n_grain_done):
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < N_test_max:
                i_test = i_test + 1
                center = np.array([random.uniform(x_min+1.1*dimension/2,x_max-1.1*dimension/2),random.uniform(y_min+1.1*dimension/2,y_max)])
                g_tempo = Grain_Tempo(id_grain-n_not_created,center,dimension,random.uniform(0,math.pi/2),Y,nu,rho_surf)
                grain_created = True
                for grain in L_g_tempo:
                    if Grains_Polyhedral_contact_f(g_tempo,grain):
                        grain_created = False
            if i_test == N_test_max and not grain_created:
                n_not_created = n_not_created + 1
                simulation_report.write('Grain '+str(id_grain)+' has not been created after '+str(i_test)+' tries\n')
            else :
                L_g_tempo.append(g_tempo)
                L_n_grain_dimension_done[i] = L_n_grain_dimension_done[i] + 1

    return L_g_tempo, L_n_grain_dimension_done

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_f(g1,g2):
  #detect contact grain-grain

  if np.linalg.norm(g1.center-g2.center) < 1.5*(g1.r_max+g2.r_max):

      #looking for the nearest nodes
      d_virtual = max(g1.r_max,g2.r_max)
      ij_min = [0,0]
      d_ij_min = 100*d_virtual #Large
      for i in range(len(g1.l_border[:-1])):
        for j in range(len(g2.l_border[:-1])):
            d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
            if d_ij < d_ij_min :
                d_ij_min = d_ij
                ij_min = [i,j]

      d_ij_min = np.dot(g2.l_border[:-1][ij_min[1]]-g1.l_border[:-1][ij_min[0]],-(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
      return d_ij_min > 0

  else:
    return False

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
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].r_max)
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

    pass

#-------------------------------------------------------------------------------

def IC(x_L,y_L,R,w,C):
  #create initial phase field, assuming a circular grain.

  pass

#-------------------------------------------------------------------------------
