# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between two grains
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact grain - grain)
            an id (a int)
            two grains (two grains)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = dict_material['mu_friction_gg']
    self.coeff_restitution = dict_material['coeff_restitution']
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def init_contact(self,L_g):
    """
    Initialize the contact.

    The grains are updated, the tangetial reaction is set to 0 and a boolean attribute is set to False (new contact grain - grain).

        Input :
            itself (a contact grain - grain)
            a list of grains (a list)
        Output :
            Nothing, but attributes are updated (two grains, two floats and a boolean)
    """
    self.g1 = L_g[self.g1.id]
    self.g2 = L_g[self.g2.id]
    self.ft = 0
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def DEM_2grains_Polyhedral_normal(self):
    """
    Compute the normal reaction of a contact grain-grain

    Here a pontual spring is considered

        Input :
            itself (a contact grain - grain)
        Output :
            Nothing, but attributes are updated
    """
    #compute angle between grains
    g1_to_g2 = self.g2.center - self.g1.center
    if g1_to_g2[1] >= 0 :
        angle_g1_to_g2 = math.acos(g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 + math.pi
    else :
        angle_g1_to_g2 = math.pi + math.acos(-g1_to_g2[0]/np.linalg.norm(g1_to_g2))
        angle_g2_to_g1 = angle_g1_to_g2 - math.pi

    #extract
    L_i_vertices_1 = extract_vertices(self.g1, angle_g1_to_g2)
    L_i_vertices_2 = extract_vertices(self.g2, angle_g2_to_g1)

    #looking for the nearest nodes
    d_virtual = max(self.g1.r_max,self.g2.r_max)
    ij_min = [0,0]
    d_ij_min = 100*d_virtual #Large
    for i in L_i_vertices_1:
        for j in L_i_vertices_2:
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
        R_eq = 1/(1/self.g1.r_mean+1/self.g2.r_mean)
        k = 4/3*Y_eq*math.sqrt(R_eq)
        F_2_1_n = -k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.Ep_n = 2/5 * k * overlap**(5/2) #-dEp/dx = F_2_1_n
        self.g1.update_f( F_2_1[0],  F_2_1[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1[0], -F_2_1[1], self.g2.l_border[:-1][ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.m*self.g2.m/(self.g1.m+self.g2.m)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][ij_min[1]])

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0
        self.Ep_n = 0

#-------------------------------------------------------------------------------

  def DEM_2grains_Polyhedral_tangential(self,dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact grain - grain)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/self.g1.r_mean+1/self.g2.r_mean)
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

        self.g1.update_f( self.ft*self.pc_tangential[0],  self.ft*self.pc_tangential[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-self.ft*self.pc_tangential[0], -self.ft*self.pc_tangential[1], self.g2.l_border[:-1][self.ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.m*self.g2.m/(self.g1.m+self.g2.m)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.ft_damp = F_2_1_damp_t
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][self.ij_min[1]])

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------

  def DEM_2grains_Polyhedral_normal_surface(self):
    """
    Compute the normal reaction of a contact grain-grain.

    Here a surface spring is considered.

        Input :
            itself (a contact grain - grain)
        Output :
            Nothing, but attributes are updated
    """
    #looking for the nearest nodes
    d_virtual = max(self.g1.r_max,self.g2.r_max) #virtual distance
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
    M1 = self.g1.l_border[:-1][ij_min[0]+1]-self.g1.l_border[:-1][ij_min[0]]
    M1 = M1/np.linalg.norm(M1)
    M3 = self.g1.l_border[:-1][ij_min[0]-1]-self.g1.l_border[:-1][ij_min[0]]
    M3 = M3/np.linalg.norm(M3)
    #reorganize the candidats
    if np.dot(M1,PB) < 0:
        Mtempo = M1.copy()
        M1 = M3.copy()
        M3 = Mtempo.copy()

    #candidats from grain 2
    M2 = self.g2.l_border[:-1][ij_min[1]+1]-self.g2.l_border[:-1][ij_min[1]]
    M2 = M2/np.linalg.norm(M2)
    M4 = self.g2.l_border[:-1][ij_min[1]-1]-self.g2.l_border[:-1][ij_min[1]]
    M4 = M4/np.linalg.norm(M4)
    #reorganize the candidats
    if np.dot(M2,PB) < 0:
      Mtempo = M2.copy()
      M2 = M4.copy()
      M4 = Mtempo.copy()

    #compute the different agnles
    theta_PB = math.pi/2
    theta_M1 =  math.acos(np.dot(M1,N))
    theta_M2 =  math.acos(np.dot(M2,N))
    theta_M3 = -math.acos(np.dot(M3,N))
    theta_M4 = -math.acos(np.dot(M4,N))

    #find the PV
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
    self.pc_normal = PC_normal
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

        #determine the caracteristics of the contact area
        self.S_inter()
        D_inter = self.d_max
        L_inter = sself.l_eq

        #Surface spring term
        Y_eq = (self.g1.y+self.g2.y)/2
        k_surf = Y_eq / L_inter
        F_2_1_n = -k_surf * D_inter * overlap #surface linear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.g1.update_f( F_2_1[0],  F_2_1[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1[0], -F_2_1[1], self.g2.l_border[:-1][ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.m*self.g2.m/(self.g1.m+self.g2.m)
        eta = 2 * gamma * math.sqrt(mass_eq*k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v,PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][ij_min[1]])

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0

#-------------------------------------------------------------------------------

  def DEM_2grains_Polyhedral_tangential_surface(self,dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a surface spring is considered.

        Input :
            itself (a contact grain - grain)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential reaction on the new tangential plane
          self.ft = self.ft*np.dot(self.tangential_old,self.pc_tangential)
        else:
          self.tangential_old_statut = True

        #determine the caracteristics of the contact area
        self.S_inter()
        D_inter = self.d_max
        L_inter = sself.l_eq

        #Surface spring term
        Y_eq = (self.g1.y+self.g2.y)/2
        kn_surf = Y_eq / L_inter #inspired by the bond formulation
        kt_surf = kn_surf

        r1 = np.linalg.norm(self.g1.l_border[:-1][self.ij_min[0]] - self.g1.center) - self.overlap_normal/2
        r2 = np.linalg.norm(self.g2.l_border[:-1][self.ij_min[1]] - self.g2.center) - self.overlap_normal/2
        Delta_Us = (np.dot(self.g1.v-self.g2.v,self.pc_tangential) + r1*self.g1.w + r2*self.g2.w)* dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - kt_surf*D_inter*Delta_Us #inspired by the bond formulation
        self.tangential_old = self.pc_tangential
        if abs(self.ft) > abs(self.mu*self.F_2_1_n) : #Coulomb criteria
            self.ft = self.mu * abs(self.F_2_1_n) * np.sign(self.ft)
        self.g1.update_f( self.ft*self.pc_tangential[0],  self.ft*self.pc_tangential[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-self.ft*self.pc_tangential[0], -self.ft*self.pc_tangential[1], self.g2.l_border[:-1][self.ij_min[1]])

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.m*self.g2.m/(self.g1.m+self.g2.m)
        eta = 2 * gamma * math.sqrt(mass_eq*kt*D_inter)
        F_2_1_damp_t = -Delta_Us/dt_DEM*eta/2
        self.ft_damp = F_2_1_damp_t
        F_2_1_damp = F_2_1_damp_t *self.pc_tangential
        self.g1.update_f( F_2_1_damp[0],  F_2_1_damp[1], self.g1.l_border[:-1][self.ij_min[0]])
        self.g2.update_f(-F_2_1_damp[0], -F_2_1_damp[1], self.g2.l_border[:-1][self.ij_min[1]])

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = 0
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------

  def S_inter(self):
      """
      Compute the intersection surface of the contact grain - grain.

        Input :
            itself (a contact grain - grain)
        Output :
            Nothing, but attributes are updated (three floats)
      """
      L_border = []
      #looking for vertices from grain 2 if they are inside of the grain 1
      for p in self.g2.l_border[:-1] :
          d = np.linalg.norm(np.array(p)-self.g1.center)
          if p[1]>self.g1.center[1]:
              theta = math.acos((p[0]-self.g1.center[0])/np.linalg.norm(self.g1.center-p))
          else :
              theta = 2*math.pi - math.acos((p[0]-self.g1.center[0])/np.linalg.norm(self.g1.center-p))
          #looking for the grain 1 radius in the same direction
          L_theta_R_i = list(abs(np.array(self.g1.l_theta_r)-theta))
          R = self.g1.l_r[L_theta_R_i.index(min(L_theta_R_i))]
          if d < R :
              L_border.append(p)
      #looking for vertices from grain 1 if they are inside of the grain 2
      for p in self.g1.l_border[:-1] :
          d = np.linalg.norm(np.array(p)-self.g2.center)
          if p[1]>self.g2.center[1]:
              theta = math.acos((p[0]-self.g2.center[0])/np.linalg.norm(self.g2.center-p))
          else :
              theta = 2*math.pi - math.acos((p[0]-self.g2.center[0])/np.linalg.norm(self.g2.center-p))
          #looking for the grain 2 radius in the same direction
          L_theta_R_i = list(abs(np.array(self.g2.l_theta_r)-theta))
          R = self.g2.l_r[L_theta_R_i.index(min(L_theta_R_i))]
          if d < R :
              L_border.append(p)

      #adaptation
      if L_border != []:
          L_id_used = [0]
          L_border_adapted = [L_border[0]]
          L_border_x_adapted = [L_border[0][0]]
          L_border_y_adapted = [L_border[0][1]]
          HighValue = 100000000 #Large

          current_node = L_border[0]
          for j in range(1,len(L_border)):
              L_d = list(np.zeros(len(L_border)))
              for i in range(0,len(L_border)):
                  node = L_border[i]
                  if  i not in L_id_used:
                      d = np.linalg.norm(np.array(node) - np.array(current_node))
                      L_d[i] = d
                  else :
                      L_d[i] = HighValue #Value need to be larger than potential distance between node

              index_nearest_node = L_d.index(min(L_d))
              nearest_node = L_border[index_nearest_node]
              current_node = nearest_node
              L_border_adapted.append(nearest_node)
              L_border_x_adapted.append(nearest_node[0])
              L_border_y_adapted.append(nearest_node[1])
              L_id_used.append(index_nearest_node)

      #compute the characteristic values of the contact
      vect_tangential = self.pc_tangential
      d_max = 0
      ij_p = (0,0)
      for i_p1 in range(len(L_border_adapted)-1):
          for i_p2 in range(i_p1+1,len(L_border_adapted)):
              d = abs(np.dot(L_border_adapted[i_p1]-L_border_adapted[i_p2],vect_tangential))
              if d > d_max :
                  d_max = d
                  ij_p = (i_p1, i_p2)
      S = DetermineSurfacePolyhedral(L_border_adapted)
      self.surface = S
      self.d_max = d_max
      self.l_eq = S/d_max

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def extract_vertices(g, angle_g_to_other_g) :
    """
    Extract a list of indices of vertices inside a angular window.

        Input :
            a grain (a grain)
            an angle (a float)
        Output :
            a list of indices (a list)
    """
    dtheta = 3*2*math.pi/len(g.l_border)
    angle_minus = angle_g_to_other_g - dtheta/2
    if angle_minus < 0 :
        angle_minus = angle_minus + 2*math.pi
    angle_plus  = angle_g_to_other_g + dtheta/2
    if 2*math.pi <= angle_plus :
        angle_plus = angle_plus - 2*math.pi
    i_minus = find_value_in_list(g.l_theta_r.copy(), angle_minus)
    i_plus = i_minus + find_value_in_list(g.l_theta_r[i_minus:].copy() + g.l_theta_r[:i_minus].copy(), angle_plus)

    L_i_vertices = []
    for i in range(i_minus, i_plus+1):
        if i < len(g.l_theta_r):
            L_i_vertices.append(i)
        else :
            L_i_vertices.append(i-len(g.l_theta_r))
    return L_i_vertices

#-------------------------------------------------------------------------------

def find_value_in_list(List_to_search, value_to_find) :
    """
    Extract the index of the nearest value from a target in a list.

        Input :
            a list of float (a list)
            a target (a float)
        Output :
            an index (an int)
    """
    L_search = list(abs(np.array(List_to_search)-value_to_find))
    return L_search.index(min(L_search))

#-------------------------------------------------------------------------------

def DetermineSurfacePolyhedral(L_border):
      """
      Use the Monte-Carlo method to compute the surface of a grain.

      Used in S_inter() to compute the surface of the intersection.

        Input :
            a list of coordinates (a list)
        Output :
            a surface (a float)
      """
      #Determine the study box
      min_max_defined = False
      for p in L_border :
          if not min_max_defined:
              box_min_x = p[0]
              box_max_x = p[0]
              box_min_y = p[1]
              box_max_y = p[1]
              min_max_defined = True
          else:
              if p[0] < box_min_x:
                  box_min_x = p[0]
              elif p[0] > box_max_x:
                  box_max_x = p[0]
              if p[1] < box_min_y:
                  box_min_y = p[1]
              elif p[1] > box_max_y:
                  box_max_y = p[1]

      #Monte Carlo method inside the box to  determine the surface
      N_MonteCarlo = 3000
      sigma = 1 #kg/µm2
      M_Mass = 0
      for i in range(N_MonteCarlo):
          P = np.array([random.uniform(box_min_x,box_max_x),random.uniform(box_min_y,box_max_y)])
          if P_is_inside(P,L_border):
              M_Mass = M_Mass + sigma
      Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Mass

      return Mass/sigma

#-------------------------------------------------------------------------------

def P_is_inside(P,L_border):
  """
  Determine if a point P is inside of a grain.

  See Franklin 1994, see Alonso-Marroquin 2009

    Input :
        a point (a 1 x 2 numpy array)
        a list of coordinates (a list)
    Output :
        a Boolean (True if the point is inside the grain)
  """
  #slide on y constant
  counter = 0
  for i_p_border in range(len(L_border)-1):
      #consider only points if the coordinates frame the y-coordinate of the point
      if (L_border[i_p_border][1]-P[1])*(L_border[i_p_border+1][1]-P[1]) < 0 :
        x_border = L_border[i_p_border][0] + (L_border[i_p_border+1][0]-L_border[i_p_border][0])*(P[1]-L_border[i_p_border][1])/(L_border[i_p_border+1][1]-L_border[i_p_border][1])
        if x_border > P[0] :
            counter = counter + 1
  if counter % 2 == 0:
    return False
  else :
    return True

#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_algorithm,dict_sample):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_Polyhedral_contact_Neighborhoods().

    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighbourood of grain_i, with i_grain < j_grain

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but attribut is updated for each grains (a list)
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(dict_sample['L_g'])):
            if np.linalg.norm(dict_sample['L_g'][i_grain].center-dict_sample['L_g'][j_grain].center) < dict_algorithm['factor_neighborhood']*(dict_sample['L_g'][i_grain].r_max+dict_sample['L_g'][j_grain].r_max):
                neighborhood.append(dict_sample['L_g'][j_grain])
        dict_sample['L_g'][i_grain].neighbourood = neighborhood

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods_bool(g1,g2):
  """
  Detect contact grain-grain.

    Input :
        two grains (two grains)
    Output :
        a Boolean (True if there is a contact)
  """
  #compute angle between grains
  g1_to_g2 = g2.center - g1.center
  if g1_to_g2[1] >= 0 :
      angle_g1_to_g2 = math.acos(g1_to_g2[0]/np.linalg.norm(g1_to_g2))
      angle_g2_to_g1 = angle_g1_to_g2 + math.pi
  else :
      angle_g1_to_g2 = math.pi + math.acos(-g1_to_g2[0]/np.linalg.norm(g1_to_g2))
      angle_g2_to_g1 = angle_g1_to_g2 - math.pi

  #extract
  L_i_vertices_1 = extract_vertices(g1, angle_g1_to_g2)
  L_i_vertices_2 = extract_vertices(g2, angle_g2_to_g1)

  #looking for the nearest nodes
  d_virtual = max(g1.r_max,g2.r_max)
  ij_min = [0,0]
  d_ij_min = 100*d_virtual #Large
  for i in L_i_vertices_1:
    for j in L_i_vertices_2:
        d_ij = np.linalg.norm(g2.l_border[:-1][j]-g1.l_border[:-1][i]+d_virtual*(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
        if d_ij < d_ij_min :
            d_ij_min = d_ij
            ij_min = [i,j]

  d_ij_min = np.dot(g2.l_border[:-1][ij_min[1]]-g1.l_border[:-1][ij_min[0]],-(g2.center-g1.center)/np.linalg.norm(g2.center-g1.center))
  return d_ij_min > 0

#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_Neighborhoods(dict_material,dict_sample):
    """
    Detect contact between a grain and grains from its neighborhood.

    The neighbourood is updated with Update_Neighborhoods_f()

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with a list of contacts grain - grain
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        for neighbour in dict_sample['L_g'][i_grain].neighbourood:
            j_grain = neighbour.id
            if Grains_Polyhedral_contact_Neighborhoods_bool(dict_sample['L_g'][i_grain],dict_sample['L_g'][j_grain]):
                if (i_grain,j_grain) not in dict_sample['L_ij_contact']:  #contact not detected previously
                   #creation of contact
                   dict_sample['L_ij_contact'].append((i_grain,j_grain))
                   dict_sample['L_contact'].append(Contact(dict_sample['id_contact'], dict_sample['L_g'][i_grain], dict_sample['L_g'][j_grain], dict_material))
                   dict_sample['id_contact'] = dict_sample['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_sample['L_ij_contact'] : #contact detected previously is not anymore
                       dict_sample['L_contact'].pop(dict_sample['L_ij_contact'].index((i_grain,j_grain)))
                       dict_sample['L_ij_contact'].remove((i_grain,j_grain))
