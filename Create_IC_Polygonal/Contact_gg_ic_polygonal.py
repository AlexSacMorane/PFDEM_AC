# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between two grains during initial condition generation with polygonal particles.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math

#Own
import Create_IC_Polygonal.Grain_ic_polygonal

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_Tempo_Polygonal:
  """
  A temporary contact grain - grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact_tempo)
            an id (a int)
            two grains (two grain_tempos)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated (a contact_tempo)
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    self.ft = 0
    self.mu = 0
    self.coeff_restitution = dict_material['coeff_restitution']
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def normal(self):
    """
    Compute the normal reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact_tempo)
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
        R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
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
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered

        Input :
            itself (a contact_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """

    #change here!!!


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
#Function
#-------------------------------------------------------------------------------

def Grains_Polyhedral_contact_f(g1,g2):
  """
  Detect the contact grain-grain.

    Input :
        two temporary grains (two grain_tempos)
    Output :
        a Boolean, True if there is contaxct between the twwo grains (a Boolean)
  """
  if np.linalg.norm(g1.center-g2.center) < 1.5*(g1.r_max+g2.r_max):
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

  else:
    return False

#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_ic):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_contact_Neighborhoods().
    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighbourood of grain_i. With i_grain < j_grain.

        Input :
            an initial condition dictionnary (a dict)
        Output :
            Nothing, but the neighborhood of the temporary grains is updated
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        neighbourood = []
        for j_grain in range(i_grain+1,len(dict_ic['L_g_tempo'])):
            if np.linalg.norm(dict_ic['L_g_tempo'][i_grain].center-dict_ic['L_g_tempo'][j_grain].center) < dict_ic['factor_neighborhood_IC']*(dict_ic['L_g_tempo'][i_grain].r_max+dict_ic['L_g_tempo'][j_grain].r_max):
                neighbourood.append(dict_ic['L_g_tempo'][j_grain])
        dict_ic['L_g_tempo'][i_grain].neighbourood = neighbourood

#-------------------------------------------------------------------------------

def Grains_contact_Neighborhoods(dict_ic,dict_material):
    """
    Detect contact between a grain and grains from its neighbourood.

    The neighbourood is updated with Update_Neighbouroods().

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with grain - grain contacts
    """
    for i_grain in range(len(dict_ic['L_g_tempo'])-1) :
        grain_i = dict_ic['L_g_tempo'][i_grain]
        for neighbour in dict_ic['L_g_tempo'][i_grain].neighbourood:
            j_grain = neighbour.id
            grain_j = neighbour
            if Grains_Polyhedral_contact_f(grain_i,grain_j):
                if (i_grain,j_grain) not in dict_ic['L_contact_ij']:  #contact not detected previously
                   #creation of contact
                   dict_ic['L_contact_ij'].append((i_grain,j_grain))
                   dict_ic['L_contact'].append(Contact_Tempo_Polygonal(dict_ic['id_contact'], grain_i, grain_j, dict_material))
                   dict_ic['id_contact'] = dict_ic['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_ic['L_contact_ij'] : #contact detected previously is not anymore
                       dict_ic['L_contact'].pop(dict_ic['L_contact_ij'].index((i_grain,j_grain)))
                       dict_ic['L_contact_ij'].remove((i_grain,j_grain))

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
