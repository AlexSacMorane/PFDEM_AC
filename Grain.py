# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the grains
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import math
from PIL import Image as im
import random
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain:

#-------------------------------------------------------------------------------

  def __init__(self, dict_ic_to_real, Id_Eta = None, V = np.array([0,0]), A = np.array([0,0])):
    #defining the grain
    #each grain is described from a tempo grain (see Create_LG_IC)

    #Id of the grain
    self.id = dict_ic_to_real['Id']
    self.id_eta = Id_Eta
    #Material property
    self.dissolved = False
    self.y = dict_ic_to_real['Y']
    self.nu = dict_ic_to_real['Nu']
    self.g = self.y /2/(1+self.nu) #shear modulus
    self.rho_surf = dict_ic_to_real['Rho_surf']
    #kinematic
    self.v = V
    self.a = A
    self.theta = 0
    self.w = 0 #dtheta/dt
    #position
    self.center = dict_ic_to_real['Center']
    self.l_border = dict_ic_to_real['L_border']
    self.l_border_x = dict_ic_to_real['L_border_x']
    self.l_border_y = dict_ic_to_real['L_border_y']
    #characteristic
    self.l_r = dict_ic_to_real['L_r']
    self.l_theta_r = dict_ic_to_real['L_theta_r']
    self.r_min = dict_ic_to_real['R_min']
    self.r_max = dict_ic_to_real['R_max']
    self.r_mean = dict_ic_to_real['R_mean']
    self.surface = dict_ic_to_real['Surface']
    self.m = dict_ic_to_real['Mass']
    self.inertia = dict_ic_to_real['Inertia']

#-------------------------------------------------------------------------------

  def __str__(self):
      pass

#-------------------------------------------------------------------------------

  def update_geometry_kinetic(self, V, A, W, DT):
    #update the acceleration and the velocity of a grain
    #update geometrical parameters as border and center nodes

    #translation
    self.v = V
    self.a = A
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + self.v*DT
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*DT
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*DT
    self.center = self.center + self.v*DT

    #rotation
    self.w = W
    self.theta = self.theta + self.w*DT

    for i_theta_r in range(len(self.l_theta_r)) :
        theta_r = self.l_theta_r[i_theta_r]
        theta_r = theta_r + self.w*DT
        while theta_r >= 2*math.pi:
            theta_r = theta_r - 2*math.pi
        while theta_r < 0 :
            theta_r = theta_r + 2*math.pi
        self.l_theta_r[i_theta_r] = theta_r

    for i in range(len(self.l_border)):
        p = self.l_border[i] - self.center
        Rot_Matrix = np.array([[math.cos(self.w*DT), -math.sin(self.w*DT)],
                               [math.sin(self.w*DT),  math.cos(self.w*DT)]])
        p = np.dot(Rot_Matrix,p)
        self.l_border[i] = p + self.center
        self.l_border_x[i] = p[0] + self.center[0]
        self.l_border_y[i] = p[1] + self.center[1]

#-------------------------------------------------------------------------------

  def init_f_control(self,dict_sollicitations):
      #initialize the force applied to the grain
      #a gravity of g is applied

      #-------------------------------------------------------------------------
      #load data needed
      g = dict_sollicitations['gravity']
      #-------------------------------------------------------------------------

      self.fx = 0
      self.fy = -g*self.m
      self.f = np.array([self.fx,self.fy])
      self.mz = 0

#-------------------------------------------------------------------------------

  def update_f(self, Fx, Fy, p_application):
    #add a force (an array [Fx,Fy]) to the grain

    self.fx = self.fx + Fx
    self.fy = self.fy + Fy
    self.f = np.array([self.fx,self.fy])

    v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
    v2 = np.array([Fx, Fy, 0])
    self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------

  def Geometricstudy(self,dict_geometry,dict_sample,simulation_report):
      #Searching limits
      #Not best method but first approach
      #We iterate on y constant, we look for a value under and over 0.5
      #If both conditions are verified, there is a limit at this y
      #Same with iteration on x constant

      #-------------------------------------------------------------------------
      #load data needed
      n = dict_geometry['grain_discretisation']
      x_L = dict_sample['x_L']
      y_L = dict_sample['y_L']
      #-------------------------------------------------------------------------

      L_border_old = []
      for y_i in range(len(y_L)):
          L_extract_x = self.etai_M[y_i][:]
          if id == 1:
              L_extract_x = list(L_extract_x)
              L_extract_x.reverse()
          if max(L_extract_x)>0.5 and min(L_extract_x)<0.5:
              y_intersect = y_L[len(y_L)-1-y_i]
              for x_i in range(len(x_L)-1):
                  if (L_extract_x[x_i]-0.5)*(L_extract_x[x_i+1]-0.5)<0:
                      x_intersect = (0.5-L_extract_x[x_i])/(L_extract_x[x_i+1]-L_extract_x[x_i])*\
                                  (x_L[x_i+1]-x_L[x_i]) + x_L[x_i]
                      L_border_old.append(np.array([x_intersect,y_intersect]))

      for x_i in range(len(x_L)):
          L_extract_y = []
          for y_i in range(len(y_L)):
              L_extract_y.append(self.etai_M[y_i][x_i])
          if max(L_extract_y)>0.5 and min(L_extract_y)<0.5:
              x_intersect = x_L[x_i]
              for y_i in range(len(y_L)-1):
                  if (L_extract_y[y_i]-0.5)*(L_extract_y[y_i+1]-0.5)<0:
                      y_intersect = (0.5-L_extract_y[y_i])/(L_extract_y[y_i+1]-L_extract_y[y_i])*\
                                  (y_L[len(y_L)-1-y_i-1]-y_L[len(y_L)-1-y_i]) + y_L[len(y_L)-1-y_i]
                      L_border_old.append(np.array([x_intersect,y_intersect]))

      #Adaptating
      L_id_used = [0]
      L_border = [L_border_old[0]]
      HighValue = 100000000 #Large

      current_node = L_border_old[0]
      for j in range(1,len(L_border_old)):
          L_d = list(np.zeros(len(L_border_old)))
          for i in range(0,len(L_border_old)):
              node = L_border_old[i]
              if  i not in L_id_used:
                  d = np.linalg.norm(node - current_node)
                  L_d[i] = d
              else :
                  L_d[i] = HighValue #Value need to be larger than potential distance between node

          index_nearest_node = L_d.index(min(L_d))
          nearest_node = L_border_old[index_nearest_node]
          current_node = nearest_node
          L_border.append(nearest_node)
          L_id_used.append(index_nearest_node)

      #Correcting
      L_d_final = []
      for i in range(len(L_border)-1):
          L_d_final.append(np.linalg.norm(L_border[i+1] - L_border[i]))

      #look for really far points, we assume the first point is accurate
      d_final_mean = np.mean(L_d_final)
      while np.max(L_d_final) > 5 * d_final_mean : #5 here is an user choixe value
          i_error = L_d_final.index(np.max(L_d_final))+1
          simulation_report.write('Point '+str(L_border[i_error])+' is deleted because it is detected as an error\n')
          L_border.pop(i_error)
          L_id_used.pop(i_error)
          L_d_final = []
          for i in range(len(L_border)-1):
              L_d_final.append(np.linalg.norm(L_border[i+1] - L_border[i]))

      #-------------------------------------------------------------------------------
      #Reduce the number of nodes for a grain
      #-------------------------------------------------------------------------------

      Perimeter = 0
      for i_p in range(len(L_border)-1):
          Perimeter = Perimeter + np.linalg.norm(L_border[i_p+1]-L_border[i_p])
      Perimeter = Perimeter + np.linalg.norm(L_border[-1]-L_border[0])
      distance_min = Perimeter/n
      L_border_adapted = [L_border[0]]
      for p in L_border[1:]:
          distance = np.linalg.norm(p-L_border_adapted[-1])
          if distance >= distance_min:
              L_border_adapted.append(p)
      L_border = L_border_adapted
      L_border.append(L_border[0])
      self.l_border = L_border

      #-------------------------------------------------------------------------------
      #Searching Surface, Center of mass and Inertia.
      #Monte Carlo Method
      #A box is defined, we take a random point and we look if it is inside or outside the grain
      #Properties are the statistic times the box properties
      #-------------------------------------------------------------------------------

      min_max_defined = False
      for p in L_border[:-1] :
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

      N_MonteCarlo = 3000 #The larger it is, the more accurate it is
      sigma = self.rho_surf #kg/µm2
      M_Mass = 0
      M_Center_Mass = np.array([0,0])
      M_Inertia = 0

      for i in range(N_MonteCarlo):
          P = np.array([random.uniform(box_min_x,box_max_x),random.uniform(box_min_y,box_max_y)])
          if self.P_is_inside(P):
              M_Mass = M_Mass + sigma
              M_Center_Mass = M_Center_Mass + sigma*P
              M_Inertia = M_Inertia + sigma*np.dot(P,P)

      Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Mass
      Center_Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Center_Mass/Mass
      Inertia = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Inertia-Mass*np.dot(Center_Mass,Center_Mass)

      #-------------------------------------------------------------------------------
      #Updating the grain geometry and properties
      #-------------------------------------------------------------------------------

      L_R = []
      L_theta_R = []
      L_border_x = []
      L_border_y = []
      for p in L_border[:-1]:
          L_R.append(np.linalg.norm(p-Center_Mass))
          L_border_x.append(p[0])
          L_border_y.append(p[1])
          if (p-Center_Mass)[1] > 0:
              theta = math.acos((p-Center_Mass)[0]/np.linalg.norm(p-Center_Mass))
          else :
              theta = 2*math.pi - math.acos((p-Center_Mass)[0]/np.linalg.norm(p-Center_Mass))
          L_theta_R.append(theta)
      L_border_x.append(L_border_x[0])
      L_border_y.append(L_border_y[0])
      #reorganize lists
      L_R.reverse()
      L_theta_R.reverse()
      i_min_theta = L_theta_R.index(min(L_theta_R))
      L_R = L_R[i_min_theta:]+L_R[:i_min_theta]
      L_theta_R = L_theta_R[i_min_theta:]+L_theta_R[:i_min_theta]

      self.r_min = np.min(L_R)
      self.r_max = np.max(L_R)
      self.r_mean = np.mean(L_R)
      self.l_r = L_R
      self.l_theta_r = L_theta_R
      self.surface = Mass/self.rho_surf
      self.m = Mass
      self.center = Center_Mass
      self.l_border_x = L_border_x
      self.l_border_y = L_border_y
      self.inertia = Inertia

#-------------------------------------------------------------------------------

  def P_is_inside(self,P):
      #Franklin 1994, see Alonso-Marroquin 2009
      #determine if a point P is inside of a grain
      #slide on y constant

      counter = 0
      for i_p_border in range(len(self.l_border)-1):
          #consider only points if the coordinates frame the y-coordinate of the point
          if (self.l_border[i_p_border][1]-P[1])*(self.l_border[i_p_border+1][1]-P[1]) < 0 :
            x_border = self.l_border[i_p_border][0] + (self.l_border[i_p_border+1][0]-self.l_border[i_p_border][0])*(P[1]-self.l_border[i_p_border][1])/(self.l_border[i_p_border+1][1]-self.l_border[i_p_border][1])
            if x_border > P[0] :
                counter = counter + 1
      if counter % 2 == 0:
        return False
      else :
        return True

#-------------------------------------------------------------------------------

  def DEMtoPF_Decons_rebuild(self,dict_material,dict_sample):
        #from the grain geometry the phase variable is rebuilt
        #the distance between the point of the mesh and the particle center determines the value of the variable
        #a cosine profile is applied inside the interface

        etai_M_new = self.etai_M.copy()
        L_R = self.l_r
        L_theta_R = self.l_theta_r

        #extract a part focused on the grain
        x_extract_min = self.center[0] - self.r_max - dict_material['w']
        x_extract_max = self.center[0] + self.r_max + dict_material['w']
        y_extract_min = self.center[1] - self.r_max - dict_material['w']
        y_extract_max = self.center[1] + self.r_max + dict_material['w']

        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_extract_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_extract_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_extract_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_extract_max)

        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))

        #---------------------------------------------------------------------------
        # Move
        #---------------------------------------------------------------------------

        for i_x in range(i_x_min,i_x_max+1):
            for i_y in range(i_y_min,i_y_max+1):
                p = np.array([dict_sample['x_L'][i_x], dict_sample['y_L'][len(dict_sample['y_L'])-1-i_y]])
                r = np.linalg.norm(self.center - p)
                if p[1]>self.center[1]:
                    theta = math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
                else :
                    theta= 2*math.pi - math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))

                L_theta_R_i = list(abs(np.array(L_theta_R)-theta))
                R = L_R[L_theta_R_i.index(min(L_theta_R_i))]
                #Cosine_Profile
                if r<R-dict_material['w']/2:
                    etai_M_new[i_y][i_x] = 1
                elif r>R+dict_material['w']/2:
                    etai_M_new[i_y][i_x] = 0
                else :
                    etai_M_new[i_y][i_x] = 0.5*(1 + np.cos(math.pi*(r-R+dict_material['w']/2)/dict_material['w']))

        self.etai_M = etai_M_new.copy()

#-------------------------------------------------------------------------------

  def Write_e_dissolution_local_txt(self,dict_algorithm,dict_sollicitations):
      #write an .txt file for MOOSE
      #this file described an homogenous dissolution field

      file_to_write = open(f"Data/e_diss_g{self.id}_ite{dict_algorithm['i_PF']}.txt",'w')
      file_to_write.write('AXIS X\n')
      line = ''
      for x in self.x_L_local:
          line = line + str(x)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('AXIS Y\n')
      line = ''
      for y in self.y_L:
        line = line + str(y)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('DATA\n')
      for l in range(len(y_L)):
          for c in range(len(y_L_local)):
              file_to_write.write(str(dict_sollicitations['Dissolution_Energy'])+'\n')

      file_to_write.close()

#-------------------------------------------------------------------------------

  def Compute_etaiM_local(self,dict_algorithm,dict_material):

      x_min_local = min(self.l_border_x)-dict_material['w']
      x_max_local = max(self.l_border_x)+dict_material['w']
      y_min_local = min(self.l_border_y)-dict_material['w']
      y_max_local = max(self.l_border_y)+dict_material['w']
      x_L_local = np.linspace(x_min_local,x_max_local,dict_algorithm['n_local'])
      y_L_local = np.linspace(x_min_local,x_max_local,dict_algorithm['n_local'])
      grain.x_L_local = x_L_local
      grain.y_L_local = y_L_local

      # compute phase field
      etai_M = np.array(np.zeros((dict_algorithm['n_local'],dict_algorithm['n_local'])))
      for i_x in range(dict_algorithm['n_local']):
          for i_y in range(dict_algorithm['n_local']):
              p = np.array([x_L_local[i_x],y_L_local[len(y_L_local)-1-i_y]])
              r = np.linalg.norm(self.center - p)
              if p[1]>self.center[1]:
                  theta = math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
              else :
                  theta= 2*math.pi - math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))

              L_theta_R_i = list(abs(np.array(self.l_theta_r)-theta))
              R = self.l_r[L_theta_R_i.index(min(L_theta_R_i))]
              #Cosine_Profile
              if r<R-dict_material['w']/2:
                  etai_M_new[i_y][i_x] = 1
              elif r>R+dict_material['w']/2:
                  etai_M_new[i_y][i_x] = 0
              else :
                  etai_M_new[i_y][i_x] = 0.5*(1 + np.cos(math.pi*(r-R+dict_material['w']/2)/dict_material['w']))
      self.etai_M = etai_M.copy()

#-------------------------------------------------------------------------------

  def Write_txt_Decons_rebuild_local(self,dict_algorithm):
      #write a .txt file
      #this file is used to define initial condition of MOOSE simulation

      file_to_write = open('Data/g'+str(self.id)+'_'+str(dict_algorithm['i_PF'])+'.txt','w')
      file_to_write.write('AXIS X\n')
      line = ''
      for x in self.x_L_local:
          line = line + str(x)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('AXIS Y\n')
      line = ''
      for y in self.y_L_local:
        line = line + str(y)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('DATA\n')
      for l in range(len(self.y_L_local)):
          for c in range(len(self.x_L_local)):
              file_to_write.write(str(self.etai_M[-l-1][c])+'\n')

      file_to_write.close()

#-------------------------------------------------------------------------------

  def DEMtoPF_Interpolation(self):

      pass

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------
