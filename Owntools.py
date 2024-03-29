# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
This tools can be : - a little function
                    - a plot to debug or print something
                    - a postprocess
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import math
import Report
import os
import shutil
import pickle
import random
import imageio
from pathlib import Path
from datetime import datetime
import Contact
import Contact_gw
import Grain

#-------------------------------------------------------------------------------

class Grain_cf:

    def __init__(self, Id, Dissolved, Center, Coordinate_x, Coordinate_y):
        """
        Defining a grain for the chain force plot.

            Input :
                itself (a Grain_cf)
                an id (a int)
                a Boolean to know if the grain is dissolvable (a Boolean)
                a center (a 1 x 2 numpy array)
                two lists of vertices coordinates x and y (two lists)
            Output :
                Nothing, but a post process grain is generated
        """
        self.id = Id
        self.dissolved = Dissolved
        self.center = Center
        self.coordinate_x = Coordinate_x
        self.coordinate_y = Coordinate_y

#-------------------------------------------------------------------------------

class Contact_cf:

    def __init__(self, Id_g1, Id_g2, L_g, Normal):
        """
        Defining a contact grain - grain for the chain force plot.

            Input :
                itself (a Contact_cf)
                the ids of the grains (two int)
                a list of post process grains (a list)
                the value of normal reaction of the contact (a float)
            Output :
                Nothing, but a post process contact grain - grain is generated
        """
        for g in L_g:
            if g.id == Id_g1:
                self.g1 = g
            elif g.id == Id_g2:
                self.g2 = g
        self.normal = Normal

    def plot(self, normal_ref):
        """
        Prepare the chain force plot.

            Input :
                itself (a Contact_cf)
                a reference value (a float)
            Output :
                Nothing, but the post process contact gets the ratio of the normal force with the reference value as a new attribut (a float)
        """
        L_x = [self.g1.center[0], self.g2.center[0]]
        L_y = [self.g1.center[1], self.g2.center[1]]
        ratio_normal = self.normal/normal_ref
        return L_x, L_y, ratio_normal

#-------------------------------------------------------------------------------

class Contact_gw_cf:

    def __init__(self, Id_g, L_g, Nature, Limit, Normal):
        """
        Defining a contact grain-wall for the chain force plot.

            Input :
                itself (a Contact_gw_cf)
                a id of the grain (a float)
                the list of the post process grain (a list)
                the nature of the wall (a string)
                the value of normal reaction of the contact (a float)
        """
        for g in L_g:
            if g.id == Id_g:
                self.g = g
        self.nature = Nature
        self.limit = Limit
        self.normal = Normal

    def plot(self, normal_ref):
        """
        Prepare the chain force plot.

            Input :
                itself (a Contact_gw_cf)
                a reference value (a float)
            Output :
                Nothing, but the post process contact gets the ratio of the normal force with the reference value as a new attribut (a float)
        """
        if self.nature == 'gwx_min' or self.nature == 'gwx_max':
            virtual_center = [self.limit, self.g.center[1]]
        elif self.nature == 'gwy_min' or self.nature == 'gwy_max':
            virtual_center = [ self.g.center[0], self.limit]
        L_x = [self.g.center[0], virtual_center[0]]
        L_y = [self.g.center[1], virtual_center[1]]
        ratio_normal = self.normal/normal_ref
        return L_x, L_y, ratio_normal

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def index_to_str(j):
      """
      Convert an integer to a float with 3 components

        Input :
            an integer (a int)
        Output :
            a string (a string)
      """
      if j < 10:
          j_str = '00'+str(j)
      elif 10 <= j and j < 100:
          j_str = '0'+str(j)
      else :
          j_str = str(j)
      return j_str

#-------------------------------------------------------------------------------

def Stop_Debug(simulation_report):
      """
      Stop simulation for debugging.

        Input :
            the simulation report (a report)
        Output :
            Nothing, but simulations stops
      """
      simulation_report.write('Stop because after it is not corrected !\n')
      simulation_report.end(datetime.now())
      raise ValueError('Stop because after it is not corrected !')

#-------------------------------------------------------------------------------

def Debug_DEM_f(dict_algorithm, dict_sample):
    """
    Plot the configuration of the grains during a DEM step.

        Input:
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    #load data needed
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_min = dict_sample['y_box_min']
    y_max = dict_sample['y_box_max']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

    fig = plt.figure(1,figsize=(16,9.12))
    for grain in dict_sample['L_g']:
        plt.plot(grain.l_border_x,grain.l_border_y,'k')
    for contact in dict_sample['L_contact'] :
        alpha = (contact.g1.r_mean*0.1 + contact.g2.r_mean*0.1)/2
        M = (contact.g1.center + contact.g2.center)/2
        plt.plot([M[0], M[0]+contact.pc_normal[0]*alpha], [M[1], M[1]+contact.pc_normal[1]*alpha],'k')
        plt.plot([M[0], M[0]+contact.pc_tangential[0]*alpha], [M[1], M[1]+contact.pc_tangential[1]*alpha],'k')
    for contact in dict_sample['L_contact_gw'] :
        alpha = contact.g.r_mean*0.1
        if contact.nature == 'gwy_min' or contact.nature == 'gwy_max':
            M = np.array([contact.g.center[0],contact.limit])
        elif contact.nature == 'gwx_min' or contact.nature == 'gwx_max':
            M = np.array([contact.limit,contact.g.center[1]])
        plt.plot([M[0], M[0]+contact.nwg[0]*alpha], [M[1], M[1]+contact.nwg[1]*alpha],'k')
        plt.plot([M[0], M[0]+contact.twg[0]*alpha], [M[1], M[1]+contact.twg[1]*alpha],'k')
    plt.plot([x_min,x_max,x_max,x_min,x_min],[y_min,y_min,y_max,y_max,y_min],'k')
    plt.axis("equal")
    fig.savefig('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/png/Config_'+str(dict_algorithm['i_DEM'])+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Debug_configuration(dict_algorithm,dict_sample):
  """
  Plot the configuration of the grains before / after the DEM step.

        Input:
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
  """
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  #load data needed
  x_min = dict_sample['x_box_min']
  x_max = dict_sample['x_box_max']
  y_min = dict_sample['y_box_min']
  y_max = dict_sample['y_box_max']
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

  #look for the name of the new plot
  template_name = 'Debug/DEM_ite/PF_ite_'
  j = 0
  plotpath = Path(template_name+str(j)+'.png')
  while plotpath.exists():
      j = j + 1
      plotpath = Path(template_name+str(j)+'.png')
  name = template_name+str(j)+'.png'

  #create the plot
  fig = plt.figure(1,figsize=(16,9.12))
  for grain in dict_sample['L_g']:
      if grain.dissolved:
          plt.plot(grain.l_border_x,grain.l_border_y,'k-.')
      else:
          plt.plot(grain.l_border_x,grain.l_border_y,'k')
  plt.plot([x_min,x_max,x_max,x_min,x_min],[y_min,y_min,y_max,y_max,y_min],'k')
  plt.axis("equal")
  fig.savefig(name)
  plt.close(1)

#-------------------------------------------------------------------------------

def Debug_Trackers(dict_tracker):
    """
    Plot the trakers used during PFDEM simulation

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but severals .png files are generated (files)
    """
    #Trackers
    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['t_L'],dict_tracker['S_grains_dissolvable_L'])
    plt.title('Evolution of the dissolvable grains surface')
    fig.savefig('Debug/Evolution_Dissolvable_Surface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'],dict_tracker['k0_xmin_L'],label='k0 with xmin')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'],dict_tracker['k0_xmax_L'],label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Percentage of dissolvable grains surface dissolved (%)')
    fig.savefig('Debug/Evolution_k0_with_percentage_dissolved.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['porosity_L'][1:])
    plt.title('Evolution of the porosity')
    plt.xlabel('Percentage of dissolvable grains surface dissolved (%)')
    fig.savefig('Debug/Evolution_porosity.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Debug_Trackers_DEM(dict_algorithm,dict_sollicitations,dict_tracker):
    """
    Plot the trakers used during DEM simulation.

        Input :
            an algorithm dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but .png file is generated (a file)
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(16,9),num=1)

    ax1.set_title('Mean kinetic energy (e-12 J)')
    ax1.plot(dict_tracker['Ecin'])
    ax1.plot([0, len(dict_tracker['Ecin'])-1],[dict_algorithm['Ecin_stop'], dict_algorithm['Ecin_stop']],'r')

    ax2.set_title('Mean force applied (µN)')
    ax2.plot(dict_tracker['Force_applied'])

    ax3.set_title('k0s (-)')
    ax3.plot(dict_tracker['k0_xmin'],label='xmin')
    ax3.plot(dict_tracker['k0_xmax'],label='xmax')
    ax3.legend()

    ax4.set_title('About the upper plate')
    ax4.plot(dict_tracker['y_box_max'], color = 'blue')
    ax4.set_ylabel('Coordinate y (µm)', color = 'blue')
    ax4.tick_params(axis ='y', labelcolor = 'blue')
    ax4a = ax4.twinx()
    ax4a.plot(range(50,len(dict_tracker['Force_on_upper_wall'])),dict_tracker['Force_on_upper_wall'][50:], color = 'orange')
    ax4a.plot([50, len(dict_tracker['Force_on_upper_wall'])-1],[dict_sollicitations['Vertical_Confinement_Force'], dict_sollicitations['Vertical_Confinement_Force']], color = 'red')
    ax4a.set_ylabel('Force applied (µN)', color = 'orange')
    ax4a.tick_params(axis ='y', labelcolor = 'orange')

    fig.savefig('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/Trackers.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Sort_Files(name_template,dict_algorithm):
     """
     Sort files generated by MOOSE to different directories.

        Input :
            a template of the simulation files (a string)
            an algorithm dictionnary (a dict)
        Output :
            Nothing, but files are sorted
     """
     #master simulation
     os.rename(name_template+'_out.e','Output/'+name_template+'_out.e')
     os.rename(name_template+'.i','Input/'+name_template+'.i')
     if Path('Output/'+name_template).exists():
         shutil.rmtree('Output/'+name_template)
     os.mkdir('Output/'+name_template)
     j = 0
     j_str = index_to_str(j)
     folderpath = Path(name_template+'_other_'+j_str+'.pvtu')
     while folderpath.exists():
         for i_proc in range(dict_algorithm['np_proc']):
            os.rename(name_template+'_other_'+j_str+'_'+str(i_proc)+'.vtu','Output/'+name_template+'/'+name_template+'_other_'+j_str+'_'+str(i_proc)+'.vtu')
         os.rename(name_template+'_other_'+j_str+'.pvtu','Output/'+name_template+'/'+name_template+'_other_'+j_str+'.pvtu')
         j = j + 1
         j_str = index_to_str(j)
         folderpath = Path(name_template+'_other_'+j_str+'.pvtu')

     return index_to_str(j-1)

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
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
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Control_y_max_NR(dict_sample,dict_sollicitations):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated, concerning the upper wall position and force applied (two floats)
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    Force_target = dict_sollicitations['Vertical_Confinement_Force']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    F = 0
    overlap_L = []
    k_L = []
    for contact in dict_sample['L_contact_gw']:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dy = 0
        ite_criteria = True
        if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy = dy - error_on_ymax_f(dy,overlap_L,k_L,Force_target)/error_on_ymax_df(dy,overlap_L,k_L)
            if i_NR > 100:
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy,overlap_L,k_L,Force_target) and error_on_ymax_f(dy,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy

    else :
        #if there is no contact with the upper wall, the wall is reset
        dict_sample['y_box_max'] = Reset_y_max(dict_sample['L_g'],Force_target)

    for contact in dict_sample['L_contact_gw']:
        if contact.nature == 'gwy_max':
            #reactualisation
            contact.limit = dict_sample['y_box_max']

    #Update dict
    dict_sollicitations['Force_on_upper_wall'] = F

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
    print('Reset of the y_max on the upper grain')
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

    k = 5*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].r_mean)
    y_max = y_max - (Force/k)**(2/3)

    return y_max

#-------------------------------------------------------------------------------

def Compute_porosity(dict_sample):
    """
    Compute the porosity = grains surface / box surface.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary gets updated value concerning the porosity
    """
    Sg = 0
    for grain in dict_sample['L_g']:
        Sg = Sg + grain.surface
    Sb = (dict_sample['x_box_max']-dict_sample['x_box_min'])*(dict_sample['y_box_max']-dict_sample['y_box_min'])
    dict_sample['porosity'] = Sg/Sb

#-------------------------------------------------------------------------------

def Compute_k0(dict_sample,dict_sollicitations):
    """
    Compute the k0 = sigma_2 / sigma_1 ratio.

        Input :
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary gets updated value concerning the k0
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #Load data needed
    L_contact_gw = dict_sample['L_contact_gw']
    x_min = dict_sample['x_box_min']
    x_max = dict_sample['x_box_max']
    y_min = dict_sample['y_box_min']
    y_max = dict_sample['y_box_max']
    F_y_max = dict_sollicitations['Force_on_upper_wall']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #compute the forces applied on left and right walls (to compute k0)
    F_x_min = 0
    F_x_max = 0
    for contact in L_contact_gw:
        if contact.nature == 'gwx_min':
            F_x_min = F_x_min + contact.Fwg_n
        elif contact.nature == 'gwx_max':
            F_x_max = F_x_max + contact.Fwg_n

    #compute k0 = (sigma_2/sigma_1)
    sigma_x_min = F_x_min / (y_max-y_min)
    sigma_x_max = F_x_max / (y_max-y_min)
    sigma_y_max = F_y_max / (x_max-x_min)

    if sigma_y_max!= 0:
        k0_xmin = abs(sigma_x_min/sigma_y_max)
        k0_xmax = abs(sigma_x_min/sigma_y_max)
    else :
        k0_xmin = 0
        k0_xmax = 0

    #update element in dicts
    dict_sample['k0_xmin'] = k0_xmin
    dict_sample['k0_xmax'] = k0_xmax

#-------------------------------------------------------------------------------

def Write_e_dissolution_txt(dict_sample,dict_sollicitations):
      """
      Write an .txt file for MOOSE. This file described an homogenous dissolution field.

        Input :
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            Nothing, but a .txt file is generated (a file)
      """
      #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
      #Load data needed
      x_L = dict_sample['x_L']
      y_L = dict_sample['y_L']
      e_dissolution = dict_sollicitations['Dissolution_Energy']
      #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

      file_to_write = open('Data/e_dissolution.txt','w')
      file_to_write.write('AXIS X\n')
      line = ''
      for x in x_L:
          line = line + str(x)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('AXIS Y\n')
      line = ''
      for y in y_L:
        line = line + str(y)+ ' '
      line = line + '\n'
      file_to_write.write(line)

      file_to_write.write('DATA\n')
      for l in range(len(y_L)):
          for c in range(len(x_L)):
              file_to_write.write(str(e_dissolution)+'\n')

      file_to_write.close()

#-------------------------------------------------------------------------------

def Plot_chain_force(i_PF,i_DEM):
    """
    Plot the chain force.

        Input :
            the iteration PFDEM (a int)
            the iteration DEM (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    normal_ref = 5*10**6 #can be changed

    file_name = 'Debug/DEM_ite/PF_'+str(i_PF)+'/txt/ite_DEM_'+str(i_DEM)+'.txt'
    file = open(file_name,'r')
    lines = file.readlines()
    file.close()

    L_g = []
    L_contact = []
    L_contact_gw = []
    Read_one_grain = False
    Read_one_contact = False
    Read_one_contact_wall = False
    for line in lines :

        if line == '<grain_c>\n':
            Read_one_grain = False
            L_g.append(Grain_cf(id,dissolved,center,coordinate_x,coordinate_y))
        elif line == '<contact_c>\n':
            Read_one_contact = False
            L_contact.append(Contact_cf(L_id_g[0], L_id_g[1], L_g, normal_reaction+normal_damping))
        elif line == '<contact_w_c>\n':
            Read_one_contact_wall = False
            L_contact_gw.append(Contact_gw_cf(id_g, L_g, type, limit, normal_reaction+normal_damping))

        if Read_one_grain:
            if line[:len('\tid : ')] == '\tid : ':
                Read_data = False
                for c in range(len('\tid : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        id = float(line[c_start:c])
                        Read_data = False
            elif line[:len('\tDissolved : ')] == '\tDissolved : ':
                Read_data = False
                for c in range(len('\tDissolved : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        if line[c_start:c] == 'True':
                            dissolved = True
                        else :
                            dissolved = False
                        Read_data = False
            elif line[:len('\tCenter : ')] == '\tCenter : ':
                Read_data = False
                center = []
                for c in range(len('\tCenter : ')-1,len(line)):
                    if (line[c]!=' ' and line[c]!='[' and line[c]!=']' and line[c]!=',') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or line[c]==',' or line[c]==']') and Read_data:
                        data = float(line[c_start:c])
                        center.append(data)
                        Read_data = False
            elif line[:len('\tCoordinate X of the border : ')] == '\tCoordinate X of the border : ':
                Read_data = False
                coordinate_x = []
                for c in range(len('\tCoordinate X of the border : ')-1,len(line)):
                    if (line[c]!=' ' and line[c]!='[' and line[c]!=']' and line[c]!=',') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or line[c]==',' or line[c]==']') and Read_data:
                        data = float(line[c_start:c])
                        coordinate_x.append(data)
                        Read_data = False
            elif line[:len('\tCoordinate Y of the border : ')] == '\tCoordinate Y of the border : ':
                Read_data = False
                coordinate_y = []
                for c in range(len('\tCoordinate Y of the border : ')-1,len(line)):
                    if (line[c]!=' ' and line[c]!='[' and line[c]!=']' and line[c]!=',') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or line[c]==',' or line[c]==']') and Read_data:
                        data = float(line[c_start:c])
                        coordinate_y.append(data)
                        Read_data = False

        if Read_one_contact:
            if line[:len('\tGrains : ')] == '\tGrains : ':
                Read_data = False
                L_id_g = []
                for c in range(len('\tGrains : ')-1,len(line)):
                    if (line[c]!=' ' and line[c]!='-') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or line[c]=='-' or c==len(line)-1) and Read_data:
                        data = float(line[c_start:c])
                        L_id_g.append(data)
                        Read_data = False
            elif line[:len('\tNormal reaction : ')] == '\tNormal reaction : ':
                Read_data = False
                for c in range(len('\tNormal reaction : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        normal_reaction = float(line[c_start:c])
                        Read_data = False
            elif line[:len('\tNormal damping : ')] == '\tNormal damping : ':
                Read_data = False
                for c in range(len('\tNormal damping : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        normal_damping = float(line[c_start:c])
                        Read_data = False

        if Read_one_contact_wall:
            if line[:len('\tType : ')] == '\tType : ':
                Read_data = False
                for c in range(len('\tType : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        type = line[c_start:c]
                        Read_data = False
            elif line[:len('\tGrain : ')] == '\tGrain : ':
                Read_data = False
                for c in range(len('\tGrain : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        id_g = float(line[c_start:c])
                        Read_data = False
            elif line[:len('\tLimit : ')] == '\tLimit : ':
                Read_data = False
                for c in range(len('\tLimit : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        limit = float(line[c_start:c])
                        Read_data = False
                if type == 'gwx_min':
                    x_min = limit
                elif type == 'gwx_max':
                    x_max = limit
                elif type == 'gwy_min':
                    y_min = limit
                elif type == 'gwy_max':
                    y_max = limit
            elif line[:len('\tNormal reaction : ')] == '\tNormal reaction : ':
                Read_data = False
                for c in range(len('\tNormal reaction : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        normal_reaction = float(line[c_start:c])
                        Read_data = False
            elif line[:len('\tNormal damping : ')] == '\tNormal damping : ':
                Read_data = False
                for c in range(len('\tNormal damping : ')-1,len(line)):
                    if (line[c]!=' ') and not Read_data:
                        c_start = c
                        Read_data = True
                    elif (line[c]==' ' or c==len(line)-1) and Read_data:
                        normal_damping = float(line[c_start:c])
                        Read_data = False

        if line == '<grain_o>\n':
            Read_one_grain = True
        elif line == '<contact_o>\n':
            Read_one_contact = True
        elif line == '<contact_w_o>\n':
            Read_one_contact_wall = True

    plt.figure(1,figsize=(16,9))
    plt.plot([x_min, x_max, x_max, x_min, x_min],[y_min, y_min, y_max, y_max, y_min],'k')
    for grain in L_g:
        if grain.dissolved:
            plt.plot(grain.coordinate_x,grain.coordinate_y,color= 'k',linestyle='-.')
        else :
            plt.plot(grain.coordinate_x,grain.coordinate_y,color= 'k')
    for contact in L_contact+L_contact_gw:
        L_x, L_y, ratio_normal = contact.plot(normal_ref)
        plt.plot(L_x,L_y,linewidth = ratio_normal,color = 'k')
    plt.axis("equal")
    plt.savefig('Debug/DEM_ite/Chain_force_'+str(i_PF)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Compute_Contact_Grain_Distribution(dict_sample):
    """
    Compute the repartition of grain-grain contacts in three categories : dissolvable-dissolvable, undissolvable-dissolvable or undissolvable-undissolvable.

    Compute the repartition of grains two categories : dissolvable or undissolvable.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing, but the dictionnary gets updated values (7 ints)
    """
    #count contacts
    n_contact = 0
    n_contact_diss_diss = 0
    n_contact_undiss_diss = 0
    n_contact_undiss_undiss = 0
    for contact in dict_sample['L_contact'] :
        n_contact = n_contact + 1
        if contact.g1.dissolved and contact.g2.dissolved :
             n_contact_diss_diss = n_contact_diss_diss + 1
        elif (contact.g1.dissolved and not contact.g2.dissolved) or (not contact.g1.dissolved and contact.g2.dissolved) :
             n_contact_undiss_diss = n_contact_undiss_diss + 1
        elif not contact.g1.dissolved and not contact.g2.dissolved :
             n_contact_undiss_undiss = n_contact_undiss_undiss + 1

    #count grains
    n_grain = 0
    n_grain_diss = 0
    n_grain_undiss = 0
    for grain in dict_sample['L_g']:
        n_grain = n_grain + 1
        if grain.dissolved :
            n_grain_diss = n_grain_diss + 1
        else :
            n_grain_undiss = n_grain_undiss + 1

    #update dict
    dict_sample['n_contact'] = n_contact
    dict_sample['n_contact_diss_diss'] = n_contact_diss_diss
    dict_sample['n_contact_undiss_diss'] = n_contact_undiss_diss
    dict_sample['n_contact_undiss_undiss'] = n_contact_undiss_undiss
    dict_sample['n_grain'] = n_grain
    dict_sample['n_grain_diss'] = n_grain_diss
    dict_sample['n_grain_undiss'] = n_grain_undiss

#-------------------------------------------------------------------------------

def Plot_YBoxMax(dict_tracker):
    """
    Plot the evolution of the upper wall position.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))

    plt.plot(dict_tracker['t_L'], dict_tracker['y_box_max_L'])
    plt.title('Upper wall position (µm)')
    plt.xlabel('Time (s)')

    plt.savefig('Debug/UpperWallPosition.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Contact_Distribution(dict_tracker):
    """
    Plot the evolution of the total contact and mean contact distributions.

    The function Compute_Contact_Grain_Distribution() needed to be runned before.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))

    plt.subplot(121)
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['n_contact_diss_diss_L'], label = 'Contact diss-diss')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['n_contact_undiss_diss_L'], label = 'Contact undiss-diss')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['n_contact_undiss_undiss_L'], label = 'Contact undiss-undiss')
    plt.legend()
    plt.title('Total contact distribution')
    plt.xlabel('Percentage of dissolvable grains surface dissolved (%)')

    plt.subplot(122)
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['mean_diss_undiss_L'], label = 'As diss with undiss')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['mean_undiss_diss_L'], label = 'As undiss with diss')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['mean_diss_diss_L'], label = 'As diss with diss')
    plt.plot(dict_tracker['S_dissolved_perc_dissolvable_L'], dict_tracker['mean_undiss_undiss_L'], label = 'As undiss with undiss')
    plt.legend()
    plt.title('Mean contact distribution, following cases')
    plt.xlabel('Percentage of dissolvable grains surface dissolved (%)')

    plt.savefig('Debug/Contact_distribution.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Dimension_Reduction(dict_sample, dict_tracker):
    """
    Plot the surface and dimension reduction.

    A post process is done to quantify the dimension reduction (comparison with initial dimension and current dimension)

        Input :
            a sample dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    Lmean_diss_L = []
    S_dL0_L = []
    S_dLi_L = []
    for i_Surface in range(len(dict_tracker['S_grains_dissolvable_L'])):
        S_g = dict_tracker['S_grains_dissolvable_L'][i_Surface] / dict_sample['n_grain_diss']
        Lmean = math.sqrt(S_g)
        Lmean_diss_L.append(Lmean)
        if i_Surface > 0 :
            S_dL0_L.append((Lmean_diss_L[0]-Lmean)/Lmean_diss_L[0]*100/i_Surface)
            S_dLi_L.append((Lmean_diss_L[i_Surface-1]-Lmean)/Lmean_diss_L[i_Surface]*100)

    #plot
    plt.figure(1,figsize = (16,9))

    plt.subplot(221)
    plt.plot(dict_tracker['S_grains_dissolvable_L'])
    plt.title('Surface')

    plt.subplot(222)
    plt.plot(Lmean_diss_L)
    plt.title('Mean dimension')

    plt.subplot(223)
    plt.plot(S_dL0_L)
    plt.title('Dimension reduction (% of L0)')

    plt.subplot(224)
    plt.plot(S_dLi_L)
    plt.title('Dimension reduction (% of Li-1)')

    plt.savefig('Debug/SizeReduction.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Radius_Reduction(dict_sample, dict_tracker):
    """
    Plot the surface and radius reduction.

    A post process is done to quantify the radius reduction (comparison with initial radius and current radius)

        Input :
            a sample dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    rmean_diss_L = []
    S_dr0_L = []
    S_dri_L = []
    for i_Surface in range(len(dict_tracker['S_grains_dissolvable_L'])):
        S_g = dict_tracker['S_grains_dissolvable_L'][i_Surface] / dict_sample['n_grain_diss']
        Rmean = math.sqrt(S_g/math.pi)
        rmean_diss_L.append(Rmean)
        if i_Surface > 0 :
            S_dr0_L.append((rmean_diss_L[0]-Rmean)/rmean_diss_L[0]*100/i_Surface)
            S_dri_L.append((rmean_diss_L[i_Surface-1]-Rmean)/rmean_diss_L[i_Surface]*100)

    #plot
    plt.figure(1,figsize = (16,9))

    plt.subplot(221)
    plt.plot(dict_tracker['S_grains_dissolvable_L'])
    plt.title('Surface')

    plt.subplot(222)
    plt.plot(rmean_diss_L)
    plt.title('Mean radius')

    plt.subplot(223)
    plt.plot(S_dr0_L)
    plt.title('Radius reduction (% of R0)')

    plt.subplot(224)
    plt.plot(S_dri_L)
    plt.title('Radius reduction (% of Ri-1)')

    plt.savefig('Debug/SizeReduction.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def make_mp4():
    """
    The goal of this function is to create a movie with all configuration pictures.

    From https://www.blog.pythonlibrary.org/2021/06/23/creating-an-animated-gif-with-python/

        Input :
            Nothing
        Output :
            Nothing, but a .mp4 file is generated (a file)
    """
    #look for the largest iteration
    template_name = 'Debug/DEM_ite/PF_ite_'
    i_f = 0
    plotpath = Path(template_name+str(i_f)+'.png')
    while plotpath.exists():
        i_f = i_f + 1
        plotpath = Path(template_name+str(i_f)+'.png')
    fileList = []
    for i in range(0,i_f):
        fileList.append(template_name+str(i)+'.png')
    duration_movie  = 10 #sec
    writer = imageio.get_writer('Debug/PF_ite.mp4', fps=int(i_f/duration_movie))
    for im in fileList:
        writer.append_data(imageio.imread(im))
    writer.close()

    #look for the largest iteration
    template_name = 'Debug/DEM_ite/Chain_force_'
    i_f = 0
    plotpath = Path(template_name+str(i_f)+'.png')
    while plotpath.exists():
        i_f = i_f + 1
        plotpath = Path(template_name+str(i_f)+'.png')
    fileList = []
    for i in range(0,i_f):
        fileList.append(template_name+str(i)+'.png')
    duration_movie  = 10 #sec
    writer = imageio.get_writer('Debug/ChainForce_ite.mp4', fps=int(i_f/duration_movie))
    for im in fileList:
        writer.append_data(imageio.imread(im))
    writer.close()

#-------------------------------------------------------------------------------

def save_DEM_tempo(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker):
    """
    Save trackers and configuration during DEM interations.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save_tempo','wb')
    dict_save = {}
    dict_save['E_cin_stop'] = dict_algorithm['Ecin_stop']
    dict_save['n_window_stop'] = dict_algorithm['n_window_stop']
    dict_save['dy_box_max_stop'] = dict_algorithm['dy_box_max_stop']
    dict_save['dk0_stop'] = dict_algorithm['dk0_stop']
    dict_save['L_g'] = dict_sample['L_g']
    dict_save['L_contact'] = dict_sample['L_contact']
    dict_save['L_ij_contact'] = dict_sample['L_ij_contact']
    dict_save['L_contact_gw'] = dict_sample['L_contact_gw']
    dict_save['L_ij_contact_gw'] = dict_sample['L_ij_contact_gw']
    dict_save['F_on_ymax'] = dict_sollicitations['Force_on_upper_wall']
    dict_save['E_cin'] = dict_tracker['Ecin']
    dict_save['Force'] = dict_tracker['Force_applied']
    dict_save['k0_xmin_tracker'] = dict_tracker['k0_xmin']
    dict_save['k0_xmax_tracker'] = dict_tracker['k0_xmax']
    dict_save['y_box_max'] = dict_tracker['y_box_max'][:-1]
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_DEM_final(dict_algorithm, dict_sample, dict_sollicitations, dict_tracker):
    """
    Save trackers and configuration at the end of DEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    os.remove('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save_tempo')
    outfile = open('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save','wb')
    dict_save = {}
    dict_save['E_cin_stop'] = dict_algorithm['Ecin_stop']
    dict_save['n_window_stop'] = dict_algorithm['n_window_stop']
    dict_save['dy_box_max_stop'] = dict_algorithm['dy_box_max_stop']
    dict_save['dk0_stop'] = dict_algorithm['dk0_stop']
    dict_save['L_g'] = dict_sample['L_g']
    dict_save['L_contact'] = dict_sample['L_contact']
    dict_save['L_ij_contact'] = dict_sample['L_ij_contact']
    dict_save['L_contact_gw'] = dict_sample['L_contact_gw']
    dict_save['L_ij_contact_gw'] = dict_sample['L_ij_contact_gw']
    dict_save['F_on_ymax'] = dict_sollicitations['Force_on_upper_wall']
    dict_save['E_cin'] = dict_tracker['Ecin']
    dict_save['Force'] = dict_tracker['Force_applied']
    dict_save['k0_xmin_tracker'] = dict_tracker['k0_xmin']
    dict_save['k0_xmax_tracker'] = dict_tracker['k0_xmax']
    dict_save['y_box_max'] = dict_tracker['y_box_max'][:-1]
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_save_dicts','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_before_pf(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_save_dicts_before_pf','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_tempo(dict_algorithm,dict_tracker):
    """
    Save trackers and configuration during  PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo','wb')
    dict_save = {}
    dict_save['t_L'] = dict_tracker['t_L']
    dict_save['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict_save['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict_save['S_grains_L'] = dict_tracker['S_grains_L']
    dict_save['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict_save['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict_save['S_grains_dissolvable_L'] = dict_tracker['S_grains_dissolvable_L']
    dict_save['S_dissolved_perc_dissolvable_L'] = dict_tracker['S_dissolved_perc_dissolvable_L']
    dict_save['porosity'] = dict_tracker['porosity_L']
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_final(dict_algorithm,dict_tracker):
    """
    Save trackers and configuration at the end of simulation.

        Input :
            an algorithm dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    os.remove('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo')
    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save','wb')
    dict_save = {}
    dict_save['t_L'] = dict_tracker['t_L']
    dict_save['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict_save['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict_save['S_grains_L'] = dict_tracker['S_grains_L']
    dict_save['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict_save['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict_save['S_grains_dissolvable_L'] = dict_tracker['S_grains_dissolvable_L']
    dict_save['S_dissolved_perc_dissolvable_L'] = dict_tracker['S_dissolved_perc_dissolvable_L']
    dict_save['porosity'] = dict_tracker['porosity_L']
    pickle.dump(dict_save,outfile)
    outfile.close()
