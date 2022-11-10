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
import math
import Report
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pickle
from pathlib import Path
from datetime import datetime
import random
import imageio
import Contact
import Contact_gw
import Grain

#-------------------------------------------------------------------------------

class Grain_pp:

    def __init__(self, Id, Dissolved, Center, Coordinate_x, Coordinate_y):
        #defining a grain for the postprocess
        #each grain is described by a id (an integer class)
        #                           a center (an array [x,y])
        #                           a list of x-coordinate of border (a list)
        #                           a list of y-coordinate of border (a list)

        self.id = Id
        self.dissolved = Dissolved
        self.center = Center
        self.coordinate_x = Coordinate_x
        self.coordinate_y = Coordinate_y

#-------------------------------------------------------------------------------

class Contact_pp:

    def __init__(self, Id_g1, Id_g2, L_g, Normal):
        #defining a contact grain-grain for the postprocess
        #each contact is described by a grain 1 id (an integer class)
        #                             a grain 2 id (an integer class)
        #                             the list of all grain (list of grain_pp)
        #                             the value of the normal reaction (a float)

        for g in L_g:
            if g.id == Id_g1:
                self.g1 = g
            elif g.id == Id_g2:
                self.g2 = g
        self.normal = Normal

    def plot(self, normal_ref):
        #prepare the chain force plot

        L_x = [self.g1.center[0], self.g2.center[0]]
        L_y = [self.g1.center[1], self.g2.center[1]]
        ratio_normal = self.normal/normal_ref
        return L_x, L_y, ratio_normal

#-------------------------------------------------------------------------------

class Contact_gw_pp:

    def __init__(self, Id_g, L_g, Nature, Limit, Normal):
        #defining a contact grain-wall for the postprocess
        #each contact is described by a grain id (an integer class)
        #                             the list of all grain (list of grain_pp)
        #                             the contact nature (a string)
        #                             the wall coordinate (a float)
        #                             the value of the normal reaction (a float)

        for g in L_g:
            if g.id == Id_g:
                self.g = g
        self.nature = Nature
        self.limit = Limit
        self.normal = Normal

    def plot(self, normal_ref):
        #prepare the chain force plot

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
      #an integer is converted to a float with 3 components

      if j < 10:
          j_str = '00'+str(j)
      elif 10 <= j and j < 100:
          j_str = '0'+str(j)
      else :
          j_str = str(j)
      return j_str

#-------------------------------------------------------------------------------

def Stop_Debug(simulation_report):
      #stop simulation

      simulation_report.write('Stop because after it is not corrected !\n')
      simulation_report.end(datetime.now())
      raise ValueError('Stop because after it is not corrected !')

#-------------------------------------------------------------------------------

def Debug_DEM_f(dict_algorithm, dict_sample):
    #plot the configuration of the grains during a DEM step
    #only for debug

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
  #plot the configuration of the grains after the DEM step

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
    #plot the trakers used during DEM simulation
    #only for debug

    #Trackers
    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['t_L'],dict_tracker['S_grains_L'])
    plt.title('Evolution of the grains surface')
    fig.savefig('Debug/Evolution_Surface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_L'][:-1],dict_tracker['k0_xmin_L'],label='k0 with xmin')
    plt.plot(dict_tracker['S_dissolved_L'][:-1],dict_tracker['k0_xmax_L'],label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Grains surface dissolved (µm2)')
    fig.savefig('Debug/Evolution_k0_with_TotalSurface.png')
    plt.close(1)

    fig = plt.figure(1,figsize=(16,9.12))
    plt.plot(dict_tracker['S_dissolved_perc_L'][:-1],dict_tracker['k0_xmin_L'],label='k0 with xmin')
    plt.plot(dict_tracker['S_dissolved_perc_L'][:-1],dict_tracker['k0_xmax_L'],label='k0 with xmax')
    plt.title('Evolution of the k0')
    plt.xlabel('Percentage of grains surface dissolved (%)')
    fig.savefig('Debug/Evolution_k0_with_percentage_dissolved.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Debug_Trackers_DEM(dict_algorithm,dict_sollicitations,dict_tracker):
    #plot the trakers used during DEM simulation
    #only for debug

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
     #sort files generated by MOOSE to different directories

     #master simulation
     os.rename(name_template+'_out.e','Output/'+name_template+'_out.e')
     os.rename(name_template+'.i','Input/'+name_template+'.i')
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

def Dissolution_Distribution(dict_sample,dict_sollicitations,simulation_report):

    i = 0
    while i < int(dict_sollicitations['frac_dissolved']*len(dict_sample['L_g'])):
        print(i,int(dict_sollicitations['frac_dissolved']*len(dict_sample['L_g'])))
        grain = random.choice(dict_sample['L_g'])
        while grain.dissolved:
            grain = random.choice(dict_sample['L_g'])
        grain.dissolved = True
        i = i + 1

    simulation_report.write_and_print(f"{int(100*i/len(dict_sample['L_g']))} % dissolvable (asked {int(100*dict_sollicitations['frac_dissolved'])})\n",f"{int(100*i/len(dict_sample['L_g']))} % dissolvable (asked {int(100*dict_sollicitations['frac_dissolved'])})")

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

def Control_y_max_NR(dict_sample,dict_sollicitations):
    #Control the upper wall to apply force
    #a Newton-Raphson method is applied

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
    #the upper wall is located as a single contact verify the target value

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

def Debug_Control_y_max_NR(L_g,y_max,mu,eta):
    #recompute the force applied on upper wall after the control step
    #debug function

    F = 0
    for grain in L_g:
        p_y_max = max(grain.l_border_y)
        #grain-wall y_max
        if p_y_max > y_max :
            overlap = p_y_max - y_max
            k = 5*4/3*grain.y/(1-grain.nu*grain.nu)*math.sqrt(grain.r_mean)
            F = F + k*overlap**(3/2)

    return F

#-------------------------------------------------------------------------------

def Compute_k0(dict_sample,dict_sollicitations):

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
      #write an .txt file for MOOSE
      #this file described an homogenous dissolution field

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
    #plot the chain force

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
            L_g.append(Grain_pp(id,dissolved,center,coordinate_x,coordinate_y))
        elif line == '<contact_c>\n':
            Read_one_contact = False
            L_contact.append(Contact_pp(L_id_g[0], L_id_g[1], L_g, normal_reaction+normal_damping))
        elif line == '<contact_w_c>\n':
            Read_one_contact_wall = False
            L_contact_gw.append(Contact_gw_pp(id_g, L_g, type, limit, normal_reaction+normal_damping))

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

def make_mp4():
    #The goal of this function is to create a movie with all configuration pictures
    #from https://www.blog.pythonlibrary.org/2021/06/23/creating-an-animated-gif-with-python/

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

#-------------------------------------------------------------------------------

def save_DEM_tempo(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker):
    #save trackers and configuration during DEM interations
    #with FollowSimulation.py, the user can follow live the simulation

    outfile = open('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save_tempo','wb')
    dict = {}
    dict['E_cin_stop'] = dict_algorithm['Ecin_stop']
    dict['n_window_stop'] = dict_algorithm['n_window_stop']
    dict['dy_box_max_stop'] = dict_algorithm['dy_box_max_stop']
    dict['dk0_stop'] = dict_algorithm['dk0_stop']
    dict['L_g'] = dict_sample['L_g']
    dict['L_contact'] = dict_sample['L_contact']
    dict['L_ij_contact'] = dict_sample['L_ij_contact']
    dict['L_contact_gw'] = dict_sample['L_contact_gw']
    dict['L_ij_contact_gw'] = dict_sample['L_ij_contact_gw']
    dict['F_on_ymax'] = dict_sollicitations['Force_on_upper_wall']
    dict['E_cin'] = dict_tracker['Ecin']
    dict['Force'] = dict_tracker['Force_applied']
    dict['k0_xmin_tracker'] = dict_tracker['k0_xmin']
    dict['k0_xmax_tracker'] = dict_tracker['k0_xmax']
    dict['y_box_max'] = dict_tracker['y_box_max'][:-1]
    pickle.dump(dict,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_DEM_final(dict_algorithm,dict_sample,dict_sollicitations,dict_tracker):
    #save trackers and configuration at the end of DEM iteration

    os.remove('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save_tempo')
    outfile = open('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/save','wb')
    dict = {}
    dict['E_cin_stop'] = dict_algorithm['Ecin_stop']
    dict['n_window_stop'] = dict_algorithm['n_window_stop']
    dict['dy_box_max_stop'] = dict_algorithm['dy_box_max_stop']
    dict['dk0_stop'] = dict_algorithm['dk0_stop']
    dict['L_g'] = dict_sample['L_g']
    dict['L_contact'] = dict_sample['L_contact']
    dict['L_ij_contact'] = dict_sample['L_ij_contact']
    dict['L_contact_gw'] = dict_sample['L_contact_gw']
    dict['L_ij_contact_gw'] = dict_sample['L_ij_contact_gw']
    dict['F_on_ymax'] = dict_sollicitations['Force_on_upper_wall']
    dict['E_cin'] = dict_tracker['Ecin']
    dict['Force'] = dict_tracker['Force_applied']
    dict['k0_xmin_tracker'] = dict_tracker['k0_xmin']
    dict['k0_xmax_tracker'] = dict_tracker['k0_xmax']
    dict['y_box_max'] = dict_tracker['y_box_max'][:-1]
    pickle.dump(dict,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_tempo(dict_algorithm,dict_tracker):
    #save trackers during PFDEM interations

    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo','wb')
    dict = {}
    dict['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict['S_grains_L'] = dict_tracker['S_grains_L']
    pickle.dump(dict,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_final(dict_algorithm,dict_tracker):
    #save trackers at the end of the simulation

    os.remove('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo')
    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save','wb')
    dict = {}
    dict['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict['S_grains_L'] = dict_tracker['S_grains_L']
    pickle.dump(dict,outfile)
    outfile.close()
