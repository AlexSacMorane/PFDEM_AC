# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the unconserved variables used in phase-field simulation.
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import Grain
import random

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Etai:

#-------------------------------------------------------------------------------

  def __init__(self,ID,L_ig):
      #defining the eta
      #each eta is described by a id (an integer class)
      #                         a list of grain (a list)

      self.id = ID
      self.l_ig = L_ig

#-------------------------------------------------------------------------------

  def update_etai_M(self,L_g):
      #sum variable field of all the grain assigned to this eta

      self.etai_M = np.array(L_g[self.l_ig[0]].etai_M.copy())
      for i in range(1,len(self.l_ig)):
          self.etai_M = self.etai_M + np.array(L_g[self.l_ig[i]].etai_M.copy())

#-------------------------------------------------------------------------------

  def __str__(self):
      pass

#-------------------------------------------------------------------------------

  def Write_txt_Decons_rebuild(self,i_PF,x_L,y_L):
      #write a .txt file
      #this file is used to define initial condition of MOOSE simulation

      file_to_write = open('Data/eta'+str(self.id+1)+'_'+str(i_PF)+'.txt','w')
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
              file_to_write.write(str(self.etai_M[l][c])+'\n')

      file_to_write.close()

#-------------------------------------------------------------------------------

  def Write_txt_Interpolation(self):
    pass

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def etai_distribution(dict_algorithm,dict_sample):
    #Assign grains to eta
    #a minimal distance between grains with same eta is computed from the factor variable

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    factor = dict_algorithm['factor_distribution_etai']
    L_g = dict_sample['L_g']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #initialisation
    for grain in L_g:
        grain.ig_near_L = []
        grain.eta_near_L = []
    #look for grains near
    for i_grain1 in range(1,len(L_g)):
        grain1 = L_g[i_grain1]
        for i_grain2 in range(0,i_grain1):
            grain2 = L_g[i_grain2]
            if np.linalg.norm(grain1.center - grain2.center) < 1.5 * (grain1.r_max+grain2.r_max):
                grain1.ig_near_L.append(i_grain2)
                grain2.ig_near_L.append(i_grain1)
    #first try
    n_etai_L = [1]
    L_etai = [0]
    L_ig_etai = [[0]]
    L_g[0].id_eta = 0
    for grain in L_g:
        if 0 in grain.ig_near_L:
            grain.eta_near_L.append(0)

    for i_grain in range(1,len(L_g)):
        grain = L_g[i_grain]
        etai_defined = False
        etai = 0
        while etai < len(L_etai) and not etai_defined:
            if etai not in grain.eta_near_L:
                grain.id_eta = etai
                etai_defined = True
                n_etai_L[etai] = n_etai_L[etai] + 1
                L_ig_etai[etai].append(i_grain)
                for grain2 in L_g:
                    if i_grain in grain2.ig_near_L:
                        grain2.eta_near_L.append(etai)
            etai = etai + 1
        if etai == len(L_etai) and not etai_defined:
            grain.id_eta = etai
            L_etai.append(etai)
            n_etai_L.append(1)
            L_ig_etai.append([i_grain])
            for grain2 in L_g:
                if i_grain in grain2.ig_near_L:
                    grain2.eta_near_L.append(etai)

    #adaptation (try to have the same number of grain assigned to all eta)
    n_etai_mean = np.mean(n_etai_L)
    adaptation_done = False
    adaptation_i = 0
    while not adaptation_done :
        adaptation_i = adaptation_i + 1

        L_ig_over =  L_ig_etai[n_etai_L.index(max(n_etai_L))]
        id_g_to_work = L_ig_over[random.randint(0,len(L_ig_over)-1)]
        etai_over = n_etai_L.index(max(n_etai_L))
        etai_under = n_etai_L.index(min(n_etai_L))

        grain = L_g[id_g_to_work]
        if etai_under not in grain.eta_near_L:
            grain.id_eta = etai_under
            n_etai_L[etai_over] = n_etai_L[etai_over] - 1
            n_etai_L[etai_under] = n_etai_L[etai_under] + 1
            L_ig_etai[etai_over].remove(id_g_to_work)
            L_ig_etai[etai_under].append(id_g_to_work)
            for grain2 in L_g:
                if id_g_to_work in grain2.ig_near_L:
                    grain2.eta_near_L.remove(etai_over)
                    grain2.eta_near_L.append(etai_under)

        #check the quality
        adaptation_done = True
        for n_etai in n_etai_L :
            if n_etai_mean-2 < n_etai and n_etai < n_etai_mean+2:
                adaptation_done = False
        if adaptation_i > len(L_g):
            adaptation_done = True

    return L_ig_etai

#-------------------------------------------------------------------------------

def etai_distribution_dissolution(dict_algorithm,dict_sample, dict_sollicitations, simulation_report):
    #Assign grains to eta
    #some etas are for dissolved grains, other for undissolved

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    L_g = dict_sample['L_g']
    frac_dissolved = dict_sollicitations['frac_dissolved']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    #run initialy etai_distribution to have a distribution of the etai
    L_ig_etai_undissolved = etai_distribution(dict_algorithm,dict_sample)

    #compute the number of grain dissolved
    n_grain_dissolved = int(frac_dissolved*len(L_g))

    #compute the number of etai needed for L_ig_etai_dissolved
    n_eta = int(frac_dissolved*len(L_ig_etai_undissolved))+1

    #prepare the selection of grain dissolved
    n_try_max = 2*len(L_ig_etai_undissolved[0])
    L_etai_dissolved = []
    L_ig_etai_dissolved = []
    for i in range(n_eta):
        L_etai_dissolved.append(len(L_ig_etai_undissolved)+i)
        L_ig_etai_dissolved.append([])

    #random selection of grain dissolved
    for i in range(n_grain_dissolved):
        etai = L_etai_dissolved[i%n_eta]

        grain_selected = False
        i_try = 0
        while not grain_selected:
            i_try = i_try + 1
            grain_i = random.choice(L_ig_etai_undissolved[i%len(L_ig_etai_undissolved)])
            grain = L_g[grain_i]
            if etai not in grain.eta_near_L:
                grain_selected = True
                grain.id_eta = etai
                grain.dissolved = True
                L_ig_etai_dissolved[i%n_eta].append(grain_i)
                L_ig_etai_undissolved[i%len(L_ig_etai_undissolved)].remove(grain_i)
                for grain2 in L_g:
                    if grain_i in grain2.ig_near_L:
                        grain2.eta_near_L.append(etai)
                        grain2.eta_near_L.remove(i%len(L_ig_etai_undissolved))

    #cntrol the algorithm
    n_dissolved = 0
    for L_ig in L_ig_etai_dissolved:
        n_dissolved = n_dissolved + len(L_ig)
    n_total = n_dissolved
    for L_ig in L_ig_etai_undissolved :
        n_total = n_total + len(L_ig)

    simulation_report.write_and_print('Real fraction dissolved : '+str(int(100*n_dissolved/n_total))+' (asked : '+str(int(100*frac_dissolved))+')\n','Real fraction dissolved : '+str(int(100*n_dissolved/n_total))+' (asked : '+str(int(100*frac_dissolved))+')')
    simulation_report.write('Number of eta for undissolved grain : '+str(len(L_ig_etai_undissolved))+'\n'+\
                            'Number of eta for dissolved grain : '+str(len(L_ig_etai_dissolved))+'\n')

    #add element in dict
    dict_sample['L_ig_etai_undissolved'] = L_ig_etai_undissolved
    dict_sample['L_ig_etai_dissolved'] = L_ig_etai_dissolved
