# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to write a new .i file.
The function use PF_base_AC.i as a template
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
from math import pi

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def Create_i_AC(dict_algorithm, dict_material, dict_sample):

  file_to_write = open('PF_'+str(dict_algorithm['i_PF'])+'.i','w')
  file_to_read = open('PF_base_AC.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j+1
    if j == 4:
      line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
    elif j == 10:
      line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
    elif j == 21:
      line = line[:-1]
      for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved'] :
          line = line + '\t[./eta'+str(etai.id+1)+']\n'+\
                             '\t\torder = FIRST\n'+\
                             '\t\tfamily = LAGRANGE\n'+\
                             '\t\t[./InitialCondition]\n'+\
                             '\t\t\ttype = FunctionIC\n'+\
                             '\t\t\tfunction = eta'+str(etai.id+1)+'_txt\n'+\
                             '\t\t[../]\n'+\
                             '\t[../]\n'
    elif j == 25:
      line = line[:-1]
      for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved'] :
          line = line + '\t[./eta_res_'+str(etai.id+1)+']\n'+\
                        '\t\ttype = AllenCahn\n'+\
                        '\t\tvariable = eta'+str(etai.id+1)+'\n'+\
                        '\t\tf_name = Fc\n'+\
                        '\t\tmob_name = L\n'+\
                        '\t[../]\n'+\
                        '\t[./eta_int_res_'+str(etai.id+1)+']\n'+\
                        '\t\ttype = ACInterface\n'+\
                        '\t\tvariable = eta'+str(etai.id+1)+'\n'+\
                        '\t\tmob_name = L\n'+\
                        '\t\tkappa_name = kappa_c\n'+\
                        '\t[../]\n'+\
                        '\t[./time_'+str(etai.id+1)+']\n'+\
                        '\t\ttype = TimeDerivative\n'+\
                        '\t\tvariable = eta'+str(etai.id+1)+'\n'+\
                        '\t[../]\n'
    elif j == 32:
        line = line[:-1] + "'"+ str(dict_material['M_pf'])+' '+ str(dict_material['kc_pf'])+"'"+'\n'
    elif j == 39:
        line = line[:-1]+"'"
        for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
            line = line +'eta'+str(etai.id+1)+' '
        line = line[:-1]+"'\n"
    elif j == 41:
      line = line[:-1] + ' ' +str(dict_material['double_well_height'])+'\n'
    elif j == 43:
        line = line[:-1]
        for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
            line = line + 'eta'+str(etai.id+1)+'^2*(1-eta'+str(etai.id+1)+')^2+'
        line = line[:-1]+')+e_dissolution*('
        for etai in dict_sample['L_etai_dissolved']:
            line = line + 'eta'+str(etai.id+1)+'^2*(3-2*eta'+str(etai.id+1)+')+'
        line = line[:-1]+')\n'
    elif j ==  56:
        line = line[:-1]
        for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
            line = line +'\t[eta'+str(etai.id+1)+'_txt]\n'+\
                         '\t\ttype = PiecewiseMultilinear\n'+\
                         '\t\tdata_file = Data/eta'+str(etai.id+1)+'_'+str(dict_algorithm['i_PF'])+'.txt\n'+\
                         '\t[]\n'
        line = line +'\t[e_dissolution_txt]\n'+\
                      '\t\ttype = PiecewiseMultilinear\n'+\
                      '\t\tdata_file = Data/e_dissolution.txt\n'+\
                      '\t[]\n'
    elif j == 80:
      line = line[:-1] + ' ' + str(dict_algorithm['n_t_PF']*dict_algorithm['dt_PF']) + '\n'
    elif j == 84:
            line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']) + '\n'
    file_to_write.write(line)

  file_to_write.close()
