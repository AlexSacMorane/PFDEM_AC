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

def Create_i_AC(i,dt,x_L,y_L,M_pf,kc_pf,n_t_PF_moose,double_well_height,L_etai_undissolved,L_etai_dissolved):

  file_to_write = open('PF_'+str(i)+'.i','w')
  file_to_read = open('PF_base_AC.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j+1
    if j == 4:
      line = line[:-1] + ' ' + str(len(x_L)-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(y_L)-1)+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(min(x_L))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(max(x_L))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(min(y_L))+'\n'
    elif j == 10:
      line = line[:-1] + ' ' + str(max(y_L))+'\n'
    elif j == 21:
      line = line[:-1]
      for etai in L_etai_undissolved+L_etai_dissolved :
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
      for etai in L_etai_undissolved+L_etai_dissolved :
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
        line = line[:-1] + "'"+ str(M_pf)+' '+ str(kc_pf)+"'"+'\n'
    elif j == 39:
        line = line[:-1]+"'"
        for etai in L_etai_undissolved+L_etai_dissolved:
            line = line +'eta'+str(etai.id+1)+' '
        line = line[:-1]+"'\n"
    elif j == 41:
      line = line[:-1] + ' ' +str(double_well_height)+'\n'
    elif j == 43:
        line = line[:-1]
        for etai in L_etai_undissolved+L_etai_dissolved:
            line = line + 'eta'+str(etai.id+1)+'^2*(1-eta'+str(etai.id+1)+')^2+'
        line = line[:-1]+')+e_dissolution*('
        for etai in L_etai_dissolved:
            line = line + 'eta'+str(etai.id+1)+'^2*(3-2*eta'+str(etai.id+1)+')+'
        line = line[:-1]+')\n'
    elif j ==  56:
        line = line[:-1]
        for etai in L_etai_undissolved+L_etai_dissolved:
            line = line +'\t[eta'+str(etai.id+1)+'_txt]\n'+\
                         '\t\ttype = PiecewiseMultilinear\n'+\
                         '\t\tdata_file = Data/eta'+str(etai.id+1)+'_'+str(i)+'.txt\n'+\
                         '\t[]\n'
        line = line +'\t[e_dissolution_txt]\n'+\
                      '\t\ttype = PiecewiseMultilinear\n'+\
                      '\t\tdata_file = Data/e_dissolution.txt\n'+\
                      '\t[]\n'
    elif j == 80:
      line = line[:-1] + ' ' + str(n_t_PF_moose*dt) + '\n'
    elif j == 84:
            line = line[:-1] + ' ' + str(dt) + '\n'
    file_to_write.write(line)

  file_to_write.close()
