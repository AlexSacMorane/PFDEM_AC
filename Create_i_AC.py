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

import Grain

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def Create_i_AC(dict_algorithm, dict_material, dict_sample, dict_sollicitations):
    """
    The goal of this file is to write a new .i file from PF_base_AC.i.

    There is a global input file for every grains.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            a file is generated (a .i file)
    """
    #Create .txt files
    Owntools.Write_txt_data(dict_algorithm, dict_material, dict_sample)

    # create .i
    file_to_write = open(f"PF_{dict_algorithm['i_PF']}.i",'w')
    file_to_read = open('PF_base_AC.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines :
        j = j+1
        if j == 4:
            line = line[:-1] + ' ' + str(len(grain.x_L_local)-1)+'\n'
        elif j == 5:
            line = line[:-1] + ' ' + str(len(grain.y_L_local)-1)+'\n'
        elif j == 7:
            line = line[:-1] + ' ' + str(min(grain.x_L_local))+'\n'
        elif j == 8:
            line = line[:-1] + ' ' + str(max(grain.x_L_local))+'\n'
        elif j == 9:
            line = line[:-1] + ' ' + str(min(grain.y_L_local))+'\n'
        elif j == 10:
            line = line[:-1] + ' ' + str(max(grain.y_L_local))+'\n'
        elif j == 21:
            line = line[:-1] + '\t[./etai]\n'+\
                              '\t\torder = FIRST\n'+\
                              '\t\tfamily = LAGRANGE\n'+\
                              '\t\t[./InitialCondition]\n'+\
                              '\t\t\ttype = FunctionIC\n'+\
                              '\t\t\tfunction = etai_txt\n'+\
                              '\t\t[../]\n'+\
                              '\t[../]\n'
        elif j == 25:
            line = line[:-1] + '\t[./eta_res]\n'+\
                               '\t\ttype = AllenCahn\n'+\
                               '\t\tvariable = etai\n'+\
                               '\t\tf_name = Fc\n'+\
                               '\t\tmob_name = L\n'+\
                               '\t[../]\n'+\
                               '\t[./eta_int_res]\n'+\
                               '\t\ttype = ACInterface\n'+\
                               '\t\tvariable = etai\n'+\
                               '\t\tmob_name = L\n'+\
                               '\t\tkappa_name = kappa_c\n'+\
                               '\t[../]\n'+\
                               '\t[./time]\n'+\
                               '\t\ttype = TimeDerivative\n'+\
                               '\t\tvariable = etai\n'+\
                               '\t[../]\n'
        elif j == 32:
            line = line[:-1] + "'"+ str(dict_material['M_pf'])+' '+ str(dict_material['kc_pf'])+"'"+'\n'
        elif j == 39:
            line = line[:-1]+"'etai'\n"
        elif j == 41:
            line = line[:-1] + ' ' +str(dict_material['double_well_height'])+'\n'
        elif j == 43:
            line = line[:-1]+ 'etai^2*(1-etai)^2)+e_dissolution*(etai^2*(3-2*etai))\n'
        elif j ==  56:
            line = line[:-1] +'\t[etai_txt]\n'+\
                              '\t\ttype = PiecewiseMultilinear\n'+\
                              '\t\tdata_file = Data/g'+str(grain.id)+'_'+str(dict_algorithm['i_PF'])+'.txt\n'+\
                              '\t[]\n'+\
                              '\t[e_dissolution_txt]\n'+\
                              '\t\ttype = PiecewiseMultilinear\n'+\
                              '\t\tdata_file = Data/e_diss_g'+str(grain.id)+'_ite'+str(dict_algorithm['i_PF'])+'.txt\n'+\
                              '\t[]\n'
        elif j == 80:
            line = line[:-1] + ' ' + str(dict_algorithm['n_t_PF']*dict_algorithm['dt_PF']) + '\n'
        elif j == 84:
            line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']) + '\n'
        file_to_write.write(line)

    file_to_write.close()
