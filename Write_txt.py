# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file allow to create a .txt.
The .txt contains data about grains, contact grain-wall and contact wall-wall.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import Grain
import Contact
import Contact_gw

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def Write_txt(dict_algorithm,dict_sample):
    """
    Write a .txt file to give information about grains and contacts.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but a .txt file is generated (a .txt)
    """
    file_to_write = open('Debug/DEM_ite/PF_'+str(dict_algorithm['i_PF'])+'/txt/ite_DEM_'+str(dict_algorithm['i_DEM'])+'.txt','w')
    file_to_write.write('Iteration PF : '+str(dict_algorithm['i_PF'])+'\n'+\
                        'Iteration DEM : '+str(dict_algorithm['i_DEM'])+'\n')
    file_to_write.write('\n')
    file_to_write.write('GRAINS LIST\n')
    file_to_write.write('\n')
    for grain in dict_sample['L_g']:
        file_to_write.write('<grain_o>\n')
        file_to_write.write('\tid : '+str(grain.id)+'\n')
        if grain.dissolved :
            file_to_write.write('\tDissolved : True\n')
        else:
            file_to_write.write('\tDissolved : False\n')
        file_to_write.write('\tCenter : ['+str(int(grain.center[0]))+', '+str(int(grain.center[1]))+']\n')
        file_to_write.write('\tSpeed : ['+str(round(grain.v[0],2))+', '+str(round(grain.v[1],2))+']\n')
        file_to_write.write('\tForce applied : ['+str(round(grain.f[0],1))+', '+str(round(grain.f[1],1))+']\n')
        file_to_write.write('\tOmega : '+str(grain.w)+'\n')
        file_to_write.write('\tSurface : '+str(int(grain.surface))+'\n')
        file_to_write.write('\tCoordinate X of the border : '+str(grain.l_border_x)+'\n')
        file_to_write.write('\tCoordinate Y of the border : '+str(grain.l_border_y)+'\n')
        file_to_write.write('<grain_c>\n')
    file_to_write.write('\n')
    file_to_write.write('CONTACTS LIST\n')
    file_to_write.write('\n')
    for contact in dict_sample['L_contact']:
        file_to_write.write('<contact_o>\n')
        file_to_write.write('\tid : '+str(contact.id)+'\n')
        file_to_write.write('\tGrains : '+str(contact.g1.id)+'-'+str(contact.g2.id)+'\n')
        file_to_write.write('\tNormal : ['+str(round(contact.pc_normal[0],2))+', '+str(round(contact.pc_normal[1],2))+']\n')
        file_to_write.write('\tNormal overlap : '+str(round(contact.overlap_normal,2))+'\n')
        file_to_write.write('\tNormal reaction : '+str(int(contact.F_2_1_n))+'\n')
        file_to_write.write('\tNormal damping : '+str(int(contact.F_2_1_damp))+'\n')
        file_to_write.write('\tTangential : ['+str(round(contact.pc_tangential[0],2))+', '+str(round(contact.pc_tangential[1],2))+']\n')
        file_to_write.write('\tTangential overlap : '+str(round(contact.overlap_tangential,2))+'\n')
        file_to_write.write('\tTangential reaction : '+str(int(contact.ft))+'\n')
        file_to_write.write('\tTangential damping : '+str(int(contact.ft_damp))+'\n')
        file_to_write.write('<contact_c>\n')
    file_to_write.write('\n')
    file_to_write.write('CONTACTS WITH WALL LIST\n')
    file_to_write.write('\n')
    for contact in dict_sample['L_contact_gw']:
        file_to_write.write('<contact_w_o>\n')
        file_to_write.write('\tid : '+str(contact.id)+'\n')
        file_to_write.write('\tType : '+str(contact.nature)+'\n')
        file_to_write.write('\tGrain : '+str(contact.g.id)+'\n')
        file_to_write.write('\tLimit : '+str(contact.limit)+'\n')
        file_to_write.write('\tNormal : ['+str(round(contact.nwg[0],2))+', '+str(round(contact.nwg[1],2))+']\n')
        file_to_write.write('\tNormal overlap : '+str(round(contact.overlap,2))+'\n')
        file_to_write.write('\tNormal reaction : '+str(int(contact.Fwg_n))+'\n')
        file_to_write.write('\tNormal damping : '+str(int(contact.Fwg_damp_n))+'\n')
        file_to_write.write('\tTangential : ['+str(round(contact.twg[0],2))+', '+str(round(contact.twg[1],2))+']\n')
        file_to_write.write('\tTangential overlap : '+str(round(contact.overlap_tangential,2))+'\n')
        file_to_write.write('\tTangential reaction : '+str(int(contact.ft))+'\n')
        file_to_write.write('<contact_w_c>\n')
    file_to_write.close()
