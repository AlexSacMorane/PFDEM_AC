# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to communicate data from MOOSE to python
A moose multiproccessor run is assumed
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def PFtoDEM_Multi(FileToRead,dict_algorithm,dict_material,dict_sample):

    #---------------------------------------------------------------------------
    #Global parameters
    #---------------------------------------------------------------------------

    L_etai_M_etai = []
    for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
        L_etai_M_etai.append(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L'])))) #etai

    id_L = None
    eta_selector_len = len('        <DataArray type="Float64" Name="etai')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    L_L_Work = [] #Debug

    for i_proc in range(dict_algorithm['np_proc']):

        id_etai = 0
        L_Work = [[], #X
                  []] #Y
        for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
            L_Work.append([]) #etai


    #---------------------------------------------------------------------------
    #Reading file
    #---------------------------------------------------------------------------

        f = open(FileToRead+'_'+str(i_proc)+'.vtu','r')
        data = f.read()
        f.close
        lines = data.splitlines()

        #iterations on line
        for line in lines:

            if line[0:eta_selector_len-1] == '        <DataArray type="Float64" Name="eta':
                if line[0:eta_selector_len] == '        <DataArray type="Float64" Name="eta'+str(int(id_etai+1)):
                    id_L = 2 + id_etai

            elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                id_L = 0

            elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey')and id_L != None:
                id_L = None
                id_etai = id_etai + 1

            elif line[0:data_jump_len] == '          ' and id_L != 0 and id_L != None: #Read etai
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))

            elif line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y, Z]
                line = line[data_jump_len:]
                XYZ_temp = []
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        XYZ_temp.append(float(line[c_start:c_end]))
                        if len(XYZ_temp)==3:
                            L_Work[id_L].append(XYZ_temp[0])
                            L_Work[id_L+1].append(XYZ_temp[1])
                            XYZ_temp = []
                        c_start = c_i+1
                XYZ_temp.append(float(line[c_start:]))
                L_Work[id_L].append(XYZ_temp[0])
                L_Work[id_L+1].append(XYZ_temp[1])

    #---------------------------------------------------------------------------
    #Adaptating data
    #---------------------------------------------------------------------------

        for etai in dict_sample['L_etai_undissolved']+dict_sample['L_etai_dissolved']:
            for i in range(len(L_Work[0])):
                #Interpolation method
                L_dy = []
                for y_i in dict_sample['y_L'] :
                    L_dy.append(abs(y_i - L_Work[1][i]))
                L_dx = []
                for x_i in dict_sample['x_L'] :
                    L_dx.append(abs(x_i - L_Work[0][i]))
                L_etai_M_etai[etai.id][L_dy.index(min(L_dy))][L_dx.index(min(L_dx))] = L_Work[2+etai.id][i]

        L_L_Work.append(L_Work)

    #---------------------------------------------------------------------------
    # Update
    #---------------------------------------------------------------------------

    for grain in dict_sample['L_g']:

        #extract a part focused on the grain
        x_extract_min = grain.center[0] - grain.r_max - dict_material['w']
        x_extract_max = grain.center[0] + grain.r_max + dict_material['w']
        y_extract_min = grain.center[1] - grain.r_max - dict_material['w']
        y_extract_max = grain.center[1] + grain.r_max + dict_material['w']

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

        for i_x in range(i_x_min,i_x_max+1):
            for i_y in range(i_y_min,i_y_max+1):
                grain.etai_M[len(dict_sample['y_L'])-1-i_y][i_x] = L_etai_M_etai[grain.id_eta][len(dict_sample['y_L'])-1-i_y][i_x]
