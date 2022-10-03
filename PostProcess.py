# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to postprocces data once the simulation is done.
Some postprocces function are already run during simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
import os
import shutil
from datetime import datetime
from pathlib import Path

#-------------------------------------------------------------------------------

class Grain_pp:

    def __init__(self, Id, Center, Coordinate_x, Coordinate_y):
        #defining a grain for the postprocess
        #each grain is described by a id (an integer class)
        #                           a center (an array [x,y])
        #                           a list of x-coordinate of border (a list)
        #                           a list of y-coordinate of border (a list)

        self.id = Id
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

def contact_g_w_pp(id_g,nature_contact,t_start,t_end,dt,i_PF):
    #plot history of overlap, normal reaction and normal damping for a grain-wall contact
    #id_g_grain is the id of the grain (a int)
    #nature_contact is the id of the wall (a str)

    L_overlap = []
    L_normal_reaction = []
    L_normal_damping = []

    t = t_start

    while t <= t_end:

        file_name = 'Debug/DEM_ite/PF_'+str(i_PF)+'/txt/ite_DEM_'+str(t)+'.txt'
        file = open(file_name,'r')
        lines = file.readlines()
        file.close()

        Read_one_contact = False
        Read_good_contact = False
        Read_good_contact_done = False
        for line in lines :

            if line == '<contact_w_c>\n':
                Read_one_contact = False
                Read_good_contact = False

            if Read_one_contact :
                if Read_good_contact and line[:len('\tNormal overlap : ')] == '\tNormal overlap : ':
                    Read_data = False
                    for c in range(len('\tNormal overlap : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    L_overlap.append(data)

                elif Read_good_contact and line[:len('\tNormal reaction : ')] == '\tNormal reaction : ':
                    Read_data = False
                    for c in range(len('\tNormal reaction : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    L_normal_reaction.append(data)

                elif Read_good_contact and line[:len('\tNormal damping : ')] == '\tNormal damping : ':
                    Read_data = False
                    for c in range(len('\tNormal damping : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    L_normal_damping.append(data)

                if line == '\tType : '+str(nature_contact)+'\n' :
                    Read_good_contact = True

                if line == '\tGrain : '+str(id_g)+'\n' and Read_good_contact:
                    Read_good_contact_done = True

            if line == '<contact_w_o>\n':
                Read_one_contact = True

        if not Read_good_contact_done :
            L_overlap.append(0)
            L_normal_reaction.append(0)
            L_normal_damping.append(0)
        t = t + dt

    plt.figure(1,figsize=(16,9))
    plt.plot(L_overlap)
    plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/Contact_g_'+str(id_g)+'_'+nature_contact+'_overlap.png')
    plt.close(1)

    plt.figure(1,figsize=(16,9))
    plt.plot(L_normal_reaction,label='Reaction')
    plt.plot(L_normal_damping,label='Damping')
    plt.legend()
    plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/Contact_g_'+str(id_g)+'_'+nature_contact+'_force.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def multiple_contact_g_w_pp(L_nature_contact,L_id_g_contact,t_start,t_end,dt,i_PF):
    #same as contact_g_w_pp but for several grain-wall contact

    for i in range(len(L_id_g_contact)):
        id_g_contact = L_id_g_contact[i]
        nature_contact = L_nature_contact[i]
        contact_g_w_pp(id_g_contact,nature_contact,t_start,t_end,dt,i_PF)

#-------------------------------------------------------------------------------

def force_wall_pp(nature_contact,t_start,t_end,dt,i_PF):
    #compute normal force, normal reaction, normal damping and limit value for a wall
    #nature_contact is the id of the wall (a str)

    L_normal_reaction = []
    L_normal_damping = []
    L_normal = []
    L_limit = []

    t = t_start

    while t <= t_end:

        file_name = 'Debug/DEM_ite/PF_'+str(i_PF)+'/txt/ite_DEM_'+str(t)+'.txt'
        file = open(file_name,'r')
        lines = file.readlines()
        file.close()

        normal_reaction = 0
        normal_damping = 0
        Read_one_contact = False
        Read_good_contact = False
        Read_good_contact_done = False
        for line in lines :

            if line == '<contact_w_c>\n':
                Read_one_contact = False
                Read_good_contact = False

            if Read_one_contact :
                if Read_good_contact and line[:len('\tNormal reaction : ')] == '\tNormal reaction : ':
                    Read_data = False
                    for c in range(len('\tNormal reaction : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    normal_reaction = normal_reaction + data

                elif Read_good_contact and line[:len('\tNormal damping : ')] == '\tNormal damping : ':
                    Read_data = False
                    for c in range(len('\tNormal damping : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    normal_damping = normal_damping + data

                elif Read_good_contact and line[:len('\tLimit : ')] == '\tLimit : ':
                    Read_data = False
                    for c in range(len('\tLimit : ')-1,len(line)):
                        if (line[c]!=' ') and not Read_data:
                            c_start = c
                            Read_data = True
                        elif (line[c]==' ' or c==len(line)-1) and Read_data:
                            data = float(line[c_start:c])
                            Read_data = False
                    limit = data

                if line == '\tType : '+str(nature_contact)+'\n' :
                    Read_good_contact = True
                    Read_good_contact_done = True

            if line == '<contact_w_o>\n':
                Read_one_contact = True

        if Read_good_contact_done :
            L_normal_reaction.append(normal_reaction)
            L_normal_damping.append(normal_damping)
            L_normal.append(normal_reaction + normal_damping)
            L_limit.append(limit)
        else :
            L_normal_reaction.append(0)
            L_normal_damping.append(0)
            L_normal.append(0)
            L_limit.append(0)
        t = t + dt

    return L_normal, L_normal_reaction, L_normal_damping, L_limit

#-------------------------------------------------------------------------------

def force_sigma_walls_pp(t_start,t_end,dt,i_PF):
    #plot history of force and pressure for all walls

    L_normal_xmin, L_normal_reaction_xmin, L_normal_damping_xmin, L_limit_xmin = force_wall_pp('gwx_min',t_start,t_end,dt,i_PF)
    L_normal_xmax, L_normal_reaction_xmax, L_normal_damping_xmax, L_limit_xmax = force_wall_pp('gwx_max',t_start,t_end,dt,i_PF)
    L_normal_ymin, L_normal_reaction_ymin, L_normal_damping_ymin, L_limit_ymin = force_wall_pp('gwy_min',t_start,t_end,dt,i_PF)
    L_normal_ymax, L_normal_reaction_ymax, L_normal_damping_ymax, L_limit_ymax = force_wall_pp('gwy_max',t_start,t_end,dt,i_PF)

    L_sigma_xmin = []
    L_sigma_xmax = []
    L_sigma_ymin = []
    L_sigma_ymax = []
    L_percentage_reaction_xmin = []
    L_percentage_reaction_xmax = []
    L_percentage_reaction_ymin = []

    for i in range(len(L_normal_xmin)):
        L_sigma_xmin.append(L_normal_xmin[i]/(L_limit_ymax[i]-L_limit_ymin[i]))
        L_sigma_xmax.append(L_normal_xmax[i]/(L_limit_ymax[i]-L_limit_ymin[i]))
        L_sigma_ymin.append(L_normal_ymin[i]/(L_limit_xmax[i]-L_limit_xmin[i]))
        L_sigma_ymax.append(L_normal_ymax[i]/(L_limit_xmax[i]-L_limit_xmin[i]))

        L_percentage_reaction_xmin.append(L_normal_reaction_xmin[i]/L_normal_xmin[i])
        L_percentage_reaction_xmax.append(L_normal_reaction_xmax[i]/L_normal_xmax[i])
        L_percentage_reaction_ymin.append(L_normal_reaction_ymin[i]/L_normal_ymin[i])

    plt.figure(1,figsize=(16,9))
    plt.plot(L_normal_xmin,label='x min')
    plt.plot(L_normal_xmax,label='x max')
    plt.plot(L_normal_ymin,label='y min')
    plt.plot(L_normal_ymax,label='y max')
    plt.legend()
    plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/Force_on_walls.png')
    plt.close(1)

    plt.figure(1,figsize=(16,9))
    plt.plot(L_sigma_xmin,label='x min')
    plt.plot(L_sigma_xmax,label='x max')
    plt.plot(L_sigma_ymin,label='y min')
    plt.plot(L_sigma_ymax,label='y max')
    plt.legend()
    plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/Pressure_on_walls.png')
    plt.close(1)

    plt.figure(1,figsize=(16,9))
    plt.plot(L_percentage_reaction_xmin,label='x min')
    plt.plot(L_percentage_reaction_xmax,label='x max')
    plt.plot(L_percentage_reaction_ymin,label='y min')
    plt.legend()
    plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/Percentage_reaction_on_walls.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def chain_force_pp(normal_ref,t_start,t_end,dt,i_PF):
    #plot chain force

    t = t_start

    while t <= t_end:

        file_name = 'Debug/DEM_ite/PF_'+str(i_PF)+'/txt/ite_DEM_'+str(t)+'.txt'
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
                L_g.append(Grain_pp(id,center,coordinate_x,coordinate_y))
            elif line == '<contact_c>\n':
                Read_one_contact = False
                L_contact.append(Contact_pp(L_id_g[0],L_id_g[1],L_g, normal_reaction+normal_damping))
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
            plt.plot(grain.coordinate_x,grain.coordinate_y,color= 'k')
        for contact in L_contact+L_contact_gw:
            L_x, L_y, ratio_normal = contact.plot(normal_ref)
            plt.plot(L_x,L_y,linewidth = ratio_normal,color = 'k')
        plt.axis("equal")
        plt.savefig('Debug/DEM_ite/PF_'+str(i_PF)+'/png/Chain_force_ite_'+str(t)+'.png')
        plt.close(1)

        t = t + dt

#-------------------------------------------------------------------------------
# User parameters
#-------------------------------------------------------------------------------

i_PF = 1
t_start = 50
t_end = 1950
dt = 50
normal_ref = 5*10**6

#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------

chain_force_pp(normal_ref,t_start,t_end,dt,i_PF)
force_sigma_walls_pp(t_start,t_end,dt,i_PF)
