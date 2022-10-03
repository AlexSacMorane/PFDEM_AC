# -*- coding: utf-8 -*-
"""
@author: https://www.blog.pythonlibrary.org/2021/06/23/creating-an-animated-gif-with-python/

The goal of this file is to create a movie with all configuration pictures
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import imageio

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def make_mp4(i_f):

    fileList = []
    for i in range(0,i_f+1):
        fileList.append('Debug/DEM_ite/PF_ite_'+str(i)+'.png')

    duration_movie  = 10 #sec
    writer = imageio.get_writer('Debug/PF_ite.mp4', fps=int(i_f/duration_movie))
    for im in fileList:
        writer.append_data(imageio.imread(im))
    writer.close()
