# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define value of phase variable.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
from math import pi

#-------------------------------------------------------------------------------
#Function Definition
#-------------------------------------------------------------------------------

def Cosine_Profile(R,w,r):
  #r is the absolute value of the distance to the current point from the center
  if r<R-w/2:
    return 1
  elif r>R+w/2:
    return 0
  else :
    return 0.5*(1 + np.cos(pi*(r-R+w/2)/w))

#-------------------------------------------------------------------------------

def IC(x_L,y_L,R,w,C):

  etai_M_IC = np.zeros((len(y_L),len(x_L)))

  for y in y_L:
    for x in x_L:
      r = np.linalg.norm(np.array([x,y])-C)
      etai_M_IC[len(y_L)-1-y_L.index(y)][x_L.index(x)] = Cosine_Profile(R,w,r)

  return etai_M_IC

#-------------------------------------------------------------------------------
