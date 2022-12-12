# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is to run several times the main file.
"""

#-------------------------------------------------------------------------------
#User parameters
#-------------------------------------------------------------------------------

N_iterations = 10

#-------------------------------------------------------------------------------
#main multiple times
#-------------------------------------------------------------------------------

for i in range(N_iterations):
    exec(open("main.py").read())
