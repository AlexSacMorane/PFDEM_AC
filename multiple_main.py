# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is to run several times the main file.
"""

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

N_iterations = 10

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

for i in range(N_iterations):
    exec(open("main.py").read())
