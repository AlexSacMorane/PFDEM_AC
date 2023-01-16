# PFDEM_AC

Phase-Field Discrete Element Modelisation applied to granular material.

## Simulation presentation
The goal is to reproduce result from <i>Shin, H., Santamarina, J.C. (2009) Mineral Dissolution and the Evolution of k0, Journal of Geotechnical and Geoenvironmental Engineering, 135, 1141-1147</i>.
This simulation is about a 2D test under oedometrical conditions with an acid injection to dissolve a fraction of the particles.

![scheme of grain-grain interaction](image/Oedometric_test_WB.png)

The chain force and the k0 (= &sigma;<sub>1</sub>/&sigma;<sub>2</sub>) evolves with the material dissolution. The influence of the shape of the dissolvable particles is discussed.

## Technical description
The PFDEM is a discussion between phase-field (PF) simulation and discrete element modelization (DEM). Even if those formulations are different, they can collaborate with data exchange.
![scheme of grain-grain interaction](image/PFDEM_Exchange_Scheme_WB.png)


A back and forth is done between phase-field (PF) simulation and discrete element modelisation (DEM). The following schemes illustrate the algorithm.


![scheme of grain-grain interaction](image/General_plan_WB.png)

 PF approach is computed with [MOOSE](https://github.com/idaholab/moose) and DEM is computed with an own script based on python 3.9.6. The launch of this project is done by executing the python script <i>main.py</i>. The user can change variables inside the python script <i>User.py</i>. A documentation explains what is the role of the different variable.

### PF
The goal of the PF step is to compute the shape of the grains as it evolves with the dissolution.
An Allen-Cahn formulation is used.

### DEM
The goal of the DEM step is to compute a steady-state configuration. Grains are assumed as polygonal. The shape is computed from an interpolation of the PF variable. At the end of the step, interactions (grain-grain and grain-wall) are computed to highlight the evolution of the chain force and k0.

## What is next ?
Some next releases are presented here.

- optimize DEM (rework on the nearest vertices detection)
- add rolling resistance
- increase the user experience (homogenize names, correct english)
- create a database of initial condition
- create AI to select algorithm parameters (as neighborhood update frequency, time step, stop criteria, etc...)n
- add some quick tests to check the script
- correct the functions P_is_inside() because if the coordinate y = the coordinate of the vertex, the point is considered outside (even  if it is not the case)

test
