# DEM documentation

## What is a grain ?

A grain is a polygonal particle defined initially by the coordinates of the vertices and some material proprieties (Young modulus, Poisson's ratio, surface mass). By apllying a [Monte Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method), the mass, the center of the mass and the inertia of the grain is determined.
A grain can know two kinds of interaction a grain-grain one or a grain-wall one.

## Grain - grain interaction

The contact between two particles is determined by apllying the method developped by <i>Nezami et al.</i>.

Once the contact is determined, Hertz laws illustrated by the following scheme are applied to the contact.
[scheme of grain-grain interaction](../image/DEM_Modele.png)

The normal behavior is based on:
  - an unlinear spring  
  - a linear dashpot

The tangential behavior is based on:
  - an unlinear spring
  - a linear dashpot
  - a Coulomb criteria to detect the sliding

The rolling behavior is not defined for the moment (release in coming...).

## Grain - wall interaction

## References

E. Nezami, Y. Hashash, Zhao D., Ghaboussi J., A fast contact detection algorithm for 3-D discrete element method (2004) Computers and Geotechnics, Vol. 31, Pages 575-587, DOI : 10.1016/j.compgeo.2004.08.002

C. O'Sullivan, Particulate Discrete Element Modelling (2011) DOI : 10.1201/9781482266498
