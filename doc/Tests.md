# Description of the test

## From tempo grain to grain (disk)

A tempo grain is created, the function From_LG_tempo_to_usable is used.
A geometric study is done to compute the center, surface, radius.
Those values are compared to parameters given to the tempo grain. Margins are considered because of the Monte-Carlo algorithm.

## From tempo grain to grain (square)

A tempo grain is created, the function From_LG_tempo_to_usable is used.
A geometric study is done to compute the center, surface, dimension.
Those values are compared to parameters given to the tempo grain. Margins are considered because of the Monte-Carlo algorithm.

## Interaction between wall and polyhedral particle (disk)

A grain from a disk is generated in a box. The box is fit at the particle size. DEM iteration to compute a settlement. The grain center is compared to the box center. Margins are considered

## Interaction between wall and polyhedral particle (square)

A grain from a square is generated in a box. The box is fit at the particle size. DEM iteration to compute a settlement. The grain center is compared to the box center. Margins are considered

## Interaction between polyhedral particles (disk)

Two grains are generated. A force is applied to confine the two particle. The motion of the particle is only on the force axis. Compare after the settlement the overlap.

## Interaction between polyhedral particles (square)

Two grains are generated. A force is applied to confine the two particle. The motion of the particle is only on the force axis. Compare after the settlement the overlap.

## Control the upper wall (i)

A grain is generated. The upper wall is not in contact. A confinement force is asked. The function Control_y_max_NR is called. The force applied on the upper wall is compared to the confinement force.

## Control the upper wall (ii)

A grain is generated. The upper wall is in contact. A confinement force is asked. The function Control_y_max_NR is called. The force applied on the upper wall is compared to the confinement force.

## Contact detection between particles

Two particles are generated (in contact, not or far). The functions Update_Neighbourhoods and Grains_Polyhedral_contact_Neighbourhoods are called. Three case mustbe seen : contact detection, no contact detection even if the grains are in the same neighborhood, particle are not in the same neighborhood.

## Contact detection between particle and wall

One particle and a wall are generated (in contact, not or far). The functions Update_wall_Neighbourhoods and Grains_Polyhedral_Wall_contact_Neighbourhoods are called. Three case mustbe seen : contact detection, no contact detection even if the grains are in the same neighborhood, particle are not in the same neighborhood.

## Create initial condition

A initial condition is asked. 

