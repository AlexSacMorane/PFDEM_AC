# PF documentation

## Base

The material is described by a variable &eta;. As explained on the website [wikipedia](https://en.wikipedia.org/wiki/Phase-field_model), the variable is equal to 1 if it is inside the grain or 0 else. A smooth interface is between those two phases.

Grains are sorted in several unconserved parameters &eta; in way that two grains at a distance less than a criteria have not the same parameter. A same variable describes several grains. The distance criteria is here to ensure there is no interferences.

## An Allen-Cahn equation

Description of the formulation in coming...

### The local free energy

The AC formulation asks for a local free energy density function. In the case of this project, a double-well function is considered : f<sub>loc</sub> = 16&times;h&times;(&Sigma;&eta;<sup>2</sup><sub>i</sub>(1-&eta;<sub>i</sub>)<sup>2</sup>).

### Include perturbations

The change of the phase-field is due to some external energy source E<sub>d</sub>. In the case of the simulation, this source term represents the acid injection (with a tendency to reduce values of the variable &eta; inside the interface).

E<sub>d</sub> = e_dissolution&times;&Sigma;(&eta;<sub>i</sub><sup>2<sup>&times;(3-2&times;&eta;<sub>i</sub>)
