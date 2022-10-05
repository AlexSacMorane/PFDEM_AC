# DEM documentation

## What is a grain ?

A grain is a polygonal particle defined initially by the coordinates of the vertices and some material proprieties (Young modulus, Poisson's ratio, surface mass). By apllying a [Monte Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method), the mass, the center of the mass and the inertia of the grain is determined.
A grain can know two kinds of interaction a grain-grain one or a grain-wall one.

## Grain - grain interaction

The contact between two particles is determined by apllying the method developped by <i>Nezami et al.</i>.

Once the contact is determined, Hertz laws illustrated by the following scheme are applied to the contact.
![scheme of grain-grain interaction](../image/DEM_Modele.png)

### Normal behavior

Some equivalent parameters must be defined :
  - an Young modulus 1/Y<sub>eq</sub> = (1-&nu;<sup>2</sup><sub>i</sub>)/Y<sub>i</sub> + (1-&nu;<sup>2</sup><sub>j</sub>)/Y<sub>j</sub>
  - a radius 1/R<sub>eq</sub> = 1/R<sub>i</sub> + 1/R<sub>j</sub>  
  - a mass m<sub>eq</sub> = (m<sub>i</sub>&times;m<sub>j</sub>) / (m<sub>i</sub>+m<sub>j</sub>)

About the spring term :
  The value of the spring is k = 4/3&times; Y<sub>eq</sub>&radic;(R<sub>eq</sub>).
  An Hertz law is applied, considering an unlinear spring : Fs<sub>n</sub> = -k&times;overlap<sub>n</sub><sup>3/2</sup>

About the damping term :
  The value of the dashpot is &eta; = 2 &times; &gamma; &radic;(m<sub>eq</sub>&times;k)
  with &gamma; = -ln(e)/&radic;(&pi;<sup>2</sup>+ln(e)<sup>2</sup>)
  e the restitution coefficient of the contact
  An Hertz law is applied, considering a linear damping : Fd<sub>n</sub> = -&eta;&times;(v<sub>i</sub>-v<sub>j</sub>).n

### Tangential behavior

Some equivalent parameters must be defined :
  - a shear modulus 1/G<sub>eq</sub> = (1-&nu;<sub>i</sub>)/G<sub>i</sub> + (1-&nu;<sub>j</sub>)/G<sub>j</sub>
    with G<sub>i</sub> = Y<sub>i</sub>/(2(1+&nu;<sub>i</sub>))

About the spring term :
  The value of the unlinear spring is kt = kt0&radic;(1-(2&times;kt0&times;overlap<sub>t</sub>)/(3&times;&mu;Fs<sub>n</sub>))
  with kt0 = 8&times;G<sub>eq</sub>&radic;(R<sub>eq</sub>overlap<sub>n</sub>)
  An incremental Hertz law is applied, considering an unlinear spring : Fs<sub>t</sub><sup>t</sup> = Fs<sub>t</sub><sup>t-1</sup> - kt&times;&Delta;overlap<sub>t</sub>
  with &Delta;overlap<sub>t</sub> = (v<sub>i</sub>-v<sub>j</sub>).t &times; dt
  An Coulomb criteria is considering to add sliding : Fs<sub>t</sub><sup>t</sup> &le; &mu; Fs<sub>t</sub><sup>n</sup>

About the damping term :
  The value of the dashpot is &eta; = 2 &times; &gamma; &radic;(m<sub>eq</sub>&times;kt)
  with &gamma; = -ln(e)/&radic;(&pi;<sup>2</sup>+ln(e)<sup>2</sup>)
  e the restitution coefficient of the contact
  An Hertz law is applied, considering a linear damping : Fd<sub>t</sub> = -&eta;&times;(v<sub>i</sub>-v<sub>j</sub>).t

### Rolling behavior

The rolling behavior is not defined for the moment (release in coming...).

## Grain - wall interaction

The contact detection is done with the coordinate of the extremum vertex of the grain :
  - minimum of the coordinate x < coordinate of the left wall -> contact with the left wall
  - maximum of the coordinate x > coordinate of the right wall -> contact with the right wall
  - minimum of the coordinate y < coordinate of the lower wall -> contact with the lower wall
  - maximum of the coordinate y > coordinate of the upper wall -> contact with the upper wall

The compute of the reaction is the same as grain - grain interaction except the following :
  - Y<sub>eq</sub> = Y/(1-&nu;<sup>2</sup>)
  - R<sub>eq</sub> = R
  - m<sub>eq</sub> = m
  - k = factor&times;4/3&times; Y<sub>eq</sub>&radic;(R<sub>eq</sub>) (factor is here to increase the stiffness)
  - there is not a normal damping term for the upper wall
  - there are not tangential damping terms for all walls

## Time integration

The Symplectic Euler method is used.

For the translation due to forces:
a<sup>t-1</sup> = &Sigma;F<sup>t-1</sup>/m
v<sup>t</sup> = v<sup>t-1</sup> + a<sup>t-1</sup> &times; dt
p<sup>t</sup> = p<sup>t-1</sup> + v<sup>t-1</sup> &times; dt

For the rotation due to moments:
d&omega;/dt<sup>t-1</sup> = &Sigma;M<sup>t-1</sup>/I
&omega;<sup>t</sup> = &omega;<sup>t-1</sup> + d&omega;/dt<sup>t-1</sup> &times; dt
&theta;<sup>t</sup> = &theta;<sup>t-1</sup> + &omega;<sup>t-1</sup> &times; dt
p<sup>t</sup> = P<sub>rot</sub> p<sup>t</sup><sub>tempo</sub> + center<sup>t</sup>
P<sub>rot</sub> is a rotation matrix with the angle &omega;&times;dt
p<sup>t</sup><sub>tempo</sub> = p<sup>t</sup> - center<sup>t</sup>

## References

E. Nezami, Y. Hashash, Zhao D., Ghaboussi J., A fast contact detection algorithm for 3-D discrete element method (2004) Computers and Geotechnics, Vol. 31, Pages 575-587, DOI : 10.1016/j.compgeo.2004.08.002

C. O'Sullivan, Particulate Discrete Element Modelling (2011) DOI : 10.1201/9781482266498
