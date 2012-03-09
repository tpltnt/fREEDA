#About
fREEDA is a multi-physics simulator that uses compact models. It can be 
used as an alternative Spice circuit simulator and will read a Spice
netlist as well as its own (better) native netlist with many more
elemenbts supported. fREEDA uses uses state variables and many new
circuit concepts supporting multi physics modeling. It minimizes energy
rather than minimizing current errors and so is much more robust than
Spice and with high dynamic range. Many phenomena that Spice can never
model can be modeled in fREEDA. Features include easy model development,
high dynamic range, and support for transient, harmonic balance and multi
resolution analysis. fREEDA's large signal noise analysis captures phase
noise and the effects of high levbels of noise. fREEDA also supports the
integration of EM models in transient analysis. There are at least 107
different models including the major transistor models, nonlinear
electro-thermal models capturing long-tail thermal effects.

fREEDA implements several types of analyses. It implements a DC,
Harmonic Balance, several Time marching transient and a unique wavelet
analysis. It also implements several device models including common three
and four terminal transistors, transmission line, Foster's canonical form
and diodes - electronic, optical and tunnel types.

#License
LGPL

#System Requirements
- gcc, g++, g77 compilers ver. 4.0+
- GNU make, flex, bison
- bash shell
