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

fREEDA(TM) can be downloaded at http://www.freeda.org

#License
LGPL

#System Requirements
- GNU Make
- GNU C compiler gcc/g++ version 4.4.x 
- GNU fortran compiler gfortran 4.4.x
- GNU flex-2.5.xx and bison-2.4.x
- The bash shell version 4.1.5
- CMake 2.8.x for the Trilinos installation

For ifreeda

- Qt 4.5.0. This includes the qmake tool, and the associated toolkit library.

It is encouraged to use Qt 4.5.0 only. Newer versions of Qt have not been
tested.

#Developers
fREEDA has evolved because of the work of many people. Below is a list of
people, many of whom are still actively involved. Please email Nikhil Kriplani
at nkriplaniHERE-COMES-THE-AT-SIGNncsu.edu to become a developer or if your
name may have been mistakenly left out.

Michael B Steer <mbsHERE-COMES-THE-AT-SIGNncsu.edu>
Carlos Christoffersen <cechristHERE-COMES-THE-AT-SIGNvision.lakeheadu.ca>
Nikhil Kriplani <nmkriplaHERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Senthil Velu <sveluHERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Sonali Luniya <srluniyaHERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Justin Lowry <jqlowryHERE-COMES-THE-AT-SIGNncsu.edu>
Houssam Kanj <houssam_kanjHERE-COMES-THE-AT-SIGNyahoo.com>
Frank Hart <fphartHERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Aaron Walker <alwalke3HERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Ramya Mohan <rmohanHERE-COMES-THE-AT-SIGNunity.ncsu.edu>
Shubha Vijaychand <svijaycHERE-COMES-THE-AT-SIGNhotmail.com>
Satish Uppathil <satish.uppathilHERE-COMES-THE-AT-SIGNericsson.com>
T. Rob Harris <trharrisHERE-COMES-THE-AT-SIGNncsu.edu>
Shivam Priyadarshi <spriyadHERE-COMES-THE-AT-SIGNncsu.edu>
Chris Snowden
Bill Batty
