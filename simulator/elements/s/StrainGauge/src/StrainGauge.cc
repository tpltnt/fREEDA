#include "StrainGauge.h"

const unsigned StrainGauge :: n_par = 6;

ItemInfo StrainGauge::einfo =
{
  "straing",
  "Nonlinear strain guage",
  "Jonathan Cantor and Richard McMunn",
  "category:mechanical"
};

ParmInfo StrainGauge :: pinfo[] =
{
  {"h1","length in y-direction (m)", TR_DOUBLE, true},
  {"h2","length in x-direction (m)", TR_DOUBLE, true},
  {"h3","length in z-direction (m)", TR_DOUBLE, true},
  {"c","Elastic constant", TR_DOUBLE, true},
  {"e","Piezoelectric coupling constant (C/m^2)", TR_DOUBLE, true},
  {"eps","dielectric constant (F/m)", TR_DOUBLE, true}
};

StrainGauge::StrainGauge(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &(h1 = 0.0);
  paramvalue[1] = &(h2 = 0.0);
  paramvalue[2] = &(h3 = 0.0);
  paramvalue[3] = &(C = 0.0);
  paramvalue[4] = &(e = 0.0);
  paramvalue[5] = &(eps = 0.0);
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  setNumberOfStates(2);
  // Set the number of terminals
  setNumTerms(3);
}

void StrainGauge::init() throw(string&)
{
  DenseIntVector var1(2);
  var1[0] = 0;
  var1[1] = 1;
  DenseIntVector var2(2);
  initializeAD(var1, var2);
}

void StrainGauge::eval(AD * x, AD * effort, AD * flow)
{
  // state variables S and E3
  // x[0] : S
  // x[1] : E3
  // x[2] : dS/dt
  // x[3] : dE3/dt
  // output will be vp and ip
  // effort[0] : S
  // flow[0] : F
  // effort[1] : vout
  // flow[1] : iout

  // force related to state variables
  flow[0] = (C * x[0] - e * x[1])* h1 * h3;
  // F = (C*S - e*E3)*h1*h3

  // Current
  flow[1] = -h1*h2*(e * x[2] + eps * x[3]);
  // I = -h1*h2*(e*(dS/dt) + eps*(dE3/dt)

  // Voltage
  effort[1] = -h3*(x[1]);
  // V = -h3*E3
  //voltage is the electric field times the z-length
}
