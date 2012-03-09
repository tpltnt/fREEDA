#include "VccsPoly.h"

// Static members
const unsigned VccsPoly::n_par = 6;

// Element information
ItemInfo VccsPoly::einfo =
{
  "vccspoly",
  "Voltage-controlled current source where Io = k*poly(Vin)",
  "Jim Hall - added time delay (NMK)",
  DEFAULT_ADDRESS"category:source",
  "2003_05_15"
} ;

// Parameter information
ParmInfo VccsPoly::pinfo[] =
{
  {"k",  "gain", TR_DOUBLE, false},
  {"ri", "Input resistance value (Ohms)", TR_DOUBLE, false},
  {"ro", "Output resistance value (Ohms)", TR_DOUBLE, false},
  {"poly_coef", "Coefficients of polynomial",TR_DOUBLE_VECTOR,false},
  {"td", "time delay", TR_DOUBLE, false},
  {"gamma", "saturation factor", TR_DOUBLE, false}
} ;

VccsPoly::VccsPoly(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(k = one);
  paramvalue[1] = &(ri = 50.0);
  paramvalue[2] = &(ro = 50.0);
  paramvalue[3] = &(poly_coef);
  paramvalue[4] = &(td = zero);
  paramvalue[5] = &(gamma = 1.5);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_FREQ_DOMAIN);

  // Set the number of states
  setNumberOfStates(2) ;
}

void VccsPoly::init() throw(string&)
{
  DenseIntVector var(2);    // 2 states, x[0] = Vin, x[1] = Vout
  var[1] = 1;
  DenseIntVector novar;
  DenseIntVector tvar(1);
  DenseDoubleVector tdelay(1);
  tdelay[0] = td;

  var[0] = 0;
  var[1] = 1;

  initializeAD(var, novar, novar, tvar, tdelay) ;
}


void VccsPoly::eval(AD * x, AD * effort, AD * flow)
{
  AD i_g;
  // Assign state variables to Vin and Vout
  effort[0] = x[0];
  effort[1] = x[1];

  // Evaluate input current
  if (ri)
  {
    flow[0] = x[0] / ri;
  }
  else
  {
    flow[0] = 0;
  }

  // Evaluate generator current
  if(!isSet(&poly_coef))
  {
    // No polynomial is specified, use Ig = k * Vin
    i_g = k * x[2] * tanh(gamma * effort[1]); // k * x[0];
  }
  else
  {
    // Evaluate ig = p0 + vin * (p1 + vin *  (p2 + .... (vin * p(N - 1)))
    i_g = poly_coef[poly_coef.length() - 1];
    for (unsigned i = (poly_coef.length() - 1); i > 0; i--)
    {
      i_g *= x[2] * tanh(gamma * effort[1]); //x[0];
      i_g += poly_coef[i - 1];
    }

    i_g *= k;
  }

  // Evaluate output current
  if (ro)
  {
    flow[1] = - (i_g - (x[1] / ro));
  }
  else
  {
    flow[1] = -i_g;
  }
}
