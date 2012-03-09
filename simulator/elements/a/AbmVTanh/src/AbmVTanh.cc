#include "AbmVTanh.h"

// Static members
const unsigned AbmVTanh::n_par = 3;

// Element information
ItemInfo AbmVTanh::einfo =
{
  "abmvtanh",
  "Tanh behavioral Model",
  "Minsheng Li",
  DEFAULT_ADDRESS"category:behavioral",
  "2005_04_10"
};

// Parameter information
ParmInfo AbmVTanh::pinfo[] =
{
  {"g", "Element parameter", TR_DOUBLE, true},
  {"l", "Element parameter", TR_DOUBLE, true},
  {"rout", "Element parameter", TR_DOUBLE, true}
};

AbmVTanh::AbmVTanh(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &g;
  paramvalue[1] = &l;
  paramvalue[2] = &rout;

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}

void AbmVTanh::init() throw(string&)
{
  DenseIntVector var(2);
  var[1] = 1;
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}

void AbmVTanh::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vin
  // x[1]: vout

  //  Assign voltages directly.
  effort[0] = x[0];
  effort[1] = x[1];

  // Use small signal model.
  // Set input current to zero.
  flow[0] = zero;

  // Equation for output current
  flow[1] = l*tanh(g/l*x[0]) / rout;
}
