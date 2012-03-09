#include "AbmVCann.h"

// Static members
const unsigned AbmVCann::n_par = 4;

// Element information
ItemInfo AbmVCann::einfo =
{
  "abmvcann",
  "Cann Behavorial Model",
  "Minsheng Li",
  DEFAULT_ADDRESS"category:behavioral",
  "2005_04_10"
};

// Parameter information
ParmInfo AbmVCann::pinfo[] =
{
  {"g", "Element parameter", TR_DOUBLE, true},
  {"l", "Element parameter", TR_DOUBLE, true},
  {"s", "Element parameter", TR_DOUBLE, true},
  {"rout", "Element parameter", TR_DOUBLE, true}
};

AbmVCann::AbmVCann(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &g;
  paramvalue[1] = &l;
  paramvalue[2] = &s;
  paramvalue[3] = &rout;

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}

void AbmVCann::init() throw(string&)
{
  DenseIntVector var(2);
  var[1] = 1;
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}


void AbmVCann::eval(AD * x, AD * effort, AD * flow)
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
  if(x[0]>=0)
  		flow[1] = (g*x[0]/pow((1+pow((g/l*x[0]),s)),1/s))/rout;
  else
 		flow[1] = (g*x[0]/pow((1+pow((-g/l*x[0]),s)),1/s))/rout;
}
