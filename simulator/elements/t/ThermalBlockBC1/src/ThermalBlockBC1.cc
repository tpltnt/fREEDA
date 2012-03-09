#include "ThermalBlockBC1.h"

// Static members
const unsigned ThermalBlockBC1::n_par = 4;
const double ThermalBlockBC1::sigma = 5.67051e-8;

// Element information
ItemInfo ThermalBlockBC1::einfo =
{
  "thermalblockbc1",
  "Thermal nonlinear boundary condition element",
  "Carlos E. Christoffersen, Bill Batty",
  DEFAULT_ADDRESS"category:thermal",
  "2001_10_10"
};

// Parameter information
ParmInfo ThermalBlockBC1::pinfo[] =
{
  {"epsilon", "Surface emissivity", TR_DOUBLE, false},
  {"hn", "Natural convection heat transfer coefficient", TR_DOUBLE, false},
  {"hf", "Forced convection heat transfer coefficient", TR_DOUBLE, false},
  {"tambient", "Ambient temperature (K)", TR_DOUBLE, false}
};

ThermalBlockBC1::ThermalBlockBC1(const string& iname)
: ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(epsilon = 0.7);
  paramvalue[1] = &(hn = zero);
  paramvalue[2] = &(hf = zero);
  paramvalue[3] = &(T0 = 300.);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(1);
}

void ThermalBlockBC1::init() throw(string&)
{
  if (!isSet(&hf))
    hf = 1.3*(1.0+3.)*4.0*T0*T0*T0*sigma*epsilon;
  DenseIntVector var(1);
  initializeAD(var);
}

void ThermalBlockBC1::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: (T/T_0)^4 - 1
  effort[0] = T0 * pow(x[0] + one, .25);
  flow[0] = T0*T0*T0*T0 * 0.0002777777 * epsilon*sigma* x[0];
  if (hn)
    flow[0] += hn * pow(effort[0]-T0,1.25);
  if (hf)
    flow[0] += hf * (effort[0]-T0);
}
