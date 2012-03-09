#include "ThermalInt.h"

// Static members
const unsigned ThermalInt::n_par = 3;

// Element information
ItemInfo ThermalInt::einfo =
{
  "thermalint",
  "Thermal Interface",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:thermal",
  "2001_10_10"
};

// Parameter information
ParmInfo ThermalInt::pinfo[] =
{
  {"ts", "Kirchhoff transformation temperature (K)", TR_DOUBLE, false},
  {"b1", "Thermal conductivity temperature exponent medium 1",
	TR_DOUBLE, false},
  {"b2", "Thermal conductivity temperature exponent medium 2",
	TR_DOUBLE, false}
};

ThermalInt::ThermalInt(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(ts = 300.);
  paramvalue[1] = &(b1 = 1.22);
  paramvalue[2] = &(b2 = 1.22);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}

void ThermalInt::init() throw(string&)
{
  bm_1 = b1 - one;
  bm_2 = b2 - one;

  DenseIntVector var(1);
  initializeAD(var);
}

void ThermalInt::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: Physical temperature - ts
  // x[1]: 100 * Power flow entering through port 1

  AD T = x[0] + ts;
  // effort[0], effort[1]: theta1, theta2
  effort[0] = ts * (b1 - pow(ts/T, bm_1)) / bm_1;
  effort[1] = ts * (b2 - pow(ts/T, bm_2)) / bm_2;
  // The input power is equal to the output power
  flow[0] = 0.01 * x[1];
  flow[1] = -flow[0];
}
