#include "CapacitorLC.h"

// Static members
const unsigned CapacitorLC::n_par = 8;

// Element information
ItemInfo CapacitorLC::einfo =
{
  "capacitorlc",
  "Liquid Cyrstal Nonlinear capacitor",
  "Sonali Luniya",
  DEFAULT_ADDRESS"category:lumped,semiconductor",
  "2003_05_15"
};

// Parameter information
ParmInfo CapacitorLC::pinfo[] =
{
  {"L", "Length (m)", TR_DOUBLE, false},
  {"W", "Width (m)", TR_DOUBLE, false},
  {"D", "Thickness of the LC cell (m)",TR_DOUBLE, false},
  {"DELTA", "Viscosity of liquid crystals (mm^2/s)",TR_DOUBLE, false},
  {"GAMMA", "Fitting Parameter(s/mm^2)",TR_DOUBLE, false},
  {"DTIME"," Delay time at each bias step (s)",TR_DOUBLE, false},
  {"VC", "Threshold voltage (V)",TR_DOUBLE, false},
  {"EPSILON_PL", "Dielectric permittivity ", TR_DOUBLE, false}
};


CapacitorLC::CapacitorLC(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &(L = 152e-6);
  paramvalue[1] = &(W = 148e-6);
  paramvalue[2] = &(D = 10.02e-6);
  paramvalue[3] = &(DELTA = 51.0);
  paramvalue[4] = &(GAMMA = 51.2e-3);
  paramvalue[5] = &(DTIME = 100e-3);
  paramvalue[6] = &(VC = 1.887);
  paramvalue[7] = &(EPSILON_PL = 3.1 );

  // Set the number of terminals
  setNumTerms(2);
  // Set number of states
  setNumberOfStates(1);
  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void CapacitorLC::init() throw(string&)
{
  EPSILON_0 = 8.854e-12;
  DenseIntVector var(1);
  var[0]=0;
  DenseIntVector dvar(1);
  dvar[0]=0;
  DenseIntVector dvar2;
  DenseIntVector tvar;
  DenseDoubleVector nodelay;
  initializeAD(var,dvar,dvar2,tvar,nodelay);
}

void CapacitorLC::eval(AD * x, AD * effort, AD * flow)
{
  //x[0] = V
  AD epsilon_ps;
  AD clc =0;
  if(x[0]<VC)
    epsilon_ps = EPSILON_PL;
  else
    epsilon_ps = EPSILON_PL + DELTA* GAMMA *exp(DTIME) *sqrt((x[0]/VC)-1.0);
  clc = EPSILON_0 * epsilon_ps * L * W / D;
  effort[0] = x[0];
  flow[0] = clc * x[1];
}

